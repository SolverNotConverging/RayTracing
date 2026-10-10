"""VTK scene inspection before tracing, and branch inspection after tracing."""
from dataclasses import dataclass
import math
from pathlib import Path
import numpy as np
from . import _core
from .scene import Box, Sphere, Cylinder, Disk, Rectangle, Triangle, faces, vector
from .beams import GaussianBeam, PlaneWave, ApertureSource, frame, propagate
from .materials import C0


@dataclass
class ViewOptions:
    off_screen: bool = False
    window_size: tuple = (1200, 850)
    beam_length: float = 1.0
    pattern_scale: float = 0.3
    geometry_opacity: float = 0.25
    beam_opacity: float = 0.25
    resolution: int = 48
    max_preview_beamlets: int = 64
    show_termination: bool = True
    frame_termination: bool = False


@dataclass
class AntennaPattern:
    """Display a native antenna's normalized angular pattern at a position.

    Used for Tx or Rx inspection; for Rx the lobe denotes outward look direction.
    The antenna's own orientation is respected. This is a display descriptor.
    """
    position: object
    antenna: object


def _mesh(pv,shape,resolution):
    if isinstance(shape,Box):
        meshes=[_mesh(pv,f,resolution) for f in faces(shape)]
        return meshes[0].merge(meshes[1:])
    if isinstance(shape,Sphere):
        return pv.Sphere(radius=shape.radius,center=shape.center,theta_resolution=resolution,phi_resolution=resolution)
    if isinstance(shape,Cylinder):
        return pv.Cylinder(center=shape.center,direction=shape.axis,radius=shape.radius,height=2*shape.half_length,
                           resolution=resolution,capping=shape.capped)
    if isinstance(shape,Disk):
        return pv.Disc(center=shape.center,normal=shape.normal,inner=0,outer=shape.radius,c_res=resolution)
    if isinstance(shape,Rectangle):
        u,v=shape.half_width*shape.u,shape.half_height*shape.v
        c=shape.center
        return pv.PolyData(np.array([c-u-v,c+u-v,c+u+v,c-u+v]),faces=[4,0,1,2,3])
    if isinstance(shape,Triangle): return pv.PolyData(np.array([shape.a,shape.b,shape.c]),faces=[3,0,1,2])
    raise TypeError('Unsupported shape')


def _beam_mesh(pv,state,material,frequency,length,resolution):
    angles=np.linspace(0,2*np.pi,resolution+1)
    ring=np.array([np.cos(angles),np.sin(angles)])
    stations=np.linspace(0,length,resolution)
    points=[]
    for z in stations:
        at=propagate(state,z,material,frequency)
        values,vectors=np.linalg.eigh(at.curvature.imag)
        if np.any(values<=0): return None
        k=2*np.pi*frequency/C0*material.index(frequency).real
        contour=vectors@np.diag(np.sqrt(2/(k*values)))@ring
        points.append(at.position[:,None]+at.basis@contour)
    points=np.asarray(points).transpose(1,2,0)
    return pv.StructuredGrid(points[0],points[1],points[2])


def _as_list(value):
    return list(value) if isinstance(value,(tuple,list)) else [] if value is None else [value]


def _add_pattern(pv,plotter,scene,pattern,frequency,options,label,color):
    if isinstance(pattern,AntennaPattern):
        position=vector(pattern.position)
        mesh=pv.Sphere(radius=1,theta_resolution=options.resolution,phi_resolution=options.resolution)
        directions=np.asarray(mesh.points,dtype=float).copy()
        medium=scene.material(scene.region_at(position))
        if medium.kind!='dielectric' or medium.index(frequency).imag!=0:
            raise ValueError('Angular antenna preview requires a lossless dielectric at the antenna')
        native_medium=_core.Medium(medium.epsilon_r.real,medium.mu_r)
        amplitude=np.array([np.linalg.norm(pattern.antenna.farfield(d/np.linalg.norm(d),frequency,native_medium)) for d in directions])
        peak=amplitude.max(initial=0)
        radius=amplitude/peak if peak>0 else amplitude
        mesh.points=position+options.pattern_scale*radius[:,None]*directions
        plotter.add_mesh(mesh,color=color,opacity=.55,name=label+' pattern',label=label+' normalized pattern')
        plotter.add_point_labels([position],[label],point_color=color,text_color=color,font_size=14,always_visible=True,name=label+' label')
        return
    if not isinstance(pattern,(GaussianBeam,PlaneWave,ApertureSource)):
        raise TypeError('Preview expects GaussianBeam, PlaneWave, ApertureSource or AntennaPattern')
    launched=pattern.launch(scene,frequency)
    states=launched if isinstance(launched,list) else [launched]
    # This only decimates rendering, never modifies a source or simulation.
    stride=max(1,math.ceil(len(states)/options.max_preview_beamlets))
    for i,state in enumerate(states[::stride]):
        medium=scene.material(state.region)
        mesh=_beam_mesh(pv,state,medium,frequency,options.beam_length,options.resolution)
        if mesh is not None:
            plotter.add_mesh(mesh,color=color,opacity=options.beam_opacity,name=f'{label} envelope {i}',
                             label=label+' beam' if i==0 else None)
        plotter.add_mesh(pv.Line(state.position,state.position+options.beam_length*state.direction),
                         color=color,line_width=3,name=f'{label} axis {i}')
    state=states[len(states)//2]
    plotter.add_point_labels([state.position],[label],point_color=color,text_color=color,font_size=14,always_visible=True,name=label+' label')
    plotter.add_arrows(state.position[None,:],state.direction[None,:],mag=options.beam_length*.15,color=color,name=label+' direction')
    if isinstance(pattern,GaussianBeam):
        plotter.add_mesh(pv.Sphere(center=pattern.waist_position,radius=options.beam_length*.012),color=color,name=label+' waist')


def _base(scene,transmitters,receivers,frequency_hz,options,show):
    import pyvista as pv
    scene.validate()
    if not math.isfinite(frequency_hz) or frequency_hz<=0: raise ValueError('Frequency must be positive')
    if not math.isfinite(options.beam_length) or options.beam_length<=0 or not math.isfinite(options.pattern_scale) or options.pattern_scale<=0:
        raise ValueError('Preview length and pattern scale must be positive')
    if not isinstance(options.resolution,int) or options.resolution<8 or options.max_preview_beamlets<1:
        raise ValueError('Invalid preview resolution or beamlet limit')
    if not 0<=options.geometry_opacity<=1 or not 0<=options.beam_opacity<=1:
        raise ValueError('Opacity must lie in [0,1]')
    p=pv.Plotter(off_screen=options.off_screen or not show,window_size=options.window_size,notebook=False)
    p.set_background('#f4f6fa')
    legend=[]
    for i,obj in enumerate(scene.objects):
        if obj.terminate:
            if not options.show_termination: continue
            color='#737b85'; label=f'{i}: terminate {obj.terminate}'
            actor=p.add_mesh(_mesh(pv,obj.shape,options.resolution),style='wireframe',color=color,line_width=1,
                            name=f'geometry {i}',label=label)
            actor.use_bounds=options.frame_termination
            legend.append(label)
        else:
            color='#5a99c9' if obj.material.kind=='dielectric' else '#d29b49'
            label=f'{i}: {obj.material.name} [{"V" if obj.volume else "S"}]'
            p.add_mesh(_mesh(pv,obj.shape,options.resolution),color=color,opacity=options.geometry_opacity,
                       show_edges=True,edge_color=color,name=f'geometry {i}',label=label)
            legend.append(label)
    for i,tx in enumerate(_as_list(transmitters)):
        _add_pattern(pv,p,scene,tx,frequency_hz,options,f'Tx {i+1}','#d64b35')
    for i,rx in enumerate(_as_list(receivers)):
        _add_pattern(pv,p,scene,rx,frequency_hz,options,f'Rx {i+1} look','#338550')
    p.add_axes()
    p.add_text(f'{frequency_hz/1e9:g} GHz | Tx red | Rx green\nUnscattered envelopes: 1/e field radius',
               font_size=10,color='#263449',position=(12,options.window_size[1]-55),name='preview caption')
    if legend:
        p.add_text('\n'.join(legend),font_size=10,color='#263449',
                   position=(options.window_size[0]-280,options.window_size[1]-25-20*len(legend)),name='materials legend')
    p.view_isometric()
    p.reset_camera()
    return p


def _finish(plotter,show,path):
    if path is not None:
        path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    if show:
        plotter.show(screenshot=str(path) if path else None,auto_close=False)
    elif path:
        plotter.show(screenshot=str(path),auto_close=False,interactive=False)
    return plotter


def inspect_scene(scene,transmitters=None,receivers=None,*,frequency_hz=77e9,options=None,show=True,path=None):
    """Inspect geometry and Tx/Rx patterns WITHOUT calling the tracing solver.

    Gaussian receivers use their direction as outward look direction; this
    preview does not define a port or integrate receive-mode coupling. Returned
    Plotter remains open for further controls/screenshots; call close() when done.
    """
    opts=options or ViewOptions()
    return _finish(_base(scene,transmitters,receivers,frequency_hz,opts,show),show,path)


def visualize(result,*,transmitters=None,receivers=None,options=None,show=True,path=None):
    """VTK view of the same geometry plus computed beam branch axes."""
    import pyvista as pv
    opts=options or ViewOptions()
    p=_base(result.scene,transmitters,receivers,result.config.frequency_hz,opts,show)
    for i,segment in enumerate(result.segments):
        kind=segment.state.history[-1][1] if segment.state.history else 'launch'
        color={'launch':'#d64b35','reflect':'#ad56a1','transmit':'#298bc0'}[kind]
        p.add_mesh(pv.Line(segment.state.position,segment.state.position+segment.length*segment.state.direction),
                   color=color,line_width=2,name=f'branch segment {i}')
    p.add_text(f'{len(result.segments)} segments | {len(result.terminals)} terminals | {len(result.warnings)} warnings',
               position=(12,12),font_size=10,color='#263449',name='result caption')
    p.reset_camera()
    return _finish(p,show,path)
