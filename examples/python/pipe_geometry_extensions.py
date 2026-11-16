"""Compute isolated specular returns in five noncanonical pipe geometries.

These are geometry examples, using PEC surfaces to identify return paths. They
do not estimate liquid-measurement accuracy. No launch direction is prescribed.
Run: python pipe_geometry_extensions.py --output geometry_extensions
Requires a build containing EllipticCylinder and Annulus, plus the viz extra.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform

import numpy as np
import raytracing as rt

FREQUENCY = 77e9
X, Y, Z = np.eye(3)


def cases():
    def make(name, title, walls, liquid, source, camera):
        scene = rt.Scene()
        for wall in walls:
            scene.add(wall, rt.Material.pec())
        liquid_index = scene.add(liquid, rt.Material.pec())
        return dict(name=name, title=title, scene=scene, liquid_index=liquid_index,
                    source=np.array(source, float), camera=camera)

    angle = np.deg2rad(15)
    axis = np.array([np.sin(angle), 0, np.cos(angle)])
    radial = np.array([np.cos(angle), 0, -np.sin(angle)])
    yield make('tilted', 'Tilted circular pipe',
               [rt.Cylinder(.75*axis, axis, .5, 1.0, False)],
               rt.Rectangle([0, 0, 0], X, Y, .65, .55),
               1.6*axis-.15*radial+.07*Y, [(3, -5, 2.7), (.2, 0, .8), Z])
    walls = [rt.Rectangle([s*.5, 0, .8], Y, Z, .5, .85) for s in [-1, 1]]
    walls += [rt.Rectangle([0, s*.5, .8], X, Z, .5, .85) for s in [-1, 1]]
    yield make('square', 'Square pipe', walls, rt.Rectangle([0, 0, 0], X, Y, .5, .5),
               [.13, .07, 1.6], [(3, -5, 2.8), (0, 0, .8), Z])
    yield make('ellipse', 'Elliptical pipe',
               [rt.EllipticCylinder([0, 0, .8], Z, X, .55, .4, .85, False)],
               rt.Rectangle([0, 0, 0], X, Y, .55, .4),
               [.13, .07, 1.6], [(3, -5, 2.8), (0, 0, .8), Z])
    yield make('bore_step', 'Upper bore transition',
               [rt.Cylinder([0, 0, .625], Z, .5, .675, False),
                rt.Cylinder([0, 0, 1.55], Z, .35, .25, False),
                rt.Annulus([0, 0, 1.3], Z, .35, .5)],
               rt.Disk([0, 0, 0], Z, .5), [.08, .06, 1.7],
               [(3, -5, 2.8), (0, 0, .85), Z])
    yield make('horizontal', 'Horizontal pipe, wall port',
               [rt.Cylinder([0, 0, 0], X, .5, 1.1, False)],
               rt.Rectangle([0, 0, -.15], X, Y, 1.1, np.sqrt(.5**2-.15**2)),
               [0, .30, .32], [(3.5, -5, 2.8), (0, 0, 0), Z])


def inspect_path(path, case):
    vertices = np.asarray(path.vertices)
    events = path.refinement.geometry.reflections
    sequence = [int(h.surface_index) for h in events]
    incoming = np.diff(vertices, axis=0)
    unit = incoming / np.linalg.norm(incoming, axis=1)[:, None]
    reflection_errors = []
    for i, event in enumerate(events):
        n = np.asarray(event.normal)
        reflection_errors.append(np.linalg.norm(unit[i]-2*np.dot(unit[i], n)*n-unit[i+1]))
    jacobian = np.asarray(path.spreading.jacobian)
    singular = np.linalg.svd(jacobian, compute_uv=False)
    # Refinement's miss is computed independently of the exported final vertex.
    last = np.asarray(events[-1].position)
    final = np.asarray(path.refinement.geometry.final_direction)
    displacement = case['source']-last
    miss = np.linalg.norm(displacement-np.dot(displacement, final)*final)
    assert path.spreading.status == rt.SpreadingStatus.VALID
    assert np.isfinite(jacobian).all() and singular[-1] > 1e-7
    assert miss < 1e-8 and max(reflection_errors, default=0) < 1e-8
    assert np.isfinite(complex(path.coefficient))
    if case['name'] == 'horizontal':
        assert np.max(np.abs(vertices[:, 0])) < 1e-8
    return dict(vertices_m=vertices.tolist(), surface_sequence=sequence,
                path_length_m=float(path.field.path_length),
                coefficient=[path.coefficient.real, path.coefficient.imag],
                launch_direction=np.asarray(path.refinement.launch_direction).tolist(),
                receiver_miss_m=float(miss), max_specular_error=float(max(reflection_errors, default=0)),
                endpoint_jacobian=jacobian.tolist(), singular_values_m=singular.tolist(),
                determinant_m2=float(np.linalg.det(jacobian)))


def solve_case(case, ray_count):
    antenna = rt.Antenna(rt.Isotropic(rt.Polarization.VERTICAL, 1+0j))
    antenna.orientation = np.array([[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]])
    result = rt.solve(case['scene'], rt.Transmitter(case['source'], antenna),
                      rt.Receiver(case['source'], antenna),
                      rt.SolverConfig(frequency_hz=FREQUENCY, ray_count=ray_count,
                                      max_reflections=8, max_distance=7,
                                      reception_radius=.15))
    paths = [p for p in result.paths if p.field is not None and p.spreading.status == rt.SpreadingStatus.VALID
             and any(h.surface_index == case['liquid_index'] for h in p.refinement.geometry.reflections)]
    paths.sort(key=lambda p: p.field.path_length)
    assert len(paths) >= 2, f"Too few liquid-dependent returns: {case['name']}"
    records = [inspect_path(p, case) for p in paths]
    # Show a small subset with distinct sequences, identifying reversed routes.
    # In the horizontal view, omit very short wall detours from the illustration
    # so that the cross-sectional return geometry is visible. Retain all data.
    shown, sequences = [], set()
    for i, record in enumerate(records):
        sequence = tuple(record['surface_sequence'])
        sequence = min(sequence, sequence[::-1])
        if case['name']=='horizontal' and i>0 and record['path_length_m'] < 1.5*records[0]['path_length_m']:
            continue
        if sequence not in sequences:
            shown.append(i)
            sequences.add(sequence)
        if len(shown) == 3:
            break
    return dict(name=case['name'], source_m=case['source'].tolist(),
                liquid_surface_index=int(case['liquid_index']),
                valid_liquid_return_count=len(records), shown_indices=shown, paths=records)


def render(case, record, output):
    import pyvista as pv
    from raytracing.view import _mesh
    plotter = pv.Plotter(off_screen=True, window_size=(700, 600))
    plotter.set_background('white')
    for i, surface in enumerate(case['scene'].surfaces):
        mesh = _mesh(pv, surface, 120)
        is_liquid = i == case['liquid_index']
        # Crop the liquid rendering to the actual cross section.
        if is_liquid and case['name'] in ('tilted', 'ellipse'):
            t = np.linspace(0, 2*np.pi, 180, endpoint=False)
            a, b = ((.5/np.cos(np.deg2rad(15)), .5) if case['name']=='tilted' else (.55, .4))
            points = np.column_stack([a*np.cos(t), b*np.sin(t), np.zeros_like(t)])
            mesh = pv.PolyData(points, faces=np.r_[len(t), np.arange(len(t))])
        plotter.add_mesh(mesh, color='#67b6da' if is_liquid else '#a6b3c1',
                         opacity=.60 if is_liquid else .18, smooth_shading=True)
        if is_liquid or isinstance(surface, rt.Annulus):
            plotter.add_mesh(mesh.extract_feature_edges(), color='#51738c', line_width=1)
    colors = ['#d55337', '#8c50a5', '#2e8054']
    for index, color in zip(record['shown_indices'], colors):
        vertices = np.array(record['paths'][index]['vertices_m'])
        plotter.add_mesh(pv.lines_from_points(vertices), color=color, line_width=2)
    plotter.add_mesh(pv.Sphere(radius=.035, center=case['source']), color='#202020')
    plotter.camera_position = case['camera']
    plotter.enable_parallel_projection()
    plotter.reset_camera()
    plotter.camera.zoom(1.12)
    plotter.screenshot(str(output/f"{case['name']}.png"))
    plotter.close()


def compose(output):
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch
    plt.rcParams.update({'font.size':7, 'axes.titlesize':7})
    fig, axes = plt.subplots(3, 2, figsize=(3.5, 4.45))
    titles=['Tilted circular pipe','Square pipe','Elliptical pipe',
            'Upper bore transition','Horizontal, wall port']
    for i, (ax, case) in enumerate(zip(axes.flat, cases())):
        ax.imshow(plt.imread(output/f"{case['name']}.png"))
        ax.set_title(f"({chr(97+i)}) {titles[i]}", pad=0)
        ax.axis('off')
    ax=axes.flat[-1]
    ax.axis('off'); ax.set_xlim(0,1); ax.set_ylim(0,1)
    ax.set_title('(f) Common sensing chain', pad=0)
    for y, text in [(.79,'Returned-ray model'), (.55,'Impulse response (CIR)'),
                    (.31,'Radar range processing'), (.07,'Distance / level estimate')]:
        ax.add_patch(FancyBboxPatch((.025,y),.95,.14,boxstyle='round,pad=0.015',
                                  facecolor='#e9f1f6',edgecolor='#628296',linewidth=.7))
        ax.text(.5,y+.07,text,ha='center',va='center',fontsize=6.5)
    for top in [.79,.55,.31]:
        ax.annotate('',xy=(.5,top-.085),xytext=(.5,top-.015),arrowprops=dict(arrowstyle='->',lw=.7))
    fig.subplots_adjust(left=.005,right=.995,bottom=.005,top=.96,wspace=.035,hspace=.13)
    fig.savefig(output/'geometry_extensions.pdf')
    fig.savefig(output/'geometry_extensions.png',dpi=300)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('geometry_extensions'))
    parser.add_argument('--ray-count',type=int,default=16000)
    parser.add_argument('--plot-only',action='store_true')
    args=parser.parse_args(); args.output.mkdir(parents=True,exist_ok=True)
    manifest=dict(frequency_hz=FREQUENCY,ray_count=args.ray_count,max_reflections=8,
                  max_distance_m=7,reception_radius_m=.15, material='PEC geometry demonstration',
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  core_sha256=hashlib.sha256(Path(rt._core.__file__).read_bytes()).hexdigest(),
                  raytracing_version=rt.__version__,python=platform.python_version(),cases=[])
    for case in cases():
        file=args.output/f"{case['name']}_paths.json"
        if args.plot_only:
            record=json.loads(file.read_text())
        else:
            record=solve_case(case,args.ray_count)
            file.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8',newline='\n')
        render(case,record,args.output)
        manifest['cases'].append({k:v for k,v in record.items() if k!='paths'})
        print(case['name'], record['valid_liquid_return_count'], 'valid returns',flush=True)
    compose(args.output)
    if not args.plot_only:
        (args.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':
    main()
