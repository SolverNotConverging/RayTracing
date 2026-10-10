"""Finite paraxial Gaussian beams with the exp(-i omega t) convention."""
from dataclasses import dataclass
import math
import numpy as np
from .materials import C0, ETA0
from .scene import vector, unit


def frame(direction):
    d = unit(direction)
    axis = np.eye(3)[np.argmin(np.abs(d))]
    u = unit(np.cross(axis,d))
    return np.column_stack((u,np.cross(d,u)))


def polarization(value, direction):
    p = np.array(value,dtype=complex)
    if p.shape != (3,) or not np.isfinite(p).all() or np.linalg.norm(p) == 0:
        raise ValueError("Polarization must be a finite nonzero complex 3-vector")
    if abs(direction@p) > 1e-10*np.linalg.norm(p):
        raise ValueError("Polarization must be transverse")
    return p/np.linalg.norm(p)


@dataclass
class GaussianBeam:
    waist_position: object
    direction: object
    waist_radius: object
    power: float = 1.0
    polarization: object = (1,0,0)
    launch_distance: float = 0.0
    phase: float = 0.0
    transverse_axis: object = None

    def beamlets(self, samples=17, radius_ratio=.5, extent=4.0):
        """Coherent waist-plane Gaussian convolution, converging with samples.

        Constituent fields must be summed, not their powers. Clipping within
        each constituent remains a local-beam approximation.
        """
        if not isinstance(samples,int) or samples<3 or not 0<radius_ratio<1 or not math.isfinite(extent) or extent<=0:
            raise ValueError("Require samples>=3, 0<radius_ratio<1 and positive extent")
        d=unit(self.direction)
        basis=frame(d) if self.transverse_axis is None else np.column_stack((unit(self.transverse_axis),np.cross(d,unit(self.transverse_axis))))
        w=np.broadcast_to(np.asarray(self.waist_radius,dtype=float),(2,))
        small=w*radius_ratio
        spread=np.sqrt(w*w-small*small)
        xs=np.linspace(-extent*spread[0],extent*spread[0],samples)
        ys=np.linspace(-extent*spread[1],extent*spread[1],samples)
        area=(xs[1]-xs[0])*(ys[1]-ys[0])
        beams=[]
        for i,x in enumerate(xs):
            for j,y in enumerate(ys):
                weight=area*(.5 if i in (0,samples-1) else 1)*(.5 if j in (0,samples-1) else 1)
                amplitude_ratio=weight*w.prod()/(math.pi*small.prod()*spread.prod())*math.exp(-(x/spread[0])**2-(y/spread[1])**2)
                beams.append(GaussianBeam(vector(self.waist_position)+basis@np.array([x,y]),d,small,
                                          self.power*amplitude_ratio**2*small.prod()/w.prod(),self.polarization,
                                          self.launch_distance,self.phase,basis[:,0]))
        return beams

    def launch(self, scene, frequency):
        d = unit(self.direction)
        origin = vector(self.waist_position)+float(self.launch_distance)*d
        basis = frame(d)
        if self.transverse_axis is not None:
            u = unit(self.transverse_axis)
            if abs(u@d) > 1e-10: raise ValueError("Transverse axis must be perpendicular to beam axis")
            basis = np.column_stack((u,np.cross(d,u)))
        w = np.broadcast_to(np.asarray(self.waist_radius,dtype=float),(2,)).copy()
        if not np.isfinite(w).all() or (w<=0).any(): raise ValueError("Waist radii must be positive")
        if not all(math.isfinite(x) for x in (self.power,self.phase,self.launch_distance)) or self.power<=0:
            raise ValueError("Require positive finite power and finite launch distance/phase")
        medium = scene.material(scene.region_at(origin))
        if medium.kind != "dielectric": raise ValueError("Cannot launch inside opaque material")
        n = medium.index(frequency)
        if n.imag != 0:
            raise ValueError("Launch Gaussian sources in a lossless region; absorption begins at subsequent interfaces")
        k = 2*math.pi*frequency/C0*n.real
        if np.max(2/(k*w)) > .2:
            raise ValueError("Gaussian divergence exceeds 0.2 rad; outside this paraxial model")
        zr = k*w*w/2
        q = self.launch_distance-1j*zr
        h = np.diag(1/q)
        admittance = (n/medium.mu_r).real/ETA0
        amplitude = math.sqrt(4*self.power/(math.pi*w.prod()*admittance))
        gouy = np.prod(np.sqrt((-1j*zr)/q))
        field = amplitude*gouy*np.exp(1j*(k*self.launch_distance+self.phase))*polarization(self.polarization,d)
        return State(origin,d,basis,h,field,scene.region_at(origin),self.power)


@dataclass
class PlaneWave:
    """Single plane-wave channel. Power is flux through unit ray-normal area.

    Intended for interface/slab calculations; spatial clipping follows the axis.
    """
    position: object
    direction: object
    polarization: object = (1,0,0)
    amplitude: complex = 1.0

    def launch(self, scene, frequency):
        origin,d = vector(self.position),unit(self.direction)
        mat = scene.material(scene.region_at(origin))
        if mat.kind != "dielectric": raise ValueError("Cannot launch inside opaque material")
        a = complex(self.amplitude)
        if not np.isfinite(a): raise ValueError("Amplitude must be finite")
        field = a*polarization(self.polarization,d)
        power = abs(a)**2*(mat.index(frequency)/mat.mu_r).real/(2*ETA0)
        return State(origin,d,frame(d),np.zeros((2,2),complex),field,scene.region_at(origin),power)


@dataclass
class ApertureSource:
    """Fit a sampled transverse complex field with Gaussian beamlets.

    field(position) returns global Cartesian E on the launch plane. Extents
    are half-widths. Width, sampling and aperture extent require convergence.
    """
    position: object
    direction: object
    half_widths: object
    samples: int
    beamlet_radius: float
    field: object
    rcond: float = 1e-10

    def launch(self,scene,frequency):
        if hasattr(self,'_frequency_hz') and frequency!=self._frequency_hz:
            raise ValueError("Sampled antenna source frequency must match solver frequency")
        if not isinstance(self.samples,int) or not 3<=self.samples<=41:
            raise ValueError("Aperture samples per axis must lie in [3,41]")
        extent=np.broadcast_to(np.asarray(self.half_widths,dtype=float),(2,))
        if not np.isfinite(extent).all() or np.any(extent<=0) or not math.isfinite(self.beamlet_radius) or self.beamlet_radius<=0:
            raise ValueError("Aperture widths and beamlet radius must be positive")
        if not math.isfinite(self.rcond) or not 0<self.rcond<1:
            raise ValueError("rcond must be in (0,1)")
        d=unit(self.direction); basis=frame(d); origin=vector(self.position)
        xy=np.array([(x,y) for x in np.linspace(-extent[0],extent[0],self.samples)
                    for y in np.linspace(-extent[1],extent[1],self.samples)])
        positions=origin+xy@basis.T
        values=np.asarray([self.field(p.copy()) for p in positions],dtype=complex)
        if values.shape!=(len(xy),3) or not np.isfinite(values).all():
            raise ValueError("Aperture field must return finite complex Cartesian vectors")
        if np.max(np.abs(values@d))>1e-8*max(1e-300,np.linalg.norm(values)):
            raise ValueError("Aperture field must be transverse to the launch axis")
        kernel=np.exp(-np.sum((xy[:,None,:]-xy[None,:,:])**2,axis=2)/self.beamlet_radius**2)
        coefficients=np.linalg.lstsq(kernel,values,rcond=self.rcond)[0]
        medium=scene.material(scene.region_at(origin))
        if medium.kind!='dielectric': raise ValueError("Cannot launch inside opaque material")
        if getattr(self,'_vacuum_only',False) and (medium.epsilon_r!=1 or medium.mu_r!=1 or medium.conductivity!=0):
            raise ValueError("Antenna aperture adapter requires a vacuum launch medium")
        n=medium.index(frequency)
        if n.imag!=0: raise ValueError("Aperture sources require a homogeneous lossless launch plane")
        states=[]
        for p,e in zip(positions,coefficients):
            if scene.material(scene.region_at(p))!=medium:
                raise ValueError("Aperture launch plane crosses a material boundary")
            amplitude=np.linalg.norm(e)
            if amplitude<=1e-30: continue
            power=math.pi*self.beamlet_radius**2*(n/medium.mu_r).real/(4*ETA0)*amplitude**2
            states.append(GaussianBeam(p,d,self.beamlet_radius,power,e/amplitude,transverse_axis=basis[:,0]).launch(scene,frequency))
        if not states: raise ValueError("Aperture field is zero")
        return states

    @classmethod
    def from_antenna(cls,antenna,antenna_position,position,direction,half_widths,samples,beamlet_radius,frequency_hz):
        """Sample a native far-field antenna on a vacuum plane in its far zone."""
        from . import _core
        origin=vector(antenna_position); d=unit(direction)
        def field(p):
            offset=p-origin; distance=np.linalg.norm(offset)
            if distance==0: raise ValueError("Antenna cannot lie on its sampled aperture")
            look=offset/distance
            if look@d<math.cos(.2): raise ValueError("Antenna aperture exceeds the 0.2 rad paraxial cone")
            value=np.asarray(antenna.farfield(look,frequency_hz,_core.Medium()))/distance
            return (value-d*(d@value))*np.exp(2j*math.pi*frequency_hz/C0*distance)
        result=cls(position,d,half_widths,samples,beamlet_radius,field)
        result._frequency_hz=frequency_hz
        result._vacuum_only=True
        return result


@dataclass
class State:
    position: np.ndarray
    direction: np.ndarray
    basis: np.ndarray
    curvature: np.ndarray
    field: np.ndarray
    region: int
    power: float
    distance: float = 0.0
    optical_path: float = 0.0
    history: tuple = ()
    branch: int = 0
    parent: int = -1


def propagate(state, distance, material, frequency):
    n = material.index(frequency)
    k0 = 2*math.pi*frequency/C0
    h = state.curvature
    transform = np.eye(2)+distance*h
    # Each eigenvalue advances in the upper half-plane; taking each square
    # root separately preserves both Gouy phase contributions across a focus.
    eigenvalues=1+distance*np.linalg.eigvals(h)
    if np.min(np.abs(eigenvalues))<1e-14:
        raise ValueError("Singular geometrical plane-wave focus; use a finite GaussianBeam")
    denominator = np.prod(np.sqrt(eigenvalues))
    field = state.field/denominator*np.exp(1j*k0*n*distance)
    return State(state.position+distance*state.direction,state.direction.copy(),state.basis.copy(),
                 h@np.linalg.inv(transform),field,state.region,
                 state.power*math.exp(-2*k0*n.imag*distance),state.distance+distance,
                 state.optical_path+n.real*distance,state.history,state.branch,state.parent)


@dataclass
class Segment:
    state: State
    length: float
    material: object
    frequency: float

    def field_at(self, point):
        offset = vector(point)-self.state.position
        s = float(offset@self.state.direction)
        # Half-open intervals avoid counting both sides of one interaction.
        if s < 0 or s >= self.length: return np.zeros(3,complex)
        at = propagate(self.state,s,self.material,self.frequency)
        x = self.state.basis.T@offset
        k = 2*math.pi*self.frequency/C0*self.material.index(self.frequency).real
        return at.field*np.exp(.5j*k*(x@at.curvature@x))

    def radii(self, distance=0.0):
        at = propagate(self.state,distance,self.material,self.frequency)
        eigen = np.linalg.eigvalsh(at.curvature.imag)
        k = 2*math.pi*self.frequency/C0*self.material.index(self.frequency).real
        return np.sqrt(2/(k*eigen)) if np.all(eigen>0) else np.full(2,np.inf)
