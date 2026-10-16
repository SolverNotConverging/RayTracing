"""Python source/state data adapters for the C++ paraxial beam solver."""
from dataclasses import dataclass
import numpy as np
from ._beam_adapter import native, record, values, pair
from .scene import vector
from .materials import Material


def frame(direction): return native.frame(vector(direction))
def polarization(value, direction): return native.polarization(np.asarray(value, dtype=complex), vector(direction))


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

    def _native(self):
        result = native.GaussianBeam()
        for name in ("waist_position", "direction", "power", "polarization", "launch_distance", "phase", "transverse_axis"):
            setattr(result, name, getattr(self, name))
        result.waist_radius = pair(self.waist_radius)
        return result

    @classmethod
    def _from_native(cls, value): return cls(**values(value, cls))

    def beamlets(self, samples=17, radius_ratio=.5, extent=4.0):
        """Coherent Gaussian convolution; fields, rather than powers, superpose."""
        if not isinstance(samples, int): raise ValueError("samples must be an integer")
        return [type(self)._from_native(b) for b in self._native().beamlets(samples, radius_ratio, extent)]

    def launch(self, scene, frequency):
        return State._from_native(self._native().launch(scene._native(), frequency))


@dataclass
class PlaneWave:
    """Plane-wave channel with power defined through unit ray-normal area."""
    position: object
    direction: object
    polarization: object = (1,0,0)
    amplitude: complex = 1.0

    def _native(self): return record(self, native.PlaneWave)

    def launch(self, scene, frequency):
        return State._from_native(self._native().launch(scene._native(), frequency))


@dataclass
class ApertureSource:
    """C++ Gaussian fit of a callable Cartesian field on a sampled plane.

    The callable supplies input samples under Python's GIL. Sampling geometry,
    least-squares fitting and beam launch all run in the native library.
    """
    position: object
    direction: object
    half_widths: object
    samples: int
    beamlet_radius: float
    field: object
    rcond: float = 1e-10

    def _native(self):
        if not isinstance(self.samples, int): raise ValueError("samples must be an integer")
        result = native.ApertureSource()
        for name in ("position", "direction", "samples", "beamlet_radius", "field", "rcond"):
            setattr(result, name, getattr(self, name))
        result.half_widths = pair(self.half_widths)
        result.frequency_hz = getattr(self, "_frequency_hz", None)
        result.vacuum_only = getattr(self, "_vacuum_only", False)
        return result

    def launch(self, scene, frequency):
        return [State._from_native(s) for s in self._native().launch(scene._native(), frequency)]

    @classmethod
    def from_antenna(cls, antenna, antenna_position, position, direction, half_widths, samples, beamlet_radius, frequency_hz):
        """Native far-field antenna sampling on a vacuum plane in its far zone."""
        if not isinstance(samples, int): raise ValueError("samples must be an integer")
        value = native.ApertureSource.from_antenna(antenna, vector(antenna_position), vector(position),
                    vector(direction), pair(half_widths), samples, beamlet_radius, frequency_hz)
        result = cls(**values(value, cls))
        result._frequency_hz, result._vacuum_only = value.frequency_hz, value.vacuum_only
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

    def _native(self): return record(self, native.State)

    @classmethod
    def _from_native(cls, value):
        data = values(value, cls)
        data["history"] = tuple(tuple(h) for h in data["history"])
        return cls(**data)


def propagate(state, distance, material, frequency):
    return State._from_native(native.propagate(state._native(), distance, material._native(), frequency))


@dataclass
class Segment:
    state: State
    length: float
    material: object
    frequency: float

    def _native(self): return native.Segment(self.state._native(), self.length, self.material._native(), self.frequency)

    @classmethod
    def _from_native(cls, value):
        return cls(State._from_native(value.state), value.length, Material._from_native(value.material), value.frequency)

    def field_at(self, point): return self._native().field_at(vector(point))
    def radii(self, distance=0.0): return self._native().radii(distance)
