"""Python data API for the shared C++ Gaussian-beam numerical solver."""
from dataclasses import dataclass, field
import numpy as np
from ._beam_adapter import native, record, values
from .scene import Scene, vector
from .beams import State, Segment


def _array_key(value):
    array = np.asarray(value)
    return array.shape, array.dtype.str, array.tobytes()


_shape_attributes = {
    "Box": ("center", "half_lengths", "axes"),
    "Rectangle": ("center", "u", "v", "half_width", "half_height"),
    "Disk": ("center", "normal", "radius"),
    "Sphere": ("center", "radius"),
    "Cylinder": ("center", "axis", "radius", "half_length", "capped"),
    "Triangle": ("a", "b", "c"),
}


@dataclass
class SolverConfig:
    frequency_hz: float = 77e9
    max_interactions: int = 16
    max_distance: float = 20.0
    min_power_fraction: float = 1e-10
    max_branches: int = 100000
    max_loss_ratio: float = 0.01

    def _native(self):
        if not isinstance(self.max_interactions, int) or self.max_interactions < 0 or not isinstance(self.max_branches, int) or self.max_branches < 1:
            raise ValueError("Invalid interaction/branch budget")
        return record(self, native.SolverConfig)

    def validate(self): self._native().validate()

    @classmethod
    def _from_native(cls, value): return cls(**values(value, cls))


@dataclass
class Terminal:
    branch: int
    reason: str
    position: np.ndarray
    power: float
    history: tuple

    def _native(self): return record(self, native.Terminal)

    @classmethod
    def _from_native(cls, value):
        data = values(value, cls)
        data["history"] = tuple(tuple(h) for h in data["history"])
        return cls(**data)


@dataclass
class Simulation:
    config: SolverConfig
    scene: object
    segments: list = field(default_factory=list)
    terminals: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    launched_power: float = 0.0
    absorbed_power: float = 0.0
    interface_flux_residual: float = 0.0

    def _native(self):
        result = native.Simulation()
        result.config, result.scene = self.config._native(), self.scene._native()
        result.segments = [s._native() for s in self.segments]
        result.terminals = [t._native() for t in self.terminals]
        result.warnings = self.warnings
        result.launched_power, result.absorbed_power = self.launched_power, self.absorbed_power
        result.interface_flux_residual = self.interface_flux_residual
        return result

    @classmethod
    def _from_native(cls, value):
        result = cls(SolverConfig._from_native(value.config), Scene._from_native(value.scene),
                   [Segment._from_native(s) for s in value.segments],
                   [Terminal._from_native(t) for t in value.terminals], list(value.warnings),
                   value.launched_power, value.absorbed_power, value.interface_flux_residual)
        result._field_native = value
        result._field_key = result._sampling_key()
        return result

    def _sampling_key(self):
        # Compare numerical input bytes so in-place NumPy edits also invalidate
        # the cache. Diagnostics and branch histories do not affect sampling.
        objects = tuple((type(o.shape).__name__,
                         tuple(_array_key(getattr(o.shape, name))
                               for name in _shape_attributes[type(o.shape).__name__]),
                         o.material, o.volume, o.terminate) for o in self.scene._objects)
        segments = []
        for segment in self.segments:
            s = segment.state
            segments.append((_array_key(s.position), _array_key(s.direction),
                             _array_key(s.basis), _array_key(s.curvature), _array_key(s.field),
                             s.power, s.distance, s.optical_path,
                             segment.length, segment.material, segment.frequency))
        return self.scene.background, self.scene.tolerance, objects, tuple(segments)

    def _sampling_native(self):
        key = self._sampling_key()
        if getattr(self, "_field_key", None) != key:
            value = native.Simulation()
            value.scene = self.scene._native()
            value.segments = [s._native() for s in self.segments]
            self._field_native, self._field_key = value, key
        return self._field_native

    def field_at(self, position): return self._sampling_native().field_at(vector(position))

    def fields_at(self, positions):
        points = np.asarray(positions, dtype=float)
        if points.size == 0: return np.empty((0, 3), dtype=complex)
        if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
            raise ValueError("Expected finite positions with shape (N,3)")
        return np.asarray(self._sampling_native().fields_at(points), dtype=complex).reshape((-1, 3))

    @property
    def power_balance(self): return self._native().power_balance()


def interact(state, event, scene, frequency):
    children, reason = native.interact(state._native(), event._native(), scene._native(), frequency)
    return [State._from_native(s) for s in children], reason


def solve(scene, source, config=None):
    """Launch sources in C++, then run native branching transport without the GIL.

    Python aperture callables supply input samples during launch. Geometry,
    propagation, interfaces, branch budgets and field evaluation share the C++ API.
    """
    cfg = config if config is not None else SolverConfig()
    settings = cfg._native()
    settings.validate()
    snapshot = Scene._from_native(scene._native())
    sources = source if isinstance(source, (tuple, list)) else [source]
    states = []
    for src in sources:
        launched = src.launch(snapshot, settings.frequency_hz)
        states.extend(s._native() for s in (launched if isinstance(launched, list) else [launched]))
    return Simulation._from_native(native.solve_states(snapshot._native(), states, settings))
