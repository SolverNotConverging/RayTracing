"""Scene data adapters; region and intersection decisions are shared C++ code."""
from dataclasses import dataclass
import numpy as np
from . import _core
from ._beam_adapter import native, record, values
from .materials import Material, VACUUM

Rectangle, Disk, Sphere = _core.Rectangle, _core.Disk, _core.Sphere
Cylinder, Triangle = _core.Cylinder, _core.Triangle


def vector(value):
    a = np.array(value, dtype=float, copy=True)
    if a.shape != (3,) or not np.isfinite(a).all():
        raise ValueError("Expected a finite 3-vector")
    return a


def unit(value):
    return native.unit(vector(value))


@dataclass
class Box:
    center: object
    half_lengths: object
    axes: object = None

    def __post_init__(self):
        self.center, self.half_lengths = vector(self.center), vector(self.half_lengths)
        self.axes = np.eye(3) if self.axes is None else np.array(self.axes, dtype=float, copy=True)
        if self.axes.shape != (3, 3):
            raise ValueError("Box axes must be 3x3")
        self._native().validate()

    def _native(self):
        return native.Box(self.center, self.half_lengths, self.axes)


def _shape_native(shape):
    return shape._native() if isinstance(shape, Box) else shape


def _shape_python(shape):
    return Box(shape.center, shape.half_lengths, shape.axes) if isinstance(shape, native.Box) else shape


def copy_shape(s):
    if isinstance(s, Box): return Box(s.center, s.half_lengths, s.axes)
    if isinstance(s, Rectangle): return Rectangle(s.center, s.u, s.v, s.half_width, s.half_height)
    if isinstance(s, Disk): return Disk(s.center, s.normal, s.radius)
    if isinstance(s, Sphere): return Sphere(s.center, s.radius)
    if isinstance(s, Cylinder): return Cylinder(s.center, s.axis, s.radius, s.half_length, s.capped)
    if isinstance(s, Triangle): return Triangle(s.a, s.b, s.c)
    raise TypeError("Unsupported geometry")


def faces(s):
    return native.faces(_shape_native(s))


def contains(s, p):
    return native.contains(_shape_native(s), vector(p))


def edge_distance(s, p):
    return native.edge_distance(s, vector(p))


def curvature(s, p, normal):
    return native.curvature(s, vector(p), vector(normal))


@dataclass(frozen=True)
class Object:
    shape: object
    material: object
    volume: bool
    terminate: str = ""


@dataclass
class Event:
    distance: float
    position: np.ndarray
    normal: np.ndarray
    shape: object
    object_id: int
    before: int
    after: int
    kind: str

    def _native(self): return record(self, native.Event)

    @classmethod
    def _from_native(cls, value): return cls(**values(value, cls))


class Scene:
    def __init__(self, background=VACUUM, *, tolerance=1e-8):
        if not isinstance(background, Material): raise ValueError("Background must be a dielectric")
        native.Scene(background._native(), tolerance)
        self.background, self.tolerance = background, tolerance
        self._objects = []

    def _native(self):
        result = native.Scene(self.background._native(), self.tolerance)
        for o in self._objects:
            result.add(_shape_native(o.shape), o.material._native() if o.material else None, o.volume, o.terminate)
        return result

    @classmethod
    def _from_native(cls, value):
        result = cls(Material._from_native(value.background), tolerance=value.tolerance)
        result._objects = [Object(_shape_python(o.shape), Material._from_native(o.material) if o.material else None,
                                  o.volume, o.terminate) for o in value.objects]
        return result

    @property
    def objects(self):
        return tuple(Object(copy_shape(o.shape), o.material, o.volume, o.terminate) for o in self._objects)

    def _add(self, shape, material, volume, terminate=""):
        s = copy_shape(shape)
        if material is not None and not isinstance(material, Material): raise TypeError("Expected Material")
        index = self._native().add(_shape_native(s), material._native() if material else None, volume, terminate)
        self._objects.append(Object(s, material, volume, terminate))
        return index

    def add_volume(self, shape, material): return self._add(shape, material, True)
    def add_surface(self, shape, material): return self._add(shape, material, False)

    def add_termination(self, shape, *, volume=False, on="entry"):
        if on not in ("entry", "exit") or not volume and on != "entry":
            raise ValueError("Use entry/exit for volumes; surfaces terminate on any hit")
        return self._add(shape, None, volume, on)

    def region_at(self, position): return self._native().region_at(vector(position))
    def material(self, region): return Material._from_native(self._native().material(region))
    def starts_terminated(self, p): return self._native().starts_terminated(vector(p))

    def validate(self):
        self._native().validate()
        return self

    def next_event(self, origin, direction, max_distance):
        event = self._native().next_event(vector(origin), vector(direction), max_distance)
        return Event._from_native(event) if event is not None else None
