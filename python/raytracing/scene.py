"""Explicit surface/volume topology and deterministic overlap resolution."""
from dataclasses import dataclass
import math
import numpy as np
from . import _core
from .materials import Material, VACUUM

Rectangle, Disk, Sphere = _core.Rectangle, _core.Disk, _core.Sphere
Cylinder, Triangle = _core.Cylinder, _core.Triangle


def vector(value):
    a = np.array(value, dtype=float, copy=True)
    if a.shape != (3,) or not np.isfinite(a).all():
        raise ValueError("Expected a finite 3-vector")
    return a


def unit(value):
    a = vector(value)
    if np.linalg.norm(a) == 0:
        raise ValueError("Direction must be nonzero")
    return a / np.linalg.norm(a)


@dataclass
class Box:
    center: object
    half_lengths: object
    axes: object = None

    def __post_init__(self):
        self.center, self.half_lengths = vector(self.center), vector(self.half_lengths)
        self.axes = np.eye(3) if self.axes is None else np.array(self.axes, dtype=float, copy=True)
        if (self.half_lengths <= 0).any() or self.axes.shape != (3, 3):
            raise ValueError("Box half lengths must be positive; axes must be 3x3")
        if not np.allclose(self.axes.T @ self.axes, np.eye(3), atol=1e-12, rtol=0) or np.linalg.det(self.axes) < 0:
            raise ValueError("Box axes must be a right-handed orthonormal frame")


def copy_shape(s):
    if isinstance(s, Box): return Box(s.center, s.half_lengths, s.axes)
    if isinstance(s, Rectangle): return Rectangle(s.center, s.u, s.v, s.half_width, s.half_height)
    if isinstance(s, Disk): return Disk(s.center, s.normal, s.radius)
    if isinstance(s, Sphere): return Sphere(s.center, s.radius)
    if isinstance(s, Cylinder): return Cylinder(s.center, s.axis, s.radius, s.half_length, s.capped)
    if isinstance(s, Triangle): return Triangle(s.a, s.b, s.c)
    raise TypeError("Unsupported geometry")


def faces(s):
    if not isinstance(s, Box): return [s]
    result = []
    for i in range(3):
        j, k = (i+1)%3, (i+2)%3
        for sign in (-1, 1):
            result.append(Rectangle(s.center + sign*s.half_lengths[i]*s.axes[:, i],
                                    s.axes[:, j], sign*s.axes[:, k], s.half_lengths[j], s.half_lengths[k]))
    return result


def contains(s, p):
    if isinstance(s, Box): return bool((np.abs(s.axes.T @ (p-s.center)) < s.half_lengths).all())
    if isinstance(s, Sphere): return np.linalg.norm(p-s.center) < s.radius
    if isinstance(s, Cylinder):
        v = p-s.center
        a = np.dot(v, s.axis)
        return abs(a) < s.half_length and np.linalg.norm(v-a*s.axis) < s.radius
    return False


def edge_distance(s, p):
    if isinstance(s, Rectangle):
        v = p-s.center
        return min(s.half_width-abs(v@s.u), s.half_height-abs(v@s.v))
    if isinstance(s, Disk): return s.radius-np.linalg.norm(p-s.center)
    if isinstance(s, Cylinder):
        v = p-s.center
        axial = v@s.axis
        radial = np.linalg.norm(v-axial*s.axis)
        return math.hypot(s.half_length-abs(axial), s.radius-radial)
    if isinstance(s, Triangle):
        return min(np.linalg.norm(np.cross(b-a, p-a))/np.linalg.norm(b-a)
                   for a,b in ((s.a,s.b), (s.b,s.c), (s.c,s.a)))
    return math.inf


def curvature(s, p, normal):
    """Derivative of the oriented unit normal with respect to position."""
    if isinstance(s, Sphere):
        outward = (p-s.center)/s.radius
        return np.sign(outward@normal)*(np.eye(3)-np.outer(outward,outward))/s.radius
    if isinstance(s, Cylinder):
        v = p-s.center
        radial = v-(v@s.axis)*s.axis
        if np.linalg.norm(radial) > 0 and abs(normal@s.axis) < .5:
            outward = radial/np.linalg.norm(radial)
            return np.sign(outward@normal)*(np.eye(3)-np.outer(s.axis,s.axis)-np.outer(outward,outward))/s.radius
    return np.zeros((3,3))


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


class Scene:
    def __init__(self, background=VACUUM, *, tolerance=1e-8):
        if not isinstance(background, Material) or background.kind != "dielectric":
            raise ValueError("Background must be a dielectric")
        if not math.isfinite(tolerance) or tolerance <= 0:
            raise ValueError("Tolerance must be positive")
        self.background, self.tolerance = background, tolerance
        self._objects, self._faces = [], []

    @property
    def objects(self):
        return tuple(Object(copy_shape(o.shape), o.material, o.volume, o.terminate) for o in self._objects)

    def _add(self, shape, material, volume, terminate=""):
        s = copy_shape(shape)
        if volume and not (isinstance(s, (Box, Sphere)) or isinstance(s, Cylinder) and s.capped):
            raise ValueError("Volumes must be boxes, spheres, or capped cylinders")
        fs = faces(s)
        for f in fs: _core._validate_surface(f)
        if not terminate:
            if not isinstance(material, Material): raise TypeError("Expected Material")
            if not volume and material.kind == "dielectric":
                raise ValueError("Dielectrics require a volume; surfaces must be PEC or metal")
        i = len(self._objects)
        self._objects.append(Object(s, material, volume, terminate))
        self._faces.extend((i, f) for f in fs)
        return i

    def add_volume(self, shape, material): return self._add(shape, material, True)
    def add_surface(self, shape, material): return self._add(shape, material, False)

    def add_termination(self, shape, *, volume=False, on="entry"):
        if on not in ("entry", "exit") or not volume and on != "entry":
            raise ValueError("Use entry/exit for volumes; surfaces terminate on any hit")
        return self._add(shape, None, volume, on)

    def region_at(self, position):
        p = vector(position)
        for i in range(len(self._objects)-1, -1, -1):
            o = self._objects[i]
            if o.volume and not o.terminate and contains(o.shape, p): return i
        return -1

    def material(self, region):
        return self.background if region == -1 else self._objects[region].material

    def starts_terminated(self, p):
        return any(o.volume and o.terminate == "entry" and contains(o.shape, p) for o in self._objects)

    def validate(self):
        """Validate shapes. Overlaps use insertion priority, junctions resolve at each hit.

        Geometry thinner than the intersection tolerance is not resolved.
        No destructive boolean operation is applied to user geometry.
        """
        for _, f in self._faces: _core._validate_surface(f)
        return self

    def next_event(self, origin, direction, max_distance):
        hits = []
        for i,f in self._faces:
            h = _core._intersect(f, origin, direction, self.tolerance, max_distance)
            if h is not None: hits.append((h[0], i, f, np.asarray(h[1])))
        # A single closed primitive may be entered without changing the winning
        # material. Requery after skipped crossings to find its exit as well.
        hits.sort(key=lambda h:h[0])
        if not hits: return None
        distance = hits[0][0]
        group = [h for h in hits if h[0]-distance <= self.tolerance]
        p = origin+distance*direction
        probe = 2*self.tolerance
        before, after = self.region_at(p-probe*direction), self.region_at(p+probe*direction)
        for _, i, f, n in group:
            o = self._objects[i]
            if o.terminate:
                was, now = contains(o.shape,p-probe*direction), contains(o.shape,p+probe*direction)
                if not o.volume or o.terminate == "entry" and not was and now or o.terminate == "exit" and was and not now:
                    return Event(distance,p,n,f,i,before,after,"terminate")
        active = [h for h in group if not self._objects[h[1]].terminate and
                  (not self._objects[h[1]].volume or h[1] in (before,after))]
        if active:
            _, i, f, n = max(active, key=lambda h:h[1])
            # An explicit opaque sheet supersedes region transitions at this hit.
            sheets = [h for h in active if not self._objects[h[1]].volume]
            if sheets: _, i, f, n = max(sheets, key=lambda h:h[1])
            incompatible = any(abs(float(n@h[3])) < 1-1e-10 for h in active)
            edge = any(edge_distance(h[2],p) <= self.tolerance for h in active)
            if incompatible or edge:
                return Event(distance,p,n,f,i,before,after,"ambiguous")
            if sheets or self.material(before) != self.material(after):
                return Event(distance,p,n,f,i,before,after,"interface")
        shifted = p+probe*direction
        remaining = max_distance-distance-probe
        if remaining <= self.tolerance: return None
        event = self.next_event(shifted,direction,remaining)
        if event is not None: event.distance += distance+probe
        return event
