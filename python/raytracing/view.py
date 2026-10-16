"""PyVista 3D views of scenes, point-ray paths and beam branches.

Each function returns the ``pyvista.Plotter``. With ``show=False`` it renders
off screen, which suits notebooks and ``screenshot`` export.
"""
import numpy as np

from . import _core
from ._core import beam as _beam

_SURFACE = "#a6b3cc"
_COLORS = {"dielectric": "#5a99c9", "opaque": "#d29b49", "termination": "#737b85",
           "launch": "#d64b35", "reflect": "#ad56a1", "transmit": "#298bc0"}


def _mesh(pv, shape, resolution=48):
    if isinstance(shape, _core.Box):
        meshes = [_mesh(pv, face, resolution) for face in _core.faces(shape)]
        return meshes[0].merge(meshes[1:])
    if isinstance(shape, _core.Sphere):
        return pv.Sphere(radius=shape.radius, center=shape.center, theta_resolution=resolution,
                         phi_resolution=resolution)
    if isinstance(shape, _core.Cylinder):
        return pv.Cylinder(center=shape.center, direction=shape.axis, radius=shape.radius,
                           height=2 * shape.half_length, resolution=resolution, capping=shape.capped)
    if isinstance(shape, _core.Disk):
        return pv.Disc(center=shape.center, normal=shape.normal, inner=0, outer=shape.radius, c_res=resolution)
    if isinstance(shape, _core.Rectangle):
        u, v, c = shape.half_width * shape.u, shape.half_height * shape.v, shape.center
        return pv.PolyData(np.array([c - u - v, c + u - v, c + u + v, c - u + v]), faces=[4, 0, 1, 2, 3])
    if isinstance(shape, _core.Triangle):
        return pv.PolyData(np.array([shape.a, shape.b, shape.c]), faces=[3, 0, 1, 2])
    raise TypeError(f"Unsupported shape {type(shape).__name__}")


def _polylines(pv, paths):
    paths = [np.asarray(p) for p in paths if len(p) >= 2]
    if not paths:
        return None
    lines, offset = [], 0
    for p in paths:
        lines.extend([len(p), *range(offset, offset + len(p))])
        offset += len(p)
    return pv.PolyData(np.concatenate(paths), lines=np.asarray(lines))


def _plotter(show, window_size, background="white"):
    import pyvista as pv

    plotter = pv.Plotter(off_screen=not show, window_size=window_size)
    plotter.set_background(background)
    plotter.add_axes()
    return pv, plotter


def _finish(plotter, show, screenshot):
    plotter.reset_camera()
    if show:
        plotter.show(screenshot=screenshot)
    elif screenshot is not None:
        plotter.screenshot(str(screenshot))
    return plotter


def show_paths(result, *, show=True, screenshot=None, window_size=(1000, 750), dynamic_range_db=60.0):
    """Point-ray result: surfaces, solved paths coloured by received level, Tx (red) and Rx (green)."""
    pv, plotter = _plotter(show, window_size)
    for surface in result.scene.surfaces:
        plotter.add_mesh(_mesh(pv, surface), color=_SURFACE, opacity=0.3)
    lines = _polylines(pv, [path.vertices for path in result.paths])
    if lines is not None:
        levels = np.array([abs(path.coefficient) if path.field is not None else np.nan for path in result.paths])
        peak = np.nanmax(levels, initial=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            db = np.maximum(20 * np.log10(levels / peak), -dynamic_range_db) if peak > 0 else levels * np.nan
        lines.cell_data["level (dB)"] = db
        plotter.add_mesh(lines, scalars="level (dB)", clim=(-dynamic_range_db, 0), cmap="viridis",
                         nan_color="grey", line_width=3)
    extent = max(np.ptp(np.array(plotter.bounds).reshape(3, 2), axis=1).max(), 1e-3)
    for endpoint, color in ((result.transmitter, "red"), (result.receiver, "green")):
        plotter.add_mesh(pv.Sphere(radius=0.02 * extent, center=endpoint.position), color=color)
    return _finish(plotter, show, screenshot)


def _add_beam_scene(pv, plotter, scene):
    for obj in scene.objects:
        mesh = _mesh(pv, obj.shape)
        if obj.material is None:
            plotter.add_mesh(mesh, style="wireframe", color=_COLORS["termination"])
        else:
            kind = "dielectric" if obj.material.kind == _core.MaterialKind.DIELECTRIC else "opaque"
            plotter.add_mesh(mesh, color=_COLORS[kind], opacity=0.25, show_edges=True)


def show_beams(simulation, *, show=True, screenshot=None, window_size=(1000, 750)):
    """Beam simulation: material objects and every branch segment, coloured by its last interaction."""
    pv, plotter = _plotter(show, window_size)
    _add_beam_scene(pv, plotter, simulation.scene)
    for segment in simulation.segments:
        history = segment.state.history
        kind = "launch" if not history else history[-1].kind.name.lower()
        start = segment.state.position
        plotter.add_mesh(pv.Line(start, start + segment.length * segment.state.direction),
                         color=_COLORS[kind], line_width=2)
    return _finish(plotter, show, screenshot)


def preview(scene, sources, *, frequency_hz=77e9, length=1.0, show=True, screenshot=None, window_size=(1000, 750)):
    """Beam scene and the unscattered 1/e envelopes of its sources, without tracing."""
    pv, plotter = _plotter(show, window_size)
    _add_beam_scene(pv, plotter, scene)
    ring = np.linspace(0, 2 * np.pi, 33)
    for source in sources if isinstance(sources, (list, tuple)) else [sources]:
        launched = source.launch(scene, frequency_hz)
        for state in launched if isinstance(launched, list) else [launched]:
            material = scene.material(state.region)
            k = 2 * np.pi * frequency_hz / _core.SPEED_OF_LIGHT * material.index(frequency_hz).real
            stations = []
            for z in np.linspace(0, length, 24):
                at = _beam.propagate(state, z, material, frequency_hz)
                values, vectors = np.linalg.eigh(at.curvature.imag)
                radii = np.sqrt(2 / (k * values))
                contour = vectors @ (radii[:, None] * np.array([np.cos(ring), np.sin(ring)]))
                stations.append(at.position[:, None] + at.basis @ contour)
            grid = np.asarray(stations).transpose(1, 2, 0)
            plotter.add_mesh(pv.StructuredGrid(*grid), color=_COLORS["launch"], opacity=0.25)
            plotter.add_mesh(pv.Line(state.position, state.position + length * state.direction),
                             color=_COLORS["launch"], line_width=3)
    return _finish(plotter, show, screenshot)
