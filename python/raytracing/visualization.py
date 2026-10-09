"""Python-only visualization. Rendering libraries are imported on demand."""
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np

from . import _core as core


def plot(response, path=None):
    """Return a Matplotlib Figure with impulse magnitude and phase; optionally save it."""
    import matplotlib.pyplot as plt

    delays = np.array([tap.delay_seconds for tap in response.taps]) * 1e9
    coefficients = np.array([tap.coefficient for tap in response.taps], dtype=complex)
    magnitudes = np.abs(coefficients)
    peak = magnitudes.max(initial=0)
    defined = magnitudes > peak * 1e-12
    figure, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True, layout="constrained")
    if len(delays):
        axes[0].stem(delays, magnitudes, linefmt="C0-", markerfmt="C0o", basefmt=" ")
    else:
        axes[0].text(0.5, 0.5, "No valid field paths", ha="center", transform=axes[0].transAxes)
    if defined.any():
        axes[1].stem(delays[defined], np.angle(coefficients[defined], deg=True),
                     linefmt="C3-", markerfmt="C3o", basefmt=" ")
    else:
        axes[1].text(0.5, 0.5, "Phase undefined: no nonzero taps", ha="center", transform=axes[1].transAxes)
    axes[0].set(title=f"Impulse response ({response.frequency_hz / 1e9:g} GHz)",
                ylabel="Magnitude (common source reference)", ylim=(0, 1.15 * peak if peak > 0 else 1))
    axes[1].set(ylabel="Wrapped phase (degrees)", xlabel="Absolute propagation delay (ns)", ylim=(-180, 180))
    axes[1].set_xlim(0, max(1.0, delays.max(initial=0) * 1.05))
    for axis in axes:
        axis.grid(True, alpha=0.3)
    if path is not None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(path, dpi=150)
    return figure


def plot_csv(path, output=None):
    """Return a Matplotlib Figure for a saved impulse CSV."""
    return plot(core.load_impulse_csv(path), output)


@dataclass
class ViewOptions:
    """PyVista view settings. Off-screen rendering supports notebooks and screenshots."""

    off_screen: bool = True
    notebook: bool | None = None
    window_size: tuple[int, int] = (1100, 800)
    pattern_scale: float = 0.0  # Automatic scene-relative sizing when zero.
    ray_indices: list[int] = field(default_factory=list)
    color_rays_by_field: bool = True
    field_dynamic_range_db: float = 60.0
    show_unresolved_candidates: bool = False
    show_polarization: bool = True


def _surface_mesh(pv, surface):
    if isinstance(surface, core.Sphere):
        return pv.Sphere(radius=surface.radius, center=surface.center, theta_resolution=64, phi_resolution=32)
    if isinstance(surface, core.Cylinder):
        return pv.Cylinder(center=surface.center, direction=surface.axis, radius=surface.radius,
                           height=2 * surface.half_length, resolution=96, capping=surface.capped)
    if isinstance(surface, core.Disk):
        return pv.Disc(center=surface.center, inner=0, outer=surface.radius, normal=surface.normal, c_res=96)
    if isinstance(surface, core.Rectangle):
        u, v = surface.half_width * surface.u, surface.half_height * surface.v
        c = surface.center
        return pv.PolyData(np.array([c - u - v, c + u - v, c + u + v, c - u + v]), faces=[4, 0, 1, 2, 3])
    if isinstance(surface, core.Triangle):
        return pv.PolyData(np.array([surface.a, surface.b, surface.c]), faces=[3, 0, 1, 2])
    raise TypeError(f"Unsupported surface: {type(surface).__name__}")


def _polylines(pv, paths):
    paths = [np.asarray(points) for points in paths if len(points) >= 2]
    if not paths:
        return None
    connectivity, offset = [], 0
    for points in paths:
        connectivity.extend([len(points), *range(offset, offset + len(points))])
        offset += len(points)
    return pv.PolyData(np.concatenate(paths), lines=np.asarray(connectivity))


def _antenna_mesh(pv, antenna, position, frequency, medium, scale):
    theta, phi = np.meshgrid(np.linspace(0, np.pi, 37), np.linspace(0, 2 * np.pi, 73), indexing="ij")
    directions = np.stack([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)], axis=-1)
    amplitudes = np.array([np.linalg.norm(antenna.farfield(d, frequency, medium))
                           for d in directions.reshape(-1, 3)]).reshape(theta.shape)
    peak = amplitudes.max(initial=0)
    if peak > 0:
        amplitudes /= peak
    points = position + scale * amplitudes[..., None] * directions
    return pv.StructuredGrid(points[..., 0], points[..., 1], points[..., 2]).extract_surface(
        algorithm="dataset_surface")


def visualize(result, *, options=None):
    """Return a customizable PyVista Plotter containing the solved scene.

    Use ``raytracing.show(result)`` for an interactive desktop window.
    Use ``plotter.show(jupyter_backend='static')`` in a notebook,
    ``plotter.screenshot(path)`` to export, and ``plotter.close()`` when finished.
    """
    import pyvista as pv

    options = options or ViewOptions()
    if not np.isfinite(options.pattern_scale) or options.pattern_scale < 0:
        raise ValueError("pattern_scale must be nonnegative; zero is automatic")
    if not np.isfinite(options.field_dynamic_range_db) or options.field_dynamic_range_db <= 0:
        raise ValueError("field_dynamic_range_db must be positive")
    rays = result.rays
    indices = options.ray_indices or list(range(len(rays)))
    if any(i < 0 or i >= len(rays) for i in indices):
        raise IndexError("ray index out of range")
    selected = [rays[i] for i in indices]
    meshes = [_surface_mesh(pv, surface) for surface in result.scene.surfaces]
    bounds = [result.transmitter.position, result.receiver.position]
    for mesh in meshes:
        b = np.asarray(mesh.bounds).reshape(3, 2)
        bounds.extend([b[:, 0], b[:, 1]])
    bounds = np.asarray(bounds)
    low, high = bounds.min(axis=0), bounds.max(axis=0)
    extent = max(float(np.max(high - low)), 1e-3)
    pattern_scale = options.pattern_scale or min(0.45, 0.15 * extent)
    plotter = pv.Plotter(off_screen=options.off_screen, notebook=options.notebook, window_size=options.window_size)
    plotter.set_background("#141923")
    for i, mesh in enumerate(meshes):
        plotter.add_mesh(mesh, color="#a6b3cc", opacity=0.22, name=f"surface-{i}")

    paths = _polylines(pv, [ray.vertices for ray in selected])
    peak = max((np.linalg.norm(ray.field.receiver_field) for ray in rays if ray.field is not None), default=0)
    if paths is not None:
        if options.color_rays_by_field:
            strengths = []
            for ray in selected:
                if ray.field is None:
                    strengths.append(np.nan)
                else:
                    amplitude = np.linalg.norm(ray.field.receiver_field)
                    strengths.append(max(-options.field_dynamic_range_db, 20 * np.log10(amplitude / peak))
                                     if amplitude > 0 and peak > 0 else -options.field_dynamic_range_db)
            paths.cell_data["Rx |E| (dB)"] = strengths
            plotter.add_mesh(paths, scalars="Rx |E| (dB)", preference="cell", line_width=2,
                             clim=(-options.field_dynamic_range_db, 0), cmap="viridis", nan_color="grey",
                             scalar_bar_args={"title": "Rx field (dB)\nrelative to peak", "vertical": True,
                                              "color": "white",
                                              "position_x": 0.84, "position_y": 0.25, "width": 0.10,
                                              "height": 0.55, "title_font_size": 14, "label_font_size": 12})
        else:
            plotter.add_mesh(paths, color="cyan", line_width=2)

    if options.show_unresolved_candidates:
        unresolved = _polylines(pv, [c.coarse_vertices for c in result.candidates
                                     if c.refinement.status in (core.RefinementStatus.UNRESOLVED_CORNER,
                                                                core.RefinementStatus.UNRESOLVED_EDGE)])
        if unresolved is not None:
            plotter.add_mesh(unresolved, color="orange", line_width=1)
    if options.show_polarization:
        phase = np.linspace(0, 2 * np.pi, 65)
        ellipses = []
        for ray in selected:
            if ray.field is None:
                continue
            e = ray.field.receiver_field
            magnitude = np.linalg.norm(e)
            if magnitude > 0 and np.isfinite(magnitude):
                ellipses.append(result.receiver.position + min(0.12, 0.1 * extent) / magnitude *
                                (np.cos(phase)[:, None] * e.real - np.sin(phase)[:, None] * e.imag))
        glyphs = _polylines(pv, ellipses)
        if glyphs is not None:
            plotter.add_mesh(glyphs, color="gold", line_width=2)

    for endpoint, color, opacity in [(result.transmitter, "tomato", 0.7), (result.receiver, "limegreen", 0.45)]:
        if endpoint.antenna.kind == core.AntennaKind.ISOTROPIC:
            mesh = pv.Sphere(radius=min(result.settings.reception_radius, 0.025 * extent), center=endpoint.position)
        else:
            mesh = _antenna_mesh(pv, endpoint.antenna, endpoint.position, result.settings.frequency_hz,
                                 result.settings.medium, pattern_scale)
        plotter.add_mesh(mesh, color=color, opacity=opacity)
    plotter.add_text(f"Tx: red; Rx: green; PEC: grey\nGold: polarization; grey rays: field unavailable\n"
                     f"{len(selected)} paths; peak Rx field: {peak:.4g} V/m", position="lower_left",
                     font_size=11, color="white")
    plotter.add_axes(xlabel="x (m)", ylabel="y (m)", zlabel="z (m)", color="white")
    center = (low + high) / 2
    plotter.camera_position = [center + extent * np.array([1.5, 1.2, 1.5]), center, [0, 0, 1]]
    plotter.reset_camera()
    return plotter


def visualize_h5(path, *, options=None):
    """Load and return a PyVista scene without rerunning the numerical solver."""
    return visualize(core.load_h5(path), options=options)


def show(result, *, options=None):
    """Open an interactive desktop scene and wait until its window is closed.

    Left-drag rotates, middle-drag or Shift+left-drag pans, and the mouse wheel
    zooms. Press R to reset the camera. Run from a standalone Python script.
    Other view settings are preserved; desktop rendering is always enabled.
    """
    desktop_options = replace(options or ViewOptions(), off_screen=False, notebook=False)
    viewer = visualize(result, options=desktop_options)
    try:
        viewer.enable_trackball_style()
        viewer.add_text("Drag: rotate   Shift+drag: pan   Wheel: zoom   R: reset",
                        position="upper_left", font_size=11, color="white")
        viewer.show(title="RayTracing — interactive scene")
    finally:
        viewer.close()
