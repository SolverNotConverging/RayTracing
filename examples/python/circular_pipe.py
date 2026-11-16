"""Manuscript case A: PEC pipe, antenna frame, returned rays and coherent CIR.

Install RayTracing with its viz extra, then run this file with --verify.
Use --no-show for reproducible off-screen figures and data exports.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import platform

import numpy as np
import raytracing as rt
from types import SimpleNamespace
from pec_circular_return import solve_pec_circular_return, fold_vertices

FREQUENCY = 77e9
RADIUS, DEPTH, OFFSET = 0.5, 5.0, 0.3
POSITION = np.array([OFFSET, 0.0, DEPTH])
ORIENTATION = np.array([[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]])
MAX_LENGTH = 2.5 * DEPTH
COLORS = {'Direct': '#50545a', 'Radial': '#2878b5', 'Regular polygon': '#39945b',
          'Star polygon': '#d27627', 'Non-periodic': '#9559aa'}


def make_case():
    scene = rt.Scene()
    # Open cylindrical wall, extended 20 mm above the colocated Tx/Rx.
    scene.add(rt.Cylinder([0, 0, (DEPTH + .02) / 2], [0, 0, 1],
                          RADIUS, (DEPTH + .02) / 2, False), rt.Material.pec())
    scene.add(rt.Disk([0, 0, 0], [0, 0, 1], RADIUS), rt.Material.pec())
    antenna = rt.Antenna(rt.Isotropic(rt.Polarization.VERTICAL, 1 + 0j))
    antenna.orientation = ORIENTATION
    return scene, antenna


def classify(path):
    hits = path.refinement.geometry.reflections
    n = sum(h.surface_index == 0 for h in hits)
    if n == 0:
        return 'Direct', n
    launch = np.array(path.refinement.launch_direction)[:2]
    launch /= np.linalg.norm(launch)
    impact = abs(POSITION[0] * launch[1] - POSITION[1] * launch[0])
    if impact < 1e-7:
        return 'Radial', n
    final = np.array(path.refinement.geometry.final_direction)[:2]
    final /= np.linalg.norm(final)
    if np.linalg.norm(final - launch) < 1e-6:
        winding = round(n * np.arccos(impact / RADIUS) / np.pi)
        primitive_winding = winding // np.gcd(n, winding)
        return ('Regular polygon' if primitive_winding == 1 else 'Star polygon'), n
    return 'Non-periodic', n


def calculate(scene, antenna, launches=8000, reflections=8):
    result = solve_pec_circular_return(rt, DEPTH, RADIUS, OFFSET, FREQUENCY,
                                      launches, reflections, MAX_LENGTH, .48*RADIUS)
    paths = [p for p in result.paths if p.field is not None
             and 2*DEPTH-1e-10 <= p.field.path_length <= MAX_LENGTH+1e-10]
    folded = [SimpleNamespace(field=p.field, coefficient=p.coefficient,
                              refinement=p.refinement, vertices=fold_vertices(p.vertices)) for p in paths]
    return sorted(folded,key=lambda p:(p.field.path_length,*p.refinement.launch_direction))


def group_paths(paths):
    groups = []
    for i, path in enumerate(paths):
        delay = path.field.delay
        if groups and abs(delay - groups[-1]['delay_s']) <= 1e-17:
            groups[-1]['indices'].append(i)
            groups[-1]['coefficient'] += complex(path.coefficient)
        else:
            groups.append(dict(delay_s=delay, indices=[i], coefficient=complex(path.coefficient)))
    return groups


def verify_paths(reference, candidate):
    assert len(reference) == len(candidate) == 54
    errors, length_errors, used = [], [], set()
    scale = max(abs(p.coefficient) for p in reference)
    directions = np.array([p.refinement.launch_direction for p in candidate])
    for path in reference:
        distance = np.linalg.norm(directions - path.refinement.launch_direction, axis=1)
        i = int(np.argmin(distance))
        assert distance[i] < 1e-5 and i not in used
        used.add(i)
        length_errors.append(abs(path.field.path_length - candidate[i].field.path_length))
        errors.append(abs(path.coefficient - candidate[i].coefficient) /
                      max(abs(path.coefficient), 1e-8 * scale))
    assert max(length_errors) < 1e-8 and max(errors) < 1e-4
    return dict(paths=len(reference), max_length_error_m=max(length_errors),
                max_relative_complex_error=max(errors))


def export_data(paths, groups, output, verification):
    delays = np.array([g['delay_s'] for g in groups])
    coefficients = np.array([g['coefficient'] for g in groups])
    assert np.isfinite(coefficients).all()
    np.savez_compressed(output / 'impulse_response.npz', delay_s=delays,
                        apparent_range_m=delays * rt.SPEED_OF_LIGHT / 2,
                        coefficients_native=coefficients,
                        coefficients_positive_time=coefficients.conj(),
                        multiplicity=np.array([len(g['indices']) for g in groups]),
                        frequency_hz=FREQUENCY)
    with (output / 'paths.csv').open('w', newline='') as file:
        writer = csv.writer(file, lineterminator="\n")
        writer.writerow(['path', 'family', 'wall_reflections', 'length_m', 'delay_ns',
                         'coefficient_real_native', 'coefficient_imag_native', 'caustic_count'])
        for i, path in enumerate(paths):
            family, n = classify(path)
            writer.writerow([i, family, n, path.field.path_length, path.field.delay * 1e9,
                             path.coefficient.real, path.coefficient.imag, path.field.caustic_count])
    (output / 'paths.json').write_text(json.dumps([
        dict(index=i, family=classify(p)[0], vertices_m=np.asarray(p.vertices).tolist())
        for i, p in enumerate(paths)], indent=2) + '\n', encoding='utf-8', newline='\n')
    manifest = dict(case='A', frequency_hz=FREQUENCY, radius_m=RADIUS, depth_m=DEPTH,
                    offset_m=OFFSET, orientation_local_to_world=ORIENTATION.tolist(),
                    ray_count=8000, max_reflections=8, total_reflection_limit=9, reception_radius_m=.48 * RADIUS,
                    path_length_window_m=[2 * DEPTH, MAX_LENGTH], directed_paths=len(paths),
                    delay_groups=len(groups), verification=verification,
                    raytracing_version=rt.__version__, python=platform.python_version(),
                    phasor_convention='native exp(-i omega t); conjugate coefficients for manuscript',
                    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    adapter_sha256=hashlib.sha256(Path(__file__).with_name('pec_circular_return.py').read_bytes()).hexdigest())
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8', newline='\n')


def pattern_mesh(antenna, center=(0, 0, 0), scale=1.):
    import pyvista as pv
    sphere = pv.Sphere(radius=1, theta_resolution=72, phi_resolution=49)
    directions = np.asarray(sphere.points, dtype=float)
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    fields = np.array([antenna.farfield(u, FREQUENCY) for u in directions])
    magnitude = np.linalg.norm(fields, axis=1)
    assert np.allclose(magnitude, 1, atol=1e-12)
    sphere.points = np.asarray(center) + scale * magnitude[:, None] * directions
    sphere['field_magnitude'] = magnitude
    sphere['polarization'] = fields.real
    return sphere


def add_pipe(plotter, paths, antenna, selected=None, reflector_label='PEC reflector'):
    import pyvista as pv
    wall = pv.Cylinder(center=[0, 0, DEPTH / 2], direction=[0, 0, 1],
                       radius=RADIUS, height=DEPTH, capping=False, resolution=120)
    plotter.add_mesh(wall, color='#a8c6d8', opacity=.12)
    plotter.add_mesh(pv.Disc(center=[0, 0, 0], normal=[0, 0, 1], outer=RADIUS,
                             inner=0, c_res=100), color='#586675', opacity=.45)
    for z in [0, DEPTH]:
        t = np.linspace(0, 2 * np.pi, 121)
        points = np.column_stack([RADIUS * np.cos(t), RADIUS * np.sin(t), t * 0 + z])
        plotter.add_mesh(pv.lines_from_points(points), color='#728896', line_width=1.5)
    for i, path in enumerate(paths):
        active = selected is None or i in selected
        family = classify(path)[0]
        plotter.add_mesh(pv.lines_from_points(path.vertices),
                         color=COLORS[family] if active else '#b6bdc5',
                         line_width=3.5 if active else .8, opacity=1 if active else .12,
                         name=f'path_{i}')
    plotter.add_mesh(pattern_mesh(antenna, POSITION, .17), color='#e7b94c', opacity=.65)
    plotter.add_mesh(pv.Sphere(radius=.025, center=POSITION), color='#1a1e26')
    plotter.add_point_labels([POSITION + [.15, .2, .25], [0, 0, -.25]],
                            ['Tx / Rx', reflector_label], font_size=28,
                            show_points=False, always_visible=True, shape=None)
    plotter.add_mesh(pv.Arrow(start=[.75, 0, DEPTH], direction=[0, 0, -1], scale=DEPTH,
                             shaft_radius=.003, tip_radius=.012, tip_length=.035), color='#586675')
    plotter.add_point_labels([[.85, 0, DEPTH / 2]], ['d'], font_size=28,
                            show_points=False, always_visible=True, shape=None)
    plotter.camera_position = [(7, -11, 5), (0, 0, 2.5), (0, 0, 1)]
    plotter.enable_parallel_projection()
    plotter.reset_camera()
    plotter.camera.zoom(.92)


def add_pattern(plotter, antenna):
    import pyvista as pv
    plotter.add_mesh(pattern_mesh(antenna), color='#e7b94c', opacity=.38,
                     smooth_shading=True)
    # Polarization vectors are computed by the same library antenna as the solver.
    for theta in [np.pi / 4, np.pi / 2, 3 * np.pi / 4]:
        for phi in np.arange(0, 2 * np.pi, np.pi / 2):
            u = ORIENTATION @ [np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)]
            field = np.real(antenna.farfield(u, FREQUENCY))
            plotter.add_mesh(pv.Arrow(start=u - .16 * field, direction=field, scale=.32),
                             color='#197f93')
    plotter.camera_position = [(4, -6, 3), (0, 0, 0), (0, 0, 1)]
    plotter.enable_parallel_projection()
    plotter.reset_camera()
    plotter.camera.zoom(.88)


def representative_paths(paths):
    # Show one actual SBR path per family in the static overview; CIR uses all paths.
    selected = []
    for family in COLORS:
        candidates = [i for i, p in enumerate(paths) if classify(p)[0] == family]
        if candidates:
            selected.append(candidates[0])
    return selected


def figures(paths, groups, antenna, output):
    import pyvista as pv
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    for name, draw, size in [('pipe_scene', lambda p: add_pipe(p, paths, antenna, representative_paths(paths)), (800, 1200)),
                             ('antenna_frame', lambda p: add_pattern(p, antenna), (850, 800))]:
        plotter = pv.Plotter(off_screen=True, window_size=size)
        plotter.set_background('white')
        draw(plotter)
        plotter.screenshot(output / f'{name}.png')
        plotter.close()
    pattern_mesh(antenna, POSITION, .17).save(output / 'antenna_pattern.vtp')
    blocks = pv.MultiBlock({f'path_{i}': pv.lines_from_points(p.vertices) for i, p in enumerate(paths)})
    blocks.save(output / 'returned_paths.vtm')
    plt.rcParams.update({'font.size': 9, 'pdf.fonttype': 42})
    fig = plt.figure(figsize=(7.16, 3.7), layout='constrained')
    grid = fig.add_gridspec(2, 2, width_ratios=[1.7, 1], height_ratios=[1.65, 1])
    for ax, name, title in [(fig.add_subplot(grid[0, 0]), 'pipe_scene', '(a) Pipe and representative returned rays'),
                             (fig.add_subplot(grid[0, 1]), 'antenna_frame', '(b) Antenna pattern and polarization')]:
        picture = plt.imread(output / f'{name}.png')
        yy, xx = np.where(np.min(picture[:, :, :3], axis=2) < .95)
        picture = picture[max(0, yy.min()-20):yy.max()+21, max(0, xx.min()-20):xx.max()+21]
        ax.imshow(picture)
        ax.axis('off')
        ax.set_title(title, fontsize=8, loc='left')
    ax = fig.add_subplot(grid[1, :])
    ranges = np.array([g['delay_s'] * rt.SPEED_OF_LIGHT / 2 for g in groups])
    amplitudes = np.array([abs(g['coefficient']) for g in groups])
    for g, x, y in zip(groups, ranges, amplitudes):
        family = classify(paths[g['indices'][0]])[0]
        ax.vlines(x, 0, y, color=COLORS[family], lw=1)
        ax.plot(x, y, 'o', color=COLORS[family], ms=3)
    ax.set(xlabel='Apparent range (m)', ylabel='Coherent field\nmagnitude', xlim=(4.97, 6.28),
           ylim=(0, amplitudes.max() * 1.35))
    ax.set_title('(c) Impulse response: coherent sum of all 54 directed returns', fontsize=8, loc='left')
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(axis='y', alpha=.2)
    handles = [Line2D([], [], color=color, label=family, marker='o', ms=3, lw=1)
               for family, color in COLORS.items()]
    ax.legend(handles=handles, ncol=5, loc='upper center', fontsize=6.4, frameon=False)
    for ext in ['png', 'pdf']:
        fig.savefig(output / f'circular_pipe.{ext}', dpi=300)
    plt.close(fig)
    # Geometry-only manuscript panel; its subsequent figure already compares CIRs.
    fig, axes = plt.subplots(1, 2, figsize=(3.5, 3.35), width_ratios=[2.1, 1], layout='constrained')
    for ax, name, title in zip(axes, ['pipe_scene', 'antenna_frame'],
                              ['(a) Returned rays', '(b) Antenna pattern']):
        picture = plt.imread(output / f'{name}.png')
        yy, xx = np.where(np.min(picture[:, :, :3], axis=2) < .95)
        ax.imshow(picture[max(0, yy.min()-20):yy.max()+21, max(0, xx.min()-20):xx.max()+21])
        ax.axis('off')
        ax.set_title(title, loc='left', fontsize=7)
    fig.legend(handles=handles, ncol=3, loc='lower center', fontsize=6.2, frameon=False,
               handlelength=1.3, columnspacing=.8)
    fig.get_layout_engine().set(rect=(0, .12, 1, .88))
    for ext in ['png', 'pdf']:
        fig.savefig(output / f'circular_pipe_geometry.{ext}', dpi=300)
    plt.close(fig)


def inspect(paths, groups, antenna, show=True):
    """Single VTK window: scene, antenna frame and CIR, linked by a group slider."""
    import pyvista as pv
    plotter = pv.Plotter(shape='1|2', window_size=(1500, 900), off_screen=not show,
                        title='Circular pipe — manuscript case A')
    plotter.set_background('white', all_renderers=True)
    plotter.subplot(0)
    add_pipe(plotter, paths, antenna, groups[0]['indices'], reflector_label='PEC')
    plotter.camera.zoom(.85)
    plotter.add_text('Case B: D = 1 m, d = 5 m, a = 0.3 m, 77 GHz', font_size=12)
    plotter.subplot(1)
    add_pattern(plotter, antenna)
    plotter.add_text('Unit field magnitude; arrows show polarization', font_size=11)
    plotter.subplot(2)
    chart = pv.Chart2D(x_label='Apparent range (m)', y_label='Coherent field magnitude')
    ranges = np.array([g['delay_s'] * rt.SPEED_OF_LIGHT / 2 for g in groups])
    magnitude = np.array([abs(g['coefficient']) for g in groups])
    for g, x, y in zip(groups, ranges, magnitude):
        color = COLORS[classify(paths[g['indices'][0]])[0]]
        chart.line([x, x], [0, y], color=color, width=2)
    marker = chart.scatter([ranges[0]], [magnitude[0]], color='black', size=12)
    chart.x_range = (4.95, 6.3)
    chart.y_range = (0, float(magnitude.max() * 1.2))
    for axis in [chart.x_axis, chart.y_axis]:
        axis.label_size = 15
        axis.tick_label_size = 12
    plotter.add_chart(chart)

    def select(value):
        i = min(len(groups) - 1, max(0, int(round(value))))
        group = groups[i]
        plotter.subplot(0)
        for j, path in enumerate(paths):
            actor = plotter.renderer.actors[f'path_{j}']
            active = j in group['indices']
            actor.prop.color = COLORS[classify(path)[0]] if active else '#b6bdc5'
            actor.prop.opacity = 1 if active else .12
            actor.prop.line_width = 4 if active else .8
        marker.update([ranges[i]], [magnitude[i]])
        value_text = (f'Group {i + 1}/{len(groups)}: {len(group["indices"])} directed rays\n'
                      f'R = {ranges[i]:.6f} m; delay = {group["delay_s"] * 1e9:.6f} ns\n'
                      f'Magnitude = {magnitude[i]:.6g}\n'
                      f'Phase = {np.angle(group["coefficient"], deg=True):.2f} deg (native)')
        plotter.add_text(value_text, position='lower_left', font_size=11, name='selection')
    plotter.subplot(0)
    plotter.add_slider_widget(select, [0, len(groups) - 1], value=0,
                              title='Delay group', pointa=(.2, .12), pointb=(.8, .12),
                              style='modern', slider_width=.02, tube_width=.004, title_height=.02,
                              interaction_event='always')
    if show:
        plotter.show()
    return plotter, select


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('circular_pipe_output'))
    parser.add_argument('--verify', action='store_true', help='Double launches and increase reflection limit')
    parser.add_argument('--no-show', action='store_true', help='Export figures and data without opening a window')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    scene, antenna = make_case()
    paths = calculate(scene, antenna)
    assert len(paths) == 54, f'Expected 54 returns in case A; found {len(paths)}'
    verification = {}
    if args.verify:
        verification['16000_launches'] = verify_paths(paths, calculate(scene, antenna, 16000, 8))
        verification['10_wall_reflections'] = verify_paths(paths, calculate(scene, antenna, 8000, 10))
    groups = group_paths(paths)
    export_data(paths, groups, args.output, verification)
    figures(paths, groups, antenna, args.output)
    print(json.dumps(dict(paths=len(paths), delay_groups=len(groups), output=str(args.output),
                          verification=verification), indent=2))
    if not args.no_show:
        inspect(paths, groups, antenna)


if __name__ == '__main__':
    main()
