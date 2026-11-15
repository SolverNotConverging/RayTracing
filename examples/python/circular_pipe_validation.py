"""Reproduce circular-pipe SBR cases against frozen analytical reference fields.

Requires the installed RayTracing package and NumPy, without visualization extras.
References are used only after the unseeded SBR solve, never to launch/refine rays.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import csv
import hashlib
import json
import platform
import time

import numpy as np
import raytracing as rt
import raytracing._core as core

REFERENCE = Path(__file__).with_name('circular_pipe_reference')
ORIENTATION = np.array([[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]])
LENGTH_TOL = 1e-8
FIELD_TOL = 1e-4


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def scene_and_ports(case):
    radius, depth, offset = (case[k] for k in ('radius_m', 'depth_m', 'offset_m'))
    scene = rt.Scene()
    # Analytical surfaces in the solver, not the polygonal display meshes.
    scene.add(rt.Cylinder([0, 0, (depth + .02) / 2], [0, 0, 1],
                          radius, (depth + .02) / 2, False), rt.Material.pec())
    scene.add(rt.Disk([0, 0, 0], [0, 0, 1], radius), rt.Material.pec())
    antenna = rt.Antenna(rt.Isotropic(rt.Polarization.VERTICAL, 1 + 0j))
    antenna.orientation = ORIENTATION
    position = [offset, 0, depth]
    return scene, rt.Transmitter(position, antenna), rt.Receiver(position, antenna)


def solve_case(case, launches, reflections):
    scene, tx, rx = scene_and_ports(case)
    radius, depth = case['radius_m'], case['depth_m']
    started = time.perf_counter()
    result = rt.solve(scene, tx, rx, rt.SolverConfig(
        frequency_hz=77e9, ray_count=launches, max_reflections=reflections,
        max_distance=2.5 * depth + .48 * radius, reception_radius=.48 * radius))
    elapsed = time.perf_counter() - started
    paths = [p for p in result.paths if p.field is not None
             and 2 * depth - 1e-10 <= p.field.path_length <= 2.5 * depth + 1e-10
             and sum(h.surface_index == 1 for h in p.refinement.geometry.reflections) == 1]
    return sorted(paths, key=lambda p: (p.field.path_length, *p.refinement.launch_direction)), elapsed


def compare_reference(paths, reference):
    if len(paths) != len(reference):
        raise ValueError(f'Expected {len(reference)} directed paths, found {len(paths)}')
    lengths = np.array([p.field.path_length for p in paths])
    # Published reference fields use exp(+i omega t); the library uses exp(-i omega t).
    coefficients = np.array([complex(p.coefficient).conjugate() for p in paths])
    directions = np.array([p.refinement.launch_direction for p in paths])
    orders = np.array([sum(h.surface_index == 0 for h in p.refinement.geometry.reflections)
                       for p in paths])
    if not all(np.isfinite(x).all() for x in (lengths, coefficients, directions)):
        raise ValueError('Non-finite SBR result')
    separation = np.linalg.norm(directions[:, None] - directions[None, :], axis=2)
    np.fill_diagonal(separation, np.inf)
    if np.min(separation) < 1e-7:
        raise ValueError('Duplicate launch direction')
    expected = np.array([complex(float(r['theory_real']), float(r['theory_imag']))
                         for r in reference])
    scale = np.max(np.abs(expected))
    used = np.zeros(len(paths), dtype=bool)
    matched, field_errors, length_errors = [], [], []
    # Frozen CSVs contain delay, wall order and complex field, not launch directions.
    # Match a multiset within each order/delay cluster; convergence below separately
    # verifies one-to-one launch-direction recovery between SBR runs.
    for row, field in zip(reference, expected):
        length = float(row['delay_ns']) * 1e-9 * rt.SPEED_OF_LIGHT
        eligible = np.flatnonzero(~used & (orders == int(row['wall_reflections']))
                                  & (np.abs(lengths - length) < LENGTH_TOL))
        if not len(eligible):
            raise ValueError(f'No unmatched SBR return for reference ray {row["ray_index"]}')
        index = int(eligible[np.argmin(np.abs(coefficients[eligible] - field))])
        used[index] = True
        field_errors.append(float(abs(coefficients[index] - field) / max(abs(field), 1e-8 * scale)))
        length_errors.append(float(abs(lengths[index] - length)))
        matched.append(dict(reference_ray=int(row['ray_index']), sbr_path=index,
                            length_m=float(lengths[index]),
                            coefficient_real_positive_time=float(coefficients[index].real),
                            coefficient_imag_positive_time=float(coefficients[index].imag),
                            relative_complex_error=field_errors[-1]))
    worst = max(field_errors)
    if worst >= FIELD_TOL:
        raise ValueError(f'Relative complex error {worst:.6g} exceeds {FIELD_TOL}')
    return dict(passed=True, directed_paths=len(paths),
                max_relative_complex_error=worst, max_length_error_m=max(length_errors)), matched


def compare_launches(reference, candidate):
    if len(reference) != len(candidate):
        raise ValueError('Path count changed during convergence check')
    directions = np.array([p.refinement.launch_direction for p in candidate])
    used = set()
    field_errors, length_errors = [], []
    scale = max(abs(p.coefficient) for p in reference)
    for path in reference:
        distances = np.linalg.norm(directions - path.refinement.launch_direction, axis=1)
        index = int(np.argmin(distances))
        if distances[index] >= 1e-5 or index in used:
            raise ValueError('Convergence check has missing or duplicate direction matches')
        used.add(index)
        length_errors.append(abs(path.field.path_length - candidate[index].field.path_length))
        field_errors.append(abs(path.coefficient - candidate[index].coefficient)
                            / max(abs(path.coefficient), 1e-8 * scale))
    if max(length_errors) >= LENGTH_TOL or max(field_errors) >= FIELD_TOL:
        raise ValueError('Convergence error exceeds tolerance')
    return dict(passed=True, max_relative_complex_error=max(field_errors),
                max_length_error_m=max(length_errors))


def export_cir(paths, target):
    delays, coefficients, multiplicities = [], [], []
    for path in paths:
        if delays and abs(path.field.delay - delays[-1]) <= 1e-17:
            coefficients[-1] += complex(path.coefficient)
            multiplicities[-1] += 1
        else:
            delays.append(path.field.delay)
            coefficients.append(complex(path.coefficient))
            multiplicities.append(1)
    np.savez_compressed(target, delay_s=delays,
                        apparent_range_m=np.asarray(delays) * rt.SPEED_OF_LIGHT / 2,
                        coefficients_native=coefficients,
                        coefficients_positive_time=np.conj(coefficients),
                        multiplicity=multiplicities, frequency_hz=77e9)


def main():
    catalog = json.loads((REFERENCE / 'cases.json').read_text())
    cases = {c['case']: c for c in catalog['cases']}
    parser = argparse.ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument('--cases', nargs='+', choices=list(cases), default=None)
    selection.add_argument('--all', action='store_true', help='Run all thirteen manuscript configurations')
    parser.add_argument('--verify', action='store_true', help='Also double launches and add two reflections')
    parser.add_argument('--output', type=Path, default=Path('circular_pipe_validation_output'))
    args = parser.parse_args()
    selected = list(cases) if args.all else (args.cases or ['A', 'B', 'E'])
    args.output.mkdir(parents=True, exist_ok=True)
    report = dict(passed=False, utc=datetime.now(timezone.utc).isoformat(),
                  python=platform.python_version(), platform=platform.platform(),
                  solver_version=rt.__version__, solver_binary_sha256=digest(core.__file__),
                  script_sha256=digest(__file__), reference_catalog_sha256=digest(REFERENCE / 'cases.json'),
                  convergence_requested=args.verify, selected_cases=selected,
                  comparison='Fresh SBR solves against frozen analytical fields; no fitted field scale or phase',
                  reference_phasor='exp(+i omega t)', results=[])
    try:
        for name in selected:
            case = cases[name]
            reference_path = REFERENCE / f'case_{name}_rays.csv'
            if digest(reference_path) != case['reference_sha256']:
                raise ValueError(f'Reference checksum mismatch for {name}')
            with reference_path.open(newline='') as file:
                reference = list(csv.DictReader(file))
            print(f'Case {name}: d={case["depth_m"]:g} m, a={case["offset_m"]:g} m; solving', flush=True)
            paths, elapsed = solve_case(case, case['launch_count'], case['max_reflections'])
            metrics, matches = compare_reference(paths, reference)
            row = dict(case=name, **metrics, solve_seconds=elapsed, convergence={})
            report['results'].append(row)
            directory = args.output / name
            directory.mkdir(exist_ok=True)
            (directory / 'matched_paths.json').write_text(json.dumps(matches, indent=2) + '\n')
            export_cir(paths, directory / 'impulse_response.npz')
            if args.verify:
                for label, launches, reflections in [
                    ('double_launches', 2 * case['launch_count'], case['max_reflections']),
                    ('extra_reflections', case['launch_count'], case['max_reflections'] + 2)]:
                    print(f'Case {name}: {label} ({launches} launches, {reflections} reflections)', flush=True)
                    candidate, seconds = solve_case(case, launches, reflections)
                    check, _ = compare_reference(candidate, reference)
                    row['convergence'][label] = dict(analytical_comparison=check,
                        launch_matching=compare_launches(paths, candidate), solve_seconds=seconds)
            print(f'Case {name}: PASS, {len(paths)} paths, complex error {metrics["max_relative_complex_error"]:.3g}', flush=True)
            (args.output / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
        report['passed'] = True
        report['directed_paths'] = sum(r['directed_paths'] for r in report['results'])
    except Exception as error:
        report['error'] = repr(error)
        raise
    finally:
        report['completed_utc'] = datetime.now(timezone.utc).isoformat()
        (args.output / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
    print(f'PASS: {len(selected)} cases; summary: {args.output / "summary.json"}', flush=True)


if __name__ == '__main__':
    main()
