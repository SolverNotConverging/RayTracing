"""Numerical parity between the Python API and a separately linked C++ program."""
import csv
import io
import os
from pathlib import Path
import subprocess
import numpy as np
import pytest
import raytracing as rt


def test_cpp_python_gaussian_slab_agree():
    executable = os.environ.get("RAYTRACING_BEAM_EXAMPLE")
    if executable is None:
        pytest.skip("Set RAYTRACING_BEAM_EXAMPLE to the gaussian_beam_example executable")
    completed = subprocess.run([str(Path(executable).resolve())], capture_output=True, text=True, check=True)
    rows = list(csv.DictReader(io.StringIO(completed.stdout)))
    points = [[float(row[k]) for k in ("x", "y", "z")] for row in rows]
    expected = np.asarray([[complex(float(row[f"E{axis}_real"]), float(row[f"E{axis}_imag"]))
                            for axis in "xyz"] for row in rows])
    scene = rt.Scene()
    scene.add_volume(rt.Box([0, 0, .5065], [2, 2, .0065]), rt.Material.dielectric(4+.02j, name="glass"))
    result = rt.solve(scene, rt.GaussianBeam([0, 0, 0], [0, 0, 1], .1),
                      rt.SolverConfig(frequency_hz=10e9, max_distance=4, max_interactions=32, min_power_fraction=0))
    np.testing.assert_allclose(result.fields_at(points), expected, rtol=1e-12, atol=1e-12)


def test_native_beam_inputs_are_snapshots():
    scene = rt.Scene()
    source = rt.GaussianBeam(np.zeros(3), [0, 0, 1], .1)
    config = rt.SolverConfig(frequency_hz=10e9, max_distance=2)
    result = rt.solve(scene, source, config)
    before = result.field_at([.02, 0, 1])
    source.waist_position[:] = 10
    config.frequency_hz = 20e9
    scene.add_volume(rt.Box([0, 0, 1], [1, 1, .1]), rt.Material.pec())
    np.testing.assert_array_equal(result.field_at([.02, 0, 1]), before)


def test_aperture_callback_failure_is_preserved():
    def field(_):
        raise ValueError("invalid measured aperture")
    source = rt.ApertureSource([0, 0, 0], [0, 0, 1], .2, 5, .1, field)
    with pytest.raises(ValueError, match="invalid measured aperture"):
        rt.solve(rt.Scene(), source, rt.SolverConfig(frequency_hz=10e9))


def test_sampling_cache_preserves_result_edits():
    result = rt.solve(rt.Scene(), rt.GaussianBeam([0, 0, 0], [0, 0, 1], .1),
                      rt.SolverConfig(frequency_hz=10e9, max_distance=2))
    point = [.02, 0, 1]
    before = result.field_at(point)
    cached = result._sampling_native()
    np.testing.assert_array_equal(result.field_at(point), before)
    assert result._sampling_native() is cached
    result.segments[0].state.field *= 2
    np.testing.assert_allclose(result.field_at(point), 2 * before)
    assert result._sampling_native() is not cached
    np.testing.assert_allclose(result.fields_at([point, point]), [2 * before, 2 * before])
    result.scene.add_volume(rt.Box(point, [.01, .01, .01]), rt.Material.pec())
    np.testing.assert_array_equal(result.field_at(point), np.zeros(3))


def test_sampling_cache_tracks_segment_list_and_geometry_edits():
    scene = rt.Scene()
    scene.add_volume(rt.Box([1, 0, 1], [.02, .02, .02]), rt.Material.pec())
    result = rt.solve(scene, rt.GaussianBeam([0, 0, 0], [0, 0, 1], .1),
                      rt.SolverConfig(frequency_hz=10e9, max_distance=2))
    point = [.02, 0, 1]
    before = result.field_at(point)
    segment = result.segments.pop()
    np.testing.assert_array_equal(result.field_at(point), np.zeros(3))
    result.segments.append(segment)
    np.testing.assert_array_equal(result.field_at(point), before)
    result.scene._objects[0].shape.center[:] = point
    np.testing.assert_array_equal(result.field_at(point), np.zeros(3))
