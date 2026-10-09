import gc
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

import raytracing as rt


@pytest.fixture(scope="module")
def patch():
    return rt.load_farfield(rt.example_pattern())


def direct_result():
    return rt.solve(rt.Scene(), rt.Transmitter([0, 0, 0], rt.Isotropic()),
                    rt.Receiver(np.array([3.0, 0, 0]), rt.Isotropic()),
                    rt.SolverConfig(ray_count=32, max_reflections=0))


def test_direct_numpy_and_complex_fields():
    result = direct_result()
    assert result.launched_rays == 33
    assert len(result.rays) == 1
    ray = result.rays[0]
    np.testing.assert_allclose(ray.vertices, [[0, 0, 0], [3, 0, 0]])
    assert ray.field.receiver_field.dtype == np.complex128
    assert abs(rt.frequency_response(result.impulse_response)) == pytest.approx(1 / 3)
    assert result.impulse_response.taps[0].delay_seconds == pytest.approx(3 / rt.SPEED_OF_LIGHT)
    # Returned nested values and arrays must outlive the temporary owning result.
    field = ray.field.receiver_field
    del result, ray
    gc.collect()
    assert np.linalg.norm(field) == pytest.approx(1 / 3)


def test_all_shapes_and_validation():
    scene = rt.Scene()
    shapes = [rt.Rectangle([0, 0, 0], [1, 0, 0], [0, 1, 0], 1, 2),
              rt.Disk([0, 0, 0], [0, 0, 1], 1), rt.Sphere([0, 0, 0], 1),
              rt.Cylinder([0, 0, 0], [0, 0, 1], 1, 2, capped=False),
              rt.Triangle([0, 0, 0], [1, 0, 0], [0, 1, 0])]
    for shape in shapes:
        scene.add(shape)
    assert len(scene) == 5
    assert [type(s) for s in scene.surfaces] == [type(s) for s in shapes]
    assert not scene.surfaces[3].capped
    with pytest.raises(ValueError):
        scene.add(rt.Sphere([0, 0, 0], -1))
    with pytest.raises(TypeError):
        rt.Transmitter([0, 0])
    with pytest.raises(ValueError):
        rt.solve(rt.Scene(), rt.Transmitter([0, 0, 0]), rt.Receiver([1, 0, 0]),
                 rt.SolverConfig(ray_count=0))


def test_antennas_rotation_and_nested_options(patch):
    assert patch.kind == rt.AntennaKind.IMPORTED
    assert len(patch.pattern.frequencies_hz) == 5
    assert patch.pattern.input_convention == rt.PhasorConvention.POSITIVE_TIME
    original = patch.farfield([0, 0, -1], 77e9)
    rotated = patch.copy().rotate([0, 1, 0], -90)
    assert np.linalg.norm(rotated.farfield([1, 0, 0], 77e9)) == pytest.approx(np.linalg.norm(original))
    np.testing.assert_allclose(patch.orientation, np.eye(3))
    for model in [rt.Isotropic(), rt.ShortDipole(), rt.ThinWireDipole(), rt.RectangularAperture()]:
        rt.Antenna(model).validate(77e9)
    tx, rx = rt.Transmitter([0, 0, 0], patch), rt.Receiver([1, 0, 0], patch)
    tx.antenna.rotate([0, 1, 0], -90)
    np.testing.assert_allclose(rx.antenna.orientation, np.eye(3))
    config = rt.SolverConfig()
    config.refinement.receiver_tolerance = 1e-11
    config.spreading.minimum_singular_value = 2e-7
    assert config.refinement.receiver_tolerance == 1e-11
    assert config.spreading.minimum_singular_value == 2e-7
    assert rt.Medium(relative_permittivity=4).refractive_index() == 2


def test_shifted_cylinder_and_pattern_round_trip(patch, tmp_path):
    scene = rt.Scene()
    scene.add(rt.Cylinder([0, 0, -0.25], [0, 0, 1], 0.127, 0.25, False))
    scene.add(rt.Disk([0, 0, -0.5], [0, 0, -1], 0.127))
    position = [0.005, 0, 0]
    result = rt.solve(scene, rt.Transmitter(position, patch), rt.Receiver(position, patch),
                      rt.SolverConfig(ray_count=1000))
    assert len(result.rays) == len(result.response_ray_indices) == 13
    for ray in result.rays:
        assert ray.refinement.launch_direction[2] < 0
        assert any(hit.surface_index == 1 for hit in ray.refinement.geometry.reflections)
        assert ray.spreading.status == rt.SpreadingStatus.VALID
    rt.save_h5(result, tmp_path / "simulation.h5")
    loaded = rt.load_h5(tmp_path / "simulation.h5")
    assert len(loaded.rays) == 13
    assert rt.frequency_response(loaded.impulse_response) == pytest.approx(rt.frequency_response(result.impulse_response))
    np.testing.assert_allclose(loaded.receiver.antenna.farfield([0, 0, -1], 77e9),
                               patch.farfield([0, 0, -1], 77e9))
    rt.save_impulse_csv(result.impulse_response, tmp_path / "response.csv")
    csv = rt.load_impulse_csv(tmp_path / "response.csv")
    assert rt.frequency_response(csv) == pytest.approx(rt.frequency_response(result.impulse_response))
    incident = loaded.rays[0].field.receiver_field.copy()
    rt.update_receiver(loaded, rt.Isotropic())
    np.testing.assert_array_equal(loaded.rays[0].field.receiver_field, incident)
    assert rt.frequency_response(loaded.impulse_response) != rt.frequency_response(result.impulse_response)


def test_matplotlib_and_pyvista_exports(tmp_path):
    import matplotlib.pyplot as plt
    import pyvista as pv
    from raytracing.visualization import _surface_mesh

    result = direct_result()
    figure = rt.plot(result.impulse_response, tmp_path / "response.png")
    assert len(figure.axes) == 2
    plt.close(figure)
    viewer = rt.visualize(result, options=rt.ViewOptions(notebook=False))
    assert isinstance(viewer, pv.Plotter)
    try:
        viewer.screenshot(tmp_path / "scene.png")
    finally:
        viewer.close()
    for name in ["response.png", "scene.png"]:
        data = (tmp_path / name).read_bytes()
        assert data.startswith(b"\x89PNG\r\n\x1a\n") and len(data) > 1000
    cylinder = _surface_mesh(pv, rt.Cylinder([0, 0, -.25], [0, 0, 1], .127, .25, False))
    np.testing.assert_allclose(cylinder.bounds, [-.127, .127, -.127, .127, -.5, 0], atol=1e-7)
    with pytest.raises(IndexError):
        rt.visualize(result, options=rt.ViewOptions(ray_indices=[1]))


@pytest.mark.skipif(os.name != "nt", reason="Windows wheel runtime check")
def test_installed_package_without_build_machine_path(tmp_path):
    env = dict(os.environ)
    env["PATH"] = str(Path(sys.executable).parent) + os.pathsep + str(Path(os.environ["SYSTEMROOT"]) / "System32")
    env["MPLBACKEND"] = "Agg"
    code = """
import sys
from pathlib import Path
import raytracing as rt
assert 'pyvista' not in sys.modules and 'matplotlib' not in sys.modules
assert rt.example_pattern().is_file()
assert not list((Path(rt.__file__).parent / '.libs').glob('*vtk*'))
assert not (Path(rt.__file__).parent / 'gnuplot').exists()
r = rt.solve(rt.Scene(), rt.Transmitter([0,0,0]), rt.Receiver([1,0,0]), rt.SolverConfig(ray_count=8,max_reflections=0))
figure = rt.plot(r.impulse_response, 'plot.png')
assert Path('plot.png').stat().st_size > 1000
viewer = rt.visualize(r)
viewer.screenshot('scene.png')
viewer.close()
assert Path('scene.png').stat().st_size > 1000
"""
    subprocess.run([sys.executable, "-c", code], cwd=tmp_path, env=env, check=True, timeout=90)
