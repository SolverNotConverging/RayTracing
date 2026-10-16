"""Point-ray bindings, shared result IO and included Python rendering."""
import numpy as np
import pytest
from raytracing import point_ray as rt


@pytest.fixture
def direct():
    return rt.solve(rt.Scene(), rt.Transmitter([0, 0, 0], rt.Isotropic()),
                    rt.Receiver([3, 0, 0], rt.Isotropic()),
                    rt.SolverConfig(ray_count=64, max_reflections=0, max_distance=4))


def test_direct_impulse_physics(direct):
    assert len(direct.rays) == len(direct.impulse_response.taps) == 1
    tap = direct.impulse_response.taps[0]
    assert tap.delay_seconds == pytest.approx(3 / 299792458.0)
    assert abs(tap.coefficient) == pytest.approx(1 / 3)
    assert rt.frequency_response(direct.impulse_response) == pytest.approx(tap.coefficient)


@pytest.mark.parametrize("distance", [.05, .12, .121])
def test_nearby_direct_receiver(distance):
    result = rt.solve(rt.Scene(), rt.Transmitter([0, 0, 0], rt.Isotropic()),
                      rt.Receiver([distance, 0, 0], rt.Isotropic()),
                      rt.SolverConfig(ray_count=16, max_reflections=0, max_distance=1))
    assert len(result.rays) == len(result.candidates) == 1
    tap = result.impulse_response.taps[0]
    assert tap.delay_seconds == pytest.approx(distance / 299792458.0)
    assert abs(tap.coefficient) == pytest.approx(1 / distance)


def test_point_ray_persistence(direct, tmp_path):
    rt.save_h5(direct, tmp_path / "simulation.h5")
    restored = rt.load_h5(tmp_path / "simulation.h5")
    rt.save_impulse_csv(direct.impulse_response, tmp_path / "impulse.csv")
    response = rt.load_impulse_csv(tmp_path / "impulse.csv")
    for actual in (restored.impulse_response, response):
        assert actual.taps[0].coefficient == pytest.approx(direct.impulse_response.taps[0].coefficient)
        assert actual.taps[0].delay_seconds == direct.impulse_response.taps[0].delay_seconds


def test_included_python_visualization(direct, tmp_path):
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    import pyvista
    from PIL import Image

    figure = rt.plot(direct.impulse_response, tmp_path / "impulse.png")
    plt.close(figure)
    viewer = rt.visualize(direct, options=rt.ViewOptions(off_screen=True, window_size=(320, 240)))
    try:
        assert isinstance(viewer, pyvista.Plotter)
        pixels = viewer.screenshot(str(tmp_path / "scene.png"))
        assert pixels.shape[:2] == (240, 320)
        assert np.std(pixels) > 0
    finally:
        viewer.close()
    for name in ("impulse.png", "scene.png"):
        with Image.open(tmp_path / name) as image:
            image.verify()
