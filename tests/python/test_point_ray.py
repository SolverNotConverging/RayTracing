"""Point-ray solver through the Python API."""
import numpy as np
import pytest
import raytracing as rt

C = rt.SPEED_OF_LIGHT


def direct_result(distance=3.0):
    return rt.solve(rt.Scene(), rt.Transmitter([0, 0, 0], rt.Isotropic()),
                    rt.Receiver([distance, 0, 0], rt.Isotropic()),
                    rt.SolverConfig(ray_count=64, max_reflections=0, max_distance=4))


@pytest.mark.parametrize("distance", [0.05, 0.12, 0.121, 3.0])
def test_direct_path(distance):
    result = direct_result(distance)
    response = result.impulse_response
    assert len(result.paths) == len(response.taps) == 1
    assert response.delays[0] == pytest.approx(distance / C)
    assert abs(response.coefficients[0]) == pytest.approx(1 / distance)
    assert response.frequency_response() == pytest.approx(response.coefficients[0])


def test_wall_reflection_and_receiver_swap():
    scene = rt.Scene()
    index = scene.add(rt.Rectangle([1, 0, 0], [0, 1, 0], [0, 0, 1], 5, 5), rt.Material.conductor(5.8e7))
    assert index == 0 and len(scene) == 1
    result = rt.solve(scene, rt.Transmitter([0, 0, 0], rt.ShortDipole()), rt.Receiver([0, 2, 0]),
                      rt.SolverConfig(ray_count=4000, max_reflections=1))
    assert len(result.paths) == 2
    def reflected():
        return next(p for p in result.paths if len(p.refinement.geometry.reflections) == 1)

    assert reflected().field.path_length == pytest.approx(2 * np.sqrt(2))
    assert reflected().vertices.shape == (3, 3)
    before = reflected().coefficient
    rt.set_receiver_antenna(result, rt.ShortDipole())
    assert reflected().coefficient != pytest.approx(before)


def test_antenna_and_material_api():
    antenna = rt.Antenna(rt.ThinWireDipole(C / 77e9 / 2))
    assert isinstance(antenna.model, rt.ThinWireDipole)
    antenna.rotate([0, 1, 0], 90)
    assert np.allclose(antenna.orientation @ [0, 0, 1], [1, 0, 0])
    glass = rt.Material.dielectric(4)
    r = rt.fresnel(rt.Material.vacuum(), glass, 1.0, 77e9)
    assert r.rs == pytest.approx(-1 / 3) and r.ts == pytest.approx(2 / 3)
    with pytest.raises(ValueError):
        rt.Scene().add(rt.Sphere([0, 0, 0], -1), glass)


def test_views_render_off_screen(tmp_path):
    pytest.importorskip("pyvista")
    import matplotlib
    matplotlib.use("Agg", force=True)
    result = direct_result()
    figure = rt.plot_impulse_response(result.impulse_response)
    figure.savefig(tmp_path / "impulse.png")
    plotter = rt.show_paths(result, show=False, screenshot=tmp_path / "paths.png", window_size=(320, 240))
    plotter.close()
    assert (tmp_path / "impulse.png").stat().st_size > 0 and (tmp_path / "paths.png").stat().st_size > 0
