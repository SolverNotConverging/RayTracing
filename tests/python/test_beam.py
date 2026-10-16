"""Gaussian-beam solver through the Python API."""
import numpy as np
import pytest
import raytracing as rt
from raytracing import beam


def test_gaussian_beam_in_vacuum():
    frequency, w = 10e9, 0.1
    k = 2 * np.pi * frequency / rt.SPEED_OF_LIGHT
    zr, e0 = k * w * w / 2, np.sqrt(4 * rt.FREE_SPACE_IMPEDANCE / (np.pi * w * w))
    result = beam.solve(beam.Scene(), beam.GaussianBeam(waist_radius=w),
                        beam.SolverConfig(frequency_hz=frequency, max_distance=4))
    z = np.array([0.5, 1.0, 2.0])
    points = np.column_stack([np.zeros_like(z), np.zeros_like(z), z])
    expected = e0 / (1 + 1j * z / zr) * np.exp(1j * k * z)
    np.testing.assert_allclose(result.fields_at(points)[:, 0], expected, rtol=1e-12)


def test_plane_wave_through_slab_matches_fabry_perot():
    frequency, thickness, epsilon = 10e9, 0.013, 4 + 0.02j
    k = 2 * np.pi * frequency / rt.SPEED_OF_LIGHT
    scene = beam.Scene()
    scene.add_volume(rt.Box([0, 0, 0.5 + thickness / 2], [2, 2, thickness / 2]), rt.Material.dielectric(epsilon))
    result = beam.solve(scene, beam.PlaneWave(), beam.SolverConfig(frequency_hz=frequency, max_distance=4,
                                                                    max_interactions=32, min_power_fraction=0))
    n = np.sqrt(epsilon)
    p, r = np.exp(1j * k * n * thickness), (n - 1) / (n + 1)
    expected = 2 / (1 + n) * 2 * n / (n + 1) * p / (1 - r * r * p * p) * np.exp(1j * k * (1 - thickness))
    assert result.field_at([0, 0, 1])[0] == pytest.approx(expected, rel=1e-12)
    balance = result.power_balance()
    assert sum(v for name, v in balance.items() if name != "launched") == pytest.approx(result.launched_power)


def test_sources_and_terminations():
    scene = beam.Scene()
    scene.add_termination(rt.Box([0, 0, 1], [1, 1, 0.1]), volume=True)
    result = beam.solve(scene, [beam.GaussianBeam(waist_radius=0.1), beam.PlaneWave([0.5, 0, 0])],
                        beam.SolverConfig(frequency_hz=10e9, max_distance=2))
    assert {t.reason for t in result.terminals} == {beam.TerminalReason.TERMINATED}
    assert result.warnings  # coherent multi-source warning


def test_gaussian_impulse_response_matches_reconstructed_field():
    result = beam.solve(beam.Scene(), beam.GaussianBeam(waist_radius=0.1, launch_distance=0.2),
                        beam.SolverConfig(frequency_hz=10e9, max_distance=3))
    point = [0.03, -0.02, 1.0]
    response = result.impulse_response(point)
    assert isinstance(response, rt.ImpulseResponse)
    assert response.frequency_hz == 10e9
    assert response.delays == pytest.approx([0.8 / rt.SPEED_OF_LIGHT])
    np.testing.assert_allclose(response.taps[0].field, result.field_at(point))
    assert response.frequency_response() == pytest.approx(result.field_at(point)[0])
    assert response.taps[0].paths == [0]
    offset = 1e6
    assert response.frequency_response(offset) == pytest.approx(
        response.coefficients[0] * np.exp(2j * np.pi * offset * response.delays[0]))


def test_slab_impulse_delays_and_echo_amplitudes():
    frequency, thickness, n = 10e9, 0.1, 2.0
    scene = beam.Scene()
    scene.add_volume(rt.Box([0, 0, 0.5], [2, 2, thickness / 2]), rt.Material.dielectric(n * n))
    result = beam.solve(scene, beam.PlaneWave(), beam.SolverConfig(
        frequency_hz=frequency, max_distance=4, max_interactions=10, min_power_fraction=0))
    response = result.impulse_response([0, 0, 1])
    order = np.arange(len(response.taps))
    assert len(order) >= 3
    optical_paths = 1 + (n - 1) * thickness + 2 * n * thickness * order
    np.testing.assert_allclose(response.delays, optical_paths / rt.SPEED_OF_LIGHT, rtol=1e-12)
    r = (n - 1) / (n + 1)
    expected = 4 * n / (1 + n)**2 * r**(2 * order) * np.exp(
        2j * np.pi * frequency / rt.SPEED_OF_LIGHT * optical_paths)
    np.testing.assert_allclose(response.coefficients, expected, rtol=1e-11)
    assert response.frequency_response() == pytest.approx(result.field_at([0, 0, 1])[0])
    for tap in response.taps:
        assert len(tap.paths) == 1
        np.testing.assert_allclose(tap.field, result.segments[tap.paths[0]].field_at([0, 0, 1]))


def test_impulse_response_coherent_grouping_and_complex_analyzer():
    analyzer = np.array([1, 1j, 0]) / np.sqrt(2)
    sources = [beam.PlaneWave(polarization=analyzer),
               beam.PlaneWave([0, 0, 0.1], polarization=analyzer, amplitude=1j)]
    result = beam.solve(beam.Scene(), sources, beam.SolverConfig(frequency_hz=10e9, max_distance=3))
    point = [0, 0, 1]
    separate = result.impulse_response(point, analyzer=analyzer, delay_tolerance=0)
    assert len(separate.taps) == 2
    merged = result.impulse_response(point, analyzer=analyzer, delay_tolerance=1e-9)
    assert len(merged.taps) == 1 and len(merged.taps[0].paths) == 2
    np.testing.assert_allclose(merged.taps[0].field, result.field_at(point))
    assert merged.frequency_response() == pytest.approx(np.vdot(analyzer, result.field_at(point)))
    assert merged.frequency_response() == pytest.approx(separate.frequency_response())


def test_impulse_response_uses_same_segment_and_region_selection_as_field():
    scene = beam.Scene()
    scene.add_volume(rt.Box([0, 0, 0.5], [2, 2, 0.1]), rt.Material.dielectric(4 + 0.02j))
    scene.add_termination(rt.Box([0, 0, 1.5], [2, 2, 0.1]), volume=True)
    result = beam.solve(scene, beam.GaussianBeam(), beam.SolverConfig(frequency_hz=10e9, max_distance=3))
    for z in [-1, 0, 0.4, 0.5, 0.6, 1, 1.5, 4]:
        point = [0, 0, z]
        response = result.impulse_response(point)
        total = sum((tap.field for tap in response.taps), np.zeros(3, dtype=complex))
        np.testing.assert_allclose(total, result.field_at(point), atol=1e-12)
    assert not result.impulse_response([0, 0, 1.5]).taps


def test_aperture_callback_errors_propagate():
    def field(_):
        raise ValueError("invalid measured aperture")

    source = beam.ApertureSource([0, 0, 0], [0, 0, 1], 0.2, 5, 0.1, field)
    with pytest.raises(ValueError, match="invalid measured aperture"):
        beam.solve(beam.Scene(), source, beam.SolverConfig(frequency_hz=10e9))


def test_beam_views_render_off_screen(tmp_path):
    pytest.importorskip("pyvista")
    import matplotlib
    matplotlib.use("Agg", force=True)
    scene = beam.Scene()
    scene.add_surface(rt.Sphere([0, 0, 0], 2), rt.Material.pec())
    source = beam.GaussianBeam(waist_radius=0.1)
    result = beam.solve(scene, source, beam.SolverConfig(frequency_hz=10e9, max_distance=3.9, max_interactions=1))
    response = result.impulse_response([0, 0, 1])
    figure = rt.plot_impulse_response(response)
    figure.savefig(tmp_path / "beam_impulse.png")
    assert (tmp_path / "beam_impulse.png").stat().st_size > 0
    for plotter, name in ((rt.show_beams(result, show=False, screenshot=tmp_path / "beams.png"), "beams.png"),
                          (rt.preview(scene, source, frequency_hz=10e9, show=False,
                                      screenshot=tmp_path / "preview.png"), "preview.png")):
        plotter.close()
        assert (tmp_path / name).stat().st_size > 0
