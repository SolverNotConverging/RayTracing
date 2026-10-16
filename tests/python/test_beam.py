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


def test_aperture_callback_errors_propagate():
    def field(_):
        raise ValueError("invalid measured aperture")

    source = beam.ApertureSource([0, 0, 0], [0, 0, 1], 0.2, 5, 0.1, field)
    with pytest.raises(ValueError, match="invalid measured aperture"):
        beam.solve(beam.Scene(), source, beam.SolverConfig(frequency_hz=10e9))


def test_beam_views_render_off_screen(tmp_path):
    pytest.importorskip("pyvista")
    scene = beam.Scene()
    scene.add_surface(rt.Sphere([0, 0, 0], 2), rt.Material.pec())
    source = beam.GaussianBeam(waist_radius=0.1)
    result = beam.solve(scene, source, beam.SolverConfig(frequency_hz=10e9, max_distance=3.9, max_interactions=1))
    for plotter, name in ((rt.show_beams(result, show=False, screenshot=tmp_path / "beams.png"), "beams.png"),
                          (rt.preview(scene, source, frequency_hz=10e9, show=False,
                                      screenshot=tmp_path / "preview.png"), "preview.png")):
        plotter.close()
        assert (tmp_path / name).stat().st_size > 0
