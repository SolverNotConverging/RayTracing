"""Gaussian beam through a glass slab; sample the field behind it."""
import raytracing as rt
from raytracing import beam

scene = beam.Scene()
scene.add_volume(rt.Box([0, 0, 0.5], [1, 1, 0.01]), rt.Material.dielectric(4 + 0.02j))

source = beam.GaussianBeam([0, 0, 0], [0, 0, 1], waist_radius=0.1)
result = beam.solve(scene, source, beam.SolverConfig(frequency_hz=10e9))

print("E(0, 0, 1) =", result.field_at([0, 0, 1]))
rt.show_beams(result)
