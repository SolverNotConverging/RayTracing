"""Multipath impulse response between two antennas in a corner reflector."""
import raytracing as rt

scene = rt.Scene()
scene.add(rt.Rectangle([2, 0, 0], [0, 1, 0], [0, 0, 1], 3, 1), rt.Material.pec())
scene.add(rt.Rectangle([0, 2, 0], [1, 0, 0], [0, 0, 1], 3, 1), rt.Material.conductor(5.8e7))

tx = rt.Transmitter([0, 0, 0], rt.ThinWireDipole(rt.SPEED_OF_LIGHT / 77e9 / 2))
rx = rt.Receiver([1, 0.5, 0], rt.Isotropic())
result = rt.solve(scene, tx, rx)

for tap in result.impulse_response.taps:
    print(f"{tap.delay * 1e9:.3f} ns  {tap.coefficient:.4g}")

rt.plot_impulse_response(result.impulse_response)
rt.show_paths(result)
