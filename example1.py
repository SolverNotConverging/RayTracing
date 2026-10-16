"""Mixed-geometry simulation with an analytical half-wave dipole transmitter.

Run with .venv/Scripts/python.exe example1.py to open the desktop viewer.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from raytracing import point_ray as rt


def main(show=True, *, ray_count=20_000, output=None):
    config = rt.SolverConfig(
        frequency_hz=77e9,
        ray_count=ray_count,
        max_reflections=8,
        max_distance=20.0,
    )

    x, y, z = np.eye(3)
    scene = rt.Scene()
    scene.add(rt.Rectangle([1.5, 0, 0], y, z, 3, 1), rt.Material.pec())
    scene.add(rt.Rectangle([-1.5, 0, 0], y, z, 3, 1), rt.Material.pec())
    scene.add(rt.Rectangle([0, -3, 0], x, z, 3, 1), rt.Material.pec())
    scene.add(rt.Rectangle([0, 3, 0], x, z, 3, 1), rt.Material.pec())
    scene.add(rt.Sphere([0.7, -1.2, 0], 0.3), rt.Material.pec())
    scene.add(rt.Disk([-0.7, -1.3, 0], y, 0.35), rt.Material.pec())
    axis = np.array([0.2, 0, 1.0])
    axis /= np.linalg.norm(axis)
    scene.add(rt.Cylinder([0.6, 2, 0], axis, 0.25, 0.7, capped=True), rt.Material.pec())
    scene.add(rt.Triangle([-1, 0, 0.8], [0, 0, 0.8], [-0.5, 1, 0.8]), rt.Material.pec())

    # Analytical half-wave dipole along z, radiating into the horizontal room.
    antenna = rt.Antenna(rt.ThinWireDipole(rt.SPEED_OF_LIGHT / config.frequency_hz / 2))
    tx = rt.Transmitter([0, 0, 0], antenna)
    rx = rt.Receiver([-0.5, 1.4, 0], rt.Isotropic())

    print("Running example 1...", flush=True)
    result = rt.solve(scene, tx, rx, config)
    print(
        f"{len(result.rays)} paths, {len(result.response_ray_indices)} paths with fields, "
        f"{len(result.impulse_response.taps)} impulse taps",
        flush=True,
    )

    output = Path(output) if output is not None else Path(__file__).resolve().parent / "results" / "example1"
    output.mkdir(parents=True, exist_ok=True)
    rt.save_h5(result, output / "simulation.h5")
    rt.save_impulse_csv(result.impulse_response, output / "impulse_response.csv")
    figure = rt.plot(result.impulse_response, output / "impulse_response.png")
    plt.close(figure)
    viewer = rt.visualize(result)
    try:
        viewer.screenshot(str(output / "scene.png"))
    finally:
        viewer.close()
    print(f"Results saved to {output}", flush=True)
    if show:
        print("Drag to rotate; Shift+drag to pan; scroll to zoom. Close the viewer to exit.", flush=True)
        rt.show(result)
    return result


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--no-show', action='store_true')
    parser.add_argument('--rays', type=int, default=20_000)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    main(show=not args.no_show, ray_count=args.rays, output=args.output)
