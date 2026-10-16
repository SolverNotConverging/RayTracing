"""Open cylinder and bottom disk with colocated analytical horizontal dipoles.

Run with .venv/Scripts/python.exe example2.py to open the desktop viewer.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from raytracing import point_ray as rt


def main(show=True, *, ray_count=200_000, output=None):
    config = rt.SolverConfig(
        frequency_hz=77e9,
        ray_count=ray_count,
        max_reflections=8,
        max_distance=20.0,
    )

    radius = 0.127
    depth = 0.5
    source_shift = 0.006  # Common lateral shift avoids the exact on-axis caustic.
    antenna_position = np.array([source_shift, 0, 0])

    scene = rt.Scene()
    scene.add(rt.Cylinder([0, 0, -depth / 2], [0, 0, 1], radius,
                          half_length=depth / 2, capped=False), rt.Material.pec())
    scene.add(rt.Disk([0, 0, -depth], [0, 0, -1], radius), rt.Material.pec())

    # A horizontal half-wave dipole has broadside radiation toward the bottom.
    antenna = rt.Antenna(rt.ThinWireDipole(rt.SPEED_OF_LIGHT / config.frequency_hz / 2))
    antenna.rotate([0, 1, 0], 90.0)
    tx = rt.Transmitter(antenna_position, antenna)
    rx = rt.Receiver(antenna_position, antenna)

    print("Running example 2...", flush=True)
    result = rt.solve(scene, tx, rx, config)
    print(
        f"{len(result.rays)} paths, {len(result.response_ray_indices)} paths with fields, "
        f"{len(result.impulse_response.taps)} impulse taps",
        flush=True,
    )

    output = Path(output) if output is not None else Path(__file__).resolve().parent / "results" / "example2"
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
    parser.add_argument('--rays', type=int, default=200_000)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    main(show=not args.no_show, ray_count=args.rays, output=args.output)
