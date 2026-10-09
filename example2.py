"""Open cylinder and bottom disk with colocated downward CST antennas.

Run with .venv/Scripts/python.exe example2.py to open the desktop viewer.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import raytracing as rt


def main():
    config = rt.SolverConfig(
        frequency_hz=77e9,
        ray_count=20_000,
        max_reflections=8,
        max_distance=20.0,
    )

    radius = 0.127
    depth = 0.5
    source_shift = 0.006  # Common lateral shift avoids the exact on-axis caustic.
    antenna_position = np.array([source_shift, 0, 0])

    scene = rt.Scene()
    scene.add(rt.Cylinder([0, 0, -depth / 2], [0, 0, 1], radius,
                          half_length=depth / 2, capped=False))
    scene.add(rt.Disk([0, 0, -depth], [0, 0, -1], radius))

    patch = rt.load_farfield(rt.example_pattern())
    # The supplied CST patch already points downward (-z). Both endpoints copy it.
    tx = rt.Transmitter(antenna_position, patch)
    rx = rt.Receiver(antenna_position, patch)

    print("Running example 2...", flush=True)
    result = rt.solve(scene, tx, rx, config)
    print(
        f"{len(result.rays)} paths, {len(result.response_ray_indices)} paths with fields, "
        f"{len(result.impulse_response.taps)} impulse taps",
        flush=True,
    )

    output = Path(__file__).resolve().parent / "results" / "example2"
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
    print("Drag to rotate; Shift+drag to pan; scroll to zoom. Close the viewer to exit.", flush=True)
    rt.show(result)


if __name__ == "__main__":
    main()
