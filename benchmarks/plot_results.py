"""Render the numerical benchmark CSVs with Matplotlib."""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, nargs="?", default=Path("benchmarks/results"))
    directory = parser.parse_args().directory
    data = np.genfromtxt(directory / "spreading.csv", delimiter=",", names=True,
                         dtype=None, encoding="utf-8")
    for name in dict.fromkeys(data["case"]):
        sample = data[data["case"] == name]
        figure, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True, layout="constrained")
        axes[0].plot(sample["distance_m"], sample["analytical_field_factor"], label="Analytical")
        axes[0].plot(sample["distance_m"], sample["numerical_field_factor"], ".", label="Jacobian")
        axes[0].set(title=name.replace("_", " "), ylabel="Field factor")
        axes[0].legend()
        axes[1].semilogy(sample["distance_m"], np.maximum(sample["max_relative_error"], 1e-16))
        axes[1].set(xlabel="Receiver distance (m)", ylabel="Maximum relative error")
        for axis in axes:
            axis.grid(True, alpha=0.3)
        figure.savefig(directory / f"{name}.png", dpi=150)
        plt.close(figure)

    data = np.genfromtxt(directory / "convergence.csv", delimiter=",", names=True)
    step, error = data["actual_angular_step"], data["relative_area_error"]
    figure, axis = plt.subplots(figsize=(9, 5), layout="constrained")
    axis.loglog(step, np.maximum(error, 1e-16), "o-", label="Jacobian area error")
    axis.loglog(step, error[0] * (step / step[0]) ** 2, "--", label="Second-order reference")
    axis.set(xlabel="Angular step (rad)", ylabel="Relative area error", title="Derivative convergence")
    axis.grid(True, which="both", alpha=0.3)
    axis.legend()
    figure.savefig(directory / "convergence.png", dpi=150)
    plt.close(figure)


if __name__ == "__main__":
    main()
