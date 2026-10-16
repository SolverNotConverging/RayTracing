"""Matplotlib plots."""
import numpy as np


def plot_impulse_response(response, ax=None):
    """Stem plot of tap magnitude and phase against delay; returns the figure."""
    import matplotlib.pyplot as plt

    if ax is None:
        _, (ax, phase_ax) = plt.subplots(2, 1, sharex=True, layout="constrained")
    else:
        phase_ax = ax.twinx()
    delays = response.delays * 1e9
    coefficients = response.coefficients
    ax.stem(delays, np.abs(coefficients), basefmt=" ")
    phase_ax.stem(delays, np.angle(coefficients, deg=True), linefmt="C3-", markerfmt="C3o", basefmt=" ")
    ax.set(title=f"Impulse response at {response.frequency_hz / 1e9:g} GHz", ylabel="|h|")
    phase_ax.set(xlabel="Delay (ns)", ylabel="Phase (deg)", ylim=(-180, 180))
    return ax.figure


def plot_field(simulation, points, ax=None):
    """|E| of a beam simulation along a polyline of sample points; returns the figure."""
    import matplotlib.pyplot as plt

    points = np.asarray(points, dtype=float)
    fields = simulation.fields_at(points)
    distance = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))])
    if ax is None:
        _, ax = plt.subplots(layout="constrained")
    ax.plot(distance, np.linalg.norm(fields, axis=1))
    ax.set(xlabel="Distance along samples (m)", ylabel="|E| (V/m)")
    return ax.figure
