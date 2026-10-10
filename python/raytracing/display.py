"""Matplotlib field-sampling plots."""
import numpy as np


def plot_field(result, positions, *, path=None, show=True):
    import matplotlib.pyplot as plt
    points = np.asarray(positions)
    values = result.fields_at(points)
    distance = np.linalg.norm(points-points[0],axis=1)
    fig,ax = plt.subplots()
    ax.plot(distance,np.linalg.norm(values,axis=1))
    ax.set(xlabel="Distance from first sample (m)",ylabel="|E| (V/m)")
    if path: fig.savefig(path,dpi=160,bbox_inches="tight")
    if show: plt.show()
    return fig
