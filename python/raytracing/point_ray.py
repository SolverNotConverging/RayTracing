"""Native point-ray solver used by the original multipath examples.

This API returns coherent receiver impulse taps. The top-level raytracing API
continues to provide Gaussian-beam transport.
"""
from ._core import *
from . import _core
from . import example_pattern
from .point_ray_visualization import ViewOptions, plot, plot_csv, visualize, visualize_h5, show

__version__ = _core.__version__
