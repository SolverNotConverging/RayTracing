"""Electromagnetic ray tracing: point-ray multipath and Gaussian-beam solvers.

The numerical types come from the native ``_core`` extension. Gaussian-beam
transport lives in :mod:`raytracing.beam`; plotting helpers import Matplotlib
or PyVista on first use.
"""
from ._core import *  # noqa: F403
from ._core import __version__
from . import beam
from .plot import plot_impulse_response, plot_field
from .view import show_paths, show_beams, preview
