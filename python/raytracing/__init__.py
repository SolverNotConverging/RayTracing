"""Electromagnetic ray tracing with a numerical C++ core and Python visualization."""
from __future__ import annotations

import os as _os
from pathlib import Path as _Path

_package_dir = _Path(__file__).resolve().parent
_dll_handles = []
if _os.name == "nt":
    for _directory in (_package_dir, _package_dir / ".libs"):
        if _directory.is_dir():
            _dll_handles.append(_os.add_dll_directory(str(_directory)))

from . import _core
from ._core import *
from .visualization import ViewOptions, plot, plot_csv, visualize, visualize_h5

__version__ = _core.__version__


def example_pattern() -> _Path:
    """Return the bundled 76–78 GHz CST patch pattern, independent of the working directory."""
    return _package_dir / "data" / "simple_patch.ffs"


__all__ = [name for name in dir(_core) if not name.startswith("_")] + [
    "example_pattern", "ViewOptions", "plot", "plot_csv", "visualize", "visualize_h5", "__version__"
]
