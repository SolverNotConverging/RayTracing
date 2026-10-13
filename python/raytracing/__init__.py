"""Coherent beam transport over explicit material regions; breaking 0.5 API."""
from __future__ import annotations

import os as _os
from pathlib import Path as _Path
from importlib import import_module as _import_module
from importlib.machinery import EXTENSION_SUFFIXES as _EXTENSION_SUFFIXES
from importlib.metadata import distribution as _distribution, PackageNotFoundError as _PackageNotFoundError

_package_dir = _Path(__file__).resolve().parent
# IDE source roots can put this directory ahead of the installed wheel. Keep
# Python code editable while using the native extension from this interpreter's
# matching installation, rather than relying on IDE-specific PYTHONPATH settings.
_native_dir = _package_dir
if not any((_package_dir / ("_core" + suffix)).is_file() for suffix in _EXTENSION_SUFFIXES):
    try:
        _installed = _distribution("raytracing")
    except _PackageNotFoundError:
        _installed = None
    if _installed is not None and _installed.version == "0.5.0":
        _candidate = _Path(_installed.locate_file("raytracing")).resolve()
        if any((_candidate / ("_core" + suffix)).is_file() for suffix in _EXTENSION_SUFFIXES):
            _native_dir = _candidate
            __path__.append(str(_native_dir))
_dll_handles = []
if _os.name == "nt":
    for _directory in (_native_dir, _native_dir / ".libs"):
        if _directory.is_dir():
            _dll_handles.append(_os.add_dll_directory(str(_directory)))

# Import explicitly before dependent modules, so a missing extension is reported
# accurately instead of as a misleading partially-initialized/circular import.
try:
    _core = _import_module("._core", __name__)
except ModuleNotFoundError as _error:
    if _error.name != __name__ + "._core":
        raise
    raise ImportError(
        "raytracing's native _core extension is missing for this Python interpreter. "
        "Build and install raytracing into this environment (on this development "
        "machine, run .\\install.ps1), then restart Python."
    ) from _error

from .materials import Material, VACUUM, fresnel, C0, EPS0, MU0, ETA0
from .scene import Scene, Box, Sphere, Cylinder, Rectangle, Disk, Triangle
from .beams import GaussianBeam, PlaneWave, ApertureSource
from ._core import Antenna, Isotropic, ShortDipole, ThinWireDipole, RectangularAperture, load_farfield, FarfieldImportOptions
from ._core import Polarization, PhasorConvention, PowerReference
from .transport import SolverConfig, Simulation, solve
from .persistence import save_result, load_result
from .display import plot_field
from .visualization import inspect_scene, visualize, ViewOptions, AntennaPattern

__version__ = "0.5.0"


def example_pattern() -> _Path:
    """Return the bundled 76–78 GHz CST patch pattern, independent of the working directory."""
    return _native_dir / "data" / "simple_patch.ffs"


__all__ = ["Material", "VACUUM", "fresnel", "C0", "EPS0", "MU0", "ETA0",
           "Scene", "Box", "Sphere", "Cylinder", "Rectangle", "Disk", "Triangle",
           "GaussianBeam", "PlaneWave", "ApertureSource", "SolverConfig", "Simulation", "solve",
           "Antenna", "Isotropic", "ShortDipole", "ThinWireDipole", "RectangularAperture", "load_farfield", "FarfieldImportOptions",
           "Polarization", "PhasorConvention", "PowerReference",
           "save_result", "load_result", "visualize", "plot_field", "inspect_scene", "ViewOptions", "AntennaPattern"]
