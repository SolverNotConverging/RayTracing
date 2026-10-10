"""Coherent beam transport over explicit material regions; breaking 0.5 API."""
from __future__ import annotations

import os as _os
from pathlib import Path as _Path

_package_dir = _Path(__file__).resolve().parent
_dll_handles = []
if _os.name == "nt":
    for _directory in (_package_dir, _package_dir / ".libs"):
        if _directory.is_dir():
            _dll_handles.append(_os.add_dll_directory(str(_directory)))

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
    return _package_dir / "data" / "simple_patch.ffs"


__all__ = ["Material", "VACUUM", "fresnel", "C0", "EPS0", "MU0", "ETA0",
           "Scene", "Box", "Sphere", "Cylinder", "Rectangle", "Disk", "Triangle",
           "GaussianBeam", "PlaneWave", "ApertureSource", "SolverConfig", "Simulation", "solve",
           "Antenna", "Isotropic", "ShortDipole", "ThinWireDipole", "RectangularAperture", "load_farfield", "FarfieldImportOptions",
           "Polarization", "PhasorConvention", "PowerReference",
           "save_result", "load_result", "visualize", "plot_field", "inspect_scene", "ViewOptions", "AntennaPattern"]
