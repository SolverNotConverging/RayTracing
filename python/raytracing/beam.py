"""Paraxial Gaussian-beam transport through material regions."""
from ._core import beam as _native

__all__ = [
    "Scene", "Object", "Event", "Termination", "EventKind",
    "GaussianBeam", "PlaneWave", "ApertureSource",
    "SolverConfig", "Simulation", "Segment", "State", "Interaction", "InteractionKind",
    "Terminal", "TerminalReason", "frame", "propagate", "solve", "solve_states",
]
globals().update({name: getattr(_native, name) for name in __all__})
