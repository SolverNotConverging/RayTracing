"""Data conversion for the shared C++ Gaussian-beam implementation."""
from dataclasses import fields
import numpy as np
from ._core import _beams as native


def record(value, cls):
    result = cls()
    for f in fields(value):
        item = getattr(value, f.name)
        setattr(result, f.name, item._native() if hasattr(item, "_native") else item)
    return result


def values(value, cls):
    result = {}
    for f in fields(cls):
        item = getattr(value, f.name)
        result[f.name] = item.copy() if isinstance(item, np.ndarray) else item
    return result


def pair(value):
    return np.broadcast_to(np.asarray(value, dtype=float), (2,)).copy()
