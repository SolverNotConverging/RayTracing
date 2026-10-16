"""Passive isotropic materials; numerical evaluation lives in C++."""
from dataclasses import dataclass
import math
from . import _core
from ._beam_adapter import native

C0 = _core.SPEED_OF_LIGHT
EPS0 = _core.VACUUM_PERMITTIVITY
MU0 = _core.VACUUM_PERMEABILITY
ETA0 = _core.FREE_SPACE_IMPEDANCE


@dataclass(frozen=True)
class Material:
    kind: str
    epsilon_r: complex = 1.0
    mu_r: float = 1.0
    conductivity: float = 0.0
    name: str = ""

    def __post_init__(self):
        value = self._native()
        object.__setattr__(self, "epsilon_r", value.epsilon_r)

    def _native(self):
        return native.Material(self.kind, complex(self.epsilon_r), self.mu_r, self.conductivity, self.name)

    @classmethod
    def _from_native(cls, value):
        return cls(value.kind, value.epsilon_r, value.mu_r, value.conductivity, value.name)

    @classmethod
    def dielectric(cls, epsilon_r, *, mu_r=1.0, conductivity=0.0, name="dielectric"):
        return cls("dielectric", epsilon_r, mu_r, conductivity, name)

    @classmethod
    def metal(cls, conductivity, *, mu_r=1.0, name="metal"):
        if not math.isfinite(conductivity) or conductivity <= 0:
            raise ValueError("Metal conductivity must be positive")
        return cls("metal", 1.0, mu_r, conductivity, name)

    @classmethod
    def pec(cls, name="PEC"):
        return cls("pec", name=name)

    def permittivity(self, frequency_hz):
        return self._native().permittivity(frequency_hz)

    def index(self, frequency_hz):
        return self._native().index(frequency_hz)


VACUUM = Material.dielectric(1.0, name="vacuum")


def fresnel(incident, outgoing, cos_incidence, frequency_hz):
    return tuple(native.fresnel(incident._native(), outgoing._native(), cos_incidence, frequency_hz))
