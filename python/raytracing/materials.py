"""Passive isotropic materials; exp(-i omega t), SI units."""
from dataclasses import dataclass
import cmath
import math

C0 = 299792458.0
EPS0 = 8.8541878128e-12
MU0 = 1.25663706212e-6
ETA0 = math.sqrt(MU0 / EPS0)


@dataclass(frozen=True)
class Material:
    kind: str
    epsilon_r: complex = 1.0
    mu_r: float = 1.0
    conductivity: float = 0.0
    name: str = ""

    def __post_init__(self):
        e = complex(self.epsilon_r)
        if self.kind not in ("dielectric", "metal", "pec"):
            raise ValueError("Unknown material kind")
        if not all(math.isfinite(x) for x in (e.real, e.imag, self.mu_r, self.conductivity)):
            raise ValueError("Material parameters must be finite")
        if e.real <= 0 or e.imag < 0 or self.mu_r <= 0 or self.conductivity < 0:
            raise ValueError("Require Re(epsilon)>0, Im(epsilon)>=0, mu>0, sigma>=0")
        object.__setattr__(self, "epsilon_r", e)

    @classmethod
    def dielectric(cls, epsilon_r, *, mu_r=1.0, conductivity=0.0, name="dielectric"):
        return cls("dielectric", epsilon_r, mu_r, conductivity, name)

    @classmethod
    def metal(cls, conductivity, *, mu_r=1.0, name="metal"):
        if conductivity <= 0:
            raise ValueError("Metal conductivity must be positive")
        return cls("metal", 1.0, mu_r, conductivity, name)

    @classmethod
    def pec(cls, name="PEC"):
        return cls("pec", name=name)

    def permittivity(self, frequency_hz):
        if not math.isfinite(frequency_hz) or frequency_hz <= 0:
            raise ValueError("Frequency must be positive and finite")
        if self.kind == "pec":
            raise ValueError("PEC has no finite permittivity")
        return self.epsilon_r + 1j * self.conductivity / (2 * math.pi * frequency_hz * EPS0)

    def index(self, frequency_hz):
        return cmath.sqrt(self.permittivity(frequency_hz) * self.mu_r)


VACUUM = Material.dielectric(1.0, name="vacuum")


def fresnel(incident, outgoing, cos_incidence, frequency_hz):
    """Return (rs, rp, ts, tp). p basis is s cross propagation direction.

    Exact for lossless incidence, and for homogeneous normal incidence with loss.
    Oblique lossy incidence uses the homogeneous-wave weak-loss approximation.
    Opaque boundaries deliberately have zero transmission.
    """
    c = float(cos_incidence)
    if not 0 <= c <= 1 or incident.kind != "dielectric":
        raise ValueError("Require dielectric incident medium and cos incidence in [0,1]")
    n1 = incident.index(frequency_hz)
    if outgoing.kind == "pec":
        return -1+0j, 1+0j, 0j, 0j
    e1, e2 = incident.permittivity(frequency_hz), outgoing.permittivity(frequency_hz)
    n2 = outgoing.index(frequency_hz)
    q1 = n1 * c
    q2 = cmath.sqrt(n2*n2 - n1*n1 * (1-c*c) + 0j)
    if q2.imag < 0 or (q2.imag == 0 and q2.real < 0):
        q2 = -q2
    def ratio(a, b):
        return (a-b)/(a+b) if abs(a+b) > 0 else 0j
    rs = ratio(outgoing.mu_r*q1, incident.mu_r*q2)
    rp = ratio(e2*q1, e1*q2)
    if outgoing.kind == "metal":
        return rs, rp, 0j, 0j
    return rs, rp, 1+rs, n1/n2*(1+rp)
