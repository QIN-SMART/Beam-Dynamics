"""
Unified relativistic reference (v0.13) — single source of γ/β/p0/velocity.

Every module derives the reference-particle kinematics through
BeamReference.from_energy_keV() instead of re-computing 1.0+E/511 etc.
Pure python (math), no numpy dependency.

Values are IDENTICAL to the previously duplicated derivations — no physics
change.
"""

import math
from dataclasses import dataclass

from shared.constants import MEC2_KEV, M_E_SI, C_SI


@dataclass(frozen=True)
class BeamReference:
    """Reference-particle relativistic kinematics at a kinetic energy."""
    energy_keV: float
    gamma: float
    beta: float
    p0: float          # [kg·m/s]
    velocity: float    # [m/s]

    @classmethod
    def from_energy_keV(cls, energy_keV):
        gamma = 1.0 + energy_keV / MEC2_KEV
        beta = math.sqrt(1.0 - 1.0 / gamma**2) if gamma > 1.0 else 0.0
        p0 = gamma * M_E_SI * beta * C_SI
        return cls(energy_keV=energy_keV, gamma=gamma, beta=beta,
                   p0=p0, velocity=beta * C_SI)

    @property
    def beta_gamma(self):
        return self.beta * self.gamma
