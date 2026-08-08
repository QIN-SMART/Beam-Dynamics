"""
Single source of physical constants (v0.12 architecture governance).

All project modules import from here instead of redefining values.
Values are IDENTICAL to the previously duplicated definitions — no physics
change.  Pure python (no numpy dependency).
"""

# ── SI constants ──────────────────────────────────────────────────────────
C_SI       = 2.99792458e8       # speed of light            [m/s]
M_E_SI     = 9.10938356e-31     # electron rest mass        [kg]
E_SI       = 1.602176634e-19    # elementary charge         [C]
EPSILON_0  = 8.854187817e-12    # vacuum permittivity       [F/m]
MEC2_KEV   = 511.0              # m_e c²                   [keV]
M_E_GEV    = 0.0005109988671734101  # m_e c²                [GeV]

# ── derived constants used across modules ─────────────────────────────────
K_BAND_HZ  = 2.856e9            # S-band RF frequency      [Hz] (default, per config)
