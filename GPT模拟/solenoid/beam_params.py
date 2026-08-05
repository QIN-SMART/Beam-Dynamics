"""
Beam Parameters — Solenoid Benchmark
=====================================
Common baseline beam + solenoid geometry and field parameters.

Physics:  k_s = e·B / (2·p)   where  p = γ·m_e·β·c

Usage:
  from beam_params import *
"""

# =====================================================================
#  Common Baseline Beam  (Phase 1 initial distribution)
# =====================================================================

spot_rms    = 85e-6           # m       RMS transverse spot size
sig_z0      = 300e-6          # m       RMS bunch length
E_keV       = 100.0           # keV     Beam kinetic energy
v_e         = 1.6435e8        # m/s     Electron velocity @ 100 keV
epsilon_n   = 0.08e-6         # m·rad   Normalized emittance
sigma_delta = 1.0e-4          # —       Relative energy spread ΔE/E

# ---- Physical constants ----
e_SI   = 1.602176634e-19      # C       Elementary charge
m_e_SI = 9.10938356e-31       # kg      Electron rest mass
c_SI   = 2.99792458e8         # m/s     Speed of light
mec2   = 511.0                # keV     m_e·c²

# ---- Relativistic (computed) ----
E_rest     = 511.0            # keV     Electron rest energy
gamma      = 1.0 + E_keV / mec2        # Lorentz factor  γ ≈ 1.1957
beta       = (1.0 - 1.0 / gamma**2)**0.5  # Relativistic beta  β ≈ 0.5482
beta_gamma = beta * gamma               # βγ ≈ 0.6556

# ---- Derived ----
epsilon_geom = epsilon_n / beta_gamma    # m·rad   Geometric emittance
sigma_xp     = epsilon_geom / spot_rms   # rad     RMS angular divergence x'
sigma_yp     = epsilon_geom / spot_rms   # rad     RMS angular divergence y'

# =====================================================================
#  Module-specific: Solenoid
# =====================================================================

# ---- Abstract solenoid (Phase 1.8) ----
k_sol    = 1.5               # —       Solenoid strength (TL1, estimated)

# ---- Physical solenoid (Phase 1.9) ----
B_sol    = 0.05               # T       On-axis magnetic field
p_SI     = gamma * m_e_SI * beta * c_SI  # kg·m/s  Electron momentum
# k_s = e·B / (2·p)  — computed at runtime

# ---- Geometry ----
L_sol       = 0.06            # m       Solenoid length
L_drift1    = 0.10            # m       Drift before solenoid
L_drift2    = 0.60            # m       Drift after solenoid

# CLI overrides: --k, --B, --epsn, --Ek
