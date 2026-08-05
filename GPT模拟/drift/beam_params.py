"""
Beam Parameters — Drift Benchmark
==================================
Common baseline beam shared across all Phase 1-2 modules.

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

# ---- Relativistic (computed) ----
E_rest     = 511.0            # keV     Electron rest energy m_e·c²
gamma      = 1.0 + E_keV / E_rest       # Lorentz factor  γ ≈ 1.1957
beta       = (1.0 - 1.0 / gamma**2)**0.5  # Relativistic beta  β ≈ 0.5482
beta_gamma = beta * gamma               # βγ ≈ 0.6556

# ---- Derived ----
epsilon_geom = epsilon_n / beta_gamma    # m·rad   Geometric emittance
sigma_xp     = epsilon_geom / spot_rms   # rad     RMS angular divergence x'
sigma_yp     = epsilon_geom / spot_rms   # rad     RMS angular divergence y'

# =====================================================================
#  Module-specific: Drift-only
# =====================================================================

drift_length = 1.0             # m       Total drift length

# CLI overrides: --epsn, --sigd
