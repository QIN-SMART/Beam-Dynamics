"""
Beam Parameters — Space Charge Drift Benchmark
===============================================
Common baseline beam + bunch charge and SC algorithm parameters.

SC algorithm (OCELOT 26.06.1):
  3D Particle-in-Cell, mesh [63, 63, 63], NGP deposition,
  spectral Poisson solver, linear field interpolation.

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
epsilon_n   = 0.08e-6         # m·rad   Normalized emittance (default)
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
#  Module-specific: Space Charge
# =====================================================================

# ---- Bunch charge ----
Q_bunch     = 100e-15          # C       Total bunch charge (100 fC default)

# ---- Drift geometry ----
drift_length = 1.0             # m       Total drift length

# ---- SC algorithm (OCELOT 26.06.1) ----
# SpaceCharge(step=1) — apply every tracking step
# nmesh_xyz = [63, 63, 63]     default 3D mesh
# method:    3D PIC, NGP deposition, spectral Poisson solver
# API:       sc.apply(p_array, dz)

# ---- Sensitivity scan ----
# Charge scan:   10, 50, 100, 500, 1000  fC  (--all / --Q)
# Emittance scan: 0.001, 0.005, 0.01, 0.08  mm·mrad  (--epsn / --epsn-scan)

# ---- Physics scaling ----
# SC defocusing:  dθ_sc/dz ∝ Q / (β²γ³ · σ_x·σ_y)
# Emittance angle: σ_x' = ε_n / (βγ · σ_x0)
# SC/ε crossover: Q / (βγ · ε_n) > threshold

# CLI overrides: --Q, --epsn, --all, --epsn-scan
