"""
Beam Parameters — RF Longitudinal Dynamics Benchmark
=====================================================
Common baseline beam + RF cavity and compression parameters.

RF kick model:
  δ_new = δ_old + (V_RF / E₀) · sin(φ + k_rf · z)
  z = -β·c·τ  (physical longitudinal position)

Drift compression:
  τ_final = τ_initial + R56 · δ    (R56 < 0 for drift)
  Bunch compresses when chirp δ(z) has correct sign.

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
e_charge = 1.602176634e-19    # C       Elementary charge
c_SI     = 2.99792458e8       # m/s     Speed of light

# ---- Relativistic (computed) ----
E_rest     = 511.0            # keV     Electron rest energy m_e·c²
gamma      = 1.0 + E_keV / E_rest       # Lorentz factor  γ ≈ 1.1957
beta       = (1.0 - 1.0 / gamma**2)**0.5  # Relativistic beta  β ≈ 0.5482
beta_gamma = beta * gamma               # βγ ≈ 0.6556
E0_eV      = E_keV * 1000     # eV      Beam energy in eV (100,000 eV)

# ---- Derived ----
epsilon_geom = epsilon_n / beta_gamma    # m·rad   Geometric emittance
sigma_xp     = epsilon_geom / spot_rms   # rad     RMS angular divergence x'
sigma_yp     = epsilon_geom / spot_rms   # rad     RMS angular divergence y'

# =====================================================================
#  Module-specific: RF Cavity
# =====================================================================

# ---- RF cavity parameters ----
f_RF       = 2.856e9          # Hz      RF frequency (S-band, standard UED)
lambda_rf  = c_SI / f_RF       # m       RF wavelength ≈ 105.0 mm
k_rf       = 2 * np.pi * f_RF / c_SI  # rad/m   RF wavenumber ≈ 59.9 rad/m
V_RF       = 30e3             # V       RF cavity voltage (30 kV default)
phi_RF     = np.pi            # rad     RF phase (π = zero-crossing, max chirp)

# ---- Chirp (analytic) ----
# dδ/dz = (V_RF / E₀) · k_rf · cos(φ)
# At φ = π:  dδ/dz ≈ -(V_RF · k_rf) / E₀ ≈ -18 m⁻¹

# ---- Drift compression geometry ----
L_drift_RF = 0.5              # m       Drift length after RF cavity

# ---- Charge (for beam generation only, SC OFF) ----
Q_bunch    = 100e-15          # C       Bunch charge (100 fC)

# ---- Bunch compression condition ----
# τ_final = τ_initial · (1 + R56 · dδ/dτ)
# Compression when:  dδ/dτ · R56 < 0
# Since R56 < 0 for drift, need dδ/dτ > 0  →  cos(φ) < 0  →  φ ∈ (π/2, 3π/2)

# CLI overrides: --V, --phi, --phi-scan
