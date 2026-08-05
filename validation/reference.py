"""
Analytic reference models used ONLY as ground truth in validation tests.

These mirror the already-validated formulas in the existing backends'
benchmark scripts (GPT模拟/solenoid/benchmark_solenoid.py* and AG):
  - drift envelope:      σ(z) = √(σ0² + (ε_geo/σ0)² z²)      (uncorrelated beam)
  - solenoid 4×4:        hard-edge Brown-Chao transfer matrix (identical to
                         OCELOT SolenoidAtom.R_main_matrix and the benchmark
                         4×4); envelope via Σ(z) = M Σ0 Mᵀ
  - RF chirp rate:       h = −e·E_rf·k_rf / (β²γ³ m_e c²)   (= AG rf_chirp_rate)

No physics is re-implemented here beyond these standard validation references.
"""

import numpy as np

C_SI = 2.99792458e8
M_E_SI = 9.10938356e-31
E_SI = 1.602176634e-19


def gamma_beta_p(energy_keV):
    gamma = 1.0 + energy_keV / 511.0
    beta = np.sqrt(1.0 - 1.0 / gamma**2)
    return gamma, beta, gamma * M_E_SI * beta * C_SI


def drift_sigma(z_arr, sigma0, sxp0):
    """Uncorrelated-beam drift envelope: σ(z)=√(σ0²+σ0'²z²)."""
    return np.sqrt(sigma0**2 + (sxp0 * z_arr)**2)


def mat_drift_4x4(s):
    return np.array([[1.0, s, 0.0, 0.0],
                     [0.0, 1.0, 0.0, 0.0],
                     [0.0, 0.0, 1.0, s],
                     [0.0, 0.0, 0.0, 1.0]])


def mat_solenoid_4x4(s, k):
    """Hard-edge solenoid matrix (Brown-Chao) — identical to OCELOT."""
    C = np.cos(k * s)
    S = np.sin(k * s)
    return np.array([[C * C,       S * C / k,   S * C,       S * S / k],
                     [-k * S * C,  C * C,      -k * S * S,   S * C],
                     [-S * C,     -S * S / k,   C * C,       S * C / k],
                     [k * S * S,  -S * C,      -k * S * C,   C * C]])


def solenoid_envelope_4x4(z_arr, L1, Lsol, k, sigma0, sxp0):
    """Projected σ_x(z) through drift-solenoid-drift using the exact 4×4 map.

    Returns (sigma_x, sigma_y, sigma_xy, eps_proj) sampled at z_arr.
    Beam assumed initially uncorrelated round: Σ0 = diag(σ0²,σ0'²,σ0²,σ0'²).
    """
    eps_geo = sigma0 * sxp0                      # ε = σ0·σ0' for round beam
    S0 = np.diag([sigma0**2, sxp0**2, sigma0**2, sxp0**2])
    z_end = max(z_arr)
    sx = np.zeros_like(z_arr)
    sy = np.zeros_like(z_arr)
    sxy = np.zeros_like(z_arr)
    epsp = np.zeros_like(z_arr)
    for i, z in enumerate(z_arr):
        if z <= L1:
            M = mat_drift_4x4(z)
        elif z <= L1 + Lsol:
            M = mat_solenoid_4x4(z - L1, k) @ mat_drift_4x4(L1)
        else:
            M = mat_drift_4x4(z - L1 - Lsol) @ mat_solenoid_4x4(Lsol, k) @ mat_drift_4x4(L1)
        Sz = M @ S0 @ M.T
        sx[i] = np.sqrt(max(Sz[0, 0], 0.0))
        sy[i] = np.sqrt(max(Sz[2, 2], 0.0))
        sxy[i] = Sz[0, 2]
        epsp[i] = np.sqrt(max(Sz[0, 0] * Sz[1, 1] - Sz[0, 1]**2, 0.0))
    return sx, sy, sxy, epsp


def rf_chirp_rate(energy_keV, E_rf, k_rf):
    """h = −e·E_rf·k_rf / (β²γ³·m_e·c²)   [m⁻¹]   (= AG rf_chirp_rate)."""
    gamma, beta, _ = gamma_beta_p(energy_keV)
    return -E_SI * E_rf * k_rf / (max(beta**2, 1e-10) * gamma**3 * M_E_SI * C_SI**2)
