#!/usr/bin/env python3
"""
6D Covariance Matrix Beam Dynamics
====================================
Phase 3: Upgrade from RMS envelope model to full 6D covariance tracking.

Phase space:  X = (x, x', y, y', z, δ)^T   where δ = Δp/p
Covariance:   Σ = ⟨X X^T⟩  (6×6 symmetric, 21 independent elements)

Evolution:    dΣ/dz = A(z) Σ + Σ A(z)^T + Q_sc(z)

where A(z) is the infinitesimal transport matrix and Q_sc is the
space-charge contribution.

Based on:  Reiser, "Theory and Design of Charged Particle Beams" (2008)
           OPTIMIZATION_PLAN.md Phase 3
"""

import numpy as np
from scipy.integrate import solve_ivp
from dataclasses import dataclass, field
from typing import Tuple, Optional, Callable, List
import warnings

# ── Constants ──
M_E       = 9.10938356e-31
M_E_EV    = 5.109989461e5
C_LIGHT   = 2.99792458e8
E_CHARGE  = 1.60217662e-19
EPSILON_0 = 8.854187817e-12
R_E       = 2.8179403262e-15
ETA       = E_CHARGE / (M_E * C_LIGHT**2)


def gamma_to_beta(g: float) -> float:
    return np.sqrt(1.0 - 1.0/g**2) if g > 1.0 else 0.0

def beamK_to_gamma(K_eV: float) -> float:
    return 1.0 + K_eV / M_E_EV

def gamma_to_beamK(g: float) -> float:
    return (g - 1.0) * M_E_EV


# ═══════════════════════════════════════════════════════════════════
# 6D Covariance Matrix
# ═══════════════════════════════════════════════════════════════════

@dataclass
class Covariance6D:
    """
    6D beam covariance matrix in the basis (x, x', y, y', z, δ).

    Σ = ⟨X X^T⟩  where X = (x, x', y, y', z, δ)^T, δ = Δp/p.

    The full 6×6 matrix has 21 independent elements (symmetric).
    Diagonal:  ⟨x²⟩=σ_x², ⟨x'²⟩=θ_x², ⟨y²⟩=σ_y², ⟨y'²⟩=θ_y²,
               ⟨z²⟩=σ_z², ⟨δ²⟩=σ_δ².
    Key off-diagonals: ⟨x x'⟩=σ_x·ν_x, ⟨z δ⟩=C_zδ, ⟨x y⟩=σ_xy, etc.
    """
    sigma: np.ndarray = field(default_factory=lambda: np.zeros((6,6)))

    # parameters
    gamma: float = 1.0
    Ne: float = 0.0

    def __post_init__(self):
        if self.sigma.shape != (6, 6):
            raise ValueError("sigma must be 6×6")

    @property
    def beta(self) -> float:
        return gamma_to_beta(self.gamma)

    @property
    def p0(self) -> float:
        return self.gamma * self.beta

    @property
    def beamK(self) -> float:
        return gamma_to_beamK(self.gamma)

    # ── convenient accessors ──
    @property
    def sigma_x(self) -> float:
        return np.sqrt(max(self.sigma[0,0], 0.0))

    @property
    def sigma_y(self) -> float:
        return np.sqrt(max(self.sigma[2,2], 0.0))

    @property
    def sigma_z(self) -> float:
        return np.sqrt(max(self.sigma[4,4], 0.0))

    @property
    def sigma_delta(self) -> float:
        return np.sqrt(max(self.sigma[5,5], 0.0))

    @property
    def nu_x(self) -> float:
        return self.sigma[0,1] / max(self.sigma_x, 1e-15)

    @property
    def nu_y(self) -> float:
        return self.sigma[2,3] / max(self.sigma_y, 1e-15)

    @property
    def nu_z(self) -> float:
        return self.sigma[4,5] / max(self.sigma_z, 1e-15)

    @property
    def C_zd(self) -> float:
        return self.sigma[4,5]

    @property
    def theta_x(self) -> float:
        return np.sqrt(max(self.sigma[1,1], 0.0))

    @property
    def theta_y(self) -> float:
        return np.sqrt(max(self.sigma[3,3], 0.0))

    @property
    def time_spread_ps(self) -> float:
        v = self.beta * C_LIGHT
        return self.sigma_z / v * 1e12 if v > 1e3 else 0.0

    @property
    def eps_nx(self) -> float:
        """Normalized transverse emittance (from determinant of 2×2 block)."""
        det = self.sigma[0,0]*self.sigma[1,1] - self.sigma[0,1]**2
        return self.p0 * np.sqrt(max(det, 0.0))

    @property
    def eps_ny(self) -> float:
        det = self.sigma[2,2]*self.sigma[3,3] - self.sigma[2,3]**2
        return self.p0 * np.sqrt(max(det, 0.0))

    @property
    def eps_nz(self) -> float:
        det = self.sigma[4,4]*self.sigma[5,5] - self.sigma[4,5]**2
        return self.p0 * np.sqrt(max(det, 0.0))

    def state_vector(self) -> np.ndarray:
        """Flatten upper-triangle of Σ into 21-element vector."""
        idx = np.triu_indices(6)
        return self.sigma[idx]

    @classmethod
    def from_state(cls, vec: np.ndarray, gamma: float, Ne: float) -> "Covariance6D":
        """Rebuild Σ from 21-element upper-triangle vector."""
        s = np.zeros((6,6))
        idx = np.triu_indices(6)
        s[idx] = vec
        # enforce symmetry
        for i in range(6):
            for j in range(i+1, 6):
                s[j,i] = s[i,j]
        return cls(sigma=s, gamma=gamma, Ne=Ne)

    def apply_transport(self, R: np.ndarray) -> "Covariance6D":
        """Apply discrete transport: Σ → R·Σ·R^T."""
        c = self.copy()
        c.sigma = R @ self.sigma @ R.T
        return c

    def copy(self) -> "Covariance6D":
        return Covariance6D(sigma=self.sigma.copy(),
                            gamma=self.gamma, Ne=self.Ne)

    def summary(self, label: str = "") -> str:
        lines = [f"══ {label} ══" if label else "──"]
        lines.append(f"  γ={self.gamma:.4f} β={self.beta:.4f} K={self.beamK:.0f}eV")
        lines.append(f"  Ne={self.Ne:.2e}")
        lines.append(f"  σ_x={self.sigma_x*1e6:7.1f}μm  σ_y={self.sigma_y*1e6:7.1f}μm  "
                     f"σ_z={self.sigma_z*1e6:7.1f}μm")
        lines.append(f"  ν_x={self.nu_x:.6f}  ν_y={self.nu_y:.6f}  ν_z={self.nu_z:.6f}")
        lines.append(f"  σ_δ={self.sigma_delta*1e3:.3f}‰  C_zδ={self.C_zd*1e6:.3f}μm")
        lines.append(f"  ε_nx={self.eps_nx*1e6:.4f}  ε_ny={self.eps_ny*1e6:.4f}  "
                     f"ε_nz={self.eps_nz*1e6:.4f} μm")
        lines.append(f"  Δt={self.time_spread_ps*1000:.0f}fs")
        return "\n".join(lines)

    def __repr__(self):
        return (f"Cov6D(σ=({self.sigma_x*1e6:.0f},{self.sigma_y*1e6:.0f},"
                f"{self.sigma_z*1e6:.0f})μm, γ={self.gamma:.3f})")


# ═══════════════════════════════════════════════════════════════════
# Transport matrices  (6×6)
# ═══════════════════════════════════════════════════════════════════

def drift_matrix(L: float, gamma: float, beta: float) -> np.ndarray:
    """
    6×6 drift transport matrix.

    x → x + L·x',  x' → x'
    y → y + L·y',  y' → y'
    z → z + R_56·δ,  δ → δ

    R_56 = L/(γ²β²)  (velocity bunching)
    """
    g2b2 = max(gamma**2 * beta**2, 1e-10)
    R56 = L / g2b2

    R = np.eye(6)
    R[0,1] = L
    R[2,3] = L
    R[4,5] = R56
    return R


def thin_lens_quadrupole(f_x: float, f_y: float) -> np.ndarray:
    """
    6×6 thin-lens quadrupole transport.

    x' → x' + x/f_x,  y' → y' + y/f_y
    (convention: 1/f > 0 → focusing)

    For a standard quadrupole: f_y = -f_x.
    """
    R = np.eye(6)
    R[1,0] = 1.0/f_x
    R[3,2] = 1.0/f_y
    return R


def thin_lens_solenoid(Bz: float, L: float, gamma: float, beta: float) -> np.ndarray:
    """
    6×6 thin-lens solenoid transport (hard-edge model, thin limit).

    Focusing + Larmor rotation.  In the thin-lens limit:
      x' → x' - k_s²L·x  (focusing)
      y' → y' - k_s²L·y  (focusing, symmetric)
    plus Larmor rotation couples x and y.

    k_s = e·Bz/(2p) = η·c·Bz/(2p0)  [1/m]
    """
    p0 = max(gamma * beta, 1e-10)
    k_s = ETA * C_LIGHT * Bz / (2.0 * p0)
    k_s2L = k_s**2 * L

    R = np.eye(6)
    # focusing
    R[1,0] = -k_s2L
    R[3,2] = -k_s2L
    # Larmor rotation (small-angle: θ = k_s·L)
    theta = k_s * L
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)

    # Rotation in (x, y) and (x', y') planes
    # x_out  =  cosθ·x + sinθ·y
    # y_out  = -sinθ·x + cosθ·y
    # (similarly for velocities)
    R = np.zeros((6,6))
    R[0,0] = cos_t;  R[0,2] = sin_t
    R[1,0] = -k_s2L*cos_t;  R[1,2] = -k_s2L*sin_t
    R[1,1] = cos_t;  R[1,3] = sin_t

    R[2,0] = -sin_t;  R[2,2] = cos_t
    R[3,0] = k_s2L*sin_t;   R[3,2] = -k_s2L*cos_t
    R[3,1] = -sin_t;  R[3,3] = cos_t

    R[4,4] = 1.0
    R[5,5] = 1.0
    return R


def thin_lens_rf_cavity(k_rf: float, L_cav: float, gamma: float, beta: float,
                         phase: float = 0.0) -> np.ndarray:
    """
    6×6 thin-lens RF cavity transport.

    Transverse RF defocusing:
      x' → x' + η·E_rf·k·L·x/(2γβ²)  (thin-lens)
    Longitudinal chirp:
      δ → δ + h·z   where h = -η·E_rf·k·L/(β²γ³)

    Parameters
    ----------
    k_rf : float [1/m]
        RF wave number 2πf/c.
    L_cav : float [m]
        Cavity length.
    gamma, beta : float
    """
    gb2 = max(gamma * beta**2, 1e-10)
    g3b2 = max(gamma**3 * beta**2, 1e-10)

    # Assume E_rf * k_rf is given indirectly through k_rf and L_cav
    # For now, use a parameterized focusing strength
    # Transverse: F_x = -η·E·k/(2γβ²) · σ_x
    # This gives a thin-lens kick: Δx' = +F_x·L = η·E·k·L/(2γβ²)·x
    k_focus = -ETA * C_LIGHT * k_rf * L_cav / (2.0 * gb2)

    # Longitudinal chirp
    h_chirp = -ETA * C_LIGHT * k_rf * L_cav / g3b2

    R = np.eye(6)
    R[1,0] = k_focus
    R[3,2] = k_focus
    R[5,4] = h_chirp
    return R


# ═══════════════════════════════════════════════════════════════════
# Space-charge contribution to dΣ/dz
# ═══════════════════════════════════════════════════════════════════

def space_charge_cov_drift(cov: Covariance6D) -> np.ndarray:
    """
    Space-charge contribution Q_sc to dΣ/dz for a Gaussian beam.

    In the covariance formalism, space charge affects the momentum
    spread terms (x', y', δ).  Using the Sacherer formalism:

    d⟨x'²⟩/dz = (2K_sc)/(σ_x(σ_x+σ_y)) · ⟨x'²⟩  approximately

    where K_sc = N_e·r_e/(√(2π)·γ³·β²·σ_z) for a Gaussian beam.

    Returns a 6×6 symmetric matrix Q_sc.
    """
    Ne = cov.Ne
    gamma = cov.gamma
    beta = cov.beta
    sx = max(cov.sigma_x, 1e-12)
    sy = max(cov.sigma_y, 1e-12)
    sz = max(cov.sigma_z, 1e-12)

    g3b2 = max(gamma**3 * beta**2, 1e-30)
    K_sc = Ne * R_E / g3b2

    # Transverse defocusing: increases ⟨x'²⟩, ⟨y'²⟩
    # The correlation ⟨x x'⟩ also gets a contribution
    dsig = np.zeros((6,6))

    # d⟨x x'⟩/dz ≈ ⟨x'²⟩ + (SC contribution)
    # The SC term: ⟨x·x''_sc⟩ = (K_sc/(sx+sy))·⟨x²⟩  (paraxial approx)
    sc_trans = K_sc / (max(sx*(sx+sy), 1e-30))

    # Contributions to the d⟨x x'⟩/dz:
    dsig[0,1] = sc_trans * cov.sigma[0,0]   # x-x' SC coupling
    dsig[1,0] = dsig[0,1]
    dsig[2,3] = sc_trans * cov.sigma[2,2]   # y-y' SC coupling
    dsig[3,2] = dsig[2,3]

    # Longitudinal SC: d⟨zδ⟩/dz ⊃ SC_z
    sc_long = K_sc / (max(sx*sy, 1e-30))
    dsig[4,5] = sc_long * cov.sigma[4,4] * (1.0/gamma**2)  # suppressed
    dsig[5,4] = dsig[4,5]

    return dsig


# ═══════════════════════════════════════════════════════════════════
# ODE:  dΣ/dz = A(z)·Σ + Σ·A(z)^T + Q_sc(z)
# ═══════════════════════════════════════════════════════════════════

def covariance_ode(z: float, y: np.ndarray, Ne: float,
                   gamma_prime: float = 0.0,
                   external_A_func: Optional[Callable] = None,
                   include_sc: bool = True
                   ) -> np.ndarray:
    """
    RHS of the 6D covariance ODE.

    State: 21-element upper-triangle of Σ, plus γ as the 22nd element.
    Returns d(state)/dz: 21 elements of dΣ/dz + dγ/dz.
    """
    gamma = y[21]
    cov = Covariance6D.from_state(y[:21], gamma, Ne)

    # Adiabatic damping: same as envelope model
    beta = cov.beta
    if gamma_prime != 0 and gamma > 1.01 and beta > 1e-6:
        damping = gamma * gamma_prime / (gamma**2 * beta**2)
    else:
        damping = 0.0

    # Transport matrix A(z): always includes baseline drift kinematics
    gamma = cov.gamma
    beta = cov.beta
    g2b2 = max(gamma**2 * beta**2, 1e-10)

    A_drift = np.zeros((6,6))
    A_drift[0,1] = 1.0
    A_drift[2,3] = 1.0
    A_drift[4,5] = 1.0 / g2b2   # R_56 contribution: dz/d(z_indep) = δ/γ²

    if external_A_func is not None:
        A = A_drift + external_A_func(z, cov)
    else:
        A = A_drift

    # dΣ/dz = A Σ + Σ A^T
    dsigma = A @ cov.sigma + cov.sigma @ A.T

    # Adiabatic damping on off-diagonal terms
    if damping != 0:
        for i in range(6):
            for j in range(6):
                if i != j:
                    dsigma[i,j] -= damping * cov.sigma[i,j]

    # Space charge
    if include_sc and Ne > 0:
        dsigma += space_charge_cov_drift(cov)

    # Flatten back
    idx = np.triu_indices(6)
    dy = np.zeros(22)
    dy[:21] = dsigma[idx]
    dy[21] = gamma_prime

    return dy


# ═══════════════════════════════════════════════════════════════════
# Integration
# ═══════════════════════════════════════════════════════════════════

def propagate_cov(cov0: Covariance6D, z_span: Tuple[float,float],
                  n_points: int = 500,
                  external_A_func: Optional[Callable] = None,
                  gamma_prime_func: Optional[Callable] = None,
                  include_sc: bool = True,
                  **solver_kwargs
                  ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Integrate the 6D covariance ODE.

    Returns (z_arr, states_arr) where states_arr has shape (n_points, 22).
    The first 21 columns are the upper-triangle of Σ, column 21 is γ.
    """
    y0 = np.zeros(22)
    idx = np.triu_indices(6)
    y0[:21] = cov0.sigma[idx]
    y0[21] = cov0.gamma

    z_start, z_end = z_span
    t_eval = np.linspace(z_start, z_end, n_points)

    def gamma_prime(z, cov=None):
        if gamma_prime_func is not None:
            # Need beam-like object; use cov0 as proxy
            return gamma_prime_func(z)
        return 0.0

    def rhs(z, y):
        gp = gamma_prime(z)
        return covariance_ode(z, y, cov0.Ne,
                              gamma_prime=gp,
                              external_A_func=external_A_func,
                              include_sc=include_sc)

    sol = solve_ivp(rhs, (z_start, z_end), y0, t_eval=t_eval,
                    method=solver_kwargs.pop('method', 'RK45'),
                    rtol=solver_kwargs.pop('rtol', 1e-8),
                    atol=solver_kwargs.pop('atol', 1e-12),
                    **solver_kwargs)

    return sol.t, sol.y.T


# ═══════════════════════════════════════════════════════════════════
# Builder
# ═══════════════════════════════════════════════════════════════════

def make_covariance_beam(
    sigma_x0_um: float = 85.0,
    sigma_y0_um: float = 85.0,
    sigma_z0_um: float = 300.0,
    sigma_delta: float = 1e-3,
    beamK_eV: float = 100_000.0,
    Ne: float = 1e5,
    eps_nx_um: float = 0.03,
    eps_ny_um: float = 0.03,
    eps_nz_um: float = 0.2,
) -> Covariance6D:
    """
    Build an initial 6D covariance matrix for a UED beam.

    Uses diagonal Gaussian initial distribution (no correlations).
    """
    gamma = beamK_to_gamma(beamK_eV)
    beta = gamma_to_beta(gamma)
    p0 = gamma * beta

    sigma_x = sigma_x0_um * 1e-6
    sigma_y = sigma_y0_um * 1e-6
    sigma_z = sigma_z0_um * 1e-6

    eps_gx = (eps_nx_um * 1e-6) / p0
    eps_gy = (eps_ny_um * 1e-6) / p0
    eps_gz = (eps_nz_um * 1e-6) / p0

    s = np.zeros((6,6))
    s[0,0] = sigma_x**2
    s[1,1] = (eps_gx / sigma_x)**2      # ⟨x'²⟩ from emittance
    s[2,2] = sigma_y**2
    s[3,3] = (eps_gy / sigma_y)**2
    s[4,4] = sigma_z**2
    s[5,5] = sigma_delta**2

    return Covariance6D(sigma=s, gamma=gamma, Ne=Ne)


# ═══════════════════════════════════════════════════════════════════
# External A-matrix builder for beamline
# ═══════════════════════════════════════════════════════════════════

@dataclass
class ARegion:
    """A region where a transport matrix A is active."""
    z_start: float
    z_end: float
    A_matrix: np.ndarray   # 6×6 infinitesimal transport matrix


def make_A_func(regions: List[ARegion]) -> Callable:
    """
    Build A(z) function that returns the active transport matrix.
    """
    def A_func(z: float, cov: Covariance6D) -> np.ndarray:
        A = np.zeros((6,6))
        for reg in regions:
            if reg.z_start <= z <= reg.z_end:
                A += reg.A_matrix
        return A
    return A_func


# ═══════════════════════════════════════════════════════════════════
# Benchmark: compare with envelope model for pure drift
# ═══════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    print("=" * 60)
    print("  Phase 3: 6D Covariance Matrix Beam Dynamics")
    print("=" * 60)

    # ── 1. Build beam ──
    cov0 = make_covariance_beam(Ne=1e5, beamK_eV=100_000)
    print("\nInitial beam:")
    print(cov0.summary("z=0"))

    # ── 2. Pure drift benchmark (no external fields) ──
    print("\n── Drift 0 → 0.5 m (pure drift, no SC) ──")
    z_drift, st_drift = propagate_cov(cov0, (0.0, 0.5), n_points=100,
                                       include_sc=False)

    cov_end_noSC = Covariance6D.from_state(st_drift[-1,:21],
                                            st_drift[-1,21], cov0.Ne)
    print(f"  σ_x: {cov0.sigma_x*1e6:.1f} → {cov_end_noSC.sigma_x*1e6:.1f} μm")

    # Analytic: σ(z) = σ₀√(1 + (ε_eff·z/σ₀²)²)
    eps_geo = cov0.eps_nx / cov0.p0
    sx_analytic = cov0.sigma_x * np.sqrt(1 + (eps_geo*0.5/cov0.sigma_x**2)**2)
    print(f"  Analytic: σ_x = {sx_analytic*1e6:.1f} μm")
    err = abs(cov_end_noSC.sigma_x - sx_analytic) / sx_analytic * 100
    print(f"  Error: {err:.2f}%")

    # ── 3. Drift with space charge ──
    print("\n── Drift 0 → 0.1 m (with SC, Ne=1e5) ──")
    z_sc, st_sc = propagate_cov(cov0, (0.0, 0.1), n_points=100,
                                 include_sc=True)
    cov_end_SC = Covariance6D.from_state(st_sc[-1,:21], st_sc[-1,21], cov0.Ne)
    growth = cov_end_SC.sigma_x / cov0.sigma_x
    print(f"  σ_x: {cov0.sigma_x*1e6:.1f} → {cov_end_SC.sigma_x*1e6:.1f} μm ({growth:.3f}x)")
    print(f"  ε_nx conserved: {cov_end_SC.eps_nx*1e6:.4f} μm (initial: {cov0.eps_nx*1e6:.4f})")

    # ── 4. Solenoid test ──
    print("\n── Solenoid: Bz=0.1T, L=0.1m at z=0.1m ──")
    # Build A matrix for solenoid: continuous focusing
    k_s = ETA * C_LIGHT * 0.1 / (2.0 * max(cov0.p0, 1e-10))
    A_sol = np.zeros((6,6))
    A_sol[0,1] = 1.0
    A_sol[1,0] = -k_s**2
    A_sol[2,3] = 1.0
    A_sol[3,2] = -k_s**2
    # Larmor coupling
    A_sol[1,3] = -2.0*k_s
    A_sol[3,1] = 2.0*k_s

    regions = [ARegion(0.1, 0.2, A_sol)]
    A_func = make_A_func(regions)

    z_sol, st_sol = propagate_cov(cov0, (0.0, 0.3), n_points=100,
                                   external_A_func=A_func,
                                   include_sc=False)
    cov_end_sol = Covariance6D.from_state(st_sol[-1,:21], st_sol[-1,21], cov0.Ne)
    print(f"  σ_x: {cov_end_sol.sigma_x*1e6:.1f} μm  σ_y: {cov_end_sol.sigma_y*1e6:.1f} μm")
    print(f"  σ_xy = ⟨xy⟩/(σ_xσ_y) = {cov_end_sol.sigma[0,2]/(cov_end_sol.sigma_x*cov_end_sol.sigma_y):.6f}")

    # ── 5. Plot comparison ──
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    z_mm_d = z_drift * 1e3
    sx_d = np.array([np.sqrt(max(Covariance6D.from_state(s[:21],s[21],cov0.Ne).sigma[0,0],0))
                      for s in st_drift]) * 1e6
    axes[0,0].plot(z_mm_d, sx_d, 'b-', lw=1.5, label='Covariance model')
    axes[0,0].axhline(sx_analytic*1e6, color='r', ls='--', label='Analytic')
    axes[0,0].set_ylabel('σ_x (μm)'); axes[0,0].legend(fontsize=8); axes[0,0].grid(alpha=0.3)
    axes[0,0].set_title('Pure drift — benchmark vs analytic')

    z_mm_sc = z_sc * 1e3
    sx_sc = np.array([np.sqrt(max(Covariance6D.from_state(s[:21],s[21],cov0.Ne).sigma[0,0],0))
                       for s in st_sc]) * 1e6
    sz_sc = np.array([np.sqrt(max(Covariance6D.from_state(s[:21],s[21],cov0.Ne).sigma[4,4],0))
                       for s in st_sc]) * 1e6
    axes[0,1].plot(z_mm_sc, sx_sc, 'b-', lw=1.5, label='σ_x')
    axes[0,1].plot(z_mm_sc, sz_sc, 'g-', lw=1.5, label='σ_z')
    axes[0,1].set_ylabel('σ (μm)'); axes[0,1].legend(fontsize=8); axes[0,1].grid(alpha=0.3)
    axes[0,1].set_title('Drift with space charge')

    z_mm_sol = z_sol * 1e3
    sx_sol = np.array([np.sqrt(max(Covariance6D.from_state(s[:21],s[21],cov0.Ne).sigma[0,0],0))
                        for s in st_sol]) * 1e6
    sy_sol = np.array([np.sqrt(max(Covariance6D.from_state(s[:21],s[21],cov0.Ne).sigma[2,2],0))
                        for s in st_sol]) * 1e6
    xy_corr = np.array([Covariance6D.from_state(s[:21],s[21],cov0.Ne).sigma[0,2]
                         / max(np.sqrt(abs(Covariance6D.from_state(s[:21],s[21],cov0.Ne).sigma[0,0]
                                         *Covariance6D.from_state(s[:21],s[21],cov0.Ne).sigma[2,2])), 1e-12)
                         for s in st_sol])

    axes[1,0].plot(z_mm_sol, sx_sol, 'b-', lw=1.5, label='σ_x')
    axes[1,0].plot(z_mm_sol, sy_sol, 'r--', lw=1.5, label='σ_y')
    axes[1,0].axvline(100, color='gray', ls=':', alpha=0.5)
    axes[1,0].axvline(200, color='gray', ls=':', alpha=0.5)
    axes[1,0].set_xlabel('z (mm)'); axes[1,0].set_ylabel('σ (μm)')
    axes[1,0].legend(fontsize=8); axes[1,0].grid(alpha=0.3)
    axes[1,0].set_title('Solenoid focusing')

    axes[1,1].plot(z_mm_sol, xy_corr, 'm-', lw=1.5)
    axes[1,1].set_xlabel('z (mm)'); axes[1,1].set_ylabel('Correlation ⟨xy⟩/(σ_xσ_y)')
    axes[1,1].grid(alpha=0.3)
    axes[1,1].set_title('x-y coupling in solenoid')

    fig.suptitle("Phase 3: 6D Covariance Matrix — Benchmark Tests", fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig("6d_covariance_benchmark.png", dpi=150)
    print("\nSaved: 6d_covariance_benchmark.png")
    print("\n=== Phase 3 Initial Benchmarks Complete ===")
