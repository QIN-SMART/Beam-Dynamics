#!/usr/bin/env python3
"""
Six-Dimensional Beam Envelope Equations
========================================
Implements the 6D coupled envelope ODE framework based on:

  Kelisani, Barzegar, Craievich & Doebert,
  "Six-Dimensional Beam-Envelope Equations: An Ultrafast Computational
   Approach for Interactive Modeling of Accelerator Structures,"
  Phys. Rev. Applied 19, 054011 (2023).

The six state variables are:
  σ_x, σ_y, σ_z   — RMS beam sizes (m)
  ν_x, ν_y, ν_z   — RMS beam slopes σ'_x, σ'_y, σ'_z (dimensionless)

Plus:
  γ               — reference particle Lorentz factor (evolves with acceleration)

Envelope equation (Eq.3):
  σ_u'' + (γ γ')/(γ² β²) σ_u' = F_u^e + F_u^s + F_u^ε
"""

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import RectBivariateSpline
from dataclasses import dataclass
from typing import Tuple, Optional, Callable, List
import warnings

# ═════════════════════════════════════════════════════════════════════════════
# Physical constants (SI)
# ═════════════════════════════════════════════════════════════════════════════
M_E       = 9.10938356e-31
M_E_EV    = 5.109989461e5
C_LIGHT   = 2.99792458e8
E_CHARGE  = 1.60217662e-19
EPSILON_0 = 8.854187817e-12
R_E       = 2.8179403262e-15    # classical electron radius (m)


# ═════════════════════════════════════════════════════════════════════════════
# Relativistic kinematics helpers
# ═════════════════════════════════════════════════════════════════════════════
def gamma_to_beta(gamma: float) -> float:
    return np.sqrt(1.0 - 1.0 / gamma**2) if gamma > 1.0 else 0.0

def beamK_to_gamma(beamK: float) -> float:
    """beamK in eV → Lorentz γ."""
    return 1.0 + beamK / M_E_EV

def gamma_to_beamK(gamma: float) -> float:
    return (gamma - 1.0) * M_E_EV


# ═════════════════════════════════════════════════════════════════════════════
# Space-charge coefficients  α_ij(k_x, k_y)
# ═════════════════════════════════════════════════════════════════════════════

def _integrate_alpha(func, kx: float, ky: float) -> float:
    """
    Compute one of the α-coefficient integrals ∫₀^∞ f(s; kx, ky) ds.
    Uses the substitution  s = t/(1-t)  to map [0,∞) → [0,1).
    """
    # Handle degenerate cases
    kx = max(kx, 1e-12)
    ky = max(ky, 1e-12)

    def integrand(t):
        s = t / (1.0 - t) if t < 1.0 else 1e15
        ds_dt = 1.0 / (1.0 - t)**2
        return func(s, kx, ky) * ds_dt

    result, _ = quad(integrand, 0.0, 0.999999, limit=200)
    return result


def compute_all_alpha(kx: float, ky: float) -> dict:
    """
    Compute all 9 space-charge coefficients α_ij for given (k_x, k_y).
    k_x = σ_x / (γ σ_z),  k_y = σ_y / (γ σ_z).

    Returns dict with keys: 'x','y','z','xx','yy','zz','xy','xz','yz'.
    """
    eps_ = 1e-12
    kx = max(kx, eps_)
    ky = max(ky, eps_)

    # α_x  (Eq.21)
    def f_ax(s, kx, ky):
        denom = np.sqrt((s**2 + kx**2)**3 * (s**2 + ky**2) * (s**2 + 1.0))
        return kx**2 * s / denom

    # α_y  (Eq.22)
    def f_ay(s, kx, ky):
        denom = np.sqrt((s**2 + kx**2) * (s**2 + ky**2)**3 * (s**2 + 1.0))
        return ky**2 * s / denom

    # α_z  (Eq.23)
    def f_az(s, kx, ky):
        denom = np.sqrt((s**2 + kx**2) * (s**2 + ky**2) * (s**2 + 1.0)**3)
        return kx * ky * s / denom

    # α_xx  (Eq.24)
    def f_axx(s, kx, ky):
        denom = np.sqrt((s**2 + kx**2)**5 * (s**2 + ky**2) * (s**2 + 1.0))
        return 3.0 * kx**4 * s / denom

    # α_yy  (Eq.25)  — symmetric: α_yy(kx, ky) = α_xx(ky, kx)
    def f_ayy(s, kx, ky):
        denom = np.sqrt((s**2 + kx**2) * (s**2 + ky**2)**5 * (s**2 + 1.0))
        return 3.0 * ky**4 * s / denom

    # α_zz  (Eq.26)
    def f_azz(s, kx, ky):
        denom = np.sqrt((s**2 + kx**2) * (s**2 + ky**2) * (s**2 + 1.0)**5)
        return 3.0 * kx * ky * s / denom

    # α_xy  (Eq.27)
    def f_axy(s, kx, ky):
        denom = np.sqrt((s**2 + kx**2)**3 * (s**2 + ky**2)**3 * (s**2 + 1.0))
        return 4.0 * kx**2 * ky**2 * s / denom

    # α_xz  (Eq.28)
    def f_axz(s, kx, ky):
        denom = np.sqrt((s**2 + kx**2)**3 * (s**2 + ky**2) * (s**2 + 1.0)**3)
        return kx**2 * s / denom

    # α_yz  (Eq.29)
    def f_ayz(s, kx, ky):
        denom = np.sqrt((s**2 + kx**2) * (s**2 + ky**2)**3 * (s**2 + 1.0)**3)
        return ky**2 * s / denom

    return {
        'x':  _integrate_alpha(f_ax, kx, ky),
        'y':  _integrate_alpha(f_ay, kx, ky),
        'z':  _integrate_alpha(f_az, kx, ky),
        'xx': _integrate_alpha(f_axx, kx, ky),
        'yy': _integrate_alpha(f_ayy, kx, ky),
        'zz': _integrate_alpha(f_azz, kx, ky),
        'xy': _integrate_alpha(f_axy, kx, ky),
        'xz': _integrate_alpha(f_axz, kx, ky),
        'yz': _integrate_alpha(f_ayz, kx, ky),
    }


# ═════════════════════════════════════════════════════════════════════════════
# Cached space-charge coefficient table (lazy interpolation)
# ═════════════════════════════════════════════════════════════════════════════

def _build_alpha_cache(
    kx_range: Tuple[float, float] = (1e-2, 1e2),
    ky_range: Tuple[float, float] = (1e-2, 1e2),
    n_points: int = 15,
) -> dict:
    """Precompute α-coefficient tables on a log-spaced grid."""
    kx_vals = np.logspace(np.log10(kx_range[0]), np.log10(kx_range[1]), n_points)
    ky_vals = np.logspace(np.log10(ky_range[0]), np.log10(ky_range[1]), n_points)
    keys = ['x', 'y', 'z', 'xx', 'yy', 'zz', 'xy', 'xz', 'yz']
    grid = {k: np.zeros((n_points, n_points)) for k in keys}

    for i, kx in enumerate(kx_vals):
        for j, ky in enumerate(ky_vals):
            coeffs = compute_all_alpha(kx, ky)
            for k in keys:
                grid[k][i, j] = coeffs[k]

    splines = {}
    for k in keys:
        splines[k] = RectBivariateSpline(kx_vals, ky_vals, grid[k])

    return splines


# Global cache — built lazily on first use
_ALPHA_CACHE: Optional[dict] = None


def get_alpha_interpolators() -> dict:
    global _ALPHA_CACHE
    if _ALPHA_CACHE is None:
        _ALPHA_CACHE = _build_alpha_cache()
    return _ALPHA_CACHE


def alpha_at(kx: float, ky: float) -> dict:
    """Interpolated α coefficients at (kx, ky).  Clips to grid range."""
    kx = np.clip(kx, 1e-2, 1e2)
    ky = np.clip(ky, 1e-2, 1e2)
    spl = get_alpha_interpolators()
    return {k: float(spl[k](kx, ky).item()) for k in spl}


# ═════════════════════════════════════════════════════════════════════════════
# Beam State
# ═════════════════════════════════════════════════════════════════════════════

@dataclass
class Beam6D:
    """
    6D beam envelope state with x-y coupling.

    Parameters
    ----------
    sigma_x, sigma_y, sigma_z : float [m]
        RMS beam sizes.
    sigma_xy : float [m²]
        x-y correlation term  ⟨δx·δy⟩.  0 for uncorrelated beam.
        Evolves in solenoids via Larmor rotation: d⟨xy⟩/dz ∝ k_L(σ_x²−σ_y²).
    nu_x, nu_y, nu_z : float [dimensionless]
        RMS beam slopes  (dσ/dz in the envelope equation).
    nu_xy : float [m]
        z-derivative of sigma_xy:  d⟨xy⟩/dz.
    gamma : float
        Lorentz factor of the reference particle.
    Ne : float
        Number of electrons in the bunch.
    eps_nx, eps_ny, eps_nz : float [m·rad]
        Normalised rms emittances.
    """
    sigma_x: float
    sigma_y: float
    sigma_z: float
    sigma_xy: float = 0.0
    sigma_delta: float = 0.0    # σ_δ = RMS relative momentum spread ⟨δ²⟩^(1/2)
    C_zd: float = 0.0           # C_zδ = ⟨z·δ⟩, z-δ correlation [m]
    nu_x: float = 0.0
    nu_y: float = 0.0
    nu_z: float = 0.0
    nu_xy: float = 0.0
    nu_delta: float = 0.0       # dσ_δ/dz
    gamma: float = 1.0
    Ne: float = 0.0
    eps_nx: float = 0.0
    eps_ny: float = 0.0
    eps_nz: float = 0.0

    def __post_init__(self):
        pass

    @property
    def beta(self) -> float:
        return gamma_to_beta(self.gamma)

    @property
    def p0(self) -> float:
        """Normalised momentum γβ."""
        return self.gamma * self.beta

    @property
    def beamK(self) -> float:
        """Kinetic energy in eV."""
        return gamma_to_beamK(self.gamma)

    @property
    def theta_x(self) -> float:
        """
        RMS divergence in x (rad).
        From emittance conservation:  ε_geo² = σ²θ² − σ²σ'²
        →  θ_x = sqrt(ε_geo² / σ_x² + ν_x²)
        """
        p0 = self.p0
        if p0 < 1e-12 or self.sigma_x < 1e-15:
            return 0.0
        eps_geo = self.eps_nx / p0
        return np.sqrt(eps_geo**2 / self.sigma_x**2 + self.nu_x**2)

    @property
    def theta_y(self) -> float:
        p0 = self.p0
        if p0 < 1e-12 or self.sigma_y < 1e-15:
            return 0.0
        eps_geo = self.eps_ny / p0
        return np.sqrt(eps_geo**2 / self.sigma_y**2 + self.nu_y**2)

    @property
    def theta_z(self) -> float:
        """Longitudinal RMS divergence (dimensionless slope)."""
        p0 = self.p0
        if p0 < 1e-12 or self.sigma_z < 1e-15:
            return 0.0
        eps_geo = self.eps_nz / p0
        return np.sqrt(eps_geo**2 / self.sigma_z**2 + self.nu_z**2)

    def state_vector(self) -> np.ndarray:
        """Return [σx,σy,σz,σxy,σδ,C_zd, νx,νy,νz,νxy,νδ]."""
        return np.array([self.sigma_x, self.sigma_y, self.sigma_z,
                         self.sigma_xy, self.sigma_delta, self.C_zd,
                         self.nu_x, self.nu_y, self.nu_z, self.nu_xy, self.nu_delta])

    @classmethod
    def from_state(cls, y: np.ndarray, gamma: float, Ne: float,
                   eps_nx: float, eps_ny: float, eps_nz: float) -> "Beam6D":
        return cls(sigma_x=y[0], sigma_y=y[1], sigma_z=y[2],
                   sigma_xy=y[3], sigma_delta=y[4], C_zd=y[5],
                   nu_x=y[6], nu_y=y[7], nu_z=y[8], nu_xy=y[9], nu_delta=y[10],
                   gamma=gamma, Ne=Ne, eps_nx=eps_nx, eps_ny=eps_ny, eps_nz=eps_nz)

    def copy(self) -> "Beam6D":
        return Beam6D(self.sigma_x, self.sigma_y, self.sigma_z, self.sigma_xy,
                      self.sigma_delta, self.C_zd,
                      self.nu_x, self.nu_y, self.nu_z, self.nu_xy, self.nu_delta,
                      self.gamma, self.Ne, self.eps_nx, self.eps_ny, self.eps_nz)

    @property
    def time_spread_ps(self) -> float:
        """RMS time spread (ps).  For UED, this is the key temporal resolution metric."""
        v = self.beta * C_LIGHT
        return self.sigma_z / v * 1e12 if v > 1e3 else 0.0

    def summary(self, label: str = "") -> str:
        lines = [f"══ {label} ══" if label else "──"]
        lines.append(f"  γ, β         = {self.gamma:.4f}, {self.beta:.4f}")
        lines.append(f"  K            = {self.beamK:.0f} eV")
        lines.append(f"  Ne           = {self.Ne:.2e}")
        lines.append(f"  σ_x, σ_y     = {self.sigma_x*1e6:.2f}, {self.sigma_y*1e6:.2f} μm")
        lines.append(f"  σ_z          = {self.sigma_z*1e6:.2f} μm")
        if abs(self.sigma_xy) > 1e-20:
            rms_xy = np.sqrt(abs(self.sigma_xy)) * 1e6
            lines.append(f"  σ_xy         = {self.sigma_xy*1e12:.4f} pm²  (RMS ~{rms_xy:.2f} μm)")
        else:
            lines.append(f"  σ_xy         = 0 (uncoupled)")
        lines.append(f"  ν_x, ν_y     = {self.nu_x:.6f}, {self.nu_y:.6f}")
        lines.append(f"  ν_z          = {self.nu_z:.6f}")
        lines.append(f"  σ_δ          = {self.sigma_delta*1e3:.4f} ‰")
        lines.append(f"  C_zδ         = {self.C_zd*1e6:.4f} μm")
        lines.append(f"  ε_nx/ny/nz   = {self.eps_nx*1e6:.4f}/{self.eps_ny*1e6:.4f}/{self.eps_nz*1e6:.4f} μm")
        lines.append(f"  ϑ_x, ϑ_y     = {self.theta_x*1e3:.4f}, {self.theta_y*1e3:.4f} mrad")
        lines.append(f"  ϑ_z          = {self.theta_z:.6f}")
        return "\n".join(lines)

    def __repr__(self):
        return (f"Beam6D(σ=({self.sigma_x*1e6:.1f},{self.sigma_y*1e6:.1f},"
                f"{self.sigma_z*1e6:.1f}) μm, γ={self.gamma:.3f})")


# ═════════════════════════════════════════════════════════════════════════════
# Force terms  (drift-only for now — no external fields)
# ═════════════════════════════════════════════════════════════════════════════

def emittance_force(sigma_u: float, nu_u: float, theta_u: float) -> float:
    """
    Emittance (thermal pressure) force  F_u^ε = (ϑ_u² − ν_u²) / σ_u.
    Eq.(6) of Kelisani 2023.
    """
    if sigma_u < 1e-15:
        return 0.0
    return (theta_u**2 - nu_u**2) / sigma_u


def _space_charge_forces(beam: Beam6D) -> Tuple[float, float, float]:
    """
    Compute space-charge forces F_x^s, F_y^s, F_z^s for a 6D Gaussian beam.

    Uses the full analytical expressions Eqs.(15)-(17) from Kelisani 2023,
    with α coefficients interpolated from a precomputed table.
    """
    gamma = beam.gamma
    beta = beam.beta
    p0 = beam.p0
    Ne = beam.Ne

    sx, sy, sz = beam.sigma_x, beam.sigma_y, beam.sigma_z
    nx, ny, nz = beam.nu_x, beam.nu_y, beam.nu_z
    tx, ty, tz = beam.theta_x, beam.theta_y, beam.theta_z

    # Aspect ratios in the beam rest frame  (clipped to table range)
    kx_ = np.clip(sx / (gamma * sz) if sz > 1e-15 else 1e3, 1e-2, 1e2)
    ky_ = np.clip(sy / (gamma * sz) if sz > 1e-15 else 1e3, 1e-2, 1e2)

    # Get α coefficients by interpolation
    alpha = alpha_at(kx_, ky_)
    ax, ay, az = alpha['x'], alpha['y'], alpha['z']
    axx, ayy, azz = alpha['xx'], alpha['yy'], alpha['zz']
    axy, axz, ayz = alpha['xy'], alpha['xz'], alpha['yz']

    # Beam factor  f_b = η q_b / (8π √π ε₀)
    # η = e/mc²,  q_b = Ne·e
    eta = E_CHARGE / (M_E * C_LIGHT**2)
    fb = eta * Ne * E_CHARGE / (8.0 * np.pi * np.sqrt(np.pi) * EPSILON_0)

    beta2 = max(beta**2, 1e-10)
    gamma3 = gamma**3
    p02 = max(p0**2, 1e-10)

    # Safety: clip denominator for very low β (avoids NaN at cathode)
    denom_x = max(beta2 * gamma3 * sx * sz, 1e-30)
    denom_y = max(beta2 * gamma3 * sy * sz, 1e-30)
    denom_z = max(beta2 * gamma3 * sx * sy, 1e-30)

    # ── F_x^s  (Eq.15) ──
    # Leading term (∝ 1/(β²γ³))
    term_x0 = fb * ax / denom_x

    # Curly-brace term A
    A_x = ((tx**2 + 2*nx**2) * ax - nx**2 * axx) / (2*sx*sz) \
        + (8*ty**2 * ax - ny**2 * axy) / (16*sx*sz) \
        + (1 - p0**2) * (2*tz**2 * ax - nz**2 * axz) / (4*sx*sz)

    # Curly-brace term B
    B_x = ((tx**2 + 2*nx**2) * ax - nx**2 * axx) / (sx*sz) \
        + nx*ny * (8*ay - axy) / (8*sz*sy) \
        + (1 - p0**2) * nx*nz * (2*gamma**2*sz**2*az - sx*sy*axz) / (2*sx*sy*gamma**2*sz**2)

    Fs_x = term_x0 - fb/gamma * (A_x + B_x)

    # ── F_y^s  (Eq.16) ── (symmetric with x↔y)
    term_y0 = fb * ay / denom_y

    A_y = ((ty**2 + 2*ny**2) * ay - ny**2 * ayy) / (2*sy*sz) \
        + (8*tx**2 * ay - nx**2 * axy) / (16*sy*sz) \
        + (1 - p0**2) * (2*tz**2 * ay - nz**2 * ayz) / (4*sy*sz)

    B_y = ((ty**2 + 2*ny**2) * ay - ny**2 * ayy) / (sy*sz) \
        + ny*nx * (8*ax - axy) / (8*sz*sx) \
        + (1 - p0**2) * ny*nz * (2*gamma**2*sz**2*az - sy*sx*ayz) / (2*sy*sx*gamma**2*sz**2)

    Fs_y = term_y0 - fb/gamma * (A_y + B_y)

    # ── F_z^s  (Eq.17) ──
    term_z0 = fb * az / denom_z

    A_z = 3*(1 - p0**2) * (2*(tz**2 + 2*nz**2)*az - nz**2*azz) / (4*sx*sy) \
        + (2*gamma**2*sz**2*tx**2*az - nx**2*sx*sy*axz) / (4*sx*sy*gamma**2*sz**2) \
        + (2*gamma**2*sz**2*ty**2*az - ny**2*sx*sy*ayz) / (4*sx*sy*gamma**2*sz**2)

    B_z = nx*nz*(2*ax - axz) / (2*sx*sz) \
        + ny*nz*(2*ay - ayz) / (2*sy*sz)

    Fs_z = term_z0 - fb/gamma * (A_z + B_z)

    return Fs_x, Fs_y, Fs_z


# ═════════════════════════════════════════════════════════════════════════════
# Space-charge model — Uniform 3D Ellipsoid  (Luiten et al., PRL 93, 2004)
# ═════════════════════════════════════════════════════════════════════════════

def _ellipsoid_form_factors(sigma_x: float, sigma_y: float, sigma_z: float
                            ) -> Tuple[float, float, float]:
    """
    Compute the geometric form factors (Mx, My, Mz) for a uniform ellipsoid.

    For a uniform ellipsoid with semi-axes A, B, C, the internal electrostatic
    field is:  E = (ρ₀/ε₀) · (Mx·x, My·y, Mz·z)
    where ρ₀ = 3Ne/(4πABC) and Mx+My+Mz = 1.

    The RMS sizes are related to the semi-axes by: σ_u = A_u/√5.
    So the aspect ratios are σ_u/σ_v = A_u/A_v.

    For the general triaxial case, Mx, My, Mz involve elliptic integrals.
    For the spheroid case (σx=σy) we use the analytical formula Eq.(4) of Luiten.

    References
    ----------
    Luiten et al., PRL 93, 094802 (2004), Eq.(4).
    Kellogg, "Foundations of Potential Theory" (1929).
    """
    # Aspect ratio ξ = √(A²/C² − 1) = √(σ_x²/σ_z² − 1)
    # For prolates (A<C), ξ becomes imaginary — handle via arctanh.
    ratio = sigma_x / max(sigma_z, 1e-15)

    if abs(sigma_x - sigma_y) / max(sigma_x, sigma_y, 1e-15) < 1e-6:
        # ── Spheroid: σ_x = σ_y ──
        if ratio >= 1.0:
            # Oblate: disk-like  (A > C)
            xi = np.sqrt(ratio**2 - 1.0)
            if xi < 1e-4:
                Mz = 1.0/3.0 + xi**2/5.0 - 3.0*xi**4/35.0   # Taylor for xi→0
            else:
                Mz = (1.0 + xi**2) / xi**3 * (xi - np.arctan(xi))
        else:
            # Prolate: cigar-like  (A < C)
            xi = np.sqrt(1.0 - ratio**2)
            if xi < 1e-4:
                Mz = 1.0/3.0 - xi**2/5.0 + 3.0*xi**4/35.0   # Taylor for xi→0
            else:
                Mz = (1.0 - xi**2) / xi**3 * (np.arctanh(xi) - xi)
        Mx = 0.5 * (1.0 - Mz)
        My = Mx
    else:
        # ── Triaxial ellipsoid — use Carlson symmetric form ──
        # This is a robust numerical integration for all aspect ratios.
        # Using R_J(x,y,z) = (3/2) ∫₀^∞ dt/√((t+x)(t+y)(t+z))
        # Mx = σ_y σ_z · R_D(σ_x², σ_y², σ_z²) etc.
        # where R_D is the Carlson elliptic integral of the second kind.
        A2 = sigma_x**2
        B2 = sigma_y**2
        C2 = sigma_z**2
        # Carlson R_D integral: R_D(x,y,z) = (3/2)∫₀^∞ dt/((t+z)√((t+x)(t+y)(t+z)))
        # Scale by √(ABC) for the form factor
        from scipy.special import elliprd
        # elliprd(x,y,z) returns R_D(x,y,z)
        scale = sigma_x * sigma_y * sigma_z
        Mx = scale * elliprd(B2, C2, A2) / 3.0
        My = scale * elliprd(A2, C2, B2) / 3.0
        # Mz from normalization: Mx + My + Mz = 1
        Mz = 1.0 - Mx - My
        # Clip for numerical noise
        Mx = max(0.0, min(1.0, Mx))
        My = max(0.0, min(1.0, My))
        Mz = max(0.0, min(1.0, Mz))

    return Mx, My, Mz


def _ellipsoid_space_charge_forces(beam: Beam6D) -> Tuple[float, float, float]:
    """
    Compute space-charge forces for a uniform 3D ellipsoidal bunch.

    Uses the analytical form factors Mx, My, Mz from Luiten et al. (2004).
    For a uniform ellipsoid the internal electric field is perfectly linear:

        E_u = (3Ne·e / (4πε₀·ABC)) · M_u · δu

    In the envelope equation (Kelisani 2023 Eq.3), the leading-order
    space-charge force is:

        F^s_u ≈ -(Ne·r_e·M_u) / (γ³·β²·σ_u·σ_v)

    where u is the transverse index (x or y) and v is the longitudinal (z).
    For the longitudinal force: F^s_z ≈ -(Ne·r_e·M_z) / (γ³·β²·σ_x·σ_y)

    Reference: Luiten et al., PRL 93, 094802 (2004), Eq.(3)-(4).
    """
    gamma = beam.gamma
    beta = beam.beta
    Ne = beam.Ne

    sx = max(beam.sigma_x, 1e-15)
    sy = max(beam.sigma_y, 1e-15)
    sz = max(beam.sigma_z, 1e-15)

    Mx, My, Mz = _ellipsoid_form_factors(sx, sy, sz)

    gamma3 = gamma**3
    beta2 = max(beta**2, 1e-10)
    common = Ne * R_E / (max(gamma3 * beta2, 1e-30))

    # Space charge is repulsive for electrons → defocusing → F_s > 0
    # (Positive drive term in σ'' = ... + F_s means σ increases)
    Fs_x = +common * Mx / (sx * sz)
    Fs_y = +common * My / (sy * sz)
    Fs_z = +common * Mz / (sx * sy)

    return Fs_x, Fs_y, Fs_z


def space_charge_forces(beam: Beam6D, model: str = 'gaussian'
                        ) -> Tuple[float, float, float]:
    """
    Compute space-charge forces using the selected model.

    Parameters
    ----------
    beam : Beam6D
    model : str
        'gaussian' — Kelisani 2023 9-coefficient model (6D correlated Gaussian).
        'ellipsoid' — Luiten 2004 uniform ellipsoid form factors.

    Returns (F^s_x, F^s_y, F^s_z).
    """
    if model == 'ellipsoid':
        return _ellipsoid_space_charge_forces(beam)
    else:
        return _space_charge_forces(beam)


# ═════════════════════════════════════════════════════════════════════════════
# Acceleration-field energy spread  (Stupakov & Huang, PRSTAB 11, 2008)
# ═════════════════════════════════════════════════════════════════════════════

def acceleration_energy_spread(beam: Beam6D, E_accel: float, gamma_i: float = 1.0,
                               gamma_f: float = None) -> float:
    """
    Estimate the RMS energy spread (eV) induced by the acceleration field.

    Stupakov & Huang (PRSTAB 11, 014401, 2008) show that during longitudinal
    acceleration, the beam's changing electromagnetic field energy creates an
    additional self-field that scales as E_accel/γ (not γ⁻² like static SC).

    For a Gaussian bunch, the RMS energy spread is (Eq.19):
      ΔE_rms ≈ (Q / 4π ε₀) · (a / γ³βc) · (1 / σ_r) · G(σ_r/σ_z)
    where a = dv/dt is the acceleration, G is a geometric factor ~0.3.

    Parameters
    ----------
    beam : Beam6D
    E_accel : float [V/m]
        Accelerating electric field.
    gamma_i, gamma_f : float
        Initial and final Lorentz factors.

    Returns
    -------
    dE_rms : float [eV]

    Notes
    -----
    For UEM (Q~0.016 pC, E~2.5 MV/m, γ~1.2): negligible (~1e-5 eV).
    For LCLS RF gun (Q~0.72 nC, E~100 MV/m, γ~20): ~0.5 keV.
    """
    if gamma_f is None:
        gamma_f = beam.gamma

    Q = beam.Ne * E_CHARGE
    sigma_r = max(beam.sigma_x, 1e-15)
    sigma_z = max(beam.sigma_z, 1e-15)

    # Acceleration: a = e·E/(m·γ³) (proper acceleration in lab frame)
    # Actually for the energy spread effect, use the lab-frame acceleration dv/dt
    gamma = beam.gamma
    beta = beam.beta
    if beta < 1e-6 or gamma < 1.001:
        return 0.0

    a_lab = E_CHARGE * E_accel / (M_E * gamma**3)   # dv/dt in lab frame

    # Geometric scale: the field varies on scale ~σ_r; Gaussian factor
    aspect = sigma_r / max(sigma_z, 1e-15)
    geom = 1.0 / sigma_r
    if aspect < 100:   # avoid exp overflow
        geom *= np.exp(-0.5 / max(aspect**2, 1e-10))

    # Energy loss scale (Stupakov Eq. 18-19)
    dE_scale = Q / (4.0 * np.pi * EPSILON_0) * a_lab / (gamma**3 * beta * C_LIGHT)

    dE_rms_joules = abs(dE_scale * geom)
    dE_rms_eV = dE_rms_joules / E_CHARGE

    # For practical purposes, clip at reasonable values
    if dE_rms_eV > 1e9:   # physically impossible for non-plasma
        dE_rms_eV = 0.0

    return dE_rms_eV

ETA = E_CHARGE / (M_E * C_LIGHT**2)          # η = e/(m_e c²)  [1/V·m]


def solenoid_force(beam: Beam6D, Bz_on_axis: float, dBz_dz: float = 0.0
                   ) -> Tuple[float, float, float, float]:
    """
    Solenoid F-type forces (lowest order, Eqs.31-33) + x-y Larmor coupling.

    Parameters
    ----------
    Bz_on_axis : float [T]
        On-axis longitudinal magnetic field μ_SM.
    dBz_dz : float [T/m]
        Axial derivative ∂μ_SM/∂z.

    Returns (Fe_x, Fe_y, Fe_z, dσ_xy/dz contribution).
        The Larmor coupling generates σ_xy ∝ σ_x²−σ_y².
    """
    p0 = beam.p0
    mu = Bz_on_axis
    mu_z = dBz_dz

    factor = -(ETA**2 * C_LIGHT**2) / (4.0 * p0**2)
    # Main focusing term ∝ μ_SM² + (μ_SM_z · σ_z)²
    ku = mu**2 + (mu_z * beam.sigma_z)**2

    Fe_x = factor * ku * beam.sigma_x
    Fe_y = factor * ku * beam.sigma_y
    Fe_z = 0.0    # solenoid does not focus longitudinally

    # Larmor wavenumber: k_L = e Bz / (2 p) = η c Bz / (2 p0)
    # (η = e/(m c²), p0 = γβ, so η c Bz / (2 p0) has units 1/m)
    k_L = ETA * C_LIGHT * mu / (2.0 * max(p0, 1e-10))
    # X-Y coupling: d⟨xy⟩/dz = 2 k_L (⟨x²⟩ − ⟨y²⟩) in solenoid
    d_sigma_xy = 2.0 * k_L * (beam.sigma_x**2 - beam.sigma_y**2)

    # The Larmor rotation also couples ν_x and ν_y:
    # dν_x ⊃ -2 k_L ν_y,  dν_y ⊃ +2 k_L ν_x
    # These are returned as additional contributions
    return Fe_x, Fe_y, Fe_z, float(d_sigma_xy)


def quadrupole_force(beam: Beam6D, k_QM: float
                     ) -> Tuple[float, float, float]:
    """
    Quadrupole F-type forces (lowest order, Eqs.45-47).

    Parameters
    ----------
    k_QM : float [1/m²]
        Quadrupole strength.  k_QM > 0 → focusing in x, defocusing in y.

    Returns (Fe_x, Fe_y, Fe_z).
    """
    Fe_x = -k_QM * beam.sigma_x
    Fe_y = +k_QM * beam.sigma_y
    Fe_z = 0.0
    return Fe_x, Fe_y, Fe_z


def electrostatic_force(beam: Beam6D, E_on_axis: float, dE_dz: float
                        ) -> Tuple[float, float, float]:
    """
    Electrostatic lens / DC acceleration gap (lowest order, Eqs.38-40).

    Parameters
    ----------
    E_on_axis : float [V/m]
        On-axis electric field ε_el (= E for uniform DC gap).
    dE_dz : float [V/m²]
        Axial derivative ∂ε_el/∂z.  For uniform gap, dE/dz=0 except at edges.

    Returns (Fe_x, Fe_y, Fe_z).
    """
    p0 = beam.p0
    gamma = beam.gamma
    beta = beam.beta
    E0 = E_on_axis
    Ez = dE_dz

    gamma_beta2 = gamma * beta**2

    if abs(gamma_beta2) < 1e-15:
        return 0.0, 0.0, 0.0

    # Transverse: defocusing proportional to dE/dz
    Fe_x = -ETA * Ez / (2.0 * gamma_beta2) * beam.sigma_x
    Fe_y = -ETA * Ez / (2.0 * gamma_beta2) * beam.sigma_y

    # Longitudinal: focusing in accelerating gap
    Fe_z = +ETA * Ez / (gamma * p0**2) * beam.sigma_z

    return Fe_x, Fe_y, Fe_z


def rf_cavity_force(beam: Beam6D, E_rf: float, dE_dz: float, dE_cdt: float
                    ) -> Tuple[float, float, float]:
    """
    RF cavity F-type forces (lowest order, Eqs.53-55).

    Parameters
    ----------
    E_rf : float [V/m]
        On-axis peak RF electric field.
    dE_dz : float [V/m²]
        Axial derivative ∂E_rf/∂z  (= k·E_rf for standing wave, ≈0 for TW).
    dE_cdt : float [V/m²]
        Time derivative (1/c)·∂E_rf/∂t  (= k·E_rf for standing wave with π/2 shift).

    Returns (Fe_x, Fe_y, Fe_z).
    """
    p0 = beam.p0
    gamma = beam.gamma
    beta = beam.beta
    sigma_z = beam.sigma_z

    gamma_beta2 = gamma * beta**2
    if abs(gamma_beta2) < 1e-15:
        return 0.0, 0.0, 0.0

    # Effective transverse "kick" from combined electric + magnetic
    Erf_z_beta = dE_dz + beta * dE_cdt

    # Transverse defocusing
    factor_trans = -ETA * Erf_z_beta / (2.0 * gamma_beta2)
    Fe_x = factor_trans * beam.sigma_x
    Fe_y = factor_trans * beam.sigma_y

    # Longitudinal focusing / defocusing
    Fe_z = +ETA * dE_dz / (gamma * p0**2) * sigma_z

    return Fe_x, Fe_y, Fe_z


# ═════════════════════════════════════════════════════════════════════════════
# Beamline element specification for ODE integration
# ═════════════════════════════════════════════════════════════════════════════

@dataclass
class ExtFieldRegion:
    """Describes a region where a specific external field is active."""
    z_start: float       # [m] region start
    z_end: float         # [m] region end
    ftype: str           # 'solenoid', 'quadrupole', 'electrostatic', 'rf'
    params: dict         # parameters passed to the force function


def make_external_force_func(regions: List[ExtFieldRegion]) -> Callable:
    """
    Build a composite external-force function from a list of field regions.

    Returns a function (z, beam) → (Fe_x,Fe_y,Fe_z, dσ_xy_sol, dν_x_sol, dν_y_sol)
    where the last three are solenoid Larmor coupling contributions.
    """
    def force_func(z: float, beam: Beam6D) -> Tuple[float,float,float,float,float,float]:
        fx, fy, fz = 0.0, 0.0, 0.0
        dsxy = 0.0          # dσ_xy/dz from solenoid Larmor rotation
        dnx_sol = 0.0       # Larmor coupling to dν_x
        dny_sol = 0.0       # Larmor coupling to dν_y
        for reg in regions:
            if reg.z_start <= z <= reg.z_end:
                p = reg.params
                if reg.ftype == 'solenoid':
                    fx1, fy1, fz1, dsxy1 = solenoid_force(beam, p['Bz'], p.get('dBz_dz', 0.0))
                    fx += fx1; fy += fy1; fz += fz1; dsxy += dsxy1
                    # Larmor coupling: dν_x ⊃ −2k_L ν_y, dν_y ⊃ +2k_L ν_x
                    k_L = ETA * C_LIGHT * p['Bz'] / (2.0 * max(beam.p0, 1e-10))
                    dnx_sol -= 2.0 * k_L * beam.nu_y
                    dny_sol += 2.0 * k_L * beam.nu_x
                elif reg.ftype == 'quadrupole':
                    fx1, fy1, fz1 = quadrupole_force(beam, p['k_QM'])
                    fx += fx1; fy += fy1; fz += fz1
                elif reg.ftype == 'electrostatic':
                    fx1, fy1, fz1 = electrostatic_force(beam, p['E0'], p['dE_dz'])
                    fx += fx1; fy += fy1; fz += fz1
                elif reg.ftype == 'rf':
                    fx1, fy1, fz1 = rf_cavity_force(beam, p['E_rf'], p['dE_dz'], p['dE_cdt'])
                    fx += fx1; fy += fy1; fz += fz1
        return fx, fy, fz, dsxy, dnx_sol, dny_sol

    return force_func


def make_rf_chirp_func(regions: List[ExtFieldRegion]) -> Callable:
    """
    Build a function that returns the RF chirp rate h = dδ/dz at position z.

    For an RF cavity operating at phase φ (w.r.t. zero-crossing),
    the energy modulation per unit length is:
      h = e·E_rf·k / (β²γ³·mc²)   (for a standing-wave cavity)
    where k = 2π/λ_rf, E_rf is the peak on-axis field.

    Returns a function (z, beam) → h [1/m].
    """
    def chirp_func(z: float, beam: Beam6D) -> float:
        h_tot = 0.0
        for reg in regions:
            if reg.z_start <= z <= reg.z_end and reg.ftype == 'rf':
                E_rf = reg.params.get('E_rf', 0.0)
                k_rf = reg.params.get('k_rf', 2.0*np.pi*2.856e9/C_LIGHT)
                gamma = beam.gamma
                beta = beam.beta
                h = -E_CHARGE * E_rf * k_rf / (max(beta**2, 1e-10) * gamma**3 * M_E * C_LIGHT**2)
                h_tot += h
        return h_tot

    return chirp_func


def apply_rf_thin_lens(beam: Beam6D, H: float) -> Beam6D:
    """
    Apply RF longitudinal thin-lens transformation.

    Thin lens: δ_out = δ_in + H·z  (H = total chirp [1/m])

    σ_δ² → σ_δ² + 2H·C_zδ + H²·σ_z²
    C_zδ → C_zδ + H·σ_z²
    ν_z   → C_zδ_new / (γ²·σ_z)   (consistent with drift kinematics)
    σ_z   unchanged by thin lens
    """
    b = beam.copy()
    sz2 = beam.sigma_z**2
    new_var_delta = beam.sigma_delta**2 + 2*H*beam.C_zd + H**2 * sz2
    b.sigma_delta = np.sqrt(max(new_var_delta, 0.0))
    b.C_zd = beam.C_zd + H * sz2
    # Update ν_z to be consistent with post-lens σ_z evolution in drift
    g2 = max(beam.gamma**2, 1.001)
    b.nu_z = b.C_zd / (g2 * max(beam.sigma_z, 1e-15))
    return b


def compute_rf_chirp_coefficient(beam: Beam6D, E_rf: float, L_cav: float,
                                  f_rf: float = 2.856e9) -> float:
    """
    Compute the total RF chirp coefficient H over a cavity of length L_cav.

    H = -e·E_rf·k·L_cav / (β²·γ³·mc²)

    where k = 2πf/c.  δ_out = δ_in + H·z after passing through cavity.
    H > 0: head gains energy (positive chirp) → compression in drift.
    H < 0: tail gains energy (negative chirp) → decompression in drift.
    """
    k_rf = 2.0 * np.pi * f_rf / C_LIGHT
    beta2 = max(beam.beta**2, 1e-10)
    gamma3 = beam.gamma**3
    return -E_CHARGE * E_rf * k_rf * L_cav / (beta2 * gamma3 * M_E * C_LIGHT**2)


def compute_R56(L_drift: float, gamma: float, beta: float) -> float:
    """
    Compute R_56 transport element for a drift of length L.

    R_56 = L / (γ²·β²)   (non-relativistic correction included)
    """
    g2b2 = max(gamma**2 * beta**2, 1e-10)
    return L_drift / g2b2


def make_acceleration_func(regions: List[ExtFieldRegion]) -> Callable:
    """
    Build a function that returns dγ/dz at position z.

    dγ/dz = e·E_on_axis / (m_e·c²) = η · E_on_axis

    This is the energy gain per unit length for a particle in an
    accelerating electric field.  Independent of β in the relativistic
    formulation (work done = F·dz, dγ = e·E·dz / mc²).
    """
    def accel_func(z: float, beam: Beam6D) -> float:
        dgamma = 0.0
        for reg in regions:
            if reg.z_start <= z <= reg.z_end:
                if reg.ftype == 'electrostatic':
                    E0 = reg.params.get('E0', 0.0)
                    dgamma += ETA * E0
                elif reg.ftype == 'rf':
                    E0 = reg.params.get('E_rf', 0.0)
                    dgamma += ETA * E0
        return dgamma

    return accel_func


# ═════════════════════════════════════════════════════════════════════════════
# ODE right-hand side  —  d/dz of state vector
# ═════════════════════════════════════════════════════════════════════════════

def envelope_ode(z: float, y: np.ndarray,
                 Ne: float,
                 eps_nx: float, eps_ny: float, eps_nz: float,
                 external_force_func: Optional[Callable] = None,
                 gamma_prime_func: Optional[Callable] = None,
                 rf_chirp_func: Optional[Callable] = None,
                 sc_model: str = 'gaussian',
                 ) -> np.ndarray:
    """
    RHS of the 11-variable envelope ODE system.

    State: [σx,σy,σz,σxy,σ_δ,C_zδ, νx,νy,νz,νxy,ν_δ, γ]
    Returns dy/dz = [νx,νy,νz,νxy,ν_δ', dC/dz, dνx/dz,..., dγ/dz]

    Parameters
    ----------
    rf_chirp_func : callable, optional
        (z, beam) → h  where h = dδ/dz at z (RF chirp rate [1/m]).
    """
    gamma = y[11]
    beam = Beam6D.from_state(y[:11], gamma, Ne, eps_nx, eps_ny, eps_nz)

    beta = beam.beta
    p0 = beam.p0

    if gamma_prime_func is not None:
        gamma_prime = gamma_prime_func(z, beam)
    else:
        gamma_prime = 0.0

    gamma2_beta2 = gamma**2 * beta**2
    if gamma2_beta2 > 1e-15 and gamma > 1.0:
        damping_coeff = gamma * gamma_prime / gamma2_beta2
    else:
        damping_coeff = 0.0

    # ── Transverse emittance forces ──
    F_eps_x = emittance_force(beam.sigma_x, beam.nu_x, beam.theta_x)
    F_eps_y = emittance_force(beam.sigma_y, beam.nu_y, beam.theta_y)
    F_eps_z = emittance_force(beam.sigma_z, beam.nu_z, beam.theta_z)

    # ── Space-charge forces ──
    Fs_x, Fs_y, Fs_z = space_charge_forces(beam, model=sc_model)

    # ── External forces + solenoid coupling ──
    if external_force_func is not None:
        Fe_x, Fe_y, Fe_z, dsxy_sol, dnx_sol, dny_sol = external_force_func(z, beam)
    else:
        Fe_x, Fe_y, Fe_z, dsxy_sol, dnx_sol, dny_sol = 0.0,0.0,0.0, 0.0,0.0,0.0

    # ── RF chirp rate at position z ──
    h_rate = 0.0
    if rf_chirp_func is not None:
        h_rate = rf_chirp_func(z, beam)

    # ── Longitudinal (z, δ) transport ──
    # Drift contributions (always active):
    #   dC_zδ/dz = σ_δ² / γ²   (drift generates z-δ correlation via R_56)
    sigma_d = max(beam.sigma_delta, 1e-12)
    gamma2 = max(gamma**2, 1.001)
    dC_zd_dz = sigma_d**2 / gamma2

    #   dσ_δ/dz = 0 in drift, but RF chirp modifies it:
    #   In an RF cavity:  δ_out = δ_in + h·z
    #     σ_δ² → σ_δ² + 2h·C_zδ + h²·σ_z²
    #   In differential form: d(σ_δ²)/dz = 2h'·C_zδ + 2h'·h·σ_z²
    #     →  dσ_δ/dz = (h'·C_zδ + h'·h·σ_z²) / σ_δ   (for small h'·dz)
    if abs(h_rate) > 1e-15 and sigma_d > 1e-12:
        dsigma_delta_dz = h_rate * (beam.C_zd + h_rate * beam.sigma_z**2) / sigma_d
        dC_zd_dz += h_rate * beam.sigma_z**2       # RF chirp adds correlation
    else:
        dsigma_delta_dz = 0.0

    # ── Transverse envelope ──
    dnu_x = -damping_coeff * beam.nu_x + Fe_x + Fs_x + F_eps_x + dnx_sol
    dnu_y = -damping_coeff * beam.nu_y + Fe_y + Fs_y + F_eps_y + dny_sol
    dnu_z = -damping_coeff * beam.nu_z + Fe_z + Fs_z + F_eps_z
    dnu_xy = dsxy_sol
    dnu_delta = 0.0  # ν_δ tracked as independent variable

    return np.array([beam.nu_x, beam.nu_y, beam.nu_z, beam.nu_xy,
                     dsigma_delta_dz, dC_zd_dz,
                     dnu_x, dnu_y, dnu_z, dnu_xy, dnu_delta,
                     gamma_prime])


# ═════════════════════════════════════════════════════════════════════════════
# Propagation / integration
# ═════════════════════════════════════════════════════════════════════════════

def propagate(beam0: Beam6D, z_span: Tuple[float, float],
              n_points: int = 500,
              external_force_func: Optional[Callable] = None,
              gamma_prime_func: Optional[Callable] = None,
              sc_model: str = 'gaussian',
              rf_chirp_func: Optional[Callable] = None,
              **solver_kwargs) -> Tuple[np.ndarray, np.ndarray]:
    """
    Integrate the 7-variable envelope ODE over [z_start, z_end].

    Parameters
    ----------
    beam0 : Beam6D
        Initial beam state at z_start.
    z_span : (z_start, z_end)
        Integration interval (m).
    n_points : int
        Number of output points.
    external_force_func : callable, optional
        (z, beam) → (F_x^e, F_y^e, F_z^e).
    gamma_prime_func : callable, optional
        (z, beam) → dγ/dz.
    sc_model : str
        'gaussian' (Kelisani 2023, default) or 'ellipsoid' (Luiten 2004).
    **solver_kwargs
        Passed to solve_ivp (e.g., rtol, atol, method).

    Returns
    -------
    z_arr : ndarray (n_points,)
    beam_states : ndarray (n_points, 12)  — [σx,σy,σz,σxy,σδ,C_zd, νx,νy,νz,νxy,νδ, γ]
    """
    y0 = np.append(beam0.state_vector(), beam0.gamma)
    z_start, z_end = z_span
    t_eval = np.linspace(z_start, z_end, n_points)

    def rhs(z, y):
        return envelope_ode(z, y, beam0.Ne,
                            beam0.eps_nx, beam0.eps_ny, beam0.eps_nz,
                            external_force_func=external_force_func,
                            gamma_prime_func=gamma_prime_func,
                            rf_chirp_func=rf_chirp_func,
                            sc_model=sc_model)

    sol = solve_ivp(rhs, (z_start, z_end), y0, t_eval=t_eval,
                    method=solver_kwargs.pop('method', 'RK45'),
                    rtol=solver_kwargs.pop('rtol', 1e-8),
                    atol=solver_kwargs.pop('atol', 1e-10),
                    **solver_kwargs)

    return sol.t, sol.y.T


# ═════════════════════════════════════════════════════════════════════════════
# Convenience builders
# ═════════════════════════════════════════════════════════════════════════════

def make_beam_100keV(Ne: float = 1e5,
                     sigma_x0_um: float = 85.0,
                     sigma_y0_um: float = 85.0,
                     sigma_z0_um: float = 300.0,
                     sigma_delta: float = 1e-3,
                     C_zd_um: float = 0.0,
                     eps_nx_um: float = 0.03,
                     eps_ny_um: float = 0.03,
                     eps_nz_um: float = 0.2,
                     beamK_eV: float = 100_000.0,
                     divergence_nu_x: float = 0.0,
                     divergence_nu_y: float = 0.0,
                     divergence_nu_z: float = 0.0,
                     ) -> Beam6D:
    """
    Create a Beam6D state matching the UEM beamline parameters.

    Default: 100 keV, Ne=10⁵, σ_x=85 μm, σ_z=300 μm,
    ε_nx=0.03 μm, ε_nz=0.2 μm — matches `build_initial_beam()`.
    """
    gamma = beamK_to_gamma(beamK_eV)
    return Beam6D(
        sigma_x=sigma_x0_um * 1e-6,
        sigma_y=sigma_y0_um * 1e-6,
        sigma_z=sigma_z0_um * 1e-6,
        sigma_xy=0.0,
        sigma_delta=sigma_delta,
        C_zd=C_zd_um * 1e-6,
        nu_x=divergence_nu_x,
        nu_y=divergence_nu_y,
        nu_z=divergence_nu_z,
        nu_xy=0.0,
        nu_delta=0.0,
        gamma=gamma,
        Ne=Ne,
        eps_nx=eps_nx_um * 1e-6,
        eps_ny=eps_ny_um * 1e-6,
        eps_nz=eps_nz_um * 1e-6,
    )


# ═════════════════════════════════════════════════════════════════════════════
# Visualization
# ═════════════════════════════════════════════════════════════════════════════

def plot_evolution(z_arr: np.ndarray, beam_states: np.ndarray,
                   title: str = "6D Envelope Evolution",
                   savepath: str = None):
    """4-panel plot of beam sizes and slopes vs z."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)

    # Panel 1: RMS beam sizes
    ax = axes[0, 0]
    ax.plot(z_arr * 1e3, beam_states[:, 0] * 1e6, 'b-', label=r'$\sigma_x$', lw=1.5)
    ax.plot(z_arr * 1e3, beam_states[:, 1] * 1e6, 'r--', label=r'$\sigma_y$', lw=1.5)
    ax.set_ylabel(r'$\sigma_{x,y}$ (μm)')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Panel 2: Longitudinal size
    ax = axes[0, 1]
    ax.plot(z_arr * 1e3, beam_states[:, 2] * 1e6, 'g-', lw=1.5)
    ax.set_ylabel(r'$\sigma_z$ (μm)')
    ax.grid(True, alpha=0.3)

    # Panel 3: Slopes ν_x, ν_y
    ax = axes[1, 0]
    ax.plot(z_arr * 1e3, beam_states[:, 3], 'b-', label=r'$\nu_x$', lw=1.5)
    ax.plot(z_arr * 1e3, beam_states[:, 4], 'r--', label=r'$\nu_y$', lw=1.5)
    ax.set_xlabel('z (mm)')
    ax.set_ylabel(r'$\nu_{x,y}$  (slope)')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Panel 4: Slope ν_z
    ax = axes[1, 1]
    ax.plot(z_arr * 1e3, beam_states[:, 5], 'g-', lw=1.5)
    ax.set_xlabel('z (mm)')
    ax.set_ylabel(r'$\nu_z$  (slope)')
    ax.grid(True, alpha=0.3)

    fig.suptitle(title, fontsize=13, fontweight='bold')
    plt.tight_layout()
    if savepath:
        plt.savefig(savepath, dpi=150)
        print(f"Saved: {savepath}")
    return fig


# ═════════════════════════════════════════════════════════════════════════════
# Test / demo
# ═════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    print("=" * 60)
    print("  Phase 2: UED Longitudinal Phase Space (Discrete RF Lens)")
    print("=" * 60)

    # ── Beam at 100 keV, 0.1% energy spread ──
    beam0 = make_beam_100keV(Ne=1e5, beamK_eV=100_000,
                             sigma_z0_um=300, sigma_delta=0.001)
    print("\nInitial beam (post-acceleration):")
    print(beam0.summary("z=0"))
    print(f"  Δt = {beam0.time_spread_ps:.1f} ps")

    _ = get_alpha_interpolators()

    z_cav  = 0.400   # RF cavity center
    z_samp = 0.777   # sample plane
    L_cav  = 0.020   # cavity effective length

    # ── Step 1: propagate 0 → z_cav (drift) ──
    print(f"\n── Step 1: Drift 0 → {z_cav*1e3:.0f} mm ──")
    z1, st1 = propagate(beam0, (0.0, z_cav), n_points=150, sc_model='ellipsoid')
    b_before = Beam6D.from_state(st1[-1,:11], st1[-1,11], beam0.Ne,
                                  beam0.eps_nx, beam0.eps_ny, beam0.eps_nz)
    print(f"  σ_z={b_before.sigma_z*1e6:.1f}μm  σ_δ={b_before.sigma_delta*1e3:.2f}‰  "
          f"C_zδ={b_before.C_zd*1e6:.2f}μm  σ_x={b_before.sigma_x*1e6:.1f}μm")

    # ── Step 2: Apply RF thin lens ──
    print(f"\n── Step 2: RF cavity at z=400mm (E_rf=10MV/m, L={L_cav*1e3:.0f}mm) ──")
    H = compute_rf_chirp_coefficient(b_before, E_rf=10e6, L_cav=L_cav)
    print(f"  Chirp coefficient H = {H:.2f} m⁻¹")
    b_after = apply_rf_thin_lens(b_before, H)
    print(f"  σ_z={b_after.sigma_z*1e6:.1f}μm (unchanged)  "
          f"σ_δ={b_after.sigma_delta*1e3:.2f}‰  "
          f"C_zδ={b_after.C_zd*1e6:.2f}μm")

    # R_56 for drift cavity → sample
    R56 = compute_R56(z_samp - z_cav, b_after.gamma, b_after.beta)
    print(f"  R_56(cav→sample) = {R56*1e3:.1f} mm")

    # Optimal chirp for compression: H_opt = -1/R_56
    H_opt = -1.0 / R56 if abs(R56) > 1e-10 else 0.0
    print(f"  Optimal H for compression = {H_opt:.2f} m⁻¹  (actual H={H:.2f})")

    # ── Step 3: propagate cavity → sample ──
    print(f"\n── Step 3: Drift {z_cav*1e3:.0f} → {z_samp*1e3:.0f} mm ──")
    # Need to update the state array: keep σ_z,ν_z same, update σ_δ,C_zδ
    st1_end = st1[-1].copy()
    st1_end[4] = b_after.sigma_delta  # σ_δ
    st1_end[5] = b_after.C_zd         # C_zδ

    # Rebuild state at cavity exit with updated ν_z
    st1_end_new = st1_end.copy()
    st1_end_new[8] = b_after.nu_z  # update ν_z in state
    beam_cav_exit = Beam6D.from_state(st1_end_new[:11], st1_end_new[11], beam0.Ne,
                                       beam0.eps_nx, beam0.eps_ny, beam0.eps_nz)

    z2, st2 = propagate(beam_cav_exit, (z_cav, z_samp), n_points=150, sc_model='ellipsoid')
    b_sample = Beam6D.from_state(st2[-1,:11], st2[-1,11], beam0.Ne,
                                  beam0.eps_nx, beam0.eps_ny, beam0.eps_nz)

    # ── Also propagate without RF for comparison ──
    z_norf, st_norf = propagate(beam0, (0.0, z_samp), n_points=200, sc_model='ellipsoid')
    b_norf = Beam6D.from_state(st_norf[-1,:11], st_norf[-1,11], beam0.Ne,
                                beam0.eps_nx, beam0.eps_ny, beam0.eps_nz)

    # ── Results ──
    print(f"\n{'='*60}")
    print(f"  UED Longitudinal Compression Results")
    print(f"{'='*60}")
    for label, b in [("No RF (drift only)", b_norf), ("With RF chirp", b_sample)]:
        dt_rms = b.time_spread_ps * 1000
        print(f"  {label}:")
        print(f"    σ_z = {b.sigma_z*1e6:.1f} μm  |  Δt = {dt_rms:.0f} fs  |  "
              f"σ_δ = {b.sigma_delta*1e3:.2f}‰")

    dt_rf_rms = b_sample.time_spread_ps * 1000
    dt_rf_fwhm = dt_rf_rms * 2.355
    print(f"\n  Time resolution (with RF): Δt_RMS = {dt_rf_rms:.0f} fs, "
          f"Δt_FWHM = {dt_rf_fwhm:.0f} fs")

    # ── Plot comparison ──
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    z_rf = np.concatenate([z1, z2])
    st_rf = np.vstack([st1[:,:11], st2[:,:11]])
    # manual gamma array
    g_rf = np.concatenate([st1[:,11], st2[:,11]])

    z_mm_rf = z_rf * 1e3
    z_mm_norf = z_norf * 1e3

    axes[0,0].plot(z_mm_rf, st_rf[:,2]*1e6, 'b-', lw=1.5, label='With RF')
    axes[0,0].plot(z_mm_norf, st_norf[:,2]*1e6, 'gray', ls='--', lw=1, label='Drift only')
    axes[0,0].axvline(z_cav*1e3, color='red', ls=':', alpha=0.5)
    axes[0,0].set_ylabel('σ_z (μm)'); axes[0,0].legend(fontsize=8); axes[0,0].grid(alpha=0.3)
    axes[0,0].set_title('Bunch length evolution')

    axes[0,1].plot(z_mm_rf, st_rf[:,2]*1e6/(b_before.beta*C_LIGHT)*1e12, 'm-', lw=1.5)
    axes[0,1].axvline(z_cav*1e3, color='red', ls=':', alpha=0.5)
    axes[0,1].set_ylabel('Δt RMS (ps)'); axes[0,1].grid(alpha=0.3)
    axes[0,1].set_title('Time spread')

    axes[1,0].plot(z_mm_rf, st_rf[:,4]*1e3, 'r-', lw=1.5, label='σ_δ')
    axes[1,0].plot(z_mm_rf, st_rf[:,5]*1e6, 'c-', lw=1.5, label='C_zδ (μm)')
    axes[1,0].axvline(z_cav*1e3, color='red', ls=':', alpha=0.5)
    axes[1,0].set_xlabel('z (mm)'); axes[1,0].set_ylabel('σ_δ (‰) / C_zδ (μm)')
    axes[1,0].legend(fontsize=8); axes[1,0].grid(alpha=0.3)
    axes[1,0].set_title('Longitudinal phase space')

    axes[1,1].plot(z_mm_rf, st_rf[:,0]*1e6, 'b-', lw=1.2, label='σ_x')
    axes[1,1].plot(z_mm_rf, st_rf[:,1]*1e6, 'r--', lw=1.2, label='σ_y')
    axes[1,1].axvline(z_cav*1e3, color='red', ls=':', alpha=0.5)
    axes[1,1].set_xlabel('z (mm)'); axes[1,1].set_ylabel('σ (μm)')
    axes[1,1].legend(fontsize=8); axes[1,1].grid(alpha=0.3)
    axes[1,1].set_title('Transverse sizes')

    fig.suptitle("Phase 2: UED RF Compression — Discrete Thin-Lens Model", fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig("6d_phase2_longitudinal.png", dpi=150)
    print("\nSaved: 6d_phase2_longitudinal.png")
    print("\n=== Phase 2 Complete ===")
