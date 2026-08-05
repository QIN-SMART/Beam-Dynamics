#!/usr/bin/env python3
"""
Phase 1.8 — Solenoid-only benchmark validation.

Validates:
  1.  Solenoid focusing produces a beam waist
  2.  Envelope matches analytical transfer-matrix prediction
  3.  x–y symmetry
  4.  Emittance conservation

No RF, no space charge.

Usage:
  python3 benchmark_solenoid.py
  python3 benchmark_solenoid.py --k 20      # stronger solenoid (visible waist)
  python3 benchmark_solenoid.py --epsn 0.04e-6
"""

import sys, os, numpy as np

print("加载 OCELOT …", flush=True)
import ocelot
from ocelot.cpbd.elements import Drift, Solenoid
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.beam import generate_parray
from ocelot.cpbd.navi import Navigator
from ocelot.cpbd.track import tracking_step

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_OUTDIR = os.path.dirname(os.path.abspath(__file__)) + "/"

# ═══════════════════════════════════════════════════════════
#  parameters
# ═══════════════════════════════════════════════════════════

spot_rms    = 85e-6
sig_z0      = 300e-6
E_keV       = 100.0
epsilon_n   = 0.08e-6
sigma_delta = 1.0e-4
k_sol       = 1.5          # solenoid strength (default: TL1)
L_sol       = 0.06         # solenoid length [m]
L_drift1    = 0.10         # drift before solenoid
L_drift2    = 0.60         # drift after solenoid

for i, a in enumerate(sys.argv):
    if a == "--k" and i + 1 < len(sys.argv):
        k_sol = float(sys.argv[i + 1])
    elif a == "--epsn" and i + 1 < len(sys.argv):
        epsilon_n = float(sys.argv[i + 1])

# relativistic
E_rest      = 511.0
gamma       = 1.0 + E_keV / E_rest
beta        = np.sqrt(1.0 - 1.0 / gamma**2)
beta_gamma  = beta * gamma
epsilon_geom = epsilon_n / beta_gamma
sigma_xp     = epsilon_geom / spot_rms

# ═══════════════════════════════════════════════════════════
#  lattice  —  Drift + Solenoid + Drift
# ═══════════════════════════════════════════════════════════

lat = MagneticLattice([
    Drift(l=L_drift1, eid="D1"),
    Solenoid(l=L_sol, k=k_sol, eid="SOL"),
    Drift(l=L_drift2, eid="D2"),
])
lat.update_transfer_maps()

# ═══════════════════════════════════════════════════════════
#  generate beam  (Phase‑1 distribution)
# ═══════════════════════════════════════════════════════════

p = generate_parray(
    sigma_x=spot_rms, sigma_y=spot_rms,
    sigma_tau=sig_z0 / beta,          # OCELOT tau = c·t [m]; σ_tau = σ_z/β
    energy=(E_keV + 511.0) * 1e-6,    # TOTAL energy in GeV (E_kin+mc²)
    charge=100_000,
)

np.random.seed(42)
N = p.rparticles.shape[1]
p.rparticles[1, :] = np.random.normal(0.0, sigma_xp, N)   # px = x'
p.rparticles[3, :] = np.random.normal(0.0, sigma_xp, N)   # py = y'
p.rparticles[5, :] = np.random.normal(0.0, sigma_delta, N)  # delta

# record initial stats
x0  = p.x().copy()
xp0 = p.px().copy()
y0  = p.y().copy()
yp0 = p.py().copy()
sigma_x0  = np.std(x0)
sigma_xp0 = np.std(xp0)
sigma_y0  = np.std(y0)
sigma_yp0 = np.std(yp0)
eps_x0    = np.sqrt(np.mean(x0**2) * np.mean(xp0**2) - np.mean(x0 * xp0)**2)
print(f"\n  Initial: σ_x={sigma_x0*1e6:.2f} μm  σ_y={sigma_y0*1e6:.2f} μm  "
      f"σ_x'={sigma_xp0*1e3:.3f} mrad  "
      f"ε_x={eps_x0*1e6:.4f} mm·mrad  k_sol={k_sol}")

# ═══════════════════════════════════════════════════════════
#  tracking  —  record at fine z resolution
# ═══════════════════════════════════════════════════════════

navi = Navigator(lat)
dz   = 0.001
total_length = L_drift1 + L_sol + L_drift2
n_steps = int(total_length / dz)

z_hist    = []
sigx_hist = []
sigy_hist = []
epsx_hist = []

for step_i in range(n_steps):
    tracking_step(lat, p, dz, navi)
    z = navi.z0
    x  = p.x()
    xp = p.px()
    y  = p.y()
    yp = p.py()

    z_hist.append(z)
    sigx_hist.append(np.std(x))
    sigy_hist.append(np.std(y))
    epsx_hist.append(np.sqrt(np.mean(x**2)*np.mean(xp**2) - np.mean(x*xp)**2))

z_arr    = np.array(z_hist)
sigx_arr = np.array(sigx_hist)
sigy_arr = np.array(sigy_hist)
epsx_arr = np.array(epsx_hist)

# final phase space
xf  = p.x().copy()
xpf = p.px().copy()

# ═══════════════════════════════════════════════════════════
#  analytical envelope from 4×4 transfer matrices
# ═══════════════════════════════════════════════════════════

def mat_drift_4x4(s):
    """4×4 transfer matrix for drift of length s (x, x', y, y')."""
    return np.array([[1.0,  s,  0.0, 0.0],
                     [0.0, 1.0, 0.0, 0.0],
                     [0.0, 0.0, 1.0,  s ],
                     [0.0, 0.0, 0.0, 1.0]])

def mat_solenoid_4x4(s, k):
    """4×4 transfer matrix for solenoid of length s, strength k."""
    C = np.cos(k * s)
    S = np.sin(k * s)
    return np.array([
        [C*C,        S*C/k,   S*C,     S*S/k],
        [-k*S*C,     C*C,    -k*S*S,    S*C ],
        [-S*C,      -S*S/k,   C*C,     S*C/k],
        [ k*S*S,    -S*C,    -k*S*C,    C*C ],
    ])

def cumulative_mat_4x4(z):
    """Cumulative 4×4 transfer matrix from z=0 to position z."""
    if z <= L_drift1:
        return mat_drift_4x4(z)
    elif z <= L_drift1 + L_sol:
        s_in_sol = z - L_drift1
        return mat_solenoid_4x4(s_in_sol, k_sol) @ mat_drift_4x4(L_drift1)
    else:
        s_after = z - L_drift1 - L_sol
        return (mat_drift_4x4(s_after) @
                mat_solenoid_4x4(L_sol, k_sol) @
                mat_drift_4x4(L_drift1))

# initial 4×4 beam matrix Σ = diag(σ_x², σ_x'², σ_y², σ_y'²)
S0 = np.diag([sigma_x0**2, sigma_xp0**2, sigma_y0**2, sigma_yp0**2])

sigx_analytic = np.zeros_like(z_arr)
sigy_analytic = np.zeros_like(z_arr)
for i, z in enumerate(z_arr):
    M = cumulative_mat_4x4(z)
    Sz = M @ S0 @ M.T
    sigx_analytic[i] = np.sqrt(max(Sz[0, 0], 0.0))
    sigy_analytic[i] = np.sqrt(max(Sz[2, 2], 0.0))

# locate beam waist in OCELOT data
idx_waist = np.argmin(sigx_arr)
z_waist   = z_arr[idx_waist]
sig_waist = sigx_arr[idx_waist]

# thin-lens focal length: 1/f = k² L_sol (paraxial)
f_thin = 1.0 / (k_sol**2 * L_sol) if k_sol > 0 else np.inf
z_waist_thin = L_drift1 + L_sol + f_thin * (L_drift1) / (f_thin - L_drift1) if f_thin != np.inf else np.inf

# ═══════════════════════════════════════════════════════════
#  summary
# ═══════════════════════════════════════════════════════════

print(f"\n  Solenoid k={k_sol}, L={L_sol} m, thin-lens f ≈ {f_thin*1e3:.1f} mm")
print(f"  Beam waist:")
print(f"    z_waist = {z_waist*1e3:.1f} mm  (OCELOT tracking)")
print(f"    σ_x at waist = {sig_waist*1e6:.1f} μm")
print(f"    σ_x at end   = {sigx_arr[-1]*1e6:.1f} μm")
print(f"  Emittance: ε_x(0)={eps_x0*1e6:.4f}, ε_x(end)={epsx_arr[-1]*1e6:.4f} mm·mrad")
print(f"  x-y symmetry: max|σ_x-σ_y| = {np.max(np.abs(sigx_arr - sigy_arr))*1e6:.3f} μm")

max_env_err_x = np.max(np.abs(sigx_arr - sigx_analytic) / sigx_analytic) * 100
max_env_err_y = np.max(np.abs(sigy_arr - sigy_analytic) / sigy_analytic) * 100
max_err_abs_x = np.max(np.abs(sigx_arr - sigx_analytic)) * 1e6
print(f"  Envelope error σ_x (max relative): {max_env_err_x:.3f} %")
print(f"  Envelope error σ_y (max relative): {max_env_err_y:.3f} %")
print(f"  Max absolute error σ_x: {max_err_abs_x:.2f} μm  (≈ waist × {max_err_abs_x/(sig_waist*1e6)*100:.1f}%)")

# ═══════════════════════════════════════════════════════════
#  figure 1 — σ_x(z)  OCELOT vs analytic
# ═══════════════════════════════════════════════════════════

fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
fig1.suptitle("Solenoid Benchmark — Phase 1.8", fontsize=13, fontweight="bold")

ax1.plot(z_arr * 1e3, sigx_arr * 1e6, "b.", markersize=1.5, alpha=0.7,
         label="OCELOT  $\\sigma_x$")
ax1.plot(z_arr * 1e3, sigy_arr * 1e6, "c.", markersize=1.5, alpha=0.7,
         label="OCELOT  $\\sigma_y$")
ax1.plot(z_arr * 1e3, sigx_analytic * 1e6, "r-", linewidth=1.5,
         label="analytic envelope $\\sigma_x$")
ax1.plot(z_arr * 1e3, sigy_analytic * 1e6, "orange", linewidth=1.0, linestyle="--",
         label="analytic envelope $\\sigma_y$")
ax1.axvspan(L_drift1 * 1e3, (L_drift1 + L_sol) * 1e3,
            alpha=0.12, color="orange", label="solenoid")
ax1.axvline(z_waist * 1e3, color="green", linestyle="--", alpha=0.6,
            label=f"waist @ {z_waist*1e3:.0f} mm")
ax1.set_ylabel(r"$\sigma_{x,y}$  [$\mu$m]")
ax1.legend(fontsize=8, loc="upper right")
ax1.grid(True, alpha=0.3)

# bottom subplot: ε_x(z)
ax2.plot(z_arr * 1e3, epsx_arr * 1e6, "b.", markersize=1.5)
ax2.axhline(eps_x0 * 1e6, color="red", linestyle="--", linewidth=1,
            label=f"initial $\\varepsilon_x$ = {eps_x0*1e6:.4f} mm·mrad")
ax2.axvspan(L_drift1 * 1e3, (L_drift1 + L_sol) * 1e3,
            alpha=0.12, color="orange")
ax2.set_xlabel("z  [mm]")
ax2.set_ylabel(r"$\varepsilon_x$  [mm$\cdot$mrad]")
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.3)

fig1.tight_layout()
fig1.savefig(os.path.join(_OUTDIR, "benchmark_solenoid_evolution.png"), dpi=150)
print("  -> benchmark_solenoid_evolution.png")

# ═══════════════════════════════════════════════════════════
#  figure 2 — x–x′ phase space  before / after solenoid
# ═══════════════════════════════════════════════════════════

n_plot = min(8000, N)
rng = np.random.default_rng(42)
idx = rng.choice(N, n_plot, replace=False)

fig2, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(11, 4.8))
fig2.suptitle("x–x′ Phase Space — Solenoid Benchmark", fontsize=12, fontweight="bold")

eps_xf = np.sqrt(np.mean(xf**2) * np.mean(xpf**2) - np.mean(xf * xpf)**2)

ax_a.scatter(x0[idx] * 1e6, xp0[idx] * 1e3, s=0.5, alpha=0.4, c="steelblue")
ax_a.set_xlabel(r"$x$  [$\mu$m]")
ax_a.set_ylabel(r"$x'$  [mrad]")
ax_a.set_title(f"z = 0   (before solenoid)\n"
               f"$\\sigma_x$={sigma_x0*1e6:.1f} μm,  "
               f"$\\varepsilon_x$={eps_x0*1e6:.4f} mm·mrad")
ax_a.grid(True, alpha=0.3)

ax_b.scatter(xf[idx] * 1e6, xpf[idx] * 1e3, s=0.5, alpha=0.4, c="darkorange")
ax_b.set_xlabel(r"$x$  [$\mu$m]")
ax_b.set_ylabel(r"$x'$  [mrad]")
ax_b.set_title(f"z = {total_length:.1f} m  (after solenoid + drift)\n"
               f"$\\sigma_x$={np.std(xf)*1e6:.1f} μm,  "
               f"$\\varepsilon_x$={eps_xf*1e6:.4f} mm·mrad")
ax_b.grid(True, alpha=0.3)

fig2.tight_layout()
fig2.savefig(os.path.join(_OUTDIR, "benchmark_solenoid_phase_space.png"), dpi=150)
print("  -> benchmark_solenoid_phase_space.png")

print("\n  Benchmark complete.")
