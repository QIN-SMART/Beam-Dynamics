#!/usr/bin/env python3
"""
Phase 1.9 — Physical Solenoid Calibration

Upgrades solenoid benchmark from abstract k_sol to physical B-field input.
Compares three models: OCELOT particle tracking, 4×4 transfer matrix, thin lens.

Physics:  k_s = e·B / (2·p)  where p = γ·m_e·β·c

Usage:
  python3 benchmark_solenoid_physical.py
  python3 benchmark_solenoid_physical.py --B 0.05      # specific field
  python3 benchmark_solenoid_physical.py --B 0.10 --Ek 150  # 150 keV
  python3 benchmark_solenoid_physical.py --all          # sweep B = 0.02, 0.05, 0.10 T
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
#  physical constants
# ═══════════════════════════════════════════════════════════

e_SI   = 1.602176634e-19    # C
m_e_SI = 9.10938356e-31     # kg
c_SI   = 2.99792458e8       # m/s
mec2   = 511.0              # keV  (m_e c²)

# ═══════════════════════════════════════════════════════════
#  beam & lattice geometry
# ═══════════════════════════════════════════════════════════

spot_rms    = 85e-6           # m
sig_z0      = 300e-6          # m
E_keV       = 100.0           # keV
epsilon_n   = 0.08e-6         # m·rad  normalized emittance
sigma_delta = 1.0e-4          # ΔE/E

L_sol       = 0.06            # solenoid length [m]
L_drift1    = 0.10            # drift before solenoid
L_drift2    = 0.60            # drift after solenoid

# ═══════════════════════════════════════════════════════════
#  CLI  —  B-field input  (replaces --k)
# ═══════════════════════════════════════════════════════════

B_sol    = 0.05              # default [T]
run_all  = False

for i, a in enumerate(sys.argv):
    if a == "--B" and i + 1 < len(sys.argv):
        B_sol = float(sys.argv[i + 1])
    elif a == "--all":
        run_all = True
    elif a == "--Ek" and i + 1 < len(sys.argv):
        E_keV = float(sys.argv[i + 1])

# ═══════════════════════════════════════════════════════════
#  relativistic + solenoid physics
# ═══════════════════════════════════════════════════════════

gamma      = 1.0 + E_keV / mec2
beta       = np.sqrt(1.0 - 1.0 / gamma**2)
beta_gamma = beta * gamma
p_SI       = gamma * m_e_SI * beta * c_SI                # kg·m/s

# solenoid Larmor focusing strength  k = eB / (2p)
def compute_k(B):
    return e_SI * B / (2.0 * p_SI)

# geometric emittance from normalized
epsilon_geom = epsilon_n / beta_gamma
sigma_xp = epsilon_geom / spot_rms
sigma_yp = epsilon_geom / spot_rms

# ═══════════════════════════════════════════════════════════
#  transfer matrix functions
# ═══════════════════════════════════════════════════════════

def mat_drift_4x4(s):
    return np.array([[1.0,  s,  0.0, 0.0],
                     [0.0, 1.0, 0.0, 0.0],
                     [0.0, 0.0, 1.0,  s ],
                     [0.0, 0.0, 0.0, 1.0]])

def mat_solenoid_4x4(s, k):
    C = np.cos(k * s)
    S = np.sin(k * s)
    return np.array([
        [C*C,        S*C/k,   S*C,     S*S/k],
        [-k*S*C,     C*C,    -k*S*S,    S*C ],
        [-S*C,      -S*S/k,   C*C,     S*C/k],
        [ k*S*S,    -S*C,    -k*S*C,    C*C ],
    ])

def mat_thin_lens(f):
    """Thin-lens 2×2 matrix (focusing)."""
    return np.array([[1.0,  0.0],
                     [-1.0/f, 1.0]])

# ═══════════════════════════════════════════════════════════
#  analytical envelope functions
# ═══════════════════════════════════════════════════════════

def envelope_4x4(z_arr, k, sx0, sxp0, sy0, syp0):
    """σ_x(z) via 4×4 transfer matrix propagation."""
    S0 = np.diag([sx0**2, sxp0**2, sy0**2, syp0**2])
    sx = np.zeros_like(z_arr)
    sy = np.zeros_like(z_arr)
    for i, z in enumerate(z_arr):
        if z <= L_drift1:
            M = mat_drift_4x4(z)
        elif z <= L_drift1 + L_sol:
            M = mat_solenoid_4x4(z - L_drift1, k) @ mat_drift_4x4(L_drift1)
        else:
            M = (mat_drift_4x4(z - L_drift1 - L_sol) @
                 mat_solenoid_4x4(L_sol, k) @
                 mat_drift_4x4(L_drift1))
        Sz = M @ S0 @ M.T
        sx[i] = np.sqrt(max(Sz[0, 0], 0.0))
        sy[i] = np.sqrt(max(Sz[2, 2], 0.0))
    return sx, sy

def envelope_thin_lens(z_arr, f, sx0, sxp0):
    """σ_x(z) via thin-lens model (2×2, uncoupled)."""
    z_c = L_drift1 + L_sol / 2.0       # lens at solenoid centre
    M_lens = mat_thin_lens(f)
    sx = np.zeros_like(z_arr)
    for i, z in enumerate(z_arr):
        if z <= z_c:
            M = np.array([[1.0, z], [0.0, 1.0]])
        else:
            M = (np.array([[1.0, z - z_c], [0.0, 1.0]]) @
                 M_lens @
                 np.array([[1.0, z_c], [0.0, 1.0]]))
        sx[i] = np.sqrt(M[0, 0]**2 * sx0**2 + M[0, 1]**2 * sxp0**2)
    return sx

# ═══════════════════════════════════════════════════════════
#  single-field benchmark
# ═══════════════════════════════════════════════════════════

def run_benchmark(B, make_plots=True, label=""):
    k_s = compute_k(B)
    f_thin = 1.0 / (k_s**2 * L_sol) if k_s > 0 else np.inf

    # ——— lattice ———
    lat = MagneticLattice([
        Drift(l=L_drift1, eid="D1"),
        Solenoid(l=L_sol, k=k_s, eid="SOL"),
        Drift(l=L_drift2, eid="D2"),
    ])
    lat.update_transfer_maps()

    # ——— beam ———
    p = generate_parray(
        sigma_x=spot_rms, sigma_y=spot_rms,
        sigma_tau=sig_z0 / beta,          # OCELOT tau = c·t [m]; σ_tau = σ_z/β
        energy=(E_keV + 511.0) * 1e-6,    # TOTAL energy in GeV (E_kin+mc²)
        charge=100_000,
    )
    np.random.seed(42)
    N = p.rparticles.shape[1]
    p.rparticles[1, :] = np.random.normal(0.0, sigma_xp, N)
    p.rparticles[3, :] = np.random.normal(0.0, sigma_xp, N)
    p.rparticles[5, :] = np.random.normal(0.0, sigma_delta, N)

    # initial stats
    x0  = p.x().copy(); xp0 = p.px().copy()
    y0  = p.y().copy(); yp0 = p.py().copy()
    sx0 = np.std(x0); sxp0 = np.std(xp0)
    sy0 = np.std(y0); syp0 = np.std(yp0)
    eps_x0 = np.sqrt(np.mean(x0**2)*np.mean(xp0**2) - np.mean(x0*xp0)**2)

    # ——— OCELOT tracking ———
    navi = Navigator(lat)
    dz   = 0.001
    total_length = L_drift1 + L_sol + L_drift2
    n_steps = int(total_length / dz)

    z_list    = []
    sigx_list = []
    sigy_list = []
    epsx_list = []

    for _ in range(n_steps):
        tracking_step(lat, p, dz, navi)
        z = navi.z0
        x  = p.x(); xp = p.px(); y = p.y(); yp = p.py()
        z_list.append(z)
        sigx_list.append(np.std(x))
        sigy_list.append(np.std(y))
        epsx_list.append(np.sqrt(np.mean(x**2)*np.mean(xp**2) - np.mean(x*xp)**2))

    z_arr    = np.array(z_list)
    sigx_oc  = np.array(sigx_list)
    sigy_oc  = np.array(sigy_list)
    epsx_oc  = np.array(epsx_list)

    # ——— analytical envelopes ———
    sigx_4x4, sigy_4x4 = envelope_4x4(z_arr, k_s, sx0, sxp0, sy0, syp0)
    sigx_tl = envelope_thin_lens(z_arr, f_thin, sx0, sxp0)

    # ——— waist statistics ———
    idx_w_oc  = np.argmin(sigx_oc)
    z_w_oc    = z_arr[idx_w_oc]
    sx_w_oc   = sigx_oc[idx_w_oc]

    idx_w_4x4 = np.argmin(sigx_4x4)
    z_w_4x4   = z_arr[idx_w_4x4]
    sx_w_4x4  = sigx_4x4[idx_w_4x4]

    idx_w_tl  = np.argmin(sigx_tl)
    z_w_tl    = z_arr[idx_w_tl]
    sx_w_tl   = sigx_tl[idx_w_tl]

    # ——— envelope errors ———
    err_4x4_x = np.max(np.abs(sigx_oc - sigx_4x4) / sigx_4x4) * 100
    err_tl_x  = np.max(np.abs(sigx_oc - sigx_tl) / sigx_tl) * 100

    # ——— print summary ———
    print(f"\n{'='*70}")
    print(f"  B = {B:.3f} T    k_s = {k_s:.4f}    thin-lens f = {f_thin*1e3:.1f} mm")
    print(f"  E_k = {E_keV:.0f} keV  γ = {gamma:.4f}  β = {beta:.4f}  p = {p_SI:.2e} kg·m/s")
    print(f"{'='*70}")
    print(f"  {'Model':<18s} {'z_waist (mm)':>14s} {'σ_x,waist (μm)':>15s} {'σ_x,end (μm)':>14s}")
    print(f"  {'-'*60}")
    print(f"  {'OCELOT tracking':<18s} {z_w_oc*1e3:14.1f} {sx_w_oc*1e6:15.1f} {sigx_oc[-1]*1e6:14.1f}")
    print(f"  {'4×4 matrix':<18s} {z_w_4x4*1e3:14.1f} {sx_w_4x4*1e6:15.1f} {sigx_4x4[-1]*1e6:14.1f}")
    print(f"  {'thin lens approx':<18s} {z_w_tl*1e3:14.1f} {sx_w_tl*1e6:15.1f} {sigx_tl[-1]*1e6:14.1f}")
    print(f"\n  ε_x: initial={eps_x0*1e6:.4f}  final={epsx_oc[-1]*1e6:.4f} mm·mrad")
    print(f"  Envelope error (4×4): {err_4x4_x:.3f} %")
    print(f"  Envelope error (thin lens): {err_tl_x:.2f} %")

    # ——— plots ———
    if make_plots:
        tag = label if label else f"B{B*1e3:.0f}mT"
        _make_plots(z_arr, sigx_oc, sigy_oc, sigx_4x4, sigy_4x4, sigx_tl,
                    epsx_oc, eps_x0, z_w_oc, k_s, B, f_thin, tag)

    return {
        "B": B, "k_s": k_s, "f_thin": f_thin,
        "z_w_oc": z_w_oc, "sx_w_oc": sx_w_oc,
        "z_w_4x4": z_w_4x4, "sx_w_4x4": sx_w_4x4,
        "z_w_tl": z_w_tl, "sx_w_tl": sx_w_tl,
        "err_4x4": err_4x4_x, "err_tl": err_tl_x,
        "sigx_end": sigx_oc[-1], "epsx_init": eps_x0, "epsx_final": epsx_oc[-1],
    }

# ═══════════════════════════════════════════════════════════
#  plotting
# ═══════════════════════════════════════════════════════════

def _make_plots(z_arr, sigx_oc, sigy_oc, sigx_4x4, sigy_4x4, sigx_tl,
                epsx_oc, eps_x0, z_w_oc, k_s, B, f_thin, tag):

    fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig1.suptitle(f"Physical Solenoid — B = {B:.3f} T  (k = {k_s:.2f})",
                  fontsize=12, fontweight="bold")

    # — top: σ_x(z)  three models —
    ax1.plot(z_arr * 1e3, sigx_oc * 1e6, "b.", markersize=1.5, alpha=0.6,
             label="OCELOT tracking  $\\sigma_x$")
    ax1.plot(z_arr * 1e3, sigy_oc * 1e6, "c.", markersize=1.5, alpha=0.6,
             label="OCELOT tracking  $\\sigma_y$")
    ax1.plot(z_arr * 1e3, sigx_4x4 * 1e6, "r-", linewidth=1.5,
             label="4×4 transfer matrix  $\\sigma_x$")
    ax1.plot(z_arr * 1e3, sigx_tl * 1e6, "g--", linewidth=1.2,
             label=f"thin lens  (f = {f_thin*1e3:.0f} mm)")

    ax1.axvspan(L_drift1 * 1e3, (L_drift1 + L_sol) * 1e3,
                alpha=0.10, color="orange")
    ax1.axvline(z_w_oc * 1e3, color="green", linestyle=":", alpha=0.5,
                label=f"OCELOT waist @ {z_w_oc*1e3:.0f} mm")
    ax1.set_ylabel(r"$\sigma_{x,y}$  [$\mu$m]")
    ax1.legend(fontsize=7, loc="upper right", ncol=2)
    ax1.grid(True, alpha=0.25)

    # — bottom: ε_x(z) —
    ax2.plot(z_arr * 1e3, epsx_oc * 1e6, "b.", markersize=1.5, alpha=0.6)
    ax2.axhline(eps_x0 * 1e6, color="red", linestyle="--", linewidth=1,
                label=f"$\\varepsilon_x$(0) = {eps_x0*1e6:.4f} mm·mrad")
    ax2.axvspan(L_drift1 * 1e3, (L_drift1 + L_sol) * 1e3,
                alpha=0.10, color="orange")
    ax2.set_xlabel("z  [mm]")
    ax2.set_ylabel(r"$\varepsilon_x$  [mm$\cdot$mrad]")
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.25)

    fig1.tight_layout()
    fname = os.path.join(_OUTDIR, f"benchmark_solenoid_physical_{tag}.png")
    fig1.savefig(fname, dpi=150)
    print(f"  -> {fname}")
    plt.close(fig1)

# ═══════════════════════════════════════════════════════════
#  main
# ═══════════════════════════════════════════════════════════

if __name__ == "__main__":
    print(f"\n  Physical constants:")
    print(f"    e = {e_SI:.4e} C")
    print(f"    m_e = {m_e_SI:.4e} kg")
    print(f"    c = {c_SI:.4e} m/s")
    print(f"    m_e c² = {mec2:.0f} keV")

    if run_all:
        B_list = [0.02, 0.05, 0.10]
        results = []
        for B in B_list:
            res = run_benchmark(B, make_plots=True, label=f"B{B*1e3:.0f}mT")
            results.append(res)

        # —— sweep summary ——
        print(f"\n{'='*70}")
        print(f"  Sweep Summary")
        print(f"{'='*70}")
        print(f"  {'B (T)':>8s}  {'k_s':>8s}  {'f_thin (mm)':>13s}  "
              f"{'z_w OC (mm)':>13s}  {'z_w 4×4 (mm)':>13s}  {'z_w TL (mm)':>13s}  "
              f"{'err TL%':>8s}")
        for r in results:
            print(f"  {r['B']:8.3f}  {r['k_s']:8.2f}  {r['f_thin']*1e3:13.1f}  "
                  f"{r['z_w_oc']*1e3:13.1f}  {r['z_w_4x4']*1e3:13.1f}  "
                  f"{r['z_w_tl']*1e3:13.1f}  {r['err_tl']:7.1f}%")

        # —— z_waist vs B plot ——
        fig_sweep, ax = plt.subplots(figsize=(7, 4.5))
        B_arr = np.array([r["B"] for r in results])
        ax.plot(B_arr * 1e3, np.array([r["z_w_oc"]*1e3 for r in results]), "bo-", label="OCELOT")
        ax.plot(B_arr * 1e3, np.array([r["z_w_4x4"]*1e3 for r in results]), "rs--", label="4×4 matrix")
        ax.plot(B_arr * 1e3, np.array([r["z_w_tl"]*1e3 for r in results]), "g^:", label="thin lens")
        ax.set_xlabel("B  [mT]")
        ax.set_ylabel("beam waist position  [mm]")
        ax.set_title("Waist position vs magnetic field")
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig_sweep.tight_layout()
        fig_sweep.savefig(os.path.join(_OUTDIR, "benchmark_solenoid_physical_sweep.png"), dpi=150)
        print("  -> benchmark_solenoid_physical_sweep.png")
        plt.close(fig_sweep)

    else:
        run_benchmark(B_sol, make_plots=True, label=f"B{B_sol*1e3:.0f}mT")

    print("\n  Benchmark complete.")
