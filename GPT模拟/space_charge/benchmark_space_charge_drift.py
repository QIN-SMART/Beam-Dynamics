#!/usr/bin/env python3
"""
Phase 2A.5 — Space Charge Sensitivity Validation

Extends Phase 2A with parameter scans to verify SC transverse effect:
  A. Charge scan:  10 fC → 1 pC  (higher Q amplifies SC)
  B. Emittance scan: ε_n = 0.001 → 0.08 mm·mrad  (lower ε reveals SC)

SC algorithm (OCELOT 26.06.1):
  3D PIC, mesh [63,63,63], NGP deposition, spectral Poisson solver.

Usage:
  python3 benchmark_space_charge_drift.py                  # default
  python3 benchmark_space_charge_drift.py --Q 500e-15      # 500 fC
  python3 benchmark_space_charge_drift.py --epsn 0.001e-6  # low emittance
  python3 benchmark_space_charge_drift.py --epsn-scan      # emittance sweep
  python3 benchmark_space_charge_drift.py --all            # full combined scan
"""

import sys, os, numpy as np

print("加载 OCELOT …", flush=True)
import ocelot
from ocelot.cpbd.elements import Drift
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.beam import generate_parray
from ocelot.cpbd.navi import Navigator
from ocelot.cpbd.track import tracking_step
from ocelot.cpbd.sc import SpaceCharge

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_OUTDIR = os.path.dirname(os.path.abspath(__file__)) + "/"

# ═══════════════════════════════════════════════════════════
#  beam parameters
# ═══════════════════════════════════════════════════════════

spot_rms    = 85e-6           # m
sig_z0      = 300e-6          # m
E_keV       = 100.0           # keV
epsilon_n   = 0.08e-6         # m·rad  normalized emittance (default)
sigma_delta = 1.0e-4          # ΔE/E  relative energy spread

Q_bunch     = 100e-15         # C  (100 fC default)

drift_length = 1.0            # m

# relativistic
E_rest      = 511.0           # keV
gamma       = 1.0 + E_keV / E_rest
beta        = np.sqrt(1.0 - 1.0 / gamma**2)
beta_gamma  = beta * gamma

# CLI
run_all    = False
run_escan  = False
for i, a in enumerate(sys.argv):
    if a == "--Q" and i + 1 < len(sys.argv):
        Q_bunch = float(sys.argv[i + 1])
    elif a == "--epsn" and i + 1 < len(sys.argv):
        epsilon_n = float(sys.argv[i + 1])
    elif a == "--epsn-scan":
        run_escan = True
    elif a == "--all":
        run_all = True

# ═══════════════════════════════════════════════════════════
#  beam generation
# ═══════════════════════════════════════════════════════════

def generate_beam(Q_C, eps_n=None):
    """Generate ParticleArray with total charge Q_C [C] and normalized emittance eps_n [m·rad]."""
    if eps_n is None:
        eps_n = epsilon_n
    eps_geom = eps_n / beta_gamma
    sxp = eps_geom / spot_rms
    syp = eps_geom / spot_rms

    p = generate_parray(
        sigma_x=spot_rms, sigma_y=spot_rms,
        sigma_tau=sig_z0 / beta,          # OCELOT tau = c·t [m]; σ_tau = σ_z/β
        energy=(E_keV + 511.0) * 1e-6,    # TOTAL energy in GeV (E_kin+mc²)
        charge=Q_C,
        nparticles=50000,
    )
    np.random.seed(42)
    N = p.rparticles.shape[1]
    p.rparticles[1, :] = np.random.normal(0.0, sxp, N)
    p.rparticles[3, :] = np.random.normal(0.0, syp, N)
    p.rparticles[5, :] = np.random.normal(0.0, sigma_delta, N)
    return p

# ═══════════════════════════════════════════════════════════
#  diagnostics
# ═══════════════════════════════════════════════════════════

def emit(x, xp):
    return np.sqrt(np.mean(x**2) * np.mean(xp**2) - np.mean(x * xp)**2)

# ═══════════════════════════════════════════════════════════
#  tracking
# ═══════════════════════════════════════════════════════════

def run_sc_sim(Q_C, eps_n=None, sc_on=True):
    """
    Track beam through 1.0 m drift.
    Returns (z_arr, sigx, sigy, sigt, epsx, epsy, sigd,
             x0, xp0, y0, yp0, t0, d0,
             xf, xpf, yf, ypf, tf, df,
             sx0, sxp0, sy0, syp0, eps_x0, eps_y0)
    """
    lat = MagneticLattice([Drift(l=drift_length, eid="DRIFT")])
    lat.update_transfer_maps()

    p = generate_beam(Q_C, eps_n)

    x0  = p.x().copy();   xp0 = p.px().copy()
    y0  = p.y().copy();   yp0 = p.py().copy()
    t0  = p.tau().copy(); d0  = p.p().copy()

    sx0  = np.std(x0);  sxp0 = np.std(xp0)
    sy0  = np.std(y0);  syp0 = np.std(yp0)
    eps_x0 = emit(x0, xp0);  eps_y0 = emit(y0, yp0)

    sc = SpaceCharge(step=1) if sc_on else None

    navi = Navigator(lat)
    dz   = 0.01
    n_steps = int(drift_length / dz)

    z_list, sx_list, sy_list, st_list = [], [], [], []
    ex_list, ey_list, sd_list = [], [], []

    for _ in range(n_steps):
        tracking_step(lat, p, dz, navi)
        if sc is not None:
            sc.apply(p, dz)

        z = navi.z0
        x = p.x(); xp = p.px(); y = p.y(); yp = p.py()
        ta = p.tau(); dd = p.p()

        z_list.append(z)
        sx_list.append(np.std(x))
        sy_list.append(np.std(y))
        st_list.append(np.std(ta))
        ex_list.append(emit(x, xp))
        ey_list.append(emit(y, yp))
        sd_list.append(np.std(dd))

    xf = p.x().copy(); xpf = p.px().copy()
    yf = p.y().copy(); ypf = p.py().copy()
    tf = p.tau().copy(); df = p.p().copy()

    return (np.array(z_list),
            np.array(sx_list), np.array(sy_list), np.array(st_list),
            np.array(ex_list), np.array(ey_list), np.array(sd_list),
            x0, xp0, y0, yp0, t0, d0,
            xf, xpf, yf, ypf, tf, df,
            sx0, sxp0, sy0, syp0, eps_x0, eps_y0)

# ═══════════════════════════════════════════════════════════
#  single-point comparison  (SC OFF vs ON)
# ═══════════════════════════════════════════════════════════

def compare_at_charge(Q_C, eps_n=None, make_plots=False):
    eps_use = eps_n if eps_n is not None else epsilon_n
    tag = f"Q{Q_C*1e15:.0f}fC_eps{eps_use*1e6:.3f}"

    print(f"  Q={Q_C*1e15:.0f}fC  ε_n={eps_use*1e6:.3f} mm·mrad  SC OFF …", end="", flush=True)
    r_off = run_sc_sim(Q_C, eps_n, sc_on=False)
    print(" done | SC ON …", end="", flush=True)
    r_on  = run_sc_sim(Q_C, eps_n, sc_on=True)
    print(" done")

    (z_off, sx_off, sy_off, st_off, ex_off, ey_off, sd_off,
     _, _, _, _, _, _, xf_off, xpf_off, _, _, tf_off, df_off,
     sx0, sxp0, _, _, eps_x0, _) = r_off

    (z_on,  sx_on,  sy_on,  st_on,  ex_on,  ey_on,  sd_on,
     _, _, _, _, _, _, xf_on, xpf_on, _, _, tf_on, df_on,
     _, _, _, _, _, _) = r_on

    dsx = (sx_on[-1] / sx_off[-1] - 1) * 100
    dst = (st_on[-1] / st_off[-1] - 1) * 100
    dex = (ex_on[-1] / ex_off[-1] - 1) * 100
    dey = (ey_on[-1] / ey_off[-1] - 1) * 100

    print(f"    Δσ_x={dsx:+.2f}%  Δσ_t={dst:+.1f}%  Δε_x={dex:+.2f}%  Δε_y={dey:+.2f}%  "
          f"σ_x'={sxp0*1e3:.3f} mrad")

    if make_plots:
        _make_plots(z_off, z_on, sx_off, sy_off, sx_on, sy_on,
                    st_off, st_on, ex_off, ex_on, ey_off, ey_on,
                    r_off[7], r_off[8],
                    xf_off, xpf_off, xf_on, xpf_on,
                    tf_off, df_off, tf_on, df_on,
                    sx0, eps_x0, Q_C, eps_use, tag)

    return {"Q": Q_C, "eps_n": eps_use,
            "sx_off": sx_off[-1], "sx_on": sx_on[-1],
            "st_off": st_off[-1], "st_on": st_on[-1],
            "ex_off": ex_off[-1], "ex_on": ex_on[-1],
            "ey_off": ey_off[-1], "ey_on": ey_on[-1],
            "sxp0": sxp0, "z_off": z_off, "z_on": z_on,
            "sx_off_arr": sx_off, "sx_on_arr": sx_on}

# ═══════════════════════════════════════════════════════════
#  plotting
# ═══════════════════════════════════════════════════════════

def _make_plots(z_off, z_on, sx_off, sy_off, sx_on, sy_on,
                st_off, st_on, ex_off, ex_on, ey_off, ey_on,
                x0, xp0, xf_off, xpf_off, xf_on, xpf_on,
                tf_off, df_off, tf_on, df_on,
                sx0, eps_x0, Q_C, eps_use, tag):

    z_off_mm = z_off * 1e3
    z_on_mm  = z_on  * 1e3
    Q_fC = Q_C * 1e15

    # ——— figure 1: σ_x,y(z) ———
    fig1, ax = plt.subplots(figsize=(9, 5))
    ax.plot(z_off_mm, sx_off * 1e6, "b-",  linewidth=1.5, alpha=0.7, label="SC OFF  $\\sigma_x$")
    ax.plot(z_off_mm, sy_off * 1e6, "b--", linewidth=1.0, alpha=0.7, label="SC OFF  $\\sigma_y$")
    ax.plot(z_on_mm,  sx_on  * 1e6, "r-",  linewidth=1.5, alpha=0.7, label="SC ON   $\\sigma_x$")
    ax.plot(z_on_mm,  sy_on  * 1e6, "r--", linewidth=1.0, alpha=0.7, label="SC ON   $\\sigma_y$")
    ax.set_xlabel("z  [mm]")
    ax.set_ylabel(r"$\sigma_{x,y}$  [$\mu$m]")
    ax.set_title(f"Transverse beam size  —  Q = {Q_fC:.0f} fC,  ε_n = {eps_use*1e6:.3f} mm·mrad",
                fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.25)
    fig1.tight_layout()
    fig1.savefig(os.path.join(_OUTDIR, f"benchmark_sc_sigmaxy_{tag}.png"), dpi=150)
    plt.close(fig1)

    # ——— figure 2: σ_t(z) ———
    fig2, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(z_off_mm, st_off * 1e12, "b-",  linewidth=1.5, label="SC OFF")
    ax.plot(z_on_mm,  st_on  * 1e12, "r-",  linewidth=1.5, label="SC ON")
    ax.set_xlabel("z  [mm]")
    ax.set_ylabel(r"$\sigma_t$  [fs]")
    ax.set_title(f"Longitudinal bunch length  —  Q = {Q_fC:.0f} fC", fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.25)
    fig2.tight_layout()
    fig2.savefig(os.path.join(_OUTDIR, f"benchmark_sc_sigmat_{tag}.png"), dpi=150)
    plt.close(fig2)

    # ——— figure 3: ε_x(z), ε_y(z) ———
    fig3, ax = plt.subplots(figsize=(9, 5))
    ax.plot(z_off_mm, ex_off * 1e6, "b-",  linewidth=1.5, alpha=0.7, label="SC OFF  $\\varepsilon_x$")
    ax.plot(z_off_mm, ey_off * 1e6, "b--", linewidth=1.0, alpha=0.7, label="SC OFF  $\\varepsilon_y$")
    ax.plot(z_on_mm,  ex_on  * 1e6, "r-",  linewidth=1.5, alpha=0.7, label="SC ON   $\\varepsilon_x$")
    ax.plot(z_on_mm,  ey_on  * 1e6, "r--", linewidth=1.0, alpha=0.7, label="SC ON   $\\varepsilon_y$")
    ax.axhline(eps_x0 * 1e6, color="gray", linestyle=":", linewidth=0.8)
    ax.set_xlabel("z  [mm]")
    ax.set_ylabel(r"$\varepsilon_{x,y}$  [mm$\cdot$mrad]")
    ax.set_title(f"Emittance evolution  —  Q = {Q_fC:.0f} fC", fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.25)
    # add emittance growth annotation
    deg_x = (ex_on[-1] / ex_off[-1] - 1) * 100
    deg_y = (ey_on[-1] / ey_off[-1] - 1) * 100
    ax.text(0.98, 0.95,
            f"$\\Delta\\varepsilon_x$ / $\\varepsilon_x$ = {deg_x:+.2f} %\n"
            f"$\\Delta\\varepsilon_y$ / $\\varepsilon_y$ = {deg_y:+.2f} %",
            transform=ax.transAxes, ha="right", va="top", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", alpha=0.8))
    fig3.tight_layout()
    fig3.savefig(os.path.join(_OUTDIR, f"benchmark_sc_emit_{tag}.png"), dpi=150)
    plt.close(fig3)

    # ——— figure 4: x–x′ phase space  before / after ———
    n_plot = min(5000, len(xf_off))
    rng = np.random.default_rng(42)
    idx = rng.choice(len(xf_off), n_plot, replace=False)
    idx_init = rng.choice(len(x0), n_plot, replace=False)

    fig4, axes = plt.subplots(2, 2, figsize=(10, 8))
    fig4.suptitle(f"x–x′ Phase Space  —  Q = {Q_fC:.0f} fC", fontsize=12, fontweight="bold")
    (ax_off_i, ax_off_f), (ax_on_i, ax_on_f) = axes

    # SC OFF — initial (same beam, before SC)
    ax_off_i.scatter(x0[idx_init] * 1e6, xp0[idx_init] * 1e3, s=0.3, alpha=0.4, c="steelblue")
    ax_off_i.set_title(f"SC OFF — z = 0\n$\\sigma_x$={sx0*1e6:.1f} μm,  $\\varepsilon_x$={eps_x0*1e6:.4f} mm·mrad")
    ax_off_i.set_xlabel(r"$x$ [$\mu$m]")
    ax_off_i.set_ylabel(r"$x'$ [mrad]")

    # SC OFF — final
    ax_off_f.scatter(xf_off[idx] * 1e6, xpf_off[idx] * 1e3, s=0.3, alpha=0.4, c="steelblue")
    eps_off_f = emit(xf_off, xpf_off)
    ax_off_f.set_title(f"SC OFF — z = 1.0 m\n$\\sigma_x$={np.std(xf_off)*1e6:.1f} μm,  $\\varepsilon_x$={eps_off_f*1e6:.4f} mm·mrad")
    ax_off_f.set_xlabel(r"$x$ [$\mu$m]")
    ax_off_f.set_ylabel(r"$x'$ [mrad]")

    # SC ON — initial (identical to SC OFF initial)
    ax_on_i.scatter(x0[idx_init] * 1e6, xp0[idx_init] * 1e3, s=0.3, alpha=0.4, c="darkorange")
    ax_on_i.set_title(f"SC ON — z = 0\n$\\sigma_x$={sx0*1e6:.1f} μm,  $\\varepsilon_x$={eps_x0*1e6:.4f} mm·mrad")
    ax_on_i.set_xlabel(r"$x$ [$\mu$m]")
    ax_on_i.set_ylabel(r"$x'$ [mrad]")

    # SC ON — final
    ax_on_f.scatter(xf_on[idx] * 1e6, xpf_on[idx] * 1e3, s=0.3, alpha=0.4, c="darkorange")
    eps_on_f = emit(xf_on, xpf_on)
    ax_on_f.set_title(f"SC ON — z = 1.0 m\n$\\sigma_x$={np.std(xf_on)*1e6:.1f} μm,  $\\varepsilon_x$={eps_on_f*1e6:.4f} mm·mrad")
    ax_on_f.set_xlabel(r"$x$ [$\mu$m]")
    ax_on_f.set_ylabel(r"$x'$ [mrad]")

    for ax in axes.flat:
        ax.grid(True, alpha=0.25)

    fig4.tight_layout()
    fig4.savefig(os.path.join(_OUTDIR, f"benchmark_sc_xxp_{tag}.png"), dpi=150)
    plt.close(fig4)

# ═══════════════════════════════════════════════════════════
#  charge scan  (fixed ε_n)
# ═══════════════════════════════════════════════════════════

def charge_sweep(eps_n=None):
    eps_use = eps_n if eps_n is not None else epsilon_n
    Q_list = [10e-15, 50e-15, 100e-15, 500e-15, 1000e-15]  # 10 fC → 1 pC
    results = []
    sc_curves = []

    print(f"\n{'='*60}")
    print(f"  Charge Sweep  (ε_n = {eps_use*1e6:.3f} mm·mrad)")
    print(f"{'='*60}")

    for Q in Q_list:
        r = compare_at_charge(Q, eps_use, make_plots=False)
        results.append(r)
    for Q in Q_list:
        r_on = run_sc_sim(Q, eps_use, sc_on=True)
        sc_curves.append((r_on[0], r_on[1]))

    # summary
    print(f"\n  {'Q':>8s}  {'σ_x OFF(μm)':>13s}  {'σ_x ON(μm)':>13s}  "
          f"{'Δσ_x(%)':>9s}  {'Δε_x(%)':>9s}  {'Δσ_t(%)':>9s}  {'σ_x′(mrad)':>10s}")
    for r in results:
        dsx = (r["sx_on"] / r["sx_off"] - 1) * 100
        dex = (r["ex_on"] / r["ex_off"] - 1) * 100
        dst = (r["st_on"] / r["st_off"] - 1) * 100
        q_str = f"{r['Q']*1e12:.0f} pC" if r['Q'] >= 1e-12 else f"{r['Q']*1e15:.0f} fC"
        print(f"  {q_str:>8s}  {r['sx_off']*1e6:13.1f}  {r['sx_on']*1e6:13.1f}  "
              f"{dsx:9.2f}  {dex:9.2f}  {dst:9.1f}  {r['sxp0']*1e3:10.3f}")

    # plot: σ_x(z) SC ON vs Q
    fig, ax = plt.subplots(figsize=(9, 5))
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(Q_list)))
    for (z, sx), Q, c in zip(sc_curves, Q_list, colors):
        q_str = f"{Q*1e12:.0f} pC" if Q >= 1e-12 else f"{Q*1e15:.0f} fC"
        ax.plot(z * 1e3, sx * 1e6, color=c, linewidth=1.5, label=q_str)
    r0 = run_sc_sim(10e-15, eps_use, sc_on=False)
    ax.plot(r0[0] * 1e3, r0[1] * 1e6, "k--", linewidth=1, alpha=0.5, label="SC OFF")
    ax.set_xlabel("z  [mm]"); ax.set_ylabel(r"$\sigma_x$  [$\mu$m]")
    ax.set_title(f"Charge scan — ε_n = {eps_use*1e6:.3f} mm·mrad", fontweight="bold")
    ax.legend(fontsize=8, ncol=2); ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fname = f"benchmark_sc_charge_sweep_eps{eps_use*1e6:.3f}.png"
    fig.savefig(os.path.join(_OUTDIR, fname), dpi=150)
    print(f"  -> {fname}")
    plt.close(fig)
    return results

# ═══════════════════════════════════════════════════════════
#  emittance scan  (fixed Q)
# ═══════════════════════════════════════════════════════════

def emittance_sweep(Q_C=None):
    Q_use = Q_C if Q_C is not None else Q_bunch
    eps_list = [0.001e-6, 0.005e-6, 0.01e-6, 0.08e-6]  # mm·mrad
    results = []
    sc_curves = []

    print(f"\n{'='*60}")
    print(f"  Emittance Sweep  (Q = {Q_use*1e15:.0f} fC)")
    print(f"{'='*60}")

    for eps in eps_list:
        r = compare_at_charge(Q_use, eps, make_plots=False)
        results.append(r)
    for eps in eps_list:
        r_on = run_sc_sim(Q_use, eps, sc_on=True)
        sc_curves.append((r_on[0], r_on[1]))

    # summary
    print(f"\n  {'ε_n(mm·mrad)':>14s}  {'σ_x OFF(μm)':>13s}  {'σ_x ON(μm)':>13s}  "
          f"{'Δσ_x(%)':>9s}  {'Δε_x(%)':>9s}  {'Δσ_t(%)':>9s}  {'σ_x′(mrad)':>10s}")
    for r in results:
        dsx = (r["sx_on"] / r["sx_off"] - 1) * 100
        dex = (r["ex_on"] / r["ex_off"] - 1) * 100
        dst = (r["st_on"] / r["st_off"] - 1) * 100
        print(f"  {r['eps_n']*1e6:14.3f}  {r['sx_off']*1e6:13.1f}  "
              f"{r['sx_on']*1e6:13.1f}  {dsx:9.2f}  {dex:9.2f}  "
              f"{dst:9.1f}  {r['sxp0']*1e3:10.3f}")

    # plot: σ_x(z) SC ON vs ε_n
    fig, ax = plt.subplots(figsize=(9, 5))
    colors = plt.cm.plasma(np.linspace(0.2, 0.9, len(eps_list)))
    for (z, sx), eps, c in zip(sc_curves, eps_list, colors):
        ax.plot(z * 1e3, sx * 1e6, color=c, linewidth=1.5,
                label=f"ε_n = {eps*1e6:.3f} mm·mrad")
    r0 = run_sc_sim(Q_use, eps_list[0], sc_on=False)
    ax.plot(r0[0] * 1e3, r0[1] * 1e6, "k--", linewidth=1, alpha=0.5, label="SC OFF")
    ax.set_xlabel("z  [mm]"); ax.set_ylabel(r"$\sigma_x$  [$\mu$m]")
    ax.set_title(f"Emittance scan — Q = {Q_use*1e15:.0f} fC  (SC ON)", fontweight="bold")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(_OUTDIR, f"benchmark_sc_emit_sweep_Q{Q_use*1e15:.0f}fC.png"), dpi=150)
    print(f"  -> benchmark_sc_emit_sweep_Q{Q_use*1e15:.0f}fC.png")
    plt.close(fig)
    return results

# ═══════════════════════════════════════════════════════════
#  combined scan:  Q × ε_n  grid
# ═══════════════════════════════════════════════════════════

def combined_sweep():
    Q_list   = [10e-15, 50e-15, 100e-15, 500e-15, 1000e-15]
    eps_list = [0.001e-6, 0.005e-6, 0.01e-6, 0.08e-6]
    grid = np.zeros((len(eps_list), len(Q_list)))  # Δσ_x (%)
    sgx = np.zeros((len(eps_list), len(Q_list)))   # σ_x ON (μm)

    print(f"\n{'='*65}")
    print(f"  Combined Sensitivity Scan:  Q × ε_n")
    print(f"{'='*65}")

    for i, eps in enumerate(eps_list):
        for j, Q in enumerate(Q_list):
            r = compare_at_charge(Q, eps, make_plots=False)
            grid[i, j] = (r["sx_on"] / r["sx_off"] - 1) * 100
            sgx[i, j] = r["sx_on"] * 1e6

    # — summary table —
    print(f"\n  Δσ_x (%)  grid  [rows = ε_n, cols = Q]")
    header = "  ε_n\\Q      " + "".join(f"{'10fC':>8s}" if Q<1e-12 else f"{Q*1e12:.0f}pC".rjust(8)
                                         for Q in Q_list)
    print(header)
    for i, eps in enumerate(eps_list):
        row = f"  {eps*1e6:.3f} mm·mrad"
        for j in range(len(Q_list)):
            row += f"  {grid[i,j]:6.2f}%"
        print(row)

    print(f"\n  σ_x(end) SC ON  grid  [μm]")
    print(header)
    for i, eps in enumerate(eps_list):
        row = f"  {eps*1e6:.3f} mm·mrad"
        for j in range(len(Q_list)):
            row += f"  {sgx[i,j]:6.1f}"
        print(row)

    # — heatmap —
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle("Space Charge Sensitivity: Q × ε_n", fontweight="bold")

    im1 = ax1.pcolormesh(grid, cmap="Reds", edgecolors="gray", linewidth=0.5)
    ax1.set_xticks(np.arange(len(Q_list)) + 0.5)
    ax1.set_xticklabels([f"{Q*1e15:.0f} fC" for Q in Q_list], rotation=30)
    ax1.set_yticks(np.arange(len(eps_list)) + 0.5)
    ax1.set_yticklabels([f"{eps*1e6:.3f}" for eps in eps_list])
    ax1.set_xlabel("Charge Q")
    ax1.set_ylabel(r"ε_n [mm·mrad]")
    ax1.set_title(r"$\Delta\sigma_x$ (%)  —  SC ON vs OFF")
    for i in range(len(eps_list)):
        for j in range(len(Q_list)):
            ax1.text(j + 0.5, i + 0.5, f"{grid[i,j]:.1f}%", ha="center", va="center",
                     fontsize=8, color="black" if grid[i,j] < 20 else "white")
    plt.colorbar(im1, ax=ax1, label=r"$\Delta\sigma_x$ (%)")

    im2 = ax2.pcolormesh(sgx, cmap="Blues", edgecolors="gray", linewidth=0.5)
    ax2.set_xticks(np.arange(len(Q_list)) + 0.5)
    ax2.set_xticklabels([f"{Q*1e15:.0f} fC" for Q in Q_list], rotation=30)
    ax2.set_yticks(np.arange(len(eps_list)) + 0.5)
    ax2.set_yticklabels([f"{eps*1e6:.3f}" for eps in eps_list])
    ax2.set_xlabel("Charge Q")
    ax2.set_ylabel(r"ε_n [mm·mrad]")
    ax2.set_title(r"$\sigma_x$ at z=1m  SC ON  [μm]")
    for i in range(len(eps_list)):
        for j in range(len(Q_list)):
            ax2.text(j + 0.5, i + 0.5, f"{sgx[i,j]:.0f}", ha="center", va="center",
                     fontsize=7)
    plt.colorbar(im2, ax=ax2, label=r"$\sigma_x$ [μm]")

    fig.tight_layout()
    fig.savefig(os.path.join(_OUTDIR, "benchmark_sc_sensitivity_heatmap.png"), dpi=150)
    print("  -> benchmark_sc_sensitivity_heatmap.png")
    plt.close(fig)
    return grid, sgx

# ═══════════════════════════════════════════════════════════
#  main
# ═══════════════════════════════════════════════════════════

if __name__ == "__main__":
    eps_geom_default = epsilon_n / beta_gamma
    sxp_default = eps_geom_default / spot_rms
    print(f"\n  OCELOT {ocelot.__version__}  |  SpaceCharge: 3D PIC, mesh [63,63,63]")
    print(f"  E_k = {E_keV:.0f} keV  βγ = {beta_gamma:.4f}")
    print(f"  ε_n = {epsilon_n*1e6:.3f} mm·mrad  →  ε_geo = {eps_geom_default*1e6:.4f} mm·mrad")
    print(f"  σ_x' = {sxp_default*1e3:.3f} mrad  (emittance angular spread)")

    if run_all:
        combined_sweep()
    elif run_escan:
        emittance_sweep()
    else:
        compare_at_charge(Q_bunch, None, make_plots=True)

    print("\n  Benchmark complete.")
