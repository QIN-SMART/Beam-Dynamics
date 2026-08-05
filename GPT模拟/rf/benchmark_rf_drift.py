#!/usr/bin/env python3
"""
Phase 2B — RF Longitudinal Dynamics Benchmark

RF kick:  δ_new = δ_old + (e·V_RF / E_0) · sin(φ + k_rf · z)
where z = -β·c·τ  (physical longitudinal position from OCELOT tau).

Drift compression:  τ_final = τ_initial + R56 · δ   (R56 < 0 for drift)
Bunch compresses when chirp δ(z) has correct sign (φ near π).

No space charge.

Usage:
  python3 benchmark_rf_drift.py                          # default
  python3 benchmark_rf_drift.py --phi 3.14               # phase scan
  python3 benchmark_rf_drift.py --V 50e3                 # 50 kV
  python3 benchmark_rf_drift.py --phi-scan               # full phase scan
"""

import sys, os, numpy as np

print("加载 OCELOT …", flush=True)
import ocelot
from ocelot.cpbd.elements import Drift
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.beam import generate_parray
from ocelot.cpbd.navi import Navigator
from ocelot.cpbd.track import tracking_step

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
epsilon_n   = 0.08e-6         # m·rad
sigma_delta = 1.0e-4          # ΔE/E
Q_bunch     = 100e-15         # C (only for beam gen, SC off)

# relativistic
E_rest      = 511.0
gamma       = 1.0 + E_keV / E_rest
beta        = np.sqrt(1.0 - 1.0 / gamma**2)
beta_gamma  = beta * gamma
E0_eV       = E_keV * 1e3     # eV
e_charge    = 1.602176634e-19

# ═══════════════════════════════════════════════════════════
#  RF parameters
# ═══════════════════════════════════════════════════════════

f_RF        = 2.856e9         # Hz  (S-band)
V_RF        = 30e3            # V   (30 kV default)
phi_RF      = np.pi           # rad (zero-crossing, max chirp)
L_drift_RF  = 0.5             # m   drift after RF cavity

k_rf        = 2 * np.pi * f_RF / c_SI if 'c_SI' in dir() else 2*np.pi*f_RF/2.99792458e8
# Compute below

# CLI
run_phi_scan = False
for i, a in enumerate(sys.argv):
    if a == "--V" and i + 1 < len(sys.argv):
        V_RF = float(sys.argv[i + 1])
    elif a == "--phi" and i + 1 < len(sys.argv):
        phi_RF = float(sys.argv[i + 1])
    elif a == "--phi-scan":
        run_phi_scan = True

c_SI        = 2.99792458e8
k_rf        = 2 * np.pi * f_RF / c_SI       # rad/m
lambda_rf   = c_SI / f_RF                    # m

# ═══════════════════════════════════════════════════════════
#  RF kick implementation
# ═══════════════════════════════════════════════════════════

def apply_rf_kick(p, V, phi):
    """
    Apply RF longitudinal kick to ParticleArray.
    δ_new = δ_old + (e·V / E₀) · sin(φ + k_rf · z)
    where z = -β·τ  (physical longitudinal position; OCELOT τ = c·t in m).
    """
    tau = p.tau()
    z_phys = -beta * tau
    d_delta = (V / (E_keV * 1000)) * np.sin(phi + k_rf * z_phys)
    p.rparticles[5, :] += d_delta

# ═══════════════════════════════════════════════════════════
#  beam generation
# ═══════════════════════════════════════════════════════════

def generate_beam():
    eps_geom = epsilon_n / beta_gamma
    sxp = eps_geom / spot_rms
    p = generate_parray(
        sigma_x=spot_rms, sigma_y=spot_rms,
        sigma_tau=sig_z0 / beta,          # OCELOT tau = c·t [m]; σ_tau = σ_z/β
        energy=(E_keV + 511.0) * 1e-6,    # TOTAL energy in GeV (E_kin+mc²)
        charge=Q_bunch,
        nparticles=50000,
    )
    np.random.seed(42)
    N = p.rparticles.shape[1]
    p.rparticles[1, :] = np.random.normal(0.0, sxp, N)
    p.rparticles[3, :] = np.random.normal(0.0, sxp, N)
    p.rparticles[5, :] = np.random.normal(0.0, sigma_delta, N)
    return p

# ═══════════════════════════════════════════════════════════
#  diagnostics
# ═══════════════════════════════════════════════════════════

def emit_long(tau, delta):
    """Longitudinal emittance: √(<τ²><δ²> - <τδ>²)."""
    return np.sqrt(np.mean(tau**2) * np.mean(delta**2) - np.mean(tau * delta)**2)

def track_and_record(p, L, dz=0.005, record_every=20):
    """Track through drift of length L, record σ_tau and σ_delta."""
    lat = MagneticLattice([Drift(l=L, eid="D")])
    lat.update_transfer_maps()
    navi = Navigator(lat)
    n_steps = int(L / dz)
    z_hist, st_hist, sd_hist, emit_hist = [], [], [], []
    for step in range(n_steps):
        tracking_step(lat, p, dz, navi)
        if step % record_every == 0:
            z_hist.append(navi.z0)
            st_hist.append(np.std(p.tau()))
            sd_hist.append(np.std(p.p()))
            emit_hist.append(emit_long(p.tau(), p.p()))
    return (np.array(z_hist), np.array(st_hist),
            np.array(sd_hist), np.array(emit_hist))

# ═══════════════════════════════════════════════════════════
#  single run
# ═══════════════════════════════════════════════════════════

def run_rf_benchmark(V, phi, make_plots=True):
    tag = f"V{V*1e-3:.0f}kV_phi{phi:.2f}"

    # —— generate beam ——
    p = generate_beam()
    tau0 = p.tau().copy()
    delta0 = p.p().copy()
    sig_t0 = np.std(tau0)
    sig_d0 = np.std(delta0)
    eps_l0 = emit_long(tau0, delta0)

    # —— RF kick ——
    apply_rf_kick(p, V, phi)
    tau_rf = p.tau().copy()
    delta_rf = p.p().copy()
    sig_t_rf = np.std(tau_rf)
    sig_d_rf = np.std(delta_rf)

    # —— chirp before/after RF ——
    z0 = -beta * tau0
    z_rf = -beta * tau_rf
    chirp0 = np.polyfit(z0, delta0, 1)[0] if len(z0) > 1 else 0
    chirp_rf = np.polyfit(z_rf, delta_rf, 1)[0] if len(z_rf) > 1 else 0

    # —— drift tracking after RF ——
    z_arr, st_arr, sd_arr, emit_arr = track_and_record(p, L_drift_RF)

    # —— results ——
    sig_t_min = np.min(st_arr)
    idx_min = np.argmin(st_arr)
    z_min = z_arr[idx_min]
    compression = sig_t_min / sig_t0

    print(f"\n  RF: V={V*1e-3:.0f} kV  φ={phi:.3f} rad ({phi*180/np.pi:.0f}°)  "
          f"f={f_RF*1e-9:.3f} GHz  λ={lambda_rf*1e3:.1f} mm")
    print(f"  Chirp: {chirp0:.2f} /m (before)  →  {chirp_rf:.2f} /m (after RF)")
    print(f"  Bunch length: σ_t0={sig_t0*1e15:.0f} fs  →  "
          f"σ_t_min={sig_t_min*1e15:.0f} fs @ z={z_min*1e3:.0f} mm")
    print(f"  Compression ratio: σ_t_min / σ_t0 = {compression:.3f}")
    print(f"  ε_long: {eps_l0*1e15:.2f} fs → {emit_arr[-1]*1e15:.2f} fs  "
          f"(Δ={(emit_arr[-1]/eps_l0-1)*100:+.1f}%)")

    if make_plots:
        _make_rf_plots(tau0, delta0, z0, tau_rf, delta_rf, z_rf,
                       z_arr, st_arr, sd_arr, emit_arr,
                       sig_t0, sig_d0, eps_l0, V, phi, chirp0, chirp_rf, tag)

    return {"V": V, "phi": phi, "sig_t0": sig_t0, "sig_t_min": sig_t_min,
            "z_min": z_min, "compression": compression,
            "chirp0": chirp0, "chirp_rf": chirp_rf, "eps_l0": eps_l0,
            "eps_l_final": emit_arr[-1]}

# ═══════════════════════════════════════════════════════════
#  plots
# ═══════════════════════════════════════════════════════════

def _make_rf_plots(tau0, delta0, z0, tau_rf, delta_rf, z_rf,
                   z_arr, st_arr, sd_arr, emit_arr,
                   sig_t0, sig_d0, eps_l0, V, phi, chirp0, chirp_rf, tag):
    n_plt = min(5000, len(tau0))
    rng = np.random.default_rng(42)
    idx = rng.choice(len(tau0), n_plt, replace=False)

    # —— figure 1: z-δ phase space  before & after RF ——
    fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
    fig1.suptitle(f"Longitudinal Phase Space  —  "
                  f"V={V*1e-3:.0f} kV  φ={phi:.2f} rad ({phi*180/np.pi:.0f}°)",
                  fontweight="bold")

    ax1.scatter(z0[idx] * 1e6, delta0[idx] * 1e3, s=0.3, alpha=0.4, c="steelblue")
    ax1.set_xlabel(r"$z$  [$\mu$m]")
    ax1.set_ylabel(r"$\delta$  [$10^{-3}$]")
    ax1.set_title(f"Before RF  (σ_τ={sig_t0*1e15:.0f} fs,  ε_L={eps_l0*1e15:.2f} fs)")
    ax1.grid(True, alpha=0.25)

    ax2.scatter(z_rf[idx] * 1e6, delta_rf[idx] * 1e3, s=0.3, alpha=0.4, c="darkorange")
    ax2.set_xlabel(r"$z$  [$\mu$m]")
    ax2.set_ylabel(r"$\delta$  [$10^{-3}$]")
    ax2.set_title(f"After RF  (chirp={chirp_rf*1e-6:.1f}/m)")
    ax2.grid(True, alpha=0.25)

    fig1.tight_layout()
    fig1.savefig(os.path.join(_OUTDIR, f"benchmark_rf_phasespace_{tag}.png"), dpi=150)
    plt.close(fig1)

    # —— figure 2: σ_τ(z) and σ_δ(z) during drift ——
    fig2, (ax_a, ax_b) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    fig2.suptitle(f"Drift Compression After RF  —  "
                  f"V={V*1e-3:.0f} kV  φ={phi:.2f} rad", fontweight="bold")

    ax_a.plot(z_arr * 1e3, st_arr * 1e15, "b-", linewidth=1.5)
    ax_a.axhline(sig_t0 * 1e15, color="gray", linestyle="--", alpha=0.5,
                 label=f"initial σ_τ = {sig_t0*1e15:.0f} fs")
    idx_min = np.argmin(st_arr)
    ax_a.axvline(z_arr[idx_min] * 1e3, color="red", linestyle=":", alpha=0.5,
                 label=f"min @ {z_arr[idx_min]*1e3:.0f} mm")
    ax_a.set_ylabel(r"$\sigma_\tau$  [fs]")
    ax_a.legend(fontsize=8); ax_a.grid(True, alpha=0.25)

    ax_b.plot(z_arr * 1e3, sd_arr * 1e3, "r-", linewidth=1.5)
    ax_b.set_xlabel("z  [mm]")
    ax_b.set_ylabel(r"$\sigma_\delta$  [$10^{-3}$]")
    ax_b.grid(True, alpha=0.25)

    fig2.tight_layout()
    fig2.savefig(os.path.join(_OUTDIR, f"benchmark_rf_compression_{tag}.png"), dpi=150)
    plt.close(fig2)

    # —— figure 3: longitudinal emittance ——
    fig3, ax = plt.subplots(figsize=(9, 4))
    ax.plot(z_arr * 1e3, emit_arr * 1e15, "g-", linewidth=1.5)
    ax.axhline(eps_l0 * 1e15, color="gray", linestyle="--",
               label=f"initial ε_L = {eps_l0*1e15:.2f} fs")
    ax.set_xlabel("z  [mm]")
    ax.set_ylabel(r"$\varepsilon_L$  [fs]")
    ax.set_title(f"Longitudinal emittance  —  V={V*1e-3:.0f} kV  φ={phi:.2f} rad",
                 fontweight="bold")
    ax.legend(); ax.grid(True, alpha=0.25)
    fig3.tight_layout()
    fig3.savefig(os.path.join(_OUTDIR, f"benchmark_rf_emit_{tag}.png"), dpi=150)
    plt.close(fig3)

    print(f"  -> benchmark_rf_phasespace_{tag}.png")
    print(f"  -> benchmark_rf_compression_{tag}.png")
    print(f"  -> benchmark_rf_emit_{tag}.png")

# ═══════════════════════════════════════════════════════════
#  φ scan
# ═══════════════════════════════════════════════════════════

def phi_scan():
    phi_list = np.linspace(0, 2 * np.pi, 25)  # 0 to 2π
    comp_list = []
    chirp_list = []
    zmin_list = []

    print(f"\n{'='*60}")
    print(f"  RF Phase Scan  (V = {V_RF*1e-3:.0f} kV)")
    print(f"{'='*60}")

    for phi in phi_list:
        r = run_rf_benchmark(V_RF, phi, make_plots=False)
        comp_list.append(r["compression"])
        chirp_list.append(r["chirp_rf"])
        zmin_list.append(r["z_min"])

    # —— plot ——
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    fig.suptitle(f"RF Phase Scan  —  V = {V_RF*1e-3:.0f} kV  f = {f_RF*1e-9:.3f} GHz",
                 fontweight="bold")

    ax1.plot(phi_list * 180 / np.pi, np.array(comp_list), "b-o", markersize=4)
    ax1.axhline(1.0, color="gray", linestyle="--", alpha=0.5, label="no compression")
    ax1.set_ylabel("σ_τ,min / σ_τ,0")
    ax1.legend(); ax1.grid(True, alpha=0.25)
    ax1.set_title("Compression ratio vs RF phase")

    ax2.plot(phi_list * 180 / np.pi, np.array(chirp_list) * 1e-6, "r-o", markersize=4)
    ax2.axhline(0, color="gray", linestyle="--", alpha=0.5)
    ax2.set_xlabel("RF phase φ  [deg]")
    ax2.set_ylabel("chirp dδ/dz  [1/m]")
    ax2.grid(True, alpha=0.25)
    ax2.set_title("Energy chirp vs RF phase")

    # annotate optimal compression
    idx_opt = np.argmin(comp_list)
    ax1.annotate(f"φ={phi_list[idx_opt]*180/np.pi:.0f}°\ncompression={comp_list[idx_opt]:.3f}",
                 xy=(phi_list[idx_opt]*180/np.pi, comp_list[idx_opt]),
                 xytext=(phi_list[idx_opt]*180/np.pi + 30, comp_list[idx_opt] + 0.1),
                 arrowprops=dict(arrowstyle="->", color="red"), fontsize=8, color="red")

    fig.tight_layout()
    fig.savefig(os.path.join(_OUTDIR, "benchmark_rf_phi_scan.png"), dpi=150)
    print("  -> benchmark_rf_phi_scan.png")
    plt.close(fig)

    return phi_list, comp_list, chirp_list

# ═══════════════════════════════════════════════════════════
#  main
# ═══════════════════════════════════════════════════════════

if __name__ == "__main__":
    eps_geom_default = epsilon_n / beta_gamma
    print(f"\n  OCELOT {ocelot.__version__}  |  RF Longitudinal Benchmark")
    print(f"  E_k = {E_keV:.0f} keV  βγ = {beta_gamma:.4f}")
    print(f"  f_RF = {f_RF*1e-9:.3f} GHz  λ_RF = {lambda_rf*1e3:.1f} mm")
    print(f"  k_rf = {k_rf:.1f} rad/m  (RF kick model: analytic sin)")

    if run_phi_scan:
        phi_scan()
    else:
        run_rf_benchmark(V_RF, phi_RF, make_plots=True)

    print("\n  Benchmark complete.")
