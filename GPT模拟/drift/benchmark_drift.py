#!/usr/bin/env python3
"""
Phase 1.5 — Drift-only benchmark validation.

Validates:
  1.  Beam divergence follows analytical formula:  σ_x(z) = √(σ₀² + σ_x'²·z²)
  2.  Emittance conservation:  ε_x(z) = constant
  3.  Phase space plots:  x-x'  and  z-δ  (tau–delta)

No solenoids, no RF, no space charge.

Usage:
  python3 benchmark_drift.py
  python3 benchmark_drift.py --epsn 0.04e-6 --sigd 2e-4
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
#  parameters  (same defaults as ocelot_beamline.py)
# ═══════════════════════════════════════════════════════════

spot_rms    = 85e-6           # m
sig_z0      = 300e-6          # m
E_keV       = 100.0           # keV
epsilon_n   = 0.08e-6         # m·rad  normalized emittance
sigma_delta = 1.0e-4          # ΔE/E   relative energy spread

# CLI overrides
for i, a in enumerate(sys.argv):
    if a == "--epsn" and i + 1 < len(sys.argv):
        epsilon_n = float(sys.argv[i + 1])
    elif a == "--sigd" and i + 1 < len(sys.argv):
        sigma_delta = float(sys.argv[i + 1])

# relativistic
E_rest     = 511.0            # keV
gamma      = 1.0 + E_keV / E_rest
beta       = np.sqrt(1.0 - 1.0 / gamma**2)
beta_gamma = beta * gamma
epsilon_geom = epsilon_n / beta_gamma
sigma_xp   = epsilon_geom / spot_rms
sigma_yp   = epsilon_geom / spot_rms

# ═══════════════════════════════════════════════════════════
#  drift-only lattice — 1.0 m pure drift
# ═══════════════════════════════════════════════════════════

drift_length = 1.0           # m

lat = MagneticLattice([Drift(l=drift_length, eid="DRIFT")])
lat.update_transfer_maps()

# ═══════════════════════════════════════════════════════════
#  generate beam with Phase‑1 distribution
# ═══════════════════════════════════════════════════════════

p = generate_parray(
    sigma_x=spot_rms, sigma_y=spot_rms,
    sigma_tau=sig_z0 / beta,          # OCELOT tau = c·t [m]; σ_tau = σ_z/β
    energy=(E_keV + 511.0) * 1e-6,    # TOTAL energy in GeV (E_kin+mc²)
    charge=100_000,
)

np.random.seed(42)
n_part = p.rparticles.shape[1]

# angular spread (transverse emittance)
p.rparticles[1, :] = np.random.normal(0.0, sigma_xp, n_part)   # px = x'
p.rparticles[3, :] = np.random.normal(0.0, sigma_yp, n_part)   # py = y'

# energy spread
p.rparticles[5, :] = np.random.normal(0.0, sigma_delta, n_part)  # delta

# record initial phase space
x0  = p.x().copy()
xp0 = p.px().copy()
y0  = p.y().copy()
yp0 = p.py().copy()
t0  = p.tau().copy()
d0  = p.p().copy()

sig_x0  = np.std(x0)
sig_xp0 = np.std(xp0)
eps_x0  = np.sqrt(np.mean(x0**2) * np.mean(xp0**2) - np.mean(x0 * xp0)**2)

print(f"\n  Initial beam:")
print(f"    σ_x  = {sig_x0*1e6:.3f} μm")
print(f"    σ_x' = {sig_xp0*1e3:.3f} mrad")
print(f"    ε_x  = {eps_x0*1e6:.4f} mm·mrad")
print(f"    σ_δ  = {np.std(d0)*1e3:.2f} e-3")
print(f"    βγ   = {beta_gamma:.4f}")

# ═══════════════════════════════════════════════════════════
#  track  —  record at fine z resolution
# ═══════════════════════════════════════════════════════════

navi = Navigator(lat)
dz   = 0.002                       # 2 mm step
n_steps = int(drift_length / dz)

z_hist      = []
sigx_hist   = []
epsx_hist   = []
sigt_hist   = []
sigd_hist   = []

for step_i in range(n_steps):
    z_before = navi.z0
    tracking_step(lat, p, dz, navi)

    z = navi.z0
    x  = p.x()
    xp = p.px()
    ta = p.tau()
    dd = p.p()

    z_hist.append(z)
    sigx_hist.append(np.std(x))
    epsx_hist.append(np.sqrt(np.mean(x**2) * np.mean(xp**2) - np.mean(x * xp)**2))
    sigt_hist.append(np.std(ta))
    sigd_hist.append(np.std(dd))

z_arr      = np.array(z_hist)
sigx_arr   = np.array(sigx_hist)
epsx_arr   = np.array(epsx_hist)
sigt_arr   = np.array(sigt_hist)
sigd_arr   = np.array(sigd_hist)

# final phase space
xf  = p.x().copy()
xpf = p.px().copy()
tf  = p.tau().copy()
df  = p.p().copy()

# ═══════════════════════════════════════════════════════════
#  analytical comparison
# ═══════════════════════════════════════════════════════════

# σ_x(z) = √(σ_x²(0) + σ_x'² · z²)   for uncorrelated initial x-x'
sigx_analytic = np.sqrt(sig_x0**2 + sig_xp0**2 * z_arr**2)

epsx_ref = eps_x0  # should be constant

# ═══════════════════════════════════════════════════════════
#  print summary
# ═══════════════════════════════════════════════════════════

print(f"\n  Final beam (z = {drift_length} m):")
print(f"    σ_x  = {sigx_arr[-1]*1e6:.1f} μm  (analytic: {sigx_analytic[-1]*1e6:.1f} μm)")
print(f"    ε_x  = {epsx_arr[-1]*1e6:.4f} mm·mrad  (initial: {epsx_ref*1e6:.4f})")
print(f"    Δε_x/ε_x = {(epsx_arr[-1] - epsx_ref) / epsx_ref * 100:.4f} %")
print(f"    σ_δ  = {sigd_arr[-1]*1e3:.2f} e-3  (initial: {np.std(d0)*1e3:.2f})")

max_sigx_err = np.max(np.abs(sigx_arr - sigx_analytic)) / sig_x0 * 100
print(f"\n  Beam-size deviation from analytic (max): {max_sigx_err:.3f} % of σ_x0")

# ═══════════════════════════════════════════════════════════
#  figure 1 — beam size & emittance vs z
# ═══════════════════════════════════════════════════════════

fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
fig1.suptitle("Drift Benchmark — Phase 1.5 Validation", fontsize=13, fontweight="bold")

# top: σ_x(z)
ax1.plot(z_arr * 1e3, sigx_arr * 1e6, "b.", markersize=2, label="OCELOT tracking")
ax1.plot(z_arr * 1e3, sigx_analytic * 1e6, "r-", linewidth=1.5,
         label=r"analytic: $\sqrt{\sigma_{x0}^2 + \sigma_{x'}^2 z^2}$")
ax1.set_ylabel(r"$\sigma_x$  [$\mu$m]")
ax1.legend(fontsize=9)
ax1.grid(True, alpha=0.3)
ax1.text(0.98, 0.05,
         rf"$\sigma_{{x0}}={sig_x0*1e6:.1f}\;\mu$m, "
         rf"$\sigma_{{x'}}={sig_xp0*1e3:.2f}\;$mrad",
         transform=ax1.transAxes, ha="right", fontsize=8)

# bottom: ε_x(z)
ax2.plot(z_arr * 1e3, epsx_arr * 1e6, "b.", markersize=2, label="OCELOT tracking")
ax2.axhline(epsx_ref * 1e6, color="r", linewidth=1.5, linestyle="--",
            label=f"initial $\\varepsilon_x$ = {epsx_ref*1e6:.4f} mm·mrad")
ax2.set_xlabel("z  [mm]")
ax2.set_ylabel(r"$\varepsilon_x$  [mm$\cdot$mrad]")
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)
ax2.text(0.98, 0.95,
         f"max relative deviation: "
         f"{np.max(np.abs(epsx_arr - epsx_ref))/epsx_ref*100:.4f} %",
         transform=ax2.transAxes, ha="right", fontsize=8, color="gray")

fig1.tight_layout()
fig1.savefig(os.path.join(_OUTDIR, "benchmark_drift_beam_evolution.png"), dpi=150)
print("  -> benchmark_drift_beam_evolution.png")

# ═══════════════════════════════════════════════════════════
#  figure 2 — phase space:  x-x'  and  tau–delta
# ═══════════════════════════════════════════════════════════

fig2, axes = plt.subplots(2, 2, figsize=(11, 9))
fig2.suptitle("Phase Space — Drift Benchmark", fontsize=13, fontweight="bold")
(ax_xxp0, ax_xxpf), (ax_td0, ax_tdf) = axes

# subsample for scatter (too many points)
n_plot = min(8000, n_part)
idx = np.random.default_rng(42).choice(n_part, n_plot, replace=False)

# -- x–x' at z=0 --
ax_xxp0.scatter(x0[idx] * 1e6, xp0[idx] * 1e3, s=0.5, alpha=0.4, c="steelblue")
ax_xxp0.set_xlabel(r"$x$  [$\mu$m]")
ax_xxp0.set_ylabel(r"$x'$  [mrad]")
ax_xxp0.set_title(f"x–x′  at z=0    $\\varepsilon_x$ = {eps_x0*1e6:.4f} mm·mrad")
ax_xxp0.grid(True, alpha=0.3)

# -- x–x' at z=final --
eps_xf = np.sqrt(np.mean(xf**2) * np.mean(xpf**2) - np.mean(xf * xpf)**2)
ax_xxpf.scatter(xf[idx] * 1e6, xpf[idx] * 1e3, s=0.5, alpha=0.4, c="darkorange")
ax_xxpf.set_xlabel(r"$x$  [$\mu$m]")
ax_xxpf.set_ylabel(r"$x'$  [mrad]")
ax_xxpf.set_title(f"x–x′  at z={drift_length:.1f} m    $\\varepsilon_x$ = {eps_xf*1e6:.4f} mm·mrad")
ax_xxpf.grid(True, alpha=0.3)

# -- tau–delta at z=0 --
sig_d0 = np.std(d0)
ax_td0.scatter(t0[idx] * 1e12, d0[idx] * 1e3, s=0.5, alpha=0.4, c="steelblue")
ax_td0.set_xlabel(r"$\tau$  [ps]")
ax_td0.set_ylabel(r"$\delta$  [$10^{-3}$]")
ax_td0.set_title(r"$\tau$–$\delta$  at z=0")
ax_td0.grid(True, alpha=0.3)

# -- tau–delta at z=final --
sig_df = np.std(df)
ax_tdf.scatter(tf[idx] * 1e12, df[idx] * 1e3, s=0.5, alpha=0.4, c="darkorange")
ax_tdf.set_xlabel(r"$\tau$  [ps]")
ax_tdf.set_ylabel(r"$\delta$  [$10^{-3}$]")
ax_tdf.set_title(r"$\tau$–$\delta$  at z=%.1f m" % drift_length)
ax_tdf.grid(True, alpha=0.3)

# correlation annotation
rho_td0 = np.corrcoef(t0, d0)[0, 1]
rho_tdf = np.corrcoef(tf, df)[0, 1]
ax_tdf.text(0.98, 0.95,
            rf"corr($\tau$,$\delta$) = {rho_tdf:.4f}" "\n"
            rf"$\sigma_\delta$ = {sig_df*1e3:.4f}e-3  (conserved)",
            transform=ax_tdf.transAxes, ha="right", va="top", fontsize=8, color="gray")

fig2.tight_layout()
fig2.savefig(os.path.join(_OUTDIR, "benchmark_drift_phase_space.png"), dpi=150)
print("  -> benchmark_drift_phase_space.png")

print("\n  Benchmark complete.")
