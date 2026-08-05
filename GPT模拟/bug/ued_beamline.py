#!/usr/bin/env python3
"""
Phase 3A — UED Beamline Assembly

Assembles validated modules (drift, solenoid, RF, SC) into a complete
UED beamline:

  Cathode → Drift → Solenoid(TL1) → Drift → RF Cavity → Drift → Sample

Outputs σ_x(z), σ_y(z), σ_z(z), ε_x(z), ε_y(z).

Usage:
  python3 ued_beamline.py
  python3 ued_beamline.py --sc        # enable space charge
"""

import sys, os, numpy as np
import yaml

print("加载 OCELOT …", flush=True)
import ocelot
from ocelot.cpbd.elements import Drift, Solenoid
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.beam import generate_parray
from ocelot.cpbd.navi import Navigator
from ocelot.cpbd.track import tracking_step
try:
    from ocelot.cpbd.sc import SpaceCharge
    _HAS_SC = True
except ImportError:
    _HAS_SC = False

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_OUTDIR = os.path.dirname(os.path.abspath(__file__)) + "/"

# ═══════════════════════════════════════════════════════════
#  load configuration
# ═══════════════════════════════════════════════════════════

with open(os.path.join(_OUTDIR, "beamline_config.yaml"), "r") as f:
    cfg = yaml.safe_load(f)

b  = cfg["beam"]
ib = cfg["initial_distribution"]
so = cfg["solenoid"]
rf = cfg["rf_cavity"]
sc_cfg = cfg["space_charge"]
lat_cfg = cfg["lattice"]
out_cfg = cfg["output"]

# CLI
use_sc = "--sc" in sys.argv
if use_sc and _HAS_SC:
    sc_cfg["enabled"] = True
elif use_sc and not _HAS_SC:
    print("  WARNING: SpaceCharge module not available, SC disabled")

# ═══════════════════════════════════════════════════════════
#  physical constants
# ═══════════════════════════════════════════════════════════

e_SI   = 1.602176634e-19
m_e_SI = 9.10938356e-31
c_SI   = 2.99792458e8
mec2   = 511.0

# ═══════════════════════════════════════════════════════════
#  relativistic parameters (from config)
# ═══════════════════════════════════════════════════════════

E_keV      = b["energy_keV"]
gamma      = 1.0 + E_keV / mec2
beta       = np.sqrt(1.0 - 1.0 / gamma**2)
beta_gamma = beta * gamma
E0_eV      = E_keV * 1e3
p_SI       = gamma * m_e_SI * beta * c_SI

# ═══════════════════════════════════════════════════════════
#  initial beam parameters (from config)
# ═══════════════════════════════════════════════════════════

spot_rms    = ib["spot_rms_um"] * 1e-6
sig_z0      = ib["bunch_length_um"] * 1e-6
epsilon_n   = ib["epsilon_n_mm_mrad"] * 1e-6
sigma_delta = ib["sigma_delta"]
Q_bunch_C   = b["charge_fC"] * 1e-15
N_part      = b["n_particles"]

epsilon_geom = epsilon_n / beta_gamma
sigma_xp     = epsilon_geom / spot_rms
sigma_yp     = epsilon_geom / spot_rms

# ═══════════════════════════════════════════════════════════
#  solenoid (from config)
# ═══════════════════════════════════════════════════════════

B_sol  = so["B_field_T"]
k_sol  = e_SI * B_sol / (2.0 * p_SI)
L_sol  = so["length_m"]
z_sol  = so["z_start_m"]

# ═══════════════════════════════════════════════════════════
#  RF cavity (from config)
# ═══════════════════════════════════════════════════════════

f_RF     = rf["frequency_GHz"] * 1e9
V_RF     = rf["voltage_kV"] * 1e3
phi_RF   = rf["phase_rad"]
L_rf_cav = rf["length_m"]
z_rf     = rf["z_start_m"]
k_rf     = 2.0 * np.pi * f_RF / c_SI

# ═══════════════════════════════════════════════════════════
#  lattice construction
# ═══════════════════════════════════════════════════════════

elements = []
z_cumulative = 0.0
for elem in lat_cfg["elements"]:
    z_start, length, etype = elem
    if etype == "cathode":
        z_cumulative = z_start
    elif etype == "drift":
        elements.append(Drift(l=length, eid=f"DRIFT_{z_start*1e3:.0f}mm"))
    elif etype == "solenoid":
        elements.append(Solenoid(l=length, k=k_sol, eid="SOL_TL1"))
    elif etype == "rf_cavity":
        # skip OCELOT Cavity; RF kick applied analytically at z_rf
        elements.append(Drift(l=length, eid=f"RFCAV_{z_start*1e3:.0f}mm"))
    elif etype == "sample":
        pass  # diagnostic plane only

lat = MagneticLattice(elements, method={"global": "SecondTM"})
lat.update_transfer_maps()

total_length = sum(e.l if hasattr(e, "l") else 0.0 for e in lat.sequence)

# ═══════════════════════════════════════════════════════════
#  diagnostics helpers
# ═══════════════════════════════════════════════════════════

def emit(x, xp):
    return np.sqrt(np.mean(x**2) * np.mean(xp**2) - np.mean(x * xp)**2)

def apply_rf_kick(parray):
    """RF longitudinal kick. delta = (p-p0)/p0 = dE/(β²·γ·m·c²)."""
    tau = parray.tau()
    z_phys = -beta * c_SI * tau
    E_total_eV = gamma * mec2 * 1e3
    d_delta = (V_RF / (beta**2 * E_total_eV)) * np.sin(phi_RF + k_rf * z_phys)
    parray.rparticles[5, :] += d_delta

def ang_momentum(x, px, y, py):
    """Canonical angular momentum Lz = x*py - y*px [normalized by P0]."""
    return x * py - y * px

# probe positions
probes = {f"z{z:.0f}mm": z * 1e-3 for z in out_cfg["z_diagnostics_mm"]}

# phase space snapshot positions
snapshot_z = {
    "cathode":     0.000,
    "after_sol":   0.160,
    "before_rf":   0.399,
    "after_rf":    0.422,
    "sample":      0.777,
}

# ═══════════════════════════════════════════════════════════
#  tracking function
# ═══════════════════════════════════════════════════════════

def run_beamline(sc_enabled=False):
    # generate fresh beam
    p = generate_parray(
        sigma_x=spot_rms, sigma_y=spot_rms,
        sigma_tau=sig_z0 / beta,   # OCELOT tau = c·t [m]; σ_tau = σ_z/β
        energy=(E_keV + 511.0) * 1e-6,   # TOTAL energy in GeV (E_kin+mc²)
        charge=Q_bunch_C,
        nparticles=N_part,
    )
    np.random.seed(42)
    N = p.rparticles.shape[1]
    p.rparticles[1, :] = np.random.normal(0.0, sigma_xp, N)
    p.rparticles[3, :] = np.random.normal(0.0, sigma_yp, N)
    p.rparticles[5, :] = np.random.normal(0.0, sigma_delta, N)

    sc = SpaceCharge(step=1) if (sc_enabled and _HAS_SC) else None

    navi = Navigator(lat, unit_step=out_cfg["step_size_m"])
    if sc is not None:
        navi.add_physics_proc(sc, lat.sequence[0], lat.sequence[-1])
    dz = out_cfg["step_size_m"]
    n_steps = int(total_length / dz)

    probes_remaining = dict(probes)
    snaps_rem = dict(snapshot_z)
    snapshots_dict = {}
    results_dict = {}

    z_hist, sigx_hist, sigy_hist, sigz_hist = [], [], [], []
    epsx_hist, epsy_hist, lz_hist = [], [], []

    rf_applied = False

    for step_i in range(n_steps):
        z_before = navi.z0

        # RF kick exactly at cavity entrance  z = 0.400 m
        if not rf_applied and z_before >= z_rf - 1e-12:
            apply_rf_kick(p)
            rf_applied = True

        # Longitudinal (tau) transport is done natively by OCELOT
        # (tau_f = tau_i + R56·δ, R56 = -L/(β²γ²)); no manual override.
        tracking_step(lat, p, dz, navi)

        # capture before_rf AFTER transport, BEFORE rf re-check
        if "before_rf" in snaps_rem and z_before < z_rf and navi.z0 >= z_rf - 1e-12:
            snaps_rem.pop("before_rf")
            snapshots_dict["before_rf"] = {
                "x": p.x().copy(), "xp": p.px().copy(),
                "y": p.y().copy(), "yp": p.py().copy(),
                "tau": p.tau().copy(), "delta": p.p().copy(),
                "z_m": navi.z0,
            }

        if not rf_applied and navi.z0 >= z_rf - 1e-12:
            apply_rf_kick(p)
            rf_applied = True

        z = navi.z0
        x = p.x(); xp = p.px(); y = p.y(); yp = p.py()

        # skip projected emittance inside solenoid (z=0.100→0.160)
        in_solenoid = z_sol < z < (z_sol + L_sol)

        z_hist.append(z * 1e3)
        sigx_hist.append(np.std(x) * 1e6)
        sigy_hist.append(np.std(y) * 1e6)
        sigz_hist.append(np.std(p.tau()) * c_SI * 1e6)
        epsx_hist.append(np.nan if in_solenoid else emit(x, xp) * 1e6)
        epsy_hist.append(np.nan if in_solenoid else emit(y, yp) * 1e6)
        lz_hist.append(np.mean(x * yp - y * xp))

        # emittance & Lz only at boundaries (before solenoid, after exit)
        for name, z_p in list(probes_remaining.items()):
            if z_before <= z_p < z:
                xa = p.x(); xpa = p.px(); ya = p.y(); ypa = p.py()
                in_sol_probe = z_sol < z_p < (z_sol + L_sol)
                results_dict[name] = {
                    "z_mm": z_p * 1e3,
                    "sigma_x_um": np.std(xa) * 1e6,
                    "sigma_y_um": np.std(ya) * 1e6,
                    "sigma_t_fs": np.std(p.tau()) * 1e15,
                    "eps_x_mm_mrad": np.nan if in_sol_probe else emit(xa, xpa) * 1e6,
                    "eps_y_mm_mrad": np.nan if in_sol_probe else emit(ya, ypa) * 1e6,
                    "sigma_delta_e3": np.std(p.p()) * 1e3,
                    "Lz_mean": np.mean(xa * ypa - ya * xpa),
                }
                del probes_remaining[name]

        for name, z_s in list(snaps_rem.items()):
            if z_before <= z_s < z:
                snapshots_dict[name] = {
                    "x": p.x().copy(), "xp": p.px().copy(),
                    "y": p.y().copy(), "yp": p.py().copy(),
                    "tau": p.tau().copy(), "delta": p.p().copy(),
                    "z_m": z,
                }
                del snaps_rem[name]

    return {
        "z_arr": np.array(z_hist),
        "sigx": np.array(sigx_hist),
        "sigy": np.array(sigy_hist),
        "sigz": np.array(sigz_hist),
        "epsx": np.array(epsx_hist),
        "epsy": np.array(epsy_hist),
        "lz": np.array(lz_hist),
        "results": results_dict,
        "snapshots": snapshots_dict,
        "sig_t_final": np.std(p.tau()),
    }

# ═══════════════════════════════════════════════════════════
#  run — SC OFF then SC ON
# ═══════════════════════════════════════════════════════════

def print_table(data, label):
    print(f"\n  [{label}]")
    print(f"  {'Probe':>10s}  {'z(mm)':>7s}  {'σ_x(μm)':>9s}  {'σ_y(μm)':>9s}  "
          f"{'σ_t(fs)':>9s}  {'ε_x(mm·mrad)':>13s}  {'σ_δ(e-3)':>9s}  {'⟨Lz⟩':>10s}")
    print(f"  {'-'*74}")
    for name in [f"z{z:.0f}mm" for z in out_cfg["z_diagnostics_mm"]]:
        if name in data["results"]:
            r = data["results"][name]
            eps_str = f"{r['eps_x_mm_mrad']:13.4f}" if not np.isnan(r['eps_x_mm_mrad']) else "      (in sol)"
            print(f"  {name:>10s}  {r['z_mm']:7.0f}  {r['sigma_x_um']:9.1f}  "
                  f"{r['sigma_y_um']:9.1f}  {r['sigma_t_fs']:9.0f}  {eps_str}  "
                  f"{r['sigma_delta_e3']:9.2f}  {r['Lz_mean']:10.2e}")

print(f"\n{'='*75}")
print(f"  UED Beamline — Phase 3B Diagnostics")
print(f"{'='*75}")
print(f"  E={E_keV:.0f} keV  Q={b['charge_fC']:.0f} fC  ε_n={epsilon_n*1e6:.3f} mm·mrad")
print(f"  Solenoid: B={B_sol:.3f}T → k={k_sol:.2f} @ z={z_sol*1e3:.0f}mm")
print(f"  RF: V={V_RF*1e-3:.0f}kV  φ={phi_RF:.2f}rad  f={f_RF*1e-9:.3f}GHz @ z={z_rf*1e3:.0f}mm")

print(f"\n  Running SC OFF …", end="", flush=True)
d_off = run_beamline(sc_enabled=False)
print(" done")

print(f"  Running SC ON  …", end="", flush=True)
d_on  = run_beamline(sc_enabled=True)
print(" done")

print_table(d_off, "SC OFF")
print(f"    Final: σ_x={d_off['sigx'][-1]:.1f} μm  σ_t={d_off['sig_t_final']*1e15:.0f} fs  "
      f"ε_x={d_off['epsx'][-1]:.4f} mm·mrad")

print_table(d_on, "SC ON")
print(f"    Final: σ_x={d_on['sigx'][-1]:.1f} μm  σ_t={d_on['sig_t_final']*1e15:.0f} fs  "
      f"ε_x={d_on['epsx'][-1]:.4f} mm·mrad")

# Phase 3B fix: verify RF z-delta chirp
print(f"\n  --- RF Chirp Verification ---")
for name, snap in d_off["snapshots"].items():
    tau = snap["tau"]; delta = snap["delta"]
    z_phys = -beta * c_SI * tau
    rho = np.corrcoef(z_phys, delta)[0, 1]
    chirp = np.polyfit(z_phys, delta, 1)[0]
    print(f"  {name:>12s}:  corr(z,δ)={rho:+.4f}  dδ/dz={chirp:+.4f} m⁻¹  σ_δ={np.std(delta)*1e3:.2f}e-3")

# ═══════════════════════════════════════════════════════════
#  Phase 3B diagnostics plots
# ═══════════════════════════════════════════════════════════

_np = min(5000, N_part)
_rng = np.random.default_rng(42)

def _subsample(arr, n=_np):
    return arr[_rng.choice(len(arr), n, replace=False)] if len(arr) > n else arr

# ---- Fig A: σ comparison SC OFF vs ON ----
figA, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
figA.suptitle("UED Beamline — SC OFF vs ON", fontsize=13, fontweight="bold")

for d, c, lab in [(d_off, "b", "SC OFF"), (d_on, "r", "SC ON")]:
    ls = "-" if lab == "SC OFF" else "--"
    ax1.plot(d["z_arr"], d["sigx"], color=c, linestyle=ls, linewidth=1.5, label=f"{lab} σ_x")
    ax2.plot(d["z_arr"], d["sigz"], color=c, linestyle=ls, linewidth=1.5, label=f"{lab} σ_z")

for ax in [ax1, ax2]:
    ax.axvspan(z_sol*1e3, (z_sol+L_sol)*1e3, alpha=0.10, color="blue")
    ax.axvline(z_rf*1e3, color="red", linestyle="--", alpha=0.4)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.25)
ax1.set_ylabel(r"$\sigma_x$  [$\mu$m]")
ax2.set_xlabel("z  [mm]"); ax2.set_ylabel(r"$\sigma_z$  [$\mu$m]")
figA.tight_layout()
figA.savefig(os.path.join(_OUTDIR, "ued_phase3b_sc_comparison.png"), dpi=150)
print("  -> ued_phase3b_sc_comparison.png")
plt.close(figA)

# ---- Fig B: emittance comparison ----
figB, ax = plt.subplots(figsize=(11, 4.5))
for d, c, lab in [(d_off, "b", "SC OFF"), (d_on, "r", "SC ON")]:
    ls = "-" if lab == "SC OFF" else "--"
    ax.plot(d["z_arr"], d["epsx"], color=c, linestyle=ls, linewidth=1.5, label=f"{lab} ε_x")
    ax.plot(d["z_arr"], d["epsy"], color=c, linestyle=ls, linewidth=1.0, alpha=0.5, label=f"{lab} ε_y")
ax.axvspan(z_sol*1e3, (z_sol+L_sol)*1e3, alpha=0.10, color="blue")
ax.axvline(z_rf*1e3, color="red", linestyle="--", alpha=0.4)
ax.set_xlabel("z  [mm]"); ax.set_ylabel(r"$\varepsilon_{x,y}$  [mm$\cdot$mrad]")
ax.set_title("Emittance evolution", fontweight="bold")
ax.legend(fontsize=7, ncol=2); ax.grid(True, alpha=0.25)
figB.tight_layout()
figB.savefig(os.path.join(_OUTDIR, "ued_phase3b_emit_comparison.png"), dpi=150)
print("  -> ued_phase3b_emit_comparison.png")
plt.close(figB)

# ---- Fig C: longitudinal phase space (z, δ) ----
snap_names = ["cathode", "before_rf", "after_rf", "sample"]
snap_labels = ["Cathode (z=0)", "Before RF (z=399 mm)", "After RF (z=422 mm)", "Sample (z=777 mm)"]
figC, axes = plt.subplots(1, 4, figsize=(16, 4))
figC.suptitle("Longitudinal Phase Space  (SC OFF)", fontsize=12, fontweight="bold")
for ax, name, lab in zip(axes, snap_names, snap_labels):
    s = d_off["snapshots"].get(name)
    if s is None: continue
    tau_s = _subsample(s["tau"])
    delta_s = _subsample(s["delta"])
    z_phys = -beta * c_SI * tau_s * 1e6
    ax.scatter(z_phys, delta_s * 1e3, s=0.3, alpha=0.4, c="steelblue")
    ax.set_title(lab, fontsize=9)
    ax.set_xlabel(r"$z$ [$\mu$m]"); ax.set_ylabel(r"$\delta$ [$10^{-3}$]")
    ax.grid(True, alpha=0.25)
figC.tight_layout()
figC.savefig(os.path.join(_OUTDIR, "ued_phase3b_longitudinal.png"), dpi=150)
print("  -> ued_phase3b_longitudinal.png")
plt.close(figC)

# ---- Fig D: transverse phase space (x-x', y-y') ----
xy_snaps = ["cathode", "after_sol", "sample"]
xy_labels = ["Cathode (z=0)", "After Solenoid (z=160 mm)", "Sample (z=777 mm)"]
figD, axes = plt.subplots(2, 3, figsize=(16, 9))
figD.suptitle("Transverse Phase Space  (SC OFF)", fontsize=12, fontweight="bold")
for col, (name, lab) in enumerate(zip(xy_snaps, xy_labels)):
    s = d_off["snapshots"].get(name)
    if s is None: continue
    x_s  = _subsample(s["x"])  * 1e6
    xp_s = _subsample(s["xp"]) * 1e3
    y_s  = _subsample(s["y"])  * 1e6
    yp_s = _subsample(s["yp"]) * 1e3
    for row, (v, vp, vl) in enumerate([(x_s, xp_s, "x"), (y_s, yp_s, "y")]):
        ax = axes[row, col]
        ax.scatter(v, vp, s=0.3, alpha=0.4, c="steelblue")
        ax.set_title(f"{lab}\nσ_{vl}={np.std(v/1e6)*1e6 if row==0 else np.std(v/1e6)*1e6:.0f} μm", fontsize=9)
        ax.set_xlabel(f"${vl}$ [$\mu$m]"); ax.set_ylabel(f"${vl}'$ [mrad]")
        ax.grid(True, alpha=0.25)
figD.tight_layout()
figD.savefig(os.path.join(_OUTDIR, "ued_phase3b_transverse.png"), dpi=150)
print("  -> ued_phase3b_transverse.png")
plt.close(figD)

print("\n  Phase 3B diagnostics complete.")
