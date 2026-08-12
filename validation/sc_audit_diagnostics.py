#!/usr/bin/env python3
"""
SC audit diagnostics (v0.14) — READ-ONLY.  No core modification.

Runs three diagnostics with the installed OCELOT SpaceCharge:
  1. smoke    : pure drift, SC OFF vs ON, identical macro particles (config seed)
  2. charge   : Q = 0..1000 fC monotonicity scan (sample σx, σz, max σx)
  3. converge : particle/mesh/SC-step resolution scan

Also reads the ACTUAL SpaceCharge attributes (nmesh_xyz, step, random_seed)
from the live object — does not trust logs.

Usage: python3 validation/sc_audit_diagnostics.py
"""

import os
import sys
import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
for p in (_THIS_DIR, os.path.dirname(_THIS_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from shared.params import load_config, derived
from shared.ocelot_coords import set_px, set_py, set_p_oc

L_DRIFT = 0.5
DZ = 0.005


def track_drift_sc(cfg, q_fC, mesh, nparticles, sc_step=1, sc_on=True):
    """Pure-drift OCELOT tracking with optional SpaceCharge (diagnostic).

    Mirrors the production beam generation: config seed → x/y/tau,
    seed+1 → px/py/δ_p.  Returns (history dict, sc_proc or None).
    """
    from ocelot.cpbd.elements import Drift
    from ocelot.cpbd.magnetic_lattice import MagneticLattice
    from ocelot.cpbd.beam import generate_parray
    from ocelot.cpbd.navi import Navigator
    from ocelot.cpbd.track import tracking_step
    from ocelot.cpbd.sc import SpaceCharge

    d = derived(cfg)
    ib = cfg["initial_distribution"]
    rng_seed = int(cfg["random"]["seed"])
    np.random.seed(rng_seed)
    p = generate_parray(
        sigma_x=ib["spot_rms_um"] * 1e-6, sigma_y=ib["spot_rms_um"] * 1e-6,
        sigma_tau=ib["bunch_length_um"] * 1e-6 / d["beta"],
        energy=(cfg["beam"]["energy_keV"] + 511.0) * 1e-6,
        charge=q_fC * 1e-15, nparticles=nparticles)
    np.random.seed(rng_seed + 1)
    N = p.rparticles.shape[1]
    set_px(p, np.random.normal(0.0, d["sigma_xp"], N))
    set_py(p, np.random.normal(0.0, d["sigma_yp"], N))
    set_p_oc(p, d["beta"] * np.random.normal(0.0, ib["sigma_delta"], N))

    lat = MagneticLattice([Drift(l=L_DRIFT)], method={"global": "SecondTM"})
    lat.update_transfer_maps()
    sc = None
    if sc_on:
        sc = SpaceCharge(step=sc_step, nmesh_xyz=list(mesh))
    navi = Navigator(lat, unit_step=DZ)
    if sc is not None:
        navi.add_physics_proc(sc, lat.sequence[0], lat.sequence[-1])

    hist = {"z_mm": [], "sigma_x_um": [], "sigma_y_um": [], "sigma_z_um": [],
            "sigma_delta_e3": [], "eps_nx_mm_mrad": [], "eps_ny_mm_mrad": []}
    n_steps = int(L_DRIFT / DZ)
    # NOTE (v0.14.1 task 1): this script deliberately KEEPS the manual
    # counter replica — it is the "manual-scheduler SC characterization"
    # reference (v0.14 smoke/charge/convergence numbers).  Production
    # (validation/backend.py, GPT模拟/ued_beamline_v2.py) now uses the
    # OCELOT NATIVE scheduler (get_next_step); SC OFF paths unchanged.
    for _ in range(n_steps):
        tracking_step(lat, p, DZ, navi)
        if sc is not None:
            sc.counter -= 1
            if sc.counter <= 0:
                sc.z0 = navi.z0
                sc.apply(p, sc.step * DZ)
                sc.counter = sc.step
        x = p.x(); xp = p.px(); y = p.y(); yp = p.py()
        hist["z_mm"].append(navi.z0 * 1e3)
        hist["sigma_x_um"].append(np.std(x) * 1e6)
        hist["sigma_y_um"].append(np.std(y) * 1e6)
        hist["sigma_z_um"].append(np.std(p.tau()) * d["beta"] * 1e6)
        hist["sigma_delta_e3"].append(np.std(p.p()) / d["beta"] * 1e3)
        ex = np.sqrt(np.mean(x**2) * np.mean(xp**2) - np.mean(x * xp)**2)
        ey = np.sqrt(np.mean(y**2) * np.mean(yp**2) - np.mean(y * yp)**2)
        hist["eps_nx_mm_mrad"].append(ex * d["beta_gamma"] * 1e6)
        hist["eps_ny_mm_mrad"].append(ey * d["beta_gamma"] * 1e6)
    return hist, sc


def main():
    cfg = load_config()
    print("== SC audit diagnostics (read-only) ==")
    print(f"config: charge_fC={cfg['beam']['charge_fC']}, "
          f"n_particles={cfg['beam']['n_particles']}, "
          f"mesh={cfg['space_charge'].get('mesh')}, "
          f"step={cfg['space_charge'].get('step')}")

    # ── smoke: SC OFF vs ON @ 500 fC (clearly visible point) ──
    q_smoke = 500.0
    off, _ = track_drift_sc(cfg, q_smoke, [63, 63, 63], 50000, sc_on=False)
    on, sc = track_drift_sc(cfg, q_smoke, [63, 63, 63], 50000, sc_on=True)
    print("\n[smoke] SC OFF vs ON @ 500 fC, identical seed, pure drift 0.5 m")
    print(f"  actual SC attrs: nmesh_xyz={sc.nmesh_xyz}, step={sc.step}, "
          f"random_seed={sc.random_seed}, random_mesh={sc.random_mesh}")
    for k in ("sigma_x_um", "sigma_z_um", "eps_nx_mm_mrad"):
        v_off = np.array(off[k]); v_on = np.array(on[k])
        d_k = v_on[-1] - v_off[-1]
        print(f"  sample {k:18s}: OFF={v_off[-1]:10.3f}  ON={v_on[-1]:10.3f}  "
              f"Δ={d_k:+.3f} ({(d_k/abs(v_off[-1]))*100:+.2f}%)")
    dsx = (np.array(on["sigma_x_um"]) - np.array(off["sigma_x_um"]))
    print(f"  max Δσx over path = {np.max(np.abs(dsx)):.3f} um "
          f"(SC effect {'PRESENT' if np.max(np.abs(dsx)) > 0.1 else 'MISSING'})")

    # ── charge scan ──
    qs = [0.0, 10.0, 50.0, 100.0, 500.0, 1000.0]
    print("\n[charge scan] Q =", qs, "fC (mesh 63, 5e4, step 1)")
    rows = []
    for q in qs:
        h, _ = track_drift_sc(cfg, q, [63, 63, 63], 50000, sc_on=(q > 0))
        rows.append((q, h["sigma_x_um"][-1], h["sigma_z_um"][-1],
                     max(h["sigma_x_um"]), h["eps_nx_mm_mrad"][-1]))
        print(f"  Q={q:7.0f} fC  sample_σx={rows[-1][1]:9.3f}  "
              f"sample_σz={rows[-1][2]:9.3f}  max_σx={rows[-1][3]:9.3f}  "
              f"εnx={rows[-1][4]:.4f}")

    # ── convergence scan ──
    print("\n[convergence] pure drift, SC ON @ 500 fC")
    qc = 500.0
    print("  A. particle number (mesh 63, step 1):")
    pa = []
    for n in (1e4, 5e4, 1e5):
        h, _ = track_drift_sc(cfg, qc, [63, 63, 63], int(n))
        pa.append((int(n), h["sigma_x_um"][-1], h["sigma_z_um"][-1],
                   h["eps_nx_mm_mrad"][-1]))
        print(f"     N={int(n):6d}  σx={pa[-1][1]:9.3f}  σz={pa[-1][2]:9.3f}  "
              f"εnx={pa[-1][3]:.4f}")
    print("  B. mesh (5e4, step 1):")
    pm = []
    for m in (33, 63, 127):
        h, _ = track_drift_sc(cfg, qc, [m, m, m], 50000)
        pm.append((m, h["sigma_x_um"][-1], h["sigma_z_um"][-1],
                   h["eps_nx_mm_mrad"][-1]))
        print(f"     mesh={m:3d}³  σx={pm[-1][1]:9.3f}  σz={pm[-1][2]:9.3f}  "
              f"εnx={pm[-1][3]:.4f}")
    print("  C. SC step (5e4, mesh 63):")
    ps = []
    for s in (1, 5, 10):
        h, _ = track_drift_sc(cfg, qc, [63, 63, 63], 50000, sc_step=s)
        ps.append((s, h["sigma_x_um"][-1], h["sigma_z_um"][-1],
                   h["eps_nx_mm_mrad"][-1]))
        print(f"     step={s:2d}  σx={ps[-1][1]:9.3f}  σz={ps[-1][2]:9.3f}  "
              f"εnx={ps[-1][3]:.4f}")

    # ── plots ──
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    ax = axes[0, 0]
    z = np.array(off["z_mm"])
    ax.plot(z, off["sigma_x_um"], "b-", lw=1.5, label="SC OFF")
    ax.plot(z, on["sigma_x_um"], "r--", lw=1.5, label="SC ON (500 fC)")
    ax.set_xlabel("z [mm]"); ax.set_ylabel(r"$\sigma_x$ [$\mu$m]")
    ax.set_title("smoke: σx(z)"); ax.legend(); ax.grid(alpha=.3)
    ax = axes[0, 1]
    ax.plot(z, off["sigma_z_um"], "b-", lw=1.5, label="SC OFF")
    ax.plot(z, on["sigma_z_um"], "r--", lw=1.5, label="SC ON")
    ax.set_xlabel("z [mm]"); ax.set_ylabel(r"$\sigma_z$ [$\mu$m]")
    ax.set_title("smoke: σz(z)"); ax.legend(); ax.grid(alpha=.3)
    ax = axes[1, 0]
    ax.plot([r[0] for r in rows], [r[1] for r in rows], "o-", label="sample σx")
    ax.plot([r[0] for r in rows], [r[3] for r in rows], "s--", label="max σx")
    ax.set_xlabel("Q [fC]"); ax.set_ylabel(r"$\sigma_x$ [$\mu$m]")
    ax.set_title("charge scan"); ax.legend(); ax.grid(alpha=.3)
    ax = axes[1, 1]
    ax.plot([r[0] for r in pa], [r[1] for r in pa], "o-", label="σx vs N")
    ax.plot([r[0] for r in pm], [r[1] for r in pm], "s--", label="σx vs mesh")
    ax.set_xlabel("resolution"); ax.set_ylabel(r"$\sigma_x$ [$\mu$m]")
    ax.set_title("convergence (particles / mesh)"); ax.legend(); ax.grid(alpha=.3)
    fig.suptitle("v0.14 SC audit diagnostics (pure drift, OCELOT)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    out = os.path.join(_THIS_DIR, "reports", "SC_diagnostics.png")
    fig.savefig(out, dpi=130)
    print(f"\nplot -> {out}")


if __name__ == "__main__":
    main()
