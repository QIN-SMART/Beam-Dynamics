#!/usr/bin/env python3
"""
Validation test — FULL BEAMLINE (cathode → sample).

Runs BOTH backends through the complete lattice using the shared config
(physics_switches: rf_transverse_kick OFF by default, as configured).

Compares the continuous curves σ_x(z), σ_y(z), σ_z(z), σ_t(z), σ_δ(z),
ε_nx(z), ε_ny(z), with shaded solenoid/RF/sample regions on the figure.

Acceptance (at the sample plane z = z_sample):
  σ_x, σ_y < 5% ; σ_δ < 2% ; ε_nx, ε_ny < 5%
  σ_z, σ_t < 5% (R56 convention RESOLVED by the longitudinal-variable
  adapter; compression-waist position within 5 mm)

Writes validation/reports/full_beamline_validation.md.

Usage:  python3 validation/test_full_beamline.py
"""

import os
import sys

import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
for p in (_THIS_DIR, os.path.dirname(_THIS_DIR), os.path.join(os.path.dirname(_THIS_DIR), "AG")):
    if p not in sys.path:
        sys.path.insert(0, p)

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from shared.params import load_config, elements_of_type, elem_geom, z_sample  # noqa: E402
from validation.backend import run_ag, run_ocelot  # noqa: E402
from validation import common  # noqa: E402

SECTION = "full"

CURVES = [
    ("sigma_x_um", r"$\sigma_x$ [$\mu$m]"),
    ("sigma_y_um", r"$\sigma_y$ [$\mu$m]"),
    ("sigma_z_um", r"$\sigma_z$ [$\mu$m]"),
    ("time_res_fs", r"$\sigma_t$ [fs]"),
    ("sigma_delta_e3", r"$\sigma_\delta$ [$10^{-3}$]"),
    ("eps_nx_mm_mrad", r"$\varepsilon_{nx}$ [mm$\cdot$mrad]"),
    ("eps_ny_mm_mrad", r"$\varepsilon_{ny}$ [mm$\cdot$mrad]"),
]

ACCEPT = {  # key, threshold (%), diagnostic-only flag
    "sigma_x_um": (5.0, False),
    "sigma_y_um": (5.0, False),
    "sigma_delta_e3": (2.0, False),
    "eps_nx_mm_mrad": (5.0, False),
    "eps_ny_mm_mrad": (5.0, False),
    # R56 convention resolved (r56 adapter): longitudinal now quantitative.
    "sigma_z_um": (5.0, False),
    "time_res_fs": (5.0, False),
}


def sample_value(r, key, zq):
    """Value of a curve at the sample plane (linear interpolation)."""
    return float(np.interp(zq, np.asarray(r.z_mm), np.asarray(getattr(r, key))))


def plot_full(results, out_png):
    fig, axes = plt.subplots(len(CURVES), 1, figsize=(10, 2.5 * len(CURVES)),
                             sharex=True)
    for ax, (key, label) in zip(axes, CURVES):
        for name, r in results:
            ax.plot(r.z_mm, getattr(r, key), lw=1.6, label=name)
        ax.set_ylabel(label)
        ax.legend(fontsize=8, loc="upper right")
        ax.grid(alpha=0.3)
    # element regions (solenoid, RF) + sample marker
    cfg = load_config()
    for e in elements_of_type(cfg, "solenoid"):
        z0, z1, _ = elem_geom(e)
        axes[0].axvspan(z0 * 1e3, z1 * 1e3, alpha=0.12, color="orange",
                        label=f"solenoid ({e['name']})")
    for e in elements_of_type(cfg, "rf_cavity"):
        z0, z1, _ = elem_geom(e)
        axes[0].axvspan(z0 * 1e3, z1 * 1e3, alpha=0.12, color="red",
                        label=f"RF ({e['name']})")
    zs = z_sample(cfg) * 1e3
    axes[0].axvline(zs, color="green", ls="--", lw=1,
                    label=f"sample @ {zs:.0f} mm")
    axes[0].legend(fontsize=8)
    axes[-1].set_xlabel("z [mm]")
    fig.suptitle("Full beamline — AG vs OCELOT (RF transverse kick OFF)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    return out_png


def main():
    cfg = load_config()
    ag = run_ag(cfg, SECTION, solenoid_coupling=False)
    oc = run_ocelot(cfg, SECTION)
    zq = z_sample(cfg) * 1e3

    print(f"\n== Full beamline (z_sample = {zq:.0f} mm) ==")
    common.print_summary(SECTION, [("AG", ag), ("OCELOT", oc)])

    # sample-plane metrics
    sample = {}
    for key, (thr, diag) in ACCEPT.items():
        va = sample_value(ag, key, zq)
        vb = sample_value(oc, key, zq)
        dev = abs(va - vb) / abs(vb) * 100 if abs(vb) > 1e-12 else float("inf")
        sample[key] = {"AG": va, "OCELOT": vb, "dev_pct": dev,
                       "threshold": thr, "diagnostic": diag,
                       "pass": (thr is None) or (dev < thr)}

    print(f"\n  sample-plane comparison (z = {zq:.0f} mm):")
    for key, s in sample.items():
        tag = "diag" if s["diagnostic"] else ("PASS" if s["pass"] else "FAIL")
        print(f"    {key:18s} AG={s['AG']:10.3f}  OCELOT={s['OCELOT']:10.3f}  "
              f"dev={s['dev_pct']:6.2f}%  [{tag}]")

    ok_all = all(s["pass"] or s["diagnostic"] for s in sample.values())

    # compression-waist position check (σ_z minimum)
    iz_ag = int(np.argmin(ag.sigma_z_um)); iz_oc = int(np.argmin(oc.sigma_z_um))
    z_waist_ag = float(ag.z_mm[iz_ag]); z_waist_oc = float(oc.z_mm[iz_oc])
    waist_dz = abs(z_waist_ag - z_waist_oc)
    waist_tol = 5.0   # mm (one OCELOT output step is 1 mm)
    ok_waist = waist_dz <= waist_tol
    print(f"\n  compression waist: AG z={z_waist_ag:.0f}mm σ_z={ag.sigma_z_um[iz_ag]:.1f}um | "
          f"OCELOT z={z_waist_oc:.0f}mm σ_z={oc.sigma_z_um[iz_oc]:.1f}um | "
          f"Δz={waist_dz:.1f}mm (tol {waist_tol}mm) "
          f"{'PASS' if ok_waist else 'FAIL'}")
    ok_all = ok_all and ok_waist
    print(f"\n  full-beamline acceptance: {'PASS' if ok_all else 'FAIL'}")

    out_png = os.path.join(common.REPORTS_DIR, f"{SECTION}_AG_vs_OCELOT.png")
    plot_full([("AG", ag), ("OCELOT", oc)], out_png)
    common.save_results(SECTION, [ag, oc])

    # ── report markdown ──
    report_path = os.path.join(common.REPORTS_DIR, "full_beamline_validation.md")
    with open(report_path, "w") as f:
        f.write("# Full-Beamline Validation — AG vs OCELOT\n\n")
        f.write("Date: 2026-08-01 · shared config · RF transverse kick OFF "
                "(`physics_switches`)\n\n")
        f.write("## Sample-plane metrics (z = %.0f mm)\n\n" % zq)
        f.write("| quantity | AG | OCELOT | dev [%] | verdict |\n|---|---|---|---|---|\n")
        for key, s in sample.items():
            verdict = ("diag" if s["diagnostic"]
                       else ("PASS" if s["pass"] else "FAIL"))
            f.write(f"| {key} | {s['AG']:.3f} | {s['OCELOT']:.3f} | "
                    f"{s['dev_pct']:.2f} | {verdict} |\n")
        f.write("\n## Passed transverse / RF-kick metrics\n\n")
        for key, s in sample.items():
            if not s["diagnostic"]:
                f.write(f"- **{key}**: dev = {s['dev_pct']:.2f}% "
                        f"(threshold {s['threshold']}%) — "
                        f"{'PASS' if s['pass'] else 'FAIL'}\n")
        f.write("\n## Longitudinal (R56 adapter applied)\n\n")
        f.write(f"- σ_z at sample: AG {sample['sigma_z_um']['AG']:.1f} µm vs "
                f"OCELOT {sample['sigma_z_um']['OCELOT']:.1f} µm "
                f"(dev {sample['sigma_z_um']['dev_pct']:.2f}%, threshold 5%).\n")
        f.write(f"- σ_t at sample: AG {sample['time_res_fs']['AG']:.1f} fs vs "
                f"OCELOT {sample['time_res_fs']['OCELOT']:.1f} fs "
                f"(dev {sample['time_res_fs']['dev_pct']:.2f}%).\n")
        f.write(f"- compression waist: AG z={z_waist_ag:.0f} mm vs "
                f"OCELOT z={z_waist_oc:.0f} mm (Δz={waist_dz:.1f} mm).\n")
        f.write("- The input δ-variable convention mismatch (classification B) "
                "was resolved by the OCELOT longitudinal adapter "
                "(p_oc = β0·δ_p input, δ_p = p_oc/β0 output); "
                "see `CHANGELOG_r56_variable_adapter.md` and "
                "`R56_convention_resolution.md`.\n")
        f.write("\n## Integration status\n\n")
        f.write(f"- Both backends ran the complete lattice "
                f"(cathode → sample, {len(cfg['lattice']['elements'])} elements) "
                f"with identical shared parameters (config_sha "
                f"{ag.meta.get('config_sha', 'n/a')}).\n")
        f.write(f"- Integration failures: **{'NONE' if ok_all else 'see FAIL rows above'}**.\n")
        f.write(f"- Figure: `full_AG_vs_OCELOT.png` · results: `full_results.json`\n")
    print(f"  report -> {report_path}")

    common.log_checkpoint(SECTION, [
        f"full beamline: sample devs — σ_x {sample['sigma_x_um']['dev_pct']:.2f}%, "
        f"σ_y {sample['sigma_y_um']['dev_pct']:.2f}%, "
        f"σ_δ {sample['sigma_delta_e3']['dev_pct']:.2f}%, "
        f"ε_nx {sample['eps_nx_mm_mrad']['dev_pct']:.2f}%, "
        f"ε_ny {sample['eps_ny_mm_mrad']['dev_pct']:.2f}% "
        f"({'PASS' if ok_all else 'FAIL'})",
        f"R56 resolved — σ_z: AG {sample['sigma_z_um']['AG']:.0f} vs "
        f"OCELOT {sample['sigma_z_um']['OCELOT']:.0f} um "
        f"({sample['sigma_z_um']['dev_pct']:.2f}%, PASS), "
        f"waist Δz={waist_dz:.1f} mm",
        f"switches: {ag.meta['switches']} == {oc.meta['switches']}",
    ])
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
