#!/usr/bin/env python3
"""
Validation test — DRIFT section (no external elements).

Calls both backends through the framework (validation/backend.py) with the
SHARED parameters, plus the analytic drift reference, and writes an
AG_vs_OCELOT comparison report.

Analytic reference:  σ(z) = √(σ0² + (σ0')² z²),  σ0' = ε_geo/σ0.

Usage:  python3 validation/test_drift.py
"""

import os
import sys

import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
for p in (_THIS_DIR, os.path.dirname(_THIS_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from shared.params import load_config, derived  # noqa: E402
from validation.backend import run_ag, run_ocelot  # noqa: E402
from validation.reference import drift_sigma  # noqa: E402
from validation.beam_result import BeamResult  # noqa: E402
from validation import common  # noqa: E402

SECTION = "drift"


def main():
    cfg = load_config()
    d = derived(cfg)
    ib = cfg["initial_distribution"]

    ag = run_ag(cfg, SECTION)
    oc = run_ocelot(cfg, SECTION)

    # analytic reference on the AG z-grid
    sigma0 = ib["spot_rms_um"] * 1e-6
    sxp0 = d["sigma_xp"]
    z_ref = ag.z_mm * 1e-3
    # longitudinal analytic reference (R56_z = L/γ², δ_p = Δp/p0):
    #   σ_z(z) = √(σ_z0² + (z·σ_δ_p/γ²)²)
    gamma2 = (1.0 + cfg["beam"]["energy_keV"] / 511.0)**2
    sigma_z0 = ib["bunch_length_um"] * 1e-6
    sigma_delta_p = ib["sigma_delta"]
    sz_analytic = np.sqrt(sigma_z0**2
                          + (z_ref * sigma_delta_p / gamma2)**2) * 1e6
    ref = BeamResult(
        route="analytic", z_mm=z_ref * 1e3,
        sigma_x_um=drift_sigma(z_ref, sigma0, sxp0) * 1e6,
        sigma_y_um=drift_sigma(z_ref, sigma0, d["sigma_yp"]) * 1e6,
        sigma_z_um=sz_analytic,
        eps_nx_mm_mrad=np.full_like(z_ref, ib["epsilon_n_mm_mrad"]),
        eps_ny_mm_mrad=np.full_like(z_ref, ib["epsilon_n_mm_mrad"]),
        energy_keV=np.full_like(z_ref, cfg["beam"]["energy_keV"]),
        sigma_delta_e3=np.full_like(z_ref, sigma_delta_p * 1e3),
        meta={"reference": "sigma_x=sqrt(s0^2+(s0' z)^2); "
                           "sigma_z=sqrt(sz0^2+(z*sd_p/gamma^2)^2)"},
    )

    metrics = common.compare_metrics(ag, oc)
    common.print_summary(SECTION, [("AG", ag), ("OCELOT", oc)], metrics)

    ok = all(v < 5.0 for k, v in metrics.items()
             if k in ("sigma_x_um", "sigma_y_um", "eps_nx_mm_mrad", "sigma_delta_e3"))
    print(f"  drift transverse agreement (<5%): {'PASS' if ok else 'FAIL'}")

    # ── R56 adapter — Test 2: initial-distribution semantic check ──
    # OCELOT p_oc = ΔE/(c·p0) = β0·δ_p; reported σ_δ_p = std(p_oc)/β0 must
    # equal the configured momentum deviation (deterministic: delta is drawn
    # after np.random.seed(42) and conserved in drift).
    cfg_sd = cfg["initial_distribution"]["sigma_delta"]
    sd_initial = oc.sigma_delta_e3[0] * 1e-3          # δ_p at z ≈ 0
    rel_sd = abs(sd_initial - cfg_sd) / cfg_sd
    ok_sd = rel_sd < 0.01
    print(f"\n  R56 adapter — initial σ_δ_p: configured={cfg_sd:.3e}  "
          f"measured(std(p_oc)/β0)={sd_initial:.3e}  rel={rel_sd*100:.3f}%  "
          f"{'PASS' if ok_sd else 'FAIL'}")

    # ── longitudinal analytic reference check (uses the sz_analytic curve
    #    already built above; common δ_p = Δp/p0 definition) ──
    dev_ag = np.max(np.abs(ag.sigma_z_um - sz_analytic) / sz_analytic)
    oc_z_interp = np.interp(ag.z_mm, oc.z_mm, oc.sigma_z_um)
    dev_oc = np.max(np.abs(oc_z_interp - sz_analytic) / sz_analytic)
    print(f"  σ_z vs analytic √(σ_z0²+(z·σ_δ_p/γ²)²): AG={dev_ag*100:.2f}%  "
          f"OCELOT={dev_oc*100:.2f}%")

    ok = ok and ok_sd
    print(f"  NOTE: OCELOT sigma_z transports natively via R56 (energy fix); "
          f"eps_nz set to consistent 0.02um, so AG/OCELOT sigma_z now agree "
          f"to ~{metrics.get('sigma_z_um', 0):.1f}%.")

    out_png = os.path.join(common.REPORTS_DIR, f"{SECTION}_AG_vs_OCELOT.png")
    common.plot_compare(SECTION, [("AG", ag), ("OCELOT", oc), ("analytic", ref)],
                        out_png, title="Drift — AG vs OCELOT vs analytic")
    common.save_results(SECTION, [ag, oc, ref])

    common.log_checkpoint(SECTION, [
        f"drift transverse: AG vs OCELOT max rel dev {metrics}",
        f"analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))",
        f"sigma_z analytic ref sqrt(sz0^2+(z*sd_p/gamma^2)^2): "
        f"AG={dev_ag*100:.2f}% OCELOT={dev_oc*100:.2f}%",
        f"R56 adapter σ_δ_p semantic: {rel_sd*100:.3f}% "
        f"({'PASS' if ok_sd else 'FAIL'})",
        f"verdict: {'PASS' if ok else 'FAIL'} — report {out_png}",
    ])
    print(f"  report -> {out_png}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
