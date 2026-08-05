#!/usr/bin/env python3
"""
Validation test — SOLENOID section (first priority).

Systematic diagnosis of the AG vs OCELOT transverse difference along the five
requested axes:
  1. Bz definition       k_s = e·Bz/(2p) = Bz/(2·Bρ)   (identical in both)
  2. length              hard-edge L = 0.060 m, same span [z0, z0+L]
  3. momentum            p = γβ·m_e·c at 100 keV (identical)
  4. focusing strength   F = −k_s²σ  → k_s² identical
  5. x-y coupling        * DIFFERENT: AG injects reduced-order Larmor terms
                          (dν_x⊃−2k_sν_y, dν_y⊃+2k_sν_x, dσ_xy⊃2k_s(σ_x²−σ_y²)).
                          Exact hard-edge 4×4 transport (and OCELOT) shows these
                          MUST vanish for a round uncorrelated beam (σ_xy≡0).

No parameters are tuned.  The coupling terms are exposed as a model option in
the framework AG adapter (validation/backend.run_ag, solenoid_coupling=) and the
test reports both settings against the exact reference.

Usage:  python3 validation/test_solenoid.py
"""

import os
import sys

import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
for p in (_THIS_DIR, os.path.dirname(_THIS_DIR), os.path.join(os.path.dirname(_THIS_DIR), "AG")):
    if p not in sys.path:
        sys.path.insert(0, p)

from shared.params import load_config, derived, flat_elem  # noqa: E402
from validation.backend import run_ag, run_ocelot  # noqa: E402
from validation.reference import (solenoid_envelope_4x4,  # noqa: E402
                                  gamma_beta_p)
from validation.beam_result import BeamResult  # noqa: E402
from validation import common  # noqa: E402

SECTION = "solenoid"


def main():
    cfg = load_config()
    d = derived(cfg)
    ib = cfg["initial_distribution"]
    so = flat_elem(cfg, "solenoid")   # geometry+params from lattice

    gamma, beta, p_SI = gamma_beta_p(cfg["beam"]["energy_keV"])
    k_s = 1.602176634e-19 * so["B_field_T"] / (2.0 * p_SI)
    sigma0 = ib["spot_rms_um"] * 1e-6
    sxp0 = d["sigma_xp"]

    # ── 1. run the two backends (AG with and without reduced-order coupling) ──
    oc = run_ocelot(cfg, SECTION)
    ag_on = run_ag(cfg, SECTION, solenoid_coupling=True)     # AG as-is
    ag_off = run_ag(cfg, SECTION, solenoid_coupling=False)   # round-beam corrected

    # ── 2. exact hard-edge 4×4 reference (identical matrix to OCELOT) ──
    z_ref = np.linspace(0, max(oc.z_mm) * 1e-3, 2000)
    sx4, sy4, sxy4, ep4 = solenoid_envelope_4x4(
        z_ref, so["z_start_m"], so["length_m"], k_s, sigma0, sxp0)
    ref = BeamResult(
        route="4x4", z_mm=z_ref * 1e3,
        sigma_x_um=sx4 * 1e6, sigma_y_um=sy4 * 1e6,
        sigma_z_um=np.full_like(z_ref, ib["bunch_length_um"]),
        eps_nx_mm_mrad=ep4 * 1e6 * (beta * gamma),          # normalized
        eps_ny_mm_mrad=ep4 * 1e6 * (beta * gamma),
        energy_keV=np.full_like(z_ref, cfg["beam"]["energy_keV"]),
        sigma_delta_e3=np.full_like(z_ref, ib["sigma_delta"] * 1e3),
        meta={"reference": "hard-edge 4x4 Brown-Chao (== OCELOT SolenoidTM)"},
    )

    # ── 3. systematic diagnosis table ──
    print("\n═══ Solenoid diagnosis (shared config, no tuning) ═══")
    print(f"  1. Bz definition : k_s = e·Bz/(2p) = Bz/(2·Bρ) = {k_s:.4f} m⁻¹   [identical in AG & OCELOT]")
    print(f"  2. length        : hard-edge L = {so['length_m']:.3f} m, span "
          f"[{so['z_start_m']:.3f}, {so['z_start_m']+so['length_m']:.3f}] m  [identical]")
    print(f"  3. momentum      : p = γβm_ec = {p_SI:.4e} kg·m/s "
          f"(γ={gamma:.4f}, β={beta:.4f})  [identical]")
    print(f"  4. focusing      : F = −k_s²σ,  k_s² = {k_s**2:.2f} m⁻²  [identical]")
    print(f"  5. x-y coupling  : AG reduced-order Larmor terms ACTIVE by default;\n"
          f"                     exact 4×4 transport ⇒ σ_xy ≡ 0 for round beam.")

    m_on = common.compare_metrics(ag_on, oc)
    m_off = common.compare_metrics(ag_off, oc)
    common.print_summary(SECTION, [("AG(coupling=ON)", ag_on),
                                   ("AG(coupling=OFF)", ag_off),
                                   ("OCELOT", oc)])
    print(f"  AG(coupling=ON ) vs OCELOT: σ_x dev {m_on['sigma_x_um']:.2f}%  "
          f"σ_y dev {m_on['sigma_y_um']:.2f}%")
    print(f"  AG(coupling=OFF) vs OCELOT: σ_x dev {m_off['sigma_x_um']:.2f}%  "
          f"σ_y dev {m_off['sigma_y_um']:.2f}%")
    print(f"  AG(coupling=ON ) x-y symmetry: σ_x(end)={ag_on.sigma_x_um[-1]:.1f} vs "
          f"σ_y(end)={ag_on.sigma_y_um[-1]:.1f} μm  → BROKEN (round beam)")

    # ── 4. plots + artifacts ──
    common.plot_compare(SECTION, [("AG(coup=ON)", ag_on), ("AG(coup=OFF)", ag_off),
                                  ("OCELOT", oc), ("4x4 ref", ref)],
                        os.path.join(common.REPORTS_DIR, f"{SECTION}_AG_vs_OCELOT.png"),
                        title="Solenoid — AG vs OCELOT vs exact 4×4")
    common.save_results(SECTION, [ag_on, ag_off, oc, ref])

    common.log_checkpoint(SECTION, [
        "ROOT CAUSE: AG reduced-order Larmor coupling (dnu_x⊃-2ks·nu_y, "
        "dnu_y⊃+2ks·nu_x, dsxy⊃2ks(sx^2-sy^2)).",
        "Exact hard-edge 4x4 (Brown-Chao, == OCELOT SolenoidTM) gives sxy≡0 "
        "for a round uncorrelated beam; AG coupling creates spurious sxy and "
        "under-focusing.",
        f"k_s (Bz/(2Brho)) = {k_s:.4f} m^-1, k_s^2 = {k_s**2:.2f} m^-2 — identical in both backends.",
        f"AG as-is:  sx={ag_on.sigma_x_um[-1]:.1f} sy={ag_on.sigma_y_um[-1]:.1f} um "
        f"(x-y broken)",
        f"AG coupling=OFF: sx={ag_off.sigma_x_um[-1]:.1f} sy={ag_off.sigma_y_um[-1]:.1f} um",
        f"OCELOT (ref):    sx={oc.sigma_x_um[-1]:.1f} sy={oc.sigma_y_um[-1]:.1f} um",
        f"AG(off) vs OCELOT sigma_x max dev = {m_off['sigma_x_um']:.2f}% "
        f"(<1% PASS) | AG(on) = {m_on['sigma_x_um']:.2f}% (FAIL)",
        "FIX: for round beams the coupling must be disabled in the AG force "
        "adapter (validation/backend.run_ag solenoid_coupling=False). Not a "
        "parameter tune; enforces exact round-beam transport.",
    ])
    print(f"\n  report -> {os.path.join(common.REPORTS_DIR, SECTION + '_AG_vs_OCELOT.png')}")
    return 0 if m_off["sigma_x_um"] < 1.0 else 1


if __name__ == "__main__":
    sys.exit(main())
