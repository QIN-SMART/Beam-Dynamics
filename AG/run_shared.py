#!/usr/bin/env python3
"""
AG-route shared driver.

Reads the SHARED beamline config (shared/beamline_config.yaml), translates it
into the route-native objects used by the EXISTING core code
(beam_dynamics_6d / external_forces / beamline_sim), runs the 6D envelope ODE
and writes the UNIFIED result file (shared/results/AG_results.json).

The core algorithm files are NOT modified. This is a thin adapter only.

Usage:  python3 AG/run_shared.py
"""

import os
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_THIS_DIR)
for p in (_REPO, _THIS_DIR):
    if p not in sys.path:
        sys.path.insert(0, p)

from shared.params import load_config, config_sha
from shared.output_schema import write_results, make_probe

from beam_dynamics_6d import (make_beam_100keV, ExtFieldRegion,
                              get_alpha_interpolators, propagate,
                              apply_rf_thin_lens, Beam6D)
from external_forces import build_all_external, build_external_force_func
from beamline_sim import compute_emittance

C_SI = 2.99792458e8
M_E_SI = 9.10938356e-31
E_SI = 1.602176634e-19

# Round-beam solenoid coupling correction (see validation/CHECKPOINTS.md):
# AG core injects reduced-order Larmor coupling (dnu_x⊃-2ks·nu_y, ...) which the
# exact hard-edge 4×4 transport (and OCELOT) shows MUST vanish for a round,
# uncorrelated beam (sigma_xy ≡ 0).  Disabled here for the round shared beam.
SOLENOID_COUPLING = False

# RF: standardized thin-lens model (same as validation/backend.py):
# longitudinal kick via apply_rf_thin_lens(H), H = eV·k·cosφ/(β²E₀);
# continuous chirp/acceleration disabled; RF transverse force kept.
def rf_thin_lens_H(cfg, rf_elem):
    p = rf_elem["parameters"]
    k_rf = 2.0 * np.pi * p["frequency_GHz"] * 1e9 / C_SI
    V = p["voltage_kV"] * 1e3
    gamma = 1.0 + cfg["beam"]["energy_keV"] / 511.0
    beta = np.sqrt(1.0 - 1.0 / gamma**2)
    return E_SI * V * k_rf * np.cos(p["phase_rad"]) / (
        beta**2 * gamma * M_E_SI * C_SI**2)


def build_regions(cfg):
    """Translate lattice.elements -> ExtFieldRegion list (lattice only).

    Only active elements produce regions: solenoid -> Bz field, rf_cavity ->
    standing-wave RF field with peak E_rf = V_RF / L_cav (transverse force
    only; the longitudinal chirp is the thin lens).  Drift / cathode / sample
    are field-free and need no region.
    """
    regions = []
    for e in cfg["lattice"]["elements"]:
        etype, z0, L = e["type"], e["z_start"], e["length"]
        if L <= 0:
            continue
        if etype == "solenoid":
            regions.append(ExtFieldRegion(z0, z0 + L, "solenoid",
                                          {"Bz": e["parameters"]["B_field_T"],
                                           "dBz_dz": 0.0}))
        elif etype == "rf_cavity":
            p = e["parameters"]
            k_rf = 2.0 * np.pi * p["frequency_GHz"] * 1e9 / C_SI
            E_rf = p["voltage_kV"] * 1e3 / L          # peak field [V/m]
            regions.append(ExtFieldRegion(z0, z0 + L, "rf",
                                          {"E_rf": E_rf, "k_rf": k_rf,
                                           "dE_dz": 0.0, "dE_cdt": k_rf * E_rf}))
    return regions


def main():
    cfg = load_config()
    b  = cfg["beam"]
    ib = cfg["initial_distribution"]
    out = cfg["output"]

    get_alpha_interpolators()

    beam0 = make_beam_100keV(
        Ne=b["n_particles"],
        beamK_eV=b["energy_keV"] * 1e3,
        sigma_x0_um=ib["spot_rms_um"],
        sigma_y0_um=ib["spot_rms_um"],
        sigma_z0_um=ib["bunch_length_um"],
        sigma_delta=ib["sigma_delta"],
        eps_nx_um=ib["epsilon_n_mm_mrad"],
        eps_ny_um=ib["epsilon_n_mm_mrad"],
        eps_nz_um=ib["epsilon_nz_mm_mrad"],
    )

    regions = build_regions(cfg)
    z_total = max(e["z_start"] + e["length"] for e in cfg["lattice"]["elements"])
    sc_enabled = bool(cfg["space_charge"]["enabled"])

    # core external fields (solenoid only); RF kept ONLY for its transverse
    # force — the continuous chirp/acceleration is replaced by the thin lens.
    regions_core = [r for r in regions if r.ftype != "rf"]
    rf_regions = [r for r in regions if r.ftype == "rf"]
    ef_core, gp_core, ch_core = build_all_external(regions_core)

    ef = ef_core
    if rf_regions:
        ef_rf = build_external_force_func(rf_regions)
        ef = (lambda z, beam, _core=ef_core, _rf=ef_rf:
              tuple(a + b for a, b in zip(_core(z, beam), _rf(z, beam))))

    if not SOLENOID_COUPLING:
        base = ef
        ef = (lambda z, beam, _base=base:
              _base(z, beam)[:3] + (0.0, 0.0, 0.0))

    beam_k = beam0.copy()
    if not sc_enabled:
        beam_k.Ne = 0.0
    sc_model = "ellipsoid" if sc_enabled else "gaussian"

    rf_elems = [e for e in cfg["lattice"]["elements"]
                if e["type"] == "rf_cavity" and e["length"] > 0]
    if rf_elems:
        # thin-lens RF (multi-instance): propagate to each cavity entrance,
        # apply H·z chirp, continue to the next one.
        z_parts, st_parts = [], []
        z_cur = 0.0
        beam_cur = beam_k
        n_per = 2000 // (len(rf_elems) + 1)
        for e in rf_elems:
            z_rf = e["z_start"]
            if z_rf > z_cur:
                z1, st1 = propagate(beam_cur, (z_cur, z_rf), n_points=n_per,
                                    external_force_func=ef, gamma_prime_func=gp_core,
                                    rf_chirp_func=ch_core, sc_model=sc_model)
                z_parts.append(z1); st_parts.append(st1)
                beam_cur = Beam6D.from_state(st1[-1, :11], st1[-1, 11], beam_k.Ne,
                                             beam_k.eps_nx, beam_k.eps_ny, beam_k.eps_nz)
            beam_cur = apply_rf_thin_lens(beam_cur, rf_thin_lens_H(cfg, e))
            z_cur = z_rf
        if z_cur < z_total:
            z2, st2 = propagate(beam_cur, (z_cur, z_total), n_points=n_per,
                                external_force_func=ef, gamma_prime_func=gp_core,
                                rf_chirp_func=ch_core, sc_model=sc_model)
            z_parts.append(z2); st_parts.append(st2)
        z_arr = np.concatenate([zz if i == 0 else zz[1:]
                                for i, zz in enumerate(z_parts)])
        st = np.concatenate([ss if i == 0 else ss[1:]
                             for i, ss in enumerate(st_parts)])
    else:
        z_arr, st = propagate(beam_k, (0.0, z_total), n_points=2000,
                              external_force_func=ef, gamma_prime_func=gp_core,
                              rf_chirp_func=ch_core, sc_model=sc_model)

    em = compute_emittance(st, beam0.eps_nx, beam0.eps_ny, beam0.eps_nz)

    def at(z_p, arr):
        return float(np.interp(z_p, z_arr, arr))

    probes = []
    for z_mm in out["z_diagnostics_mm"]:
        z_p = z_mm * 1e-3
        probes.append(make_probe(
            z_mm=z_mm,
            sigma_x_um=at(z_p, st[:, 0]) * 1e6,
            sigma_y_um=at(z_p, st[:, 1]) * 1e6,
            sigma_z_um=at(z_p, st[:, 2]) * 1e6,
            sigma_delta_e3=at(z_p, st[:, 4]) * 1e3,
            eps_nx_mm_mrad=at(z_p, em["eps_n_x"]) * 1e6,
            eps_ny_mm_mrad=at(z_p, em["eps_n_y"]) * 1e6,
        ))

    history = {
        "z_mm": (z_arr * 1e3).tolist(),
        "sigma_x_um": (st[:, 0] * 1e6).tolist(),
        "sigma_y_um": (st[:, 1] * 1e6).tolist(),
        "sigma_z_um": (st[:, 2] * 1e6).tolist(),
        "sigma_delta_e3": (st[:, 4] * 1e3).tolist(),
        "eps_nx_mm_mrad": (em["eps_n_x"] * 1e6).tolist(),
        "eps_ny_mm_mrad": (em["eps_n_y"] * 1e6).tolist(),
    }

    meta = {
        "model": "Kelisani 6D envelope ODE",
        "sc_model": "ellipsoid" if sc_enabled else "none (Ne=0)",
        "rf": "thin-lens (apply_rf_thin_lens H) + transverse kick",
        "solenoid_coupling": SOLENOID_COUPLING,
        "solenoid_coupling_note": "round-beam correction, see validation/CHECKPOINTS.md",
        "lattice_regions": [(r.z_start, r.z_end, r.ftype) for r in regions],
    }
    path = write_results("AG", probes, history, config_sha(cfg),
                         sc_enabled, meta=meta)
    print(f"AG done -> {path}")
    print(f"  final: sigma_x={history['sigma_x_um'][-1]:.1f} um, "
          f"sigma_z={history['sigma_z_um'][-1]:.1f} um, "
          f"eps_nx={history['eps_nx_mm_mrad'][-1]:.4f} mm.mrad")
    return path


if __name__ == "__main__":
    main()
