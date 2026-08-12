#!/usr/bin/env python3
"""
UED Beamline Integration — GPT main route (single-source lattice).

The ONLY geometry/parameter source is shared/beamline_config.yaml →
lattice.elements (same as the validation OCELOT route).  Step semantics:

  step 1: only drifts active; all other length-bearing elements are kept as
          plain Drifts so the total length and sample position are unchanged
  step 2: drift + ALL solenoids
  step 3: drift + ALL solenoids + ALL RF cavities
  step 4: same lattice as step 3 + SpaceCharge process (per CLI/SC flags)

RF kicks: one per rf_cavity instance, at its own z_start, in lattice order;
the shared physics_switches.rf_transverse_kick gates the transverse part.
OCELOT native tau/R56 unmodified; p_oc = β0·δ_p conversion retained.

Usage:
  python3 ued_beamline_v2.py --step 1    # drift only (inactive as drifts)
  python3 ued_beamline_v2.py --step 2    # + solenoids
  python3 ued_beamline_v2.py --step 3    # + RF
  python3 ued_beamline_v2.py --step 4    # + space charge (full)
"""

import sys
import os

import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_THIS_DIR)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)
from shared.params import load_config, config_sha, _lattice_elements  # noqa: E402
from shared.constants import MEC2_KEV, E_SI, M_E_SI, C_SI  # noqa: E402
from shared.beam_physics import BeamReference  # noqa: E402
from shared.ocelot_coords import (add_p_oc, add_px, add_py,  # noqa: E402
                                  set_px, set_py, set_p_oc)
from shared.output_schema import write_results, make_probe  # noqa: E402

print("加载 OCELOT …", flush=True)
import ocelot  # noqa: E402
from ocelot.cpbd.elements import Drift, Solenoid  # noqa: E402
from ocelot.cpbd.magnetic_lattice import MagneticLattice  # noqa: E402
from ocelot.cpbd.beam import generate_parray  # noqa: E402
from ocelot.cpbd.navi import Navigator  # noqa: E402
from ocelot.cpbd.track import tracking_step  # noqa: E402
try:
    from ocelot.cpbd.sc import SpaceCharge  # noqa: E402
    _HAS_SC = True
except ImportError:
    _HAS_SC = False

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# ═══════════════════════════════════════════════════════════
#  load ALL parameters from SHARED YAML (single source of truth)
# ═══════════════════════════════════════════════════════════

cfg = load_config()

b   = cfg["beam"]
ib  = cfg["initial_distribution"]
sc  = cfg["space_charge"]
out = cfg["output"]
switches = cfg.get("physics_switches", {})

# beam
E_keV       = b["energy_keV"]
Q_C         = b["charge_fC"] * 1e-15
N_part      = b["n_particles"]

# initial distribution
spot_rms    = ib["spot_rms_um"] * 1e-6
sig_z0      = ib["bunch_length_um"] * 1e-6
epsilon_n   = ib["epsilon_n_mm_mrad"] * 1e-6
sigma_delta = ib["sigma_delta"]

# output
z_probes_mm = out["z_diagnostics_mm"]
dz_track    = out["step_size_m"]

# ═══════════════════════════════════════════════════════════
#  relativistic (from config energies only)
# ═══════════════════════════════════════════════════════════

_br = BeamReference.from_energy_keV(E_keV)   # single γ/β/p0 source (v0.13)
gamma      = _br.gamma
beta       = _br.beta
beta_gamma = _br.beta_gamma
p_SI       = _br.p0

epsilon_geom = epsilon_n / beta_gamma
sigma_xp     = epsilon_geom / spot_rms
sigma_yp     = epsilon_geom / spot_rms

E_total_eV = gamma * MEC2_KEV * 1e3

# ═══════════════════════════════════════════════════════════
#  generic lattice builder — lattice.elements is the ONLY geometry source
# ═══════════════════════════════════════════════════════════

def build_lattice_from_shared(cfg, active_types, keep_zero_markers=False):
    """Build an OCELOT lattice from lattice.elements (single source).

    active_types : set of element types that keep their real physics
      ({"drift"}, {"drift","solenoid"}, {"drift","solenoid","rf_cavity"}).
      Every other length-bearing element is kept as a plain Drift of the
      SAME length, so the total length and sample position are preserved.

    keep_zero_markers : with SC ON (v0.14.1 task 1) the zero-length
      cathode/sample markers are kept in the runtime sequence so the
      SpaceCharge PhysProc can be attached as cathode → sample (coverage
      ends at the sample position, not at the start of the last non-zero
      element).  SC OFF keeps the previous sequence (markers skipped).

    Returns (lattice, rf_elems) where rf_elems is the ordered list of ACTIVE
    rf_cavity elements (one analytic kick per instance, in lattice order).
    """
    elems = []
    rf_elems = []
    for e in _lattice_elements(cfg):
        etype, L = e["type"], e["length"]
        if L <= 0:
            if keep_zero_markers and etype in ("cathode", "sample"):
                elems.append(Drift(l=0.0, eid=e["name"] + "_MARKER"))
            continue                       # cathode / sample markers
        if etype == "solenoid" and "solenoid" in active_types:
            B = e["parameters"]["B_field_T"]
            k = E_SI * B / (2.0 * p_SI)
            elems.append(Solenoid(l=L, k=k, eid=e["name"]))
        elif etype == "rf_cavity" and "rf_cavity" in active_types:
            elems.append(Drift(l=L, eid=e["name"] + "_BODY"))
            rf_elems.append(e)             # analytic kick applied at z_start
        else:
            elems.append(Drift(l=L, eid=e["name"]))   # drift or inactive
    lat = MagneticLattice(elems, method={"global": "SecondTM"})
    lat.update_transfer_maps()
    return lat, rf_elems


STEP_ACTIVE = {
    1: {"drift"},
    2: {"drift", "solenoid"},
    3: {"drift", "solenoid", "rf_cavity"},
    4: {"drift", "solenoid", "rf_cavity"},
}


def apply_rf_kick(parray, rf_elem, rf_transverse=False):
    """Standardized thin-lens RF kick for ONE rf_cavity instance.

    Longitudinal: dδ_p = (V/(β²E_total))·sin(φ + k·z_phys), stored as
    p_oc = β0·δ_p (OCELOT p() = ΔE/(c·p0), NOT Δp/p0 — R56 audit: B).
    Transverse (switch-gated, unchanged equation): x' += K_trans·x.
    """
    pv = rf_elem["parameters"]
    V = pv["voltage_kV"] * 1e3
    phi = pv["phase_rad"]
    k = 2.0 * np.pi * pv["frequency_GHz"] * 1e9 / C_SI
    tau = parray.tau()
    z_phys = -beta * tau                    # physical z = β·c·t = β·tau  [m]
    d_delta = (V / (beta**2 * E_total_eV)) * np.sin(phi + k * z_phys)
    add_p_oc(parray, beta * d_delta)               # p_oc += β0·δ_p (adapter boundary)
    if rf_transverse:
        K_trans = -E_SI * k * V / (2.0 * gamma * beta * M_E_SI * C_SI**2)
        x = parray.x(); y = parray.y()
        add_px(parray, K_trans * x)
        add_py(parray, K_trans * y)


# ═══════════════════════════════════════════════════════════
#  diagnostics
# ═══════════════════════════════════════════════════════════

def emit(x, xp):
    return np.sqrt(np.mean(x**2) * np.mean(xp**2) - np.mean(x * xp)**2)


# ═══════════════════════════════════════════════════════════
#  tracking (no manual tau/z override — trust OCELOT)
# ═══════════════════════════════════════════════════════════

def run_beamline(lat, rf_elems, sc_enabled=False, nparticles=None):
    """Track the beam through the given lattice, kicking every rf instance."""
    nparticles = nparticles or N_part
    # Random policy (v0.13): configured seed → x/y/tau; seed+1 → px/py/δ_p.
    rng_seed = int(cfg["random"]["seed"])
    np.random.seed(rng_seed)
    p = generate_parray(
        sigma_x=spot_rms, sigma_y=spot_rms,
        sigma_tau=sig_z0 / beta,   # OCELOT tau = c·t [m]; σ_tau = σ_z/β
        energy=(E_keV + 511.0) * 1e-6,   # TOTAL energy in GeV (E_kin+mc²)
        charge=Q_C,
        nparticles=nparticles,
    )
    np.random.seed(rng_seed + 1)
    N = p.rparticles.shape[1]
    set_px(p, np.random.normal(0.0, sigma_xp, N))
    set_py(p, np.random.normal(0.0, sigma_yp, N))
    # OCELOT p() = ΔE/(c·p0), NOT Δp/p0 → p_oc = β0·δ_p  (R56 audit: B)
    set_p_oc(p, beta * np.random.normal(0.0, sigma_delta, N))

    sc_proc = (SpaceCharge(step=sc.get("step", 1),
                           nmesh_xyz=list(sc.get("mesh", [63, 63, 63])))
               if (sc_enabled and _HAS_SC) else None)

    navi = Navigator(lat, unit_step=dz_track)
    if sc_proc is not None and len(lat.sequence) >= 2:
        navi.add_physics_proc(sc_proc, lat.sequence[0], lat.sequence[-1])

    total_length = sum(e.l for e in lat.sequence)

    probes_rem = {f"z{z:.0f}mm": z * 1e-3 for z in z_probes_mm}
    results = {}
    hist = {"z_mm": [], "sigma_x_um": [], "sigma_y_um": [], "sigma_z_um": [],
            "sigma_delta_e3": [], "eps_nx_mm_mrad": [], "eps_ny_mm_mrad": []}
    rf_done = set()
    rf_transverse = bool(switches.get("rf_transverse_kick", False))

    def _kick_rf_if_due(z_start):
        # one analytic kick per rf instance, at its own z_start, once
        for rf_elem in rf_elems:
            z_rf = rf_elem["z_start"]
            if z_rf not in rf_done and z_start >= z_rf - 1e-12:
                apply_rf_kick(p, rf_elem, rf_transverse)
                rf_done.add(z_rf)

    def _sample(z, z_before):
        # full z-history for the unified output (shared schema)
        xa = p.x(); xpa = p.px(); ya = p.y(); ypa = p.py()
        hist["z_mm"].append(z * 1e3)
        hist["sigma_x_um"].append(np.std(xa) * 1e6)
        hist["sigma_y_um"].append(np.std(ya) * 1e6)
        hist["sigma_z_um"].append(np.std(p.tau()) * beta * 1e6)
        hist["sigma_delta_e3"].append(np.std(p.p()) / beta * 1e3)   # δ_p = p_oc/β0
        hist["eps_nx_mm_mrad"].append(emit(xa, xpa) * beta_gamma * 1e6)
        hist["eps_ny_mm_mrad"].append(emit(ya, ypa) * beta_gamma * 1e6)
        for name, z_p in list(probes_rem.items()):
            if z_before <= z_p < z:
                xa = p.x(); xpa = p.px(); ya = p.y(); ypa = p.py()
                results[name] = {
                    "z_mm": z_p * 1e3,
                    "sigma_x_um": np.std(xa) * 1e6,
                    "sigma_y_um": np.std(ya) * 1e6,
                    "sigma_t_fs": np.std(p.tau()) / C_SI * 1e15,
                    "eps_x_mm_mrad": emit(xa, xpa) * 1e6,
                    "eps_y_mm_mrad": emit(ya, ypa) * 1e6,
                    "sigma_delta_e3": np.std(p.p()) / beta * 1e3,   # δ_p = p_oc/β0
                }
                del probes_rem[name]

    sc_apply_count = 0
    sc_events = []
    if sc_proc is not None:
        # ── SC ON: OCELOT NATIVE scheduler (v0.14.1 task 1) ──
        # get_next_step() is the exact mechanism of ocelot.track(); the
        # SpaceCharge PhysProc is triggered by the Navigator process counter
        # with coverage [s_start, s_stop) set by the add_physics_proc anchors
        # (cathode → sample).  The manual counter clone is retired.
        for t_maps, dz_step, proc_list, phys_steps in navi.get_next_step():
            z_start = navi.z0 - dz_step
            _kick_rf_if_due(z_start)
            for tm in t_maps:
                tm.apply(p)
            for proc, zstep in zip(proc_list, phys_steps):
                proc.z0 = navi.z0
                proc.apply(p, zstep)
                if proc is sc_proc:
                    sc_apply_count += 1
                    sc_events.append((float(navi.z0), float(zstep)))
            _sample(navi.z0, z_start)
    else:
        # ── SC OFF: unchanged loop (tracking_step; no PhysProc attached) ──
        n_steps = int(total_length / dz_track) if total_length > 0 else 1
        for step_i in range(n_steps):
            z_before = navi.z0
            _kick_rf_if_due(z_before)
            tracking_step(lat, p, dz_track, navi)
            _sample(navi.z0, z_before)

    results["history"] = hist
    results["rf_kicks_applied"] = len(rf_done)
    if sc_proc is not None:
        results["sc_scheduler"] = "ocelot_native"
        results["sc_apply_count"] = sc_apply_count
        results["sc_coverage_start_m"] = float(sc_proc.s_start)
        results["sc_coverage_stop_m"] = float(sc_proc.s_stop)
        results["sc_events"] = sc_events
    return results


# ═══════════════════════════════════════════════════════════
#  main
# ═══════════════════════════════════════════════════════════

def main():
    step = 4  # default: full beamline
    for i, a in enumerate(sys.argv):
        if a == "--step" and i + 1 < len(sys.argv):
            step = int(sys.argv[i + 1])
    if step not in STEP_ACTIVE:
        print(f"  unknown step {step}; must be 1..4")
        sys.exit(1)

    sc_on = step >= 4 and _HAS_SC
    # keep_zero_markers: SC ON keeps the cathode/sample zero-length markers
    # in the runtime sequence so SpaceCharge attaches as cathode → sample
    # (v0.14.1 task 1; SC OFF sequence unchanged).
    lat, rf_elems = build_lattice_from_shared(cfg, STEP_ACTIVE[step],
                                              keep_zero_markers=sc_on)
    total_length = sum(e.l for e in lat.sequence)
    n_sol = sum(1 for e in _lattice_elements(cfg)
                if e["type"] == "solenoid" and e["length"] > 0)
    n_rf = len([e for e in _lattice_elements(cfg)
                if e["type"] == "rf_cavity" and e["length"] > 0])

    print(f"\n{'='*65}")
    print(f"  UED Beamline Integration — Step {step}")
    print(f"{'='*65}")
    print(f"  E={E_keV:.0f} keV  ε_n={epsilon_n*1e6:.3f} mm·mrad  "
          f"Q={b['charge_fC']:.0f} fC")
    print(f"  lattice: {len(lat.sequence)} elements, total {total_length*1e3:.0f} mm "
          f"(solenoids={n_sol}, rf={n_rf})")
    print(f"  RF kicks: {len(rf_elems)} instance(s) "
          f"({'OFF' if not rf_elems else [round(e['z_start']*1e3) for e in rf_elems]})")
    print(f"  RF transverse kick: "
          f"{'ON' if switches.get('rf_transverse_kick', False) else 'OFF'} (shared switch)")
    print(f"  SC: {'ON' if (step >= 4 and _HAS_SC) else 'OFF'}  "
          f"(Navigator physics process, mesh={sc.get('mesh')})" if step >= 4 else "")

    r = run_beamline(lat, rf_elems, sc_enabled=(step >= 4))

    print(f"  {'Probe':>8s}  {'σ_x(μm)':>9s}  {'σ_y(μm)':>9s}  {'σ_t(fs)':>9s}  "
          f"{'ε_x(mm·mrad)':>13s}  {'σ_δ(e-3)':>9s}")
    for name in [f"z{z:.0f}mm" for z in z_probes_mm]:
        if name in r:
            d = r[name]
            print(f"  {name:>8s}  {d['sigma_x_um']:9.1f}  {d['sigma_y_um']:9.1f}  "
                  f"{d['sigma_t_fs']:9.0f}  {d['eps_x_mm_mrad']:13.4f}  "
                  f"{d['sigma_delta_e3']:9.2f}")

    # ═══════════════════════════════════════════════════════
    #  Step-specific validation
    # ═══════════════════════════════════════════════════════
    if step == 1:
        print(f"\n  --- Step 1 Validation: Drift vs Analytic ---")
        for name in [f"z{z:.0f}mm" for z in z_probes_mm]:
            if name in r:
                z_m = r[name]["z_mm"] * 1e-3
                sx_an = np.sqrt(spot_rms**2 + sigma_xp**2 * z_m**2) * 1e6
                sx_oc = r[name]["sigma_x_um"]
                err = abs(sx_oc - sx_an) / sx_an * 100
                flag = "✓" if err < 5 else "✗ FAIL"
                print(f"    z={z_m*1e3:5.0f}mm  σ_x OC={sx_oc:.1f}  "
                      f"analytic={sx_an:.1f}μm  err={err:.2f}%  {flag}")

    elif step == 2:
        print(f"\n  --- Step 2 Validation: Solenoid Focusing ---")
        probe_names = [f"z{z:.0f}mm" for z in z_probes_mm]
        sx_vals = [r[n]["sigma_x_um"] for n in probe_names if n in r]
        waist = min(sx_vals)
        print(f"    σ_x before solenoid (z=100mm): "
              f"{r.get('z100mm',{}).get('sigma_x_um','?'):.1f} μm")
        print(f"    σ_x after  solenoid (z=160mm): "
              f"{r.get('z160mm',{}).get('sigma_x_um','?'):.1f} μm")
        print(f"    Min σ_x in drift after: {waist:.1f} μm")
        print(f"    ε_x conserved: {r.get('z0mm',{}).get('eps_x_mm_mrad',0):.4f} → "
              f"{r.get('z777mm',{}).get('eps_x_mm_mrad',0):.4f} mm·mrad")

    elif step >= 3:
        print(f"\n  --- Step 3 Validation: RF Chirp (analytic kick) ---")
        rf0 = rf_elems[0]["parameters"]
        V0 = rf0["voltage_kV"] * 1e3
        phi0 = rf0["phase_rad"]
        k0 = 2.0 * np.pi * rf0["frequency_GHz"] * 1e9 / C_SI
        test_tau = np.linspace(-2 * sig_z0, 2 * sig_z0, 1000)   # tau = c·t [m]
        z_phys = -beta * test_tau
        d_delta = (V0 / (beta**2 * E_total_eV)) * np.sin(phi0 + k0 * z_phys)
        chirp = np.polyfit(z_phys, d_delta, 1)[0]
        rho = np.corrcoef(z_phys, d_delta)[0, 1]
        expected_chirp = -(V0 * k0 * np.cos(phi0)) / (beta**2 * E_total_eV)
        print(f"    δ = (V/(β²·E_total))·sin(φ + kz)")
        print(f"    V={V0*1e-3:.0f} kV  φ={phi0:.3f} rad  "
              f"E_total={E_total_eV*1e-3:.0f} keV")
        print(f"    dδ/dz = {chirp:.2f} m⁻¹  (expected: {expected_chirp:.2f} m⁻¹)")
        print(f"    corr(z,δ) = {rho:.4f}  σ_δ = {np.std(d_delta)*1e3:.2f}e-3")
        print(f"    Linear chirp: {'✓' if abs(rho) > 0.99 else '✗ FAIL'}")

        if step >= 4:
            print(f"\n  --- Step 4 Validation: Space Charge ---")
            r_sc = run_beamline(lat, rf_elems, sc_enabled=False)
            probe_names = [f"z{z:.0f}mm" for z in z_probes_mm]
            dsx = [(r[n]["sigma_x_um"] / r_sc[n]["sigma_x_um"] - 1) * 100
                   for n in probe_names if n in r and n in r_sc]
            dex = [(r[n]["eps_x_mm_mrad"] / r_sc[n]["eps_x_mm_mrad"] - 1) * 100
                   for n in probe_names if n in r and n in r_sc]
            print(f"    Max Δσ_x (SC ON vs OFF): {max(dsx):.2f}%")
            print(f"    Max Δε_x (SC ON vs OFF): {max(dex):.2f}%")
            print(f"    SC effect: {'DETECTED' if max(dsx) > 0.01 else 'NEGLIGIBLE (expected at 100 fC)'}")

    print("\n  Validation complete.")

    # ═══════════════════════════════════════════════════════
    #  unified output (shared schema) — SC flag follows the config
    # ═══════════════════════════════════════════════════════
    sc_flag = sc["enabled"]
    r_uni = r if (sc_flag == (step >= 4)) else run_beamline(lat, rf_elems,
                                                            sc_enabled=sc_flag)

    probes = []
    for z_mm in z_probes_mm:
        name = f"z{z_mm:.0f}mm"
        if name in r_uni:
            d = r_uni[name]
            probes.append(make_probe(
                z_mm=d["z_mm"],
                sigma_x_um=d["sigma_x_um"],
                sigma_y_um=d["sigma_y_um"],
                sigma_z_um=d["sigma_t_fs"] * C_SI * beta * 1e-9,
                sigma_delta_e3=d["sigma_delta_e3"],
                eps_nx_mm_mrad=d["eps_x_mm_mrad"] * beta_gamma,
                eps_ny_mm_mrad=d["eps_y_mm_mrad"] * beta_gamma,
            ))

    meta = {"model": "OCELOT macroparticle tracking (SecondTM)",
            "rf": "analytic thin kick per rf instance",
            "rf_kicks_applied": len(rf_elems),
            "step_flag": step,
            "switches": switches,
            "longitudinal_native_coordinate": "p_oc = dE/(c*p0)",
            "reported_delta": "delta_p = dp/p0",
            "conversion_beta": float(beta)}
    if r_uni.get("sc_scheduler") == "ocelot_native":
        meta["sc_scheduler"] = r_uni["sc_scheduler"]
        meta["sc_apply_count"] = r_uni["sc_apply_count"]
        meta["sc_coverage_start_m"] = r_uni["sc_coverage_start_m"]
        meta["sc_coverage_stop_m"] = r_uni["sc_coverage_stop_m"]
    path = write_results("GPT", probes, r_uni["history"], config_sha(cfg),
                         sc_flag, meta=meta)
    print(f"\n  Unified output -> {path}")


if __name__ == "__main__":
    main()
