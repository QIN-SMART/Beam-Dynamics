#!/usr/bin/env python3
"""
SC scheduler equivalence characterization (v0.14.1 task 1) — READ-ONLY w.r.t. production.

Compares TWO OCELOT PhysProc (SpaceCharge) scheduling implementations on the
IDENTICAL initial ParticleArray:

  A. NATIVE : Navigator.get_next_step() generator loop — the exact mechanism
              used by ocelot.cpbd.track.track() (track.py:470-478):
                  for t_maps, dz, proc_list, phys_steps in navi.get_next_step():
                      for tm in t_maps: tm.apply(p_array)
                      for proc, z_step in zip(proc_list, phys_steps):
                          proc.z0 = navi.z0
                          proc.apply(p_array, z_step)
  B. MANUAL : the project's current production replica (validation/backend.py
              run_ocelot loop, v0.14 P0 fix): fixed tracking_step(unit_step)
              + manual counter ("counter -= 1; if counter <= 0: apply;
              counter = step").

Every SpaceCharge apply event is logged with:
  event index, navi.z0 (apply position), proc.z0, zstep passed to apply,
  counter before, counter after.

Comparison levels:
  L1 apply count | L2 per-event z | L3 per-event zstep
  L4 final rparticles (bitwise + max abs diff) | L5 sigma stats

Scenarios (independent custom lattices, production config beam):
  T1 single drift, step=1                       (baseline, aligned)
  T2 single drift, step=5                       (multi-step counter)
  T3 multi-element lattice (drift+solenoid+drift)
  T4 element boundary NOT a multiple of unit_step
  T5 total length NOT a multiple of unit_step   (tail-segment case)
  T6 SC process covering only PART of the lattice

This file does NOT modify any production module.  It only reports.
Usage: /opt/anaconda3/bin/python3 validation/test_sc_scheduler_equivalence.py
"""

import os
import sys
import copy

import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
for p in (_THIS_DIR, os.path.dirname(_THIS_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from shared.params import load_config, derived  # noqa: E402
from shared.constants import C_SI, E_SI, M_E_SI  # noqa: E402
from shared.ocelot_coords import set_px, set_py, set_p_oc  # noqa: E402

# ── test parameters ─────────────────────────────────────────────────────────
Q_TEST_FEMTO = 500.0       # fC — SC diagnostics smoke charge
N_TEST = 10000             # macroparticles (scheduler test, not physics conv.)
# mesh 33³ for the SCHEDULING tests: event semantics (count/z/zstep) do not
# depend on the mesh; 63³ (production) would make pyfftw re-planning ~0.7 s per
# apply and the suite ~50 min.  Production-mesh numerical validation is a
# v0.15 concern (SC convergence hardening).
MESH = [33, 33, 33]
UNIT = 0.005               # Navigator.unit_step [m]
B_SOL = 0.05               # T  (solenoid strength, same as production config)

# ═════════════════════════════════════════════════════════════════════════
# beam generation — mirrors validation/backend.py run_ocelot (seed/seed+1)
# ═════════════════════════════════════════════════════════════════════════

def make_beam(cfg, nparticles):
    """Identical to the production beam generation (v0.13 random policy):
    seed -> generate_parray (x/y/tau); seed+1 -> px/py/delta_p (p_oc=beta*dp)."""
    from ocelot.cpbd.beam import generate_parray
    d = derived(cfg)
    ib = cfg["initial_distribution"]
    rng_seed = int(cfg["random"]["seed"])
    np.random.seed(rng_seed)
    p = generate_parray(
        sigma_x=ib["spot_rms_um"] * 1e-6, sigma_y=ib["spot_rms_um"] * 1e-6,
        sigma_tau=ib["bunch_length_um"] * 1e-6 / d["beta"],
        energy=(cfg["beam"]["energy_keV"] + 511.0) * 1e-6,   # TOTAL GeV
        charge=Q_TEST_FEMTO * 1e-15,
        nparticles=nparticles,
    )
    np.random.seed(rng_seed + 1)
    N = p.rparticles.shape[1]
    set_px(p, np.random.normal(0.0, d["sigma_xp"], N))
    set_py(p, np.random.normal(0.0, d["sigma_yp"], N))
    set_p_oc(p, d["beta"] * np.random.normal(0.0, ib["sigma_delta"], N))
    return p


def build_lattice(elems_spec):
    """Custom lattice from specs [(type, length), ...]; solenoid k from shared
    kinematics.  Returns (lat, total_len)."""
    from ocelot.cpbd.elements import Drift, Solenoid
    from ocelot.cpbd.magnetic_lattice import MagneticLattice
    cfg = load_config()
    p_SI = derived(cfg)["p_SI"]
    elems = []
    for etype, L in elems_spec:
        if etype == "drift":
            elems.append(Drift(l=L))
        elif etype == "solenoid":
            elems.append(Solenoid(l=L, k=E_SI * B_SOL / (2.0 * p_SI)))
        else:
            raise ValueError(f"unknown element type {etype}")
    lat = MagneticLattice(elems, method={"global": "SecondTM"})
    lat.update_transfer_maps()
    return lat, sum(e.l for e in elems)


# ═════════════════════════════════════════════════════════════════════════
# the two scheduling implementations (characterization only)
# ═════════════════════════════════════════════════════════════════════════

def run_native(lat, p, sc_step, unit, sc_span=None):
    """A. Native: Navigator.get_next_step() loop == ocelot track() core."""
    from ocelot.cpbd.navi import Navigator
    from ocelot.cpbd.sc import SpaceCharge
    navi = Navigator(lat, unit_step=unit)
    sc = SpaceCharge(step=sc_step, nmesh_xyz=list(MESH))
    if sc_span is None:
        navi.add_physics_proc(sc, lat.sequence[0], lat.sequence[-1])
    else:
        i0, i1 = sc_span
        navi.add_physics_proc(sc, lat.sequence[i0], lat.sequence[i1])
    events = []
    for t_maps, dz, proc_list, phys_steps in navi.get_next_step():
        for tm in t_maps:
            tm.apply(p)
        for proc, zstep in zip(proc_list, phys_steps):
            cb = proc.counter
            proc.z0 = navi.z0
            proc.apply(p, zstep)
            events.append({
                "idx": len(events), "z": float(navi.z0), "proc_z0": float(proc.z0),
                "zstep": float(zstep), "counter_before": cb,
                "counter_after": float(proc.counter),
                "proc": proc.__class__.__name__})
    return events, p


def run_manual(lat, p, sc_step, unit, sc_span=None):
    """B. Manual: production replica (validation/backend.py v0.14 loop):
    fixed tracking_step(unit) + manual counter, zstep = step*unit."""
    from ocelot.cpbd.navi import Navigator
    from ocelot.cpbd.sc import SpaceCharge
    from ocelot.cpbd.track import tracking_step
    navi = Navigator(lat, unit_step=unit)
    sc = SpaceCharge(step=sc_step, nmesh_xyz=list(MESH))
    if sc_span is None:
        navi.add_physics_proc(sc, lat.sequence[0], lat.sequence[-1])
    else:
        i0, i1 = sc_span
        navi.add_physics_proc(sc, lat.sequence[i0], lat.sequence[i1])
    total = sum(e.l for e in lat.sequence)
    n_steps = int(total / unit) if total > 0 else 1
    events = []
    for _ in range(n_steps):
        tracking_step(lat, p, unit, navi)
        if sc is not None:
            cb = sc.counter
            sc.counter -= 1
            if sc.counter <= 0:
                sc.z0 = navi.z0
                sc.apply(p, sc.step * unit)          # nominal zstep (production)
                sc.counter = sc.step
                events.append({
                    "idx": len(events), "z": float(navi.z0), "proc_z0": float(sc.z0),
                    "zstep": float(sc.step * unit), "counter_before": cb,
                    "counter_after": float(sc.counter),
                    "proc": sc.__class__.__name__})
    return events, p


# ═════════════════════════════════════════════════════════════════════════
# comparison
# ═════════════════════════════════════════════════════════════════════════

def emit(x, xp):
    return np.sqrt(np.mean(x**2) * np.mean(xp**2) - np.mean(x * xp)**2)


def stats(p, beta, beta_gamma):
    return {
        "sigma_x_um": float(np.std(p.x()) * 1e6),
        "sigma_y_um": float(np.std(p.y()) * 1e6),
        "sigma_z_um": float(np.std(p.tau()) * beta * 1e6),
        "sigma_delta_e3": float(np.std(p.p()) / beta * 1e3),
        "eps_nx_mm_mrad": float(emit(p.x(), p.px()) * beta_gamma * 1e6),
        "eps_ny_mm_mrad": float(emit(p.y(), p.py()) * beta_gamma * 1e6),
    }


def compare(name, ev_n, p_n, ev_m, p_m, beta, beta_gamma):
    print(f"\n--- {name} ---")
    d = derived(load_config())
    L = {}
    # L1: apply count
    L["L1_count"] = len(ev_n) == len(ev_m)
    n = min(len(ev_n), len(ev_m))
    # L2/L3: per-event z and zstep
    L2_diff = L3_diff = None
    for i in range(n):
        if abs(ev_n[i]["z"] - ev_m[i]["z"]) > 1e-12 and L2_diff is None:
            L2_diff = (i, ev_n[i]["z"], ev_m[i]["z"])
        if abs(ev_n[i]["zstep"] - ev_m[i]["zstep"]) > 1e-12 and L3_diff is None:
            L3_diff = (i, ev_n[i]["zstep"], ev_m[i]["zstep"])
    L["L2_z"] = L2_diff is None
    L["L3_zstep"] = L3_diff is None
    # L4: final rparticles
    bit = np.array_equal(p_n.rparticles, p_m.rparticles)
    md = float(np.max(np.abs(p_n.rparticles - p_m.rparticles)))
    L["L4_bitwise"] = bit
    # L5: sigma stats
    sn, sm = stats(p_n, beta, beta_gamma), stats(p_m, beta, beta_gamma)
    dev = {k: abs(sn[k] - sm[k]) / abs(sm[k]) if abs(sm[k]) > 1e-12 else 0.0
           for k in sn}
    L["L5_max_dev"] = max(dev.values())

    print(f"  events: native={len(ev_n)}  manual={len(ev_m)}   "
          f"[L1 {'PASS' if L['L1_count'] else 'FAIL'}]")
    if L2_diff:
        print(f"  first z diff  @event {L2_diff[0]}: native={L2_diff[1]:.6f} "
              f"manual={L2_diff[2]:.6f}  [L2 FAIL]")
    else:
        print(f"  z per event identical up to {n} compared events  [L2 {'PASS' if L['L2_z'] else 'FAIL'}]")
    if L3_diff:
        print(f"  first zstep diff @event {L3_diff[0]}: native={L3_diff[1]:.6f} "
              f"manual={L3_diff[2]:.6f}  [L3 FAIL]")
    else:
        print(f"  zstep per event identical up to {n} compared events  [L3 {'PASS' if L['L3_zstep'] else 'FAIL'}]")
    print(f"  rparticles: bitwise={bit}  max_abs_diff={md:.3e}  [L4]")
    print(f"  sigma stats max rel dev = {max(dev.values()):.3e}  [L5]")
    if L["L1_count"] and L["L2_z"] and L["L3_zstep"]:
        verdict = "STRICTLY EQUIVALENT" if bit else \
            "EVENTS EQUIVALENT (rparticles differs only by FP rounding)"
    else:
        verdict = "SCHEDULING DIFFERENCE FOUND"
    print(f"  verdict: {verdict}")
    return L, (L2_diff, L3_diff), dev


# ═════════════════════════════════════════════════════════════════════════
# scenarios
# ═════════════════════════════════════════════════════════════════════════

SCENARIOS = [
    # full-coverage scenarios end with a ZERO-LENGTH element so that
    # s_stop = totalLen (otherwise s_stop = start of the LAST element,
    # i.e. the SC interval misses the final element — see T7).
    dict(name="T1_baseline_fullcov_step1",
         elems=[("drift", 0.500), ("drift", 0.000)], step=1, span=None,
         expect_equivalent=True),
    dict(name="T2_fullcov_step5",
         elems=[("drift", 0.500), ("drift", 0.000)], step=5, span=None,
         expect_equivalent=True),
    dict(name="T3_multi_elem_drift_solenoid",
         elems=[("drift", 0.100), ("solenoid", 0.060), ("drift", 0.100),
                ("drift", 0.000)],
         step=1, span=None, expect_equivalent=True),
    dict(name="T4_elem_boundary_not_unit_multiple",
         elems=[("drift", 0.050), ("solenoid", 0.022), ("drift", 0.033),
                ("drift", 0.000)],
         step=1, span=None, expect_equivalent=True),
    dict(name="T5_total_len_not_unit_multiple",
         elems=[("drift", 0.453), ("drift", 0.050)], step=1, span=None,
         expect_equivalent=False),
    dict(name="T6_partial_coverage_sc_mid_lattice",
         elems=[("drift", 0.100), ("drift", 0.200), ("drift", 0.100)],
         step=1, span=(1, 2), expect_equivalent=False),
]


def production_lattice(keep_markers=False):
    """Production full-beamline lattice exactly as validation/backend.py
    run_ocelot builds it (lattice.elements -> Drift/Solenoid; rf_cavity body
    is a Drift + analytic kick outside the scheduler).  keep_markers=True
    reproduces the SC-ON sequence (cathode/sample zero-length markers kept)."""
    from ocelot.cpbd.elements import Drift, Solenoid
    from ocelot.cpbd.magnetic_lattice import MagneticLattice
    from shared.params import _lattice_elements
    cfg = load_config()
    p_SI = derived(cfg)["p_SI"]
    elems = []
    for e in _lattice_elements(cfg):
        etype, L = e["type"], e["length"]
        if L <= 0:
            if keep_markers and etype in ("cathode", "sample"):
                elems.append(Drift(l=0.0, eid=f"{e['name']}_MARKER"))
            continue                       # cathode / sample markers
        if etype == "solenoid":
            elems.append(Solenoid(l=L, k=E_SI * e["parameters"]["B_field_T"]
                                  / (2.0 * p_SI), eid=e["name"]))
        elif etype == "rf_cavity":
            elems.append(Drift(l=L, eid=f"{e['name']}_BODY"))
        else:
            elems.append(Drift(l=L, eid=e["name"]))
    lat = MagneticLattice(elems, method={"global": "SecondTM"})
    lat.update_transfer_maps()
    return lat, sum(e.l for e in elems)


def production_acceptance(cfg, beta, beta_gamma):
    """Production-path acceptance (v0.14.1 task 1 requirements A–F).

    A. configured SC coverage end == shared sample position (lattice.elements)
    B. runtime execution: sc_apply_count > 0
    C. step=1 / step=5 event count / z / zstep == native reference
    D. non-integer tail: no lost tail in the native scheduler
    E. partial coverage: no apply outside [s_start, s_stop)
    F. SC OFF bitwise: production run_ocelot output array hash (baseline)
    """
    import hashlib
    from validation.backend import run_ocelot
    from shared.params import _lattice_elements
    print("\n" + "=" * 66)
    print("  PRODUCTION-PATH ACCEPTANCE (run_ocelot, SC ON = native scheduler)")
    print("=" * 66)
    ok = True
    N_AC = 3000                   # event-semantics checks do not need 1e4
    # mesh 33³ (not production 63³): scheduler event semantics do not depend
    # on the mesh; production-mesh numerical results are v0.15 (see MESH note).
    MESH_AC = [33, 33, 33]

    # ── A+B+C(step=1): one production run covers coverage / count / events ──
    sample = [e for e in _lattice_elements(cfg) if e["type"] == "sample"][0]
    z_sample = sample["z_start"]
    r = run_ocelot(cfg, "full", sc_enabled=True, n_particles=N_AC,
                   sc_mesh=MESH_AC)
    m = r.meta
    a_ok = (m["sc_scheduler"] == "ocelot_native"
            and abs(m["sc_coverage_start_m"]) < 1e-12
            and abs(m["sc_coverage_stop_m"] - z_sample) < 1e-12)
    print(f"  A. scheduler={m['sc_scheduler']} coverage=[{m['sc_coverage_start_m']:.3f}, "
          f"{m['sc_coverage_stop_m']:.3f}] sample(z_start)={z_sample:.3f}  "
          f"{'PASS' if a_ok else 'FAIL'}")
    ok &= a_ok

    b_ok = m["sc_apply_count"] > 0
    print(f"  B. sc_apply_count = {m['sc_apply_count']}  "
          f"{'PASS' if b_ok else 'FAIL'}")
    ok &= b_ok

    lat_p, _ = production_lattice(keep_markers=True)
    p_ref = make_beam(cfg, N_AC)
    ev_ref, _ = run_native(lat_p, p_ref, 1, 0.001)
    ev = m["sc_events"]
    c1_ok = len(ev) == len(ev_ref) and all(
        abs(ev[i][0] - ev_ref[i]["z"]) < 1e-12
        and abs(ev[i][1] - ev_ref[i]["zstep"]) < 1e-12
        for i in range(len(ev)))
    print(f"  C. step=1: production events={len(ev)} "
          f"native reference={len(ev_ref)}  "
          f"{'PASS' if c1_ok else 'FAIL'}")
    ok &= c1_ok

    # ── C(step=5) ──
    rp5 = run_ocelot(cfg, "full", sc_enabled=True, n_particles=N_AC, sc_step=5,
                     sc_mesh=MESH_AC)
    p_ref5 = make_beam(cfg, N_AC)
    ev_ref5, _ = run_native(lat_p, p_ref5, 5, 0.001)
    ev5 = rp5.meta["sc_events"]
    c5_ok = len(ev5) == len(ev_ref5) and all(
        abs(ev5[i][0] - ev_ref5[i]["z"]) < 1e-12
        and abs(ev5[i][1] - ev_ref5[i]["zstep"]) < 1e-12
        for i in range(len(ev5)))
    print(f"  C. step=5: production events={len(ev5)} "
          f"native reference={len(ev_ref5)}  "
          f"{'PASS' if c5_ok else 'FAIL'}")
    ok &= c5_ok

    # ── D: non-integer tail — native scheduler must not lose the tail ──
    lat5, total5 = build_lattice([("drift", 0.453), ("drift", 0.050)])
    p5 = make_beam(cfg, N_AC)
    ev5t, _ = run_native(lat5, p5, 1, UNIT)
    d_ok = ev5t[-1]["z"] == 0.453 and abs(ev5t[-1]["zstep"] - 0.003) < 1e-9
    print(f"  D. tail (0.503 m, unit 0.005): last event z={ev5t[-1]['z']:.3f} "
          f"zstep={ev5t[-1]['zstep']:.4f} (expect 0.453 / 0.003)  "
          f"{'PASS' if d_ok else 'FAIL'}")
    ok &= d_ok

    # ── E: partial coverage — no apply outside [s_start, s_stop) ──
    lat6, total6 = build_lattice([("drift", 0.100), ("drift", 0.200),
                                  ("drift", 0.100)])
    p6 = make_beam(cfg, N_AC)
    ev6, _ = run_native(lat6, p6, 1, UNIT, sc_span=(1, 2))
    e_ok = all(0.1 - 1e-9 <= ev["z"] <= 0.3 + 1e-9 for ev in ev6)
    print(f"  E. partial coverage [0.1, 0.3): events={len(ev6)} all inside "
          f"({'PASS' if e_ok else 'FAIL'})")
    ok &= e_ok

    # ── F: SC OFF bitwise — production output array hash ──
    # CANONICAL no-SC regression hash definition (v0.14.1):
    #   SHA1( contiguous bytes of
    #         [z_mm, sigma_x_um, sigma_y_um, sigma_z_um,
    #          eps_nx_mm_mrad, eps_ny_mm_mrad, sigma_delta_e3] )[:12]
    #   canonical value = 7790fd9c2a2b  (run_ocelot SC OFF, config seed 42,
    #   N=5e4, dz=0.001, full beamline).
    # The old e041d6ae9fb7a0d2 is NOT tracked: its input object/algorithm
    # were never recorded, so it is not strictly reproducible.
    r_off = run_ocelot(cfg, "full")               # SC OFF (default, N=5e4)
    arr = np.array([r_off.z_mm, r_off.sigma_x_um, r_off.sigma_y_um,
                    r_off.sigma_z_um, r_off.eps_nx_mm_mrad,
                    r_off.eps_ny_mm_mrad, r_off.sigma_delta_e3])
    h = hashlib.sha1(np.ascontiguousarray(arr).tobytes()).hexdigest()[:12]
    f_ok = h == "7790fd9c2a2b"
    print(f"  F. SC OFF array hash = {h} (pre-fix baseline 7790fd9c2a2b)  "
          f"{'PASS' if f_ok else 'FAIL'}")
    print(f"     sample sigma_x={np.interp(777, r_off.z_mm, r_off.sigma_x_um):.3f} "
          f"sigma_z={np.interp(777, r_off.z_mm, r_off.sigma_z_um):.3f} "
          f"(ref 1996.205/474.022)")
    ok &= f_ok

    print("  production-path acceptance: " + ("ALL PASS" if ok else "FAILURES — see above"))
    return ok


def check_expectation(name, L, expect_equivalent):
    """Check a characterization scenario against its design expectation.

    T1–T4 (expect_equivalent=True) must show event-level equivalence;
    T5–T7 (expect_equivalent=False) must show a scheduling DIFFERENCE
    (they exist to demonstrate the manual counter is not generally
    equivalent).  A T5–T7 scenario that unexpectedly becomes equivalent is
    reported as FAIL — it may mean the characterization was broken, not
    that the manual clone got better.
    """
    events_equiv = bool(L["L1_count"] and L["L2_z"] and L["L3_zstep"])
    if expect_equivalent:
        ok = events_equiv
        msg = ("PASS (equivalent as expected)" if ok else
               "FAIL (expected equivalence, found a difference)")
    else:
        ok = not events_equiv
        msg = ("PASS (difference as expected)" if ok else
               "FAIL (expected difference, found equivalence — "
               "scenario may be broken)")
    print(f"  expectation[{name}]: {msg}")
    return ok


def main():
    cfg = load_config()
    d = derived(cfg)
    beta, beta_gamma = d["beta"], d["beta_gamma"]
    print("=" * 66)
    print("  SC scheduler equivalence — native get_next_step vs manual counter")
    print("=" * 66)
    print(f"  beam: 100 keV, Q={Q_TEST_FEMTO:.0f} fC, N={N_TEST}, mesh={MESH}")
    print(f"  unit_step = {UNIT} m, solenoid B = {B_SOL} T")
    print("  (READ-ONLY characterization; no production module modified)")

    results = {}
    exp_ok_all = True
    for sc in SCENARIOS:
        lat, total = build_lattice(sc["elems"])
        # identical initial beam for both paths (deterministic same seed)
        p_n = make_beam(cfg, N_TEST)
        p_m = make_beam(cfg, N_TEST)
        ev_n, _ = run_native(lat, p_n, sc["step"], UNIT, sc["span"])
        ev_m, _ = run_manual(lat, p_m, sc["step"], UNIT, sc["span"])
        L, diffs, dev = compare(sc["name"], ev_n, p_n, ev_m, p_m,
                                beta, beta_gamma)
        exp_ok = check_expectation(sc["name"], L, sc["expect_equivalent"])
        exp_ok_all &= exp_ok
        results[sc["name"]] = (L, diffs, total, exp_ok)
        # show the first few events side by side (for the record)
        print("  first events (idx, z, zstep, cb, ca):")
        for k in range(min(4, len(ev_n), len(ev_m))):
            e_n, e_m = ev_n[k], ev_m[k]
            print(f"    native  #{k}: z={e_n['z']:.6f} zstep={e_n['zstep']:.6f} "
                  f"cb={e_n['counter_before']} ca={e_n['counter_after']}")
            print(f"    manual  #{k}: z={e_m['z']:.6f} zstep={e_m['zstep']:.6f} "
                  f"cb={e_m['counter_before']} ca={e_m['counter_after']}")

    # ── T7: production usage characterization ──
    # (old anchor seq[0]..seq[-1] → SC ends at the last non-zero element
    #  start; manual has no coverage bound at all → expected difference)
    print("\n--- T7_production_lattice_usage (seq[0]..seq[-1], dz=0.005) ---")
    lat, total = production_lattice()
    p_n = make_beam(cfg, N_TEST)
    p_m = make_beam(cfg, N_TEST)
    ev_n, _ = run_native(lat, p_n, 1, UNIT)
    ev_m, _ = run_manual(lat, p_m, 1, UNIT)
    L7, _, dev7 = compare("T7_production_lattice_usage", ev_n, p_n, ev_m, p_m,
                          beta, beta_gamma)
    exp7 = check_expectation("T7_production_lattice_usage", L7,
                             expect_equivalent=False)
    exp_ok_all &= exp7
    from shared.params import _lattice_elements
    last_L = [e["length"] for e in _lattice_elements(cfg) if e["length"] > 0][-1]
    print(f"  note: s_stop = totalLen - L(last element) = "
          f"{total - last_L:.3f} m (last element = drift3, {last_L} m);")
    print(f"        native applies SC on [0, {total - last_L:.3f}) only; "
          f"manual keeps applying through the last drift to the sample.")
    results["T7_production_lattice_usage"] = (L7, None, total, exp7)

    # ── production-path acceptance (requirements A–F) ──
    prod_ok = production_acceptance(cfg, beta, beta_gamma)

    # ── summary ──
    print("\n" + "=" * 66)
    print("  SUMMARY (L1 count / L2 z / L3 zstep / L4 bitwise / expectation)")
    print("=" * 66)
    for name, (L, diffs, total, exp_ok) in results.items():
        flags = "".join("P" if v else "F" for v in
                        (L["L1_count"], L["L2_z"], L["L3_zstep"], L["L4_bitwise"]))
        verdict = ("EQUIV" if flags.startswith("PPP")
                   else "DIFFER")
        print(f"  {name:<38s} L={flags}  total={total:.3f}m  {verdict}  "
              f"expect={'PASS' if exp_ok else 'FAIL'}")
    print("\n  characterization expectations: "
          + ("ALL PASS" if exp_ok_all else "FAILURES — see above"))
    print("  production-path acceptance: " + ("ALL PASS" if prod_ok else "FAIL"))
    print("  (T1–T4 equivalence; T5–T7 documented differences = manual clone")
    print("   retired; production now uses the native scheduler)")
    overall = exp_ok_all and prod_ok
    print("  OVERALL: " + ("PASS" if overall else "FAIL"))
    return 0 if overall else 1


if __name__ == "__main__":
    sys.exit(main())
