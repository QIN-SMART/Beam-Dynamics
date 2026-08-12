#!/usr/bin/env python3
"""
SC runtime state contract validation (v0.14.1 task 3).

Verifies the unified runtime state machine and HARD FAIL rules:

  A. requested=False → sc_effective=False, no-SC canonical hash unchanged
  B. requested=True, normal run → available/configured/attached=True,
     apply_count>0, effective=True, coverage = cathode → sample
  C. requested=True + SpaceCharge import unavailable → HARD FAIL (raise)
  D. requested=True + final apply_count==0 → HARD FAIL (raise)
  E. GPT route state consistency: step defines capability only; config
     defines requested; step 1-3 never auto-runs SC
  F. saved JSON top-level sc_enabled == runtime sc_effective

Usage: /opt/anaconda3/bin/python3 validation/test_sc_runtime_state.py
"""

import os
import sys
import copy
import tempfile
import hashlib

import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
for p in (_THIS_DIR, os.path.dirname(_THIS_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from shared.params import load_config, _lattice_elements  # noqa: E402
from shared.output_schema import write_results, load_results  # noqa: E402
from validation.backend import run_ocelot  # noqa: E402

CANONICAL_HASH = "7790fd9c2a2b"
MESH_AC = [33, 33, 33]          # scheduler/state tests: mesh-independent
N_AC = 3000


def sc_off_hash(r):
    arr = np.array([r.z_mm, r.sigma_x_um, r.sigma_y_um, r.sigma_z_um,
                    r.eps_nx_mm_mrad, r.eps_ny_mm_mrad, r.sigma_delta_e3])
    return hashlib.sha1(np.ascontiguousarray(arr).tobytes()).hexdigest()[:12]


def cfg_with(cfg, **kw):
    c = copy.deepcopy(cfg)
    for k, v in kw.items():
        c["beam"][k] = v
    return c


def test_a_requested_false(cfg):
    print("== A. requested=False → effective=False + canonical hash ==")
    r = run_ocelot(cfg, "full")                    # SC OFF (default)
    ok = (r.meta["sc_requested"] is False
          and r.meta["sc_effective"] is False)
    h = sc_off_hash(r)
    ok &= h == CANONICAL_HASH
    print(f"  sc_requested={r.meta['sc_requested']} "
          f"sc_effective={r.meta['sc_effective']}  hash={h} "
          f"(canonical {CANONICAL_HASH})  "
          f"{'PASS' if ok else 'FAIL'}")
    return ok


def test_b_requested_normal(cfg):
    print("\n== B. requested=True, normal run ==")
    r = run_ocelot(cfg, "full", sc_enabled=True, n_particles=N_AC,
                   sc_mesh=MESH_AC)
    m = r.meta
    sample_elem = [e for e in _lattice_elements(cfg) if e["type"] == "sample"]
    z_sample = sample_elem[0]["z_start"]
    ok = (m["sc_requested"] is True
          and m["sc_available"] is True
          and m["sc_configured"] is True
          and m["sc_attached"] is True
          and m["sc_apply_count"] > 0
          and m["sc_effective"] is True
          and abs(m["sc_coverage_start_m"]) < 1e-9
          and abs(m["sc_coverage_stop_m"] - z_sample) < 1e-9
          and m["sc_scheduler"] == "ocelot_native")
    print(f"  requested={m['sc_requested']} available={m['sc_available']} "
          f"configured={m['sc_configured']} attached={m['sc_attached']} "
          f"apply_count={m['sc_apply_count']} effective={m['sc_effective']}")
    print(f"  coverage=[{m['sc_coverage_start_m']:.3f}, "
          f"{m['sc_coverage_stop_m']:.3f}] sample={z_sample:.3f}  "
          f"{'PASS' if ok else 'FAIL'}")
    return ok


def test_c_import_unavailable(cfg):
    print("\n== C. requested=True + import unavailable → HARD FAIL ==")
    import sys as _sys
    saved = _sys.modules.get("ocelot.cpbd.sc")
    _sys.modules["ocelot.cpbd.sc"] = None         # "halted import" trick
    failed = False
    try:
        run_ocelot(cfg, "full", sc_enabled=True, n_particles=N_AC,
                   sc_mesh=MESH_AC)
    except Exception as e:
        failed = True
        print(f"  raised {type(e).__name__}: {str(e)[:80]}")
    finally:
        if saved is not None:
            _sys.modules["ocelot.cpbd.sc"] = saved
        else:
            _sys.modules.pop("ocelot.cpbd.sc", None)
    print(f"  HARD FAIL on import: {'PASS' if failed else 'FAIL (no raise!)'}")
    return failed


def test_d_apply_count_zero(cfg):
    print("\n== D. requested=True + apply_count==0 → HARD FAIL ==")
    # zero-total-length lattice (only the zero-length cathode/sample markers
    # kept) → no tracking steps → no SC apply possible; parse() stays valid.
    c = copy.deepcopy(cfg)
    c["lattice"]["elements"] = [e for e in c["lattice"]["elements"]
                                if e["type"] in ("cathode", "sample")]
    failed = False
    msg = ""
    try:
        run_ocelot(c, "full", sc_enabled=True, n_particles=N_AC,
                   sc_mesh=MESH_AC)
    except RuntimeError as e:
        failed = True
        msg = str(e)
        print(f"  raised RuntimeError: {msg[:100]}")
    except Exception as e:
        failed = True
        msg = str(e)
        print(f"  raised {type(e).__name__}: {msg[:100]}")
    # the HARD FAIL must come from the state contract (apply_count==0)
    ok = failed and "apply_count" in msg
    print(f"  HARD FAIL on apply_count==0: "
          f"{'PASS' if ok else 'FAIL (no raise or wrong cause!)'}")
    return ok


def test_e_gpt_route_consistency(cfg):
    print("\n== E. GPT route state consistency (step = capability) ==")
    from GPT模拟 import ued_beamline_v2 as gpt
    cfg_on = copy.deepcopy(cfg)
    cfg_on["space_charge"]["enabled"] = True
    cfg_on["space_charge"]["mesh"] = list(MESH_AC)
    ok = True
    # decision function: config AND step>=4
    cases = [
        (cfg, 4, False),        # step4 + config OFF → no SC
        (cfg_on, 4, True),      # step4 + config ON  → SC
        (cfg_on, 3, False),     # step1-3 must NOT auto-run SC
        (cfg_on, 2, False),
        (cfg_on, 1, False),
        (cfg, 3, False),        # step3 + config OFF
    ]
    for c, step, want in cases:
        got = gpt.sc_requested_from(c, step)
        ok &= (got is want)
        print(f"  step={step} config_enabled="
              f"{c['space_charge']['enabled']} → requested={got} "
              f"(expect {want})  {'PASS' if got is want else 'FAIL'}")

    # runtime: step4+config OFF → effective False; step4+config ON → True
    saved_cfg = gpt.cfg
    try:
        gpt.cfg = cfg                       # patch module-level config (SC OFF)
        lat, rf = gpt.build_lattice_from_shared(cfg, gpt.STEP_ACTIVE[4],
                                                keep_zero_markers=False)
        r_off = gpt.run_beamline(lat, rf, sc_enabled=False, nparticles=N_AC)
        off_ok = (r_off.get("sc_effective") is False)
        print(f"  runtime step4+config OFF: effective="
              f"{r_off.get('sc_effective')}  {'PASS' if off_ok else 'FAIL'}")
        ok &= off_ok

        gpt.cfg = cfg_on                    # config SC ON (mesh 33³)
        lat, rf = gpt.build_lattice_from_shared(cfg_on, gpt.STEP_ACTIVE[4],
                                                keep_zero_markers=True)
        r_on = gpt.run_beamline(lat, rf, sc_enabled=True, nparticles=N_AC)
        on_ok = (r_on.get("sc_effective") is True
                 and r_on.get("sc_apply_count", 0) > 0)
        print(f"  runtime step4+config ON : effective="
              f"{r_on.get('sc_effective')} apply_count="
              f"{r_on.get('sc_apply_count')}  {'PASS' if on_ok else 'FAIL'}")
        ok &= on_ok
    finally:
        gpt.cfg = saved_cfg
    return ok


def test_f_saved_output_consistency(cfg):
    print("\n== F. saved JSON sc_enabled == runtime sc_effective ==")
    # main() now writes the single run's effective state (requested AND
    # step>=4); verify the write_results mapping end-to-end semantics.
    ok = True
    with tempfile.TemporaryDirectory() as td:
        for eff in (False, True):
            meta = {"sc_requested": eff, "sc_effective": eff}
            path = write_results("GPT", [], {"z_mm": []}, "test-sha",
                                 sc_enabled=eff, out_dir=td, meta=meta)
            d = load_results(path)
            same = d["sc_enabled"] is eff and d["meta"]["sc_effective"] is eff
            ok &= same
            print(f"  saved sc_enabled={d['sc_enabled']} "
                  f"meta.sc_effective={d['meta']['sc_effective']}  "
                  f"{'PASS' if same else 'FAIL'}")
    # main() wiring: sc_flag = r['sc_effective'] (single run, no second run)
    from GPT模拟 import ued_beamline_v2 as gpt
    src = open(gpt.__file__).read()
    wiring_ok = ("sc_effective = bool(r.get(" in src
                 and "r_uni = r" in src
                 and "run_beamline(lat, rf_elems, sc_enabled=sc_flag)" not in src)
    ok &= wiring_ok
    print(f"  main() single-run wiring (sc_flag = effective): "
          f"{'PASS' if wiring_ok else 'FAIL'}")
    return ok


def main():
    cfg = load_config()
    print("=" * 66)
    print("  SC runtime state contract validation (v0.14.1 task 3)")
    print("=" * 66)
    ok = True
    ok &= test_a_requested_false(cfg)
    ok &= test_b_requested_normal(cfg)
    ok &= test_c_import_unavailable(cfg)
    ok &= test_d_apply_count_zero(cfg)
    ok &= test_e_gpt_route_consistency(cfg)
    ok &= test_f_saved_output_consistency(cfg)
    print("\n  OVERALL: " + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
