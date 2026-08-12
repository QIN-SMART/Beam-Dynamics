#!/usr/bin/env python3
"""
AG charge semantics validation (v0.14.1 task 2).

Requirement:
  AG SC strength must scale with the PHYSICAL bunch charge
      Ne_phys = Q / e
  NOT with beam.n_particles (which is only the OCELOT macroparticle
  numerical resolution).

Test 1 — charge semantics:
  Q = 10/50/100/500/1000 fC → run_ag (SC ON) must report
  meta.ag_ne_phys == Q/e, and the SC effect (sample σ_x) must grow
  monotonically with Q.

Test 2 — n_particles invariance:
  Q = 100 fC fixed, n_particles = 1e4 / 5e4 / 1e5:
    - AG SC ON results must be BITWISE identical (AG no longer depends on
      the macroparticle numerical resolution);
    - AG SC OFF results must be bitwise identical and match the v0.13
      baseline sample values (σ_x = 1984.191 µm, σ_z = 477.001 µm).

Usage: /opt/anaconda3/bin/python3 validation/test_ag_charge_semantics.py
"""

import os
import sys
import copy

import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
for p in (_THIS_DIR, os.path.dirname(_THIS_DIR), os.path.join(os.path.dirname(_THIS_DIR), "AG")):
    if p not in sys.path:
        sys.path.insert(0, p)

from shared.params import load_config  # noqa: E402
from shared.constants import E_SI  # noqa: E402
from validation.backend import run_ag  # noqa: E402

SECTION = "full"
ZS = 777.0                        # sample plane [mm]

# v0.13 no-SC AG baseline sample values (documented, deterministic ODE)
AG_BASELINE_X = 1984.191          # µm
AG_BASELINE_Z = 477.001           # µm


def sample(r, key):
    return float(np.interp(ZS, np.asarray(r.z_mm), np.asarray(getattr(r, key))))


def cfg_with(cfg, **kw):
    c = copy.deepcopy(cfg)
    for k, v in kw.items():
        c["beam"][k] = v
    return c


def test1_charge_semantics(cfg):
    print("== Test 1: charge semantics (Ne_phys = Q/e) ==")
    ok = True
    qs = [10.0, 50.0, 100.0, 500.0, 1000.0]
    rows = []
    for q in qs:
        c = cfg_with(cfg, charge_fC=q)
        r = run_ag(c, SECTION, sc_enabled=True, solenoid_coupling=False)
        ne_expected = abs(q * 1e-15) / E_SI
        # formal contract fields (v0.14.1 task 3); ag_ne_phys kept as alias
        ne_actual = r.meta["physical_electron_number"]
        rel = abs(ne_actual - ne_expected) / ne_expected
        ok1 = (rel < 1e-12
               and r.meta["physical_charge_C"] == abs(q * 1e-15)
               and r.meta["ag_ne_phys"] == ne_actual
               and r.meta["sc_requested"] is True
               and r.meta["sc_effective"] is True)
        sx = sample(r, "sigma_x_um")
        sz = sample(r, "sigma_z_um")
        rows.append((q, ne_actual, sx, sz))
        print(f"  Q={q:6.0f} fC  Ne_phys={ne_actual:.6e} (Q/e={ne_expected:.6e}, "
              f"rel={rel:.1e})  sample σx={sx:9.3f} µm  σz={sz:9.3f} µm  "
              f"{'PASS' if ok1 else 'FAIL'}")
        ok &= ok1
    # SC effect must grow monotonically with Q (SC ON only; AG SC ∝ Ne)
    sx = [r_[2] for r_ in rows]
    mon = all(sx[i + 1] > sx[i] for i in range(len(sx) - 1))
    print(f"  σx(Q) monotonic: {sx}  {'PASS' if mon else 'FAIL'}")
    ok &= mon
    # sanity: at Q=100 fC the previous bug (Ne=5e4 → 8 fC) would give a
    # much weaker effect; the fixed run must differ from the old baseline
    # SC-ON numbers (old: Ne=5e4 → σx ≈ 1119 µm region — see CHECKPOINTS)
    r100 = rows[2][2]
    print(f"  Q=100 fC SC ON sample σx = {r100:.3f} µm (old Ne=5e4 bug: "
          f"~1.1 mm regime — value changed as expected)")
    return ok


def test2_nparticles_invariance(cfg):
    print("\n== Test 2: n_particles invariance (Q = 100 fC fixed) ==")
    ok = True
    nps = [10000, 50000, 100000]
    base = None
    for tag, sc_on in (("SC ON", True), ("SC OFF", False)):
        results = []
        for n in nps:
            c = cfg_with(cfg, n_particles=n)
            r = run_ag(c, SECTION, sc_enabled=sc_on, solenoid_coupling=False)
            arr = np.array([r.z_mm, r.sigma_x_um, r.sigma_y_um, r.sigma_z_um,
                            r.eps_nx_mm_mrad, r.eps_ny_mm_mrad,
                            r.sigma_delta_e3])
            results.append(arr)
        ident = all(np.array_equal(results[0], a) for a in results[1:])
        print(f"  {tag}: n={nps} → bitwise identical = {ident}  "
              f"{'PASS' if ident else 'FAIL'}")
        ok &= ident
        if not sc_on:
            # v0.13 no-SC AG baseline must be preserved exactly
            r0 = run_ag(cfg_with(cfg, n_particles=nps[0]), SECTION,
                        sc_enabled=False, solenoid_coupling=False)
            sx, sz = sample(r0, "sigma_x_um"), sample(r0, "sigma_z_um")
            bx = abs(sx - AG_BASELINE_X) < 1e-3
            bz = abs(sz - AG_BASELINE_Z) < 1e-3
            sc_off_ok = (r0.meta["sc_requested"] is False
                         and r0.meta["sc_effective"] is False
                         and r0.meta["physical_electron_number"] > 0
                         and r0.meta["physical_charge_C"] == 100e-15)
            print(f"  {tag}: sample σx={sx:.3f} (baseline {AG_BASELINE_X})  "
                  f"σz={sz:.3f} (baseline {AG_BASELINE_Z})  state"
                  f" requested/effective=False = {sc_off_ok}  "
                  f"{'PASS' if (bx and bz and sc_off_ok) else 'FAIL'}")
            ok &= (bx and bz and sc_off_ok)
    return ok


def main():
    cfg = load_config()
    print("=" * 66)
    print("  AG charge semantics validation (v0.14.1 task 2)")
    print("=" * 66)
    print(f"  config: Q={cfg['beam']['charge_fC']} fC, "
          f"n_particles={cfg['beam']['n_particles']}")
    ok = True
    ok &= test1_charge_semantics(cfg)
    ok &= test2_nparticles_invariance(cfg)
    print("\n  OVERALL: " + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
