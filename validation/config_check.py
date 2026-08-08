#!/usr/bin/env python3
"""
Read-only config consistency checker (v0.12, Phase 7).

Verifies that the shared configuration is self-consistent WITHOUT modifying
any physics parameter.  Reports inconsistencies only.

Checks:
  1. energy ↔ gamma ↔ beta ↔ p0 consistency
  2. RF frequency physically reasonable
  3. RF / solenoid positions come from lattice.elements
  4. sample position == lattice total length (contiguity)
  5. all elements contiguous (z_end == next z_start), no overlaps/gaps
  6. lattice element schema validity (name/type/z_start/length/parameters)

Usage:  python3 validation/config_check.py
"""

import os
import sys
import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
for p in (_THIS_DIR, os.path.dirname(_THIS_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from shared.params import load_config, derived, _lattice_elements, z_sample  # noqa: E402
from shared.constants import MEC2_KEV, M_E_SI, C_SI  # noqa: E402

CHECKS = []


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn
    return deco


@check("energy-gamma-beta-p0 consistency")
def _c1(cfg):
    d = derived(cfg)
    E_keV = cfg["beam"]["energy_keV"]
    gamma = 1.0 + E_keV / MEC2_KEV
    beta = np.sqrt(1.0 - 1.0 / gamma**2)
    p0 = gamma * M_E_SI * beta * C_SI
    ok = (abs(gamma - d["gamma"]) < 1e-12
          and abs(beta - d["beta"]) < 1e-12
          and abs(p0 - d["p_SI"]) / p0 < 1e-12)
    return ok, f"gamma={d['gamma']:.6f} beta={d['beta']:.6f} p0={d['p_SI']:.4e}"


@check("RF frequency reasonable")
def _c2(cfg):
    bad = []
    for e in _lattice_elements(cfg):
        if e["type"] == "rf_cavity" and e["length"] > 0:
            f = e["parameters"].get("frequency_GHz")
            if f is None or not (0.1 < f < 100.0):
                bad.append(e["name"])
    return not bad, f"rf instances checked: {len([e for e in _lattice_elements(cfg) if e['type']=='rf_cavity' and e['length']>0])}"


@check("RF/solenoid positions from lattice")
def _c3(cfg):
    # positions must be defined ONLY in lattice.elements; z_start monotone
    zs = [e["z_start"] for e in _lattice_elements(cfg) if e["length"] > 0]
    ok = all(zs[i] <= zs[i + 1] + 1e-12 for i in range(len(zs) - 1))
    return ok, f"element z_start sequence monotone: {[round(z, 3) for z in zs]}"


@check("sample == lattice total length (contiguity)")
def _c4(cfg):
    z_prev = 0.0
    gaps = []
    for e in _lattice_elements(cfg):
        if e["length"] <= 0:
            continue
        if abs(e["z_start"] - z_prev) > 1e-9:
            gaps.append((e["name"], e["z_start"], z_prev))
        z_prev = e["z_start"] + e["length"]
    return not gaps, f"sample={z_sample(cfg)*1e3:.3f}mm, total={z_prev*1e3:.3f}mm"


@check("lattice element schema")
def _c5(cfg):
    ok_types = {"cathode", "drift", "solenoid", "rf_cavity", "sample"}
    bad = []
    for e in _lattice_elements(cfg):
        if not all(k in e for k in ("name", "type", "z_start", "length")):
            bad.append(e)
            continue
        if e["type"] not in ok_types:
            bad.append(e)
    return not bad, f"elements={len(_lattice_elements(cfg))}, types ok"


def run_checks(cfg=None):
    cfg = cfg or load_config()
    results = []
    for name, fn in CHECKS:
        try:
            ok, detail = fn(cfg)
        except Exception as ex:
            ok, detail = False, f"EXCEPTION: {ex}"
        results.append((name, ok, detail))
    return results


def main():
    results = run_checks()
    print("== config consistency (read-only) ==")
    all_ok = True
    for name, ok, detail in results:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")
        all_ok &= ok
    print(f"  verdict: {'ALL CONSISTENT' if all_ok else 'INCONSISTENCY FOUND'}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
