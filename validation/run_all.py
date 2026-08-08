#!/usr/bin/env python3
"""
Run the full validation suite (drift, solenoid, rf) and print a summary.

Usage:  python3 validation/run_all.py
"""

import os
import sys

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
for p in (_THIS_DIR, os.path.dirname(_THIS_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from validation import (test_config_consistency, test_drift, test_solenoid,  # noqa: E402
                        test_rf, test_full_beamline, test_gpt_route_equivalence)


def main():
    print("=" * 66)
    print("  UED framework validation suite  (AG vs OCELOT)")
    print("=" * 66)

    results = {}
    for name, mod in (("config_schema", test_config_consistency),
                      ("drift", test_drift), ("solenoid", test_solenoid),
                      ("rf", test_rf), ("full_beamline", test_full_beamline),
                      ("gpt_route", test_gpt_route_equivalence)):
        print("\n" + "#" * 66)
        print(f"  TEST: {name}")
        print("#" * 66)
        results[name] = mod.main()

    print("\n" + "=" * 66)
    print("  SUMMARY")
    print("=" * 66)
    labels = {"config_schema": "Level-1 config consistency",
              "drift": "transverse+σ_z",
              "solenoid": "coupling resolved",
              "rf": "thin-lens; σ_z R56 resolved",
              "full_beamline": "all 7 metrics quantitative (R56 adapter)",
              "gpt_route": "lattice single-source equivalence"}
    for name, rc in results.items():
        verdict = "PASS" if rc == 0 else "FAIL"
        print(f"  {name:<13s} {verdict}  ({labels.get(name, '')})")
    print(f"\n  reports -> {os.path.join(_THIS_DIR, 'reports')}")
    print(f"  checkpoints -> {os.path.join(_THIS_DIR, 'CHECKPOINTS.md')}")


if __name__ == "__main__":
    main()
