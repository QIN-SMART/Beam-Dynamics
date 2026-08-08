#!/usr/bin/env python3
"""
Level-1 test: shared config schema & consistency (v0.12, Phase 7/10).

Read-only — never modifies config values.  Runs every check in
validation/config_check.py and gates the suite on them.

Usage:  python3 validation/test_config_consistency.py
"""

import os
import sys

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
for p in (_THIS_DIR, os.path.dirname(_THIS_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from validation.config_check import run_checks  # noqa: E402


def main():
    results = run_checks()
    print("\n== Level-1: config schema & consistency ==")
    ok_all = True
    for name, ok, detail in results:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
        ok_all &= ok
    print(f"  config consistency: {'PASS' if ok_all else 'FAIL'}")
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
