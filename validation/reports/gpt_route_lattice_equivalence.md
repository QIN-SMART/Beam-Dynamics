# GPT Route — Lattice Single-Source Equivalence

Date: 2026-08-05 · structural refactor (no physics change)

## What changed
`GPT模拟/ued_beamline_v2.py` no longer hardcodes drift lengths or element
positions.  It now builds its OCELOT lattice exclusively from
`shared/beamline_config.yaml → lattice.elements` via the new generic builder:

```python
lat, rf_elems = build_lattice_from_shared(cfg, active_types)
```

- `active_types` implements the step semantics:
  step1 `{drift}` · step2 `{drift,solenoid}` · step3/4
  `{drift,solenoid,rf_cavity}`.
- Inactive length-bearing elements are kept as plain Drifts of the SAME
  length → total length and sample position are preserved for every step.
- Multi-instance support: any number of solenoids (per-instance B) and RF
  cavities (kick at each instance's own z_start, in lattice order, once each).
- RF kick: longitudinal p_oc = β0·δ_p conversion retained; transverse part
  gated by the shared `physics_switches.rf_transverse_kick`.
- SC unchanged: process added only for step 4 / sc_enabled, mesh/step from
  config.
- The module is now importable (`main()` guard) for tests.

## Verification
`validation/test_gpt_route_equivalence.py` — **PASS**:
- A. Geometry: element names/types/lengths/order, solenoid/RF counts, RF kick
  positions, total length, sample position — exact match to lattice.elements.
- B. Step routing: step1/2 zero RF kicks and zero active solenoids (runtime
  σ_δ unchanged ≈ 0.1e-3); step3/4 all solenoids + all RF kicks (σ_δ ≈ 2.94e-3);
  every step total length = 777 mm.
- C. Sample-plane regression vs validation OCELOT full route (all < 2 %):
  σ_x 0.65 % · σ_y 0.47 % · σ_z 0.43 % · σ_δ_p 0.43 % · ε_nx 0.07 % · ε_ny 0.02 %.
- D. Baseline regression (run_all + test_r56_convention): all four acceptance
  tests PASS with numbers identical to the frozen v0.10 baseline; AG arrays
  unchanged; config SHA unchanged.

## Data flow (new)
```
shared/beamline_config.yaml
   └─ lattice.elements ──→ GPT main route build_lattice_from_shared()
                          → validation route run_ocelot("full")
                          (both consume the same single geometry source)
```

## Files
- `GPT模拟/ued_beamline_v2.py` — refactored (only file with behavior change)
- `validation/test_gpt_route_equivalence.py` — new equivalence test
- `validation/reports/gpt_route_geometry.json`, `gpt_route_before_after.json`
- `validation/reports/review_summary_gpt_lattice.png`
- `CHANGELOG_gpt_lattice_single_source.md`

## Residual risks
- OCELOT first per-process beam (x/y/tau) is unseeded → sample regression has
  MC noise ~0.5–1 %; thresholds set at 2 %.
- GPT route step-1/2 unified output files still exercise the full lattice
  geometry (as drifts) — intended.
- No SC verification performed in this task (unchanged by design).
