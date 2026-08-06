# CHANGELOG — GPT Main Route: Lattice Single Source

Date: 2026-08-05
Type: structural refactor, no physics change.
Baseline: v0.10-noSC-longitudinal-validated (config SHA dd8ada3d4cb2).

## Problem
`GPT模拟/ued_beamline_v2.py::build_lattice()` hardcoded drift lengths and
positions (0.100/0.240/0.355/0.777 m), duplicating `lattice.elements`.
Risks: GUI/YAML changes not reaching the main route; the two routes could
simulate different beamlines; no multi-solenoid/multi-RF support.

## Change
- New `build_lattice_from_shared(cfg, active_types)` builds the OCELOT
  lattice exclusively from `lattice.elements`; inactive length-bearing
  elements are kept as equal-length Drifts (total length/sample preserved).
- Step semantics via `STEP_ACTIVE`:
  1 `{drift}`, 2 `{drift,solenoid}`, 3/4 `{drift,solenoid,rf_cavity}`.
- Multi-instance: per-element B (solenoid), per-instance RF kick at its own
  z_start (once, in order).  RF transverse kick gated by the shared switch.
- Longitudinal conversion retained: p_oc = β0·δ_p in, δ_p = p_oc/β0 out.
- SC unchanged: process added only for step 4 / sc_enabled.
- Module made importable (`main()` guard) for the equivalence test.

## Verification
`validation/test_gpt_route_equivalence.py` PASS:
geometry exact match; step routing (0/0/N_rf/N_rf kicks, runtime σ_δ gating);
sample-plane regression vs validation route all < 2 %
(σ_x 0.65 %, σ_z 0.43 %, σ_δ_p 0.43 %); baseline run_all + r56 characterization
unchanged; AG arrays bit-identical; config SHA unchanged.

## Rollback
```
git checkout v0.10-noSC-longitudinal-validated -- GPT模拟/ued_beamline_v2.py
```
(or reset to the pre-refactor commit).
