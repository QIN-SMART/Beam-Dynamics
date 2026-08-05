# CHANGELOG — OCELOT Longitudinal-Variable Adapter (R56 convention fix)

Date: 2026-08-01
Reference audit: `validation/reports/R56_convention_resolution.md`
(classification B — input δ-variable convention mismatch).

## What changed

OCELOT's native sixth coordinate is `p_oc = ΔE/(c·p0)`, NOT Δp/p0 (documented
in the installed OCELOT source; its R56 = −L/(β²γ²) is exact for that pair).
The framework fed the shared momentum deviation directly into that coordinate,
making OCELOT's physical δ_p = input/β0 (the 1/β discrepancy).  This change
converts at the adapter boundary only — AG core, OCELOT core, transfer maps,
and config values are untouched.

| step | before | after |
|---|---|---|
| initial spread | `rparticles[5] = N(0, σ_δ)` | `rparticles[5] = β0 · N(0, σ_δ)` |
| RF kick | `rparticles[5] += dδ_p` | `rparticles[5] += β0 · dδ_p` |
| reported σ_δ | `std(p_oc)` | `std(p_oc)/β0` (= δ_p) |
| τ / σ_z | `σ_z = β0·std(tau)` | unchanged (verified) |

## Files modified
- `validation/backend.py` — `run_ocelot`: initial feed, RF kick extracted to
  `_ocelot_rf_kick` (×β0 + comment), output ÷β0, meta gains
  `longitudinal_native_coordinate`, `conversion_beta`, `rf_kicks_applied`.
- `GPT模拟/ued_beamline_v2.py` — same three conversions (initial feed,
  `apply_rf_kick` ×β0, σ_δ reporting ÷β0).
- `shared/output_schema.py` — `sigma_delta_e3` documented as δ_p = Δp/p0.
- `validation/test_drift.py` — Test 2 (initial σ_δ_p semantic, deterministic)
  + longitudinal analytic reference √(σ_z0²+(z·σ_δ_p/γ²)²).
- `validation/test_rf.py` — Test 3 (RF kick semantics on a controlled set) +
  structural section-routing test (drift/solenoid 0 kicks, rf/full N_rf).
- `validation/test_full_beamline.py` — σ_z/σ_t promoted to quantitative
  (5% threshold) + compression-waist position check (5 mm tolerance).
- `validation/run_all.py` — labels updated.
- `validation/reports/r56_adapter_before_after.json` — before/after metrics.

## Validation results

| test | before | after |
|---|---|---|
| test_r56_convention | raw slope = R56_tm (unchanged, independent) | **unchanged** ✓ |
| test_drift | σ_z dev 1.99% | 1.35%; σ_δ_p semantic 0.043% PASS |
| test_rf | σ_z sample 1112.8 vs AG 477.0 (2.3×) | **476.3 vs 477.0 (0.15%)**; kick residual 2.7e-14; routing PASS |
| test_full_beamline | σ_z 57% diag | **0.68% PASS** (σ_z, σ_t, waist Δz=0.6 mm) |

Invariants: AG arrays bit-identical (0.0 rel diff), config SHA unchanged,
transverse within MC tolerance, number/order of RF kicks unchanged
(routing verified).

## Remaining residual (explained)
σ_z curves agree to <1 % at the sample; the largest relative deviation (~44 %)
occurs only at the compression waist where σ_z → 10–15 µm (absolute difference
a few µm): the τ-space (OCELOT) and z-space (AG) RMS evolutions differ slightly
in the deep-compression nonlinear regime.  This is a higher-order/numerical
residual, not the previous 1/β convention gap.

## Rollback
Revert the three conversion points in `validation/backend.py::run_ocelot`
(initial feed, `_ocelot_rf_kick`, `sd` extraction) and the matching three in
`GPT模拟/ued_beamline_v2.py`; tests and schema docs can remain (they document
the corrected semantics).
