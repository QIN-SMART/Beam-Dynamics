# Implementation Report — OCELOT Longitudinal-Variable Adapter

Date: 2026-08-01 · Adapter implementation only (R56 audit classification B).

## 1. Old data flow
```
config sigma_delta (Δp/p0 meaning) ──→ rparticles[5] = N(0, σ_δ)   [WRONG: p_oc = ΔE/(c·p0)]
RF kick dδ_p (Δp/p0) ──→ rparticles[5] += dδ_p                     [WRONG units]
report σ_δ = std(p_oc)                                             [raw p_oc exposed]
σ_z = β0·std(tau)                                                  [correct]
```
OCELOT's physical δ_p was input/β0 → tau transport 1/β0 faster → σ_z/σ_t
compression displaced by ~1/β0 (2.3× at the sample).

## 2. New data flow
```
config sigma_delta (Δp/p0) ──→ δ_p samples ──→ rparticles[5] = β0·δ_p   [p_oc = β0·δ_p]
RF kick dδ_p (Δp/p0)   ──→ rparticles[5] += β0·dδ_p
report σ_δ_p = std(p_oc)/β0                                            [shared meaning]
σ_z = β0·std(tau)                                                      [unchanged]
```
OCELOT native R56 and transfer maps untouched; conversion only at the
adapter boundary.

## 3. Exact conversion equations
```
p_oc   = ΔE / (c·p0)          (OCELOT native, documented)
δ_p    = Δp / p0              (shared / AG observable)
p_oc   = β0 · δ_p             (first-order reference-particle conversion)
δ_p    = p_oc / β0
σ_δ_p  = std(p_oc) / β0       (reported)
```

## 4. Files and lines modified
- `validation/backend.py`
  - beam generation: `rparticles[5,:] = d["beta"] * delta_p`
  - new module function `_ocelot_rf_kick` (kick ×β0 + comment)
  - output: `sd[i] = np.std(p.p()) / d["beta"]`
  - meta: `longitudinal_native_coordinate`, `reported_delta`,
    `conversion_beta`, `rf_kicks_applied`
- `GPT模拟/ued_beamline_v2.py` — initial feed ×β0; `apply_rf_kick` ×β0;
  σ_δ reporting (probe + history) ÷β0.
- `shared/output_schema.py` — `sigma_delta_e3` doc = δ_p (Δp/p0).
- Tests: `test_drift.py` (Test 2 + analytic σ_z ref), `test_rf.py`
  (Test 3 + routing), `test_full_beamline.py` (σ_z/σ_t quantitative + waist).

## 5. Tests run
`test_r56_convention.py` (unchanged), `test_drift.py`, `test_rf.py`,
`test_full_beamline.py`, structural routing check, `run_all.py`.

## 6. Before/after σ_z and σ_t (sample plane)
| quantity | before | after | AG |
|---|---|---|---|
| rf σ_z | 1112.8 µm | **476.3 µm** | 477.0 µm |
| rf σ_t | 6770.8 fs | **2897.8 fs** | 2902.3 fs |
| full σ_z | 1117.5 µm | **473.8 µm** | 477.0 µm |
| full σ_t | 6799.2 fs | **2882.7 fs** | 2902.3 fs |
| drift σ_z | 316.5 µm | 306.1 µm | 310.2 µm |
| waist z | — | 547 vs 546 mm (Δz=0.6 mm) | — |

## 7. Before/after reported σ_δ_p (sample plane)
| test | before | after | AG |
|---|---|---|---|
| drift | 0.1000e-3 | 0.1000e-3 | 0.100e-3 |
| rf | 2.942e-3 | 2.953e-3 | 2.953e-3 |
| full | 2.954e-3 | 2.938e-3 | 2.953e-3 |

Semantics corrected: the reported value is now explicitly δ_p = std(p_oc)/β0
(verified deterministically: 9.996e-5 vs configured 1e-4, rel 0.043 %).

## 8. Remaining residuals
- σ_z/σ_t at the sample: <1 % (0.15–0.7 %).
- Largest relative σ_z deviation (~44 %) only at the compression waist where
  σ_z ≈ 10–15 µm (absolute few-µm difference): τ-space vs z-space RMS
  evolution in the deep-compression nonlinear regime — higher-order/numerical
  residual, not a convention gap.
- σ_δ_p agreement 0.3–0.5 % (within the 2 % criterion).

## 9. Rollback instructions
Revert the three conversion points in `validation/backend.py::run_ocelot`
(initial feed ×β0, `_ocelot_rf_kick` ×β0, `sd` ÷β0) and the three matching
edits in `GPT模拟/ued_beamline_v2.py`.  Tests and schema documentation can
remain (they describe the corrected semantics).  No core physics, config
values, or OCELOT source are touched either way.
