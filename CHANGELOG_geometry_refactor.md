# CHANGELOG — Beamline Geometry Refactor (lattice = single source of truth)

Date: 2026-08-01
Scope: shared config schema + framework/route adapters.
**No physics change. No parameter tuning. No model modification.**

---

## 1. Data-flow change

**Before:** geometry was duplicated in three places —
`lattice.elements` (list format `[z, L, type]`), `solenoid:` (B, length,
z_start) and `rf_cavity:` (f, V, φ, length, z_start).  Different backends
read different sources (e.g. OCELOT kick timing read `rf_cavity.z_start_m`
while the AG adapter read the lattice) — moving an element in the lattice
did not move it for OCELOT.

**After:** `lattice.elements` is the ONLY geometry definition.  Each element
is a dict `{name, type, z_start, length, parameters}` and carries its own
physical parameters:

```yaml
- name: solenoid1
  type: solenoid
  z_start: 0.100
  length: 0.060
  parameters: {B_field_T: 0.05}
- name: rf1
  type: rf_cavity
  z_start: 0.400
  length: 0.022
  parameters: {frequency_GHz: 2.856, voltage_kV: 30.0, phase_rad: 3.1416}
```

The top-level `solenoid:` / `rf_cavity:` sections are **removed**.

**Multi-instance support** (removed `next(element)` single-element
assumptions): multiple solenoids and multiple RF cavities are allowed.
- AG adapter (`validation/backend.py::run_ag`): builds one ExtFieldRegion per
  active element; applies `apply_rf_thin_lens(H)` per RF cavity in z order.
- OCELOT adapter (`validation/backend.py::run_ocelot`): builds one ocelot
  element per lattice element; applies the thin-kick + transverse kick at each
  RF cavity entrance; inactive elements keep their length as plain drifts so
  the total beamline length is preserved.
- `shared/params.py`: new helpers `elements_of_type`, `first_of_type`,
  `elem_geom`, `elem_params`, `flat_elem` (backward-compatible flat view for
  the route scripts); `SolenoidParams`/`RFParams` now read from the lattice.
- Route scripts (`AG/run_shared.py`, `GPT模拟/ued_beamline_v2.py`) and the
  section tests parse the lattice only.

## 2. No physics change

Parameter key names are preserved (`B_field_T`, `frequency_GHz`,
`voltage_kV`, `phase_rad`); all physics kernels (AG core, OCELOT tracking)
are untouched.

## 3. Validation results (before/after on test_drift/test_solenoid/test_rf)

| section | AG | OCELOT | verdict |
|---|---|---|---|
| drift | max rel diff **2.8e-16** (bit-identical) | ≤0.66% (pre-existing unseeded-MC noise of the first `generate_parray` per process; not caused by refactor) | PASS |
| solenoid | **9.0e-16** | **1.1e-12** (bit-identical) | PASS |
| rf | 0.00e+00 in smooth regions (σ_δ, ε, σ_x); σ_z 3e-4 sampling artifact | **2.0e-14** (bit-identical) | PASS |

All three tests still PASS with unchanged acceptance numbers
(drift σ_x ≈ 0.3%, solenoid 0.40%, rf σ_δ 0.24% / σ_x 2.8%).

**Note:** an earlier attempt to make the OCELOT beam fully reproducible by
seeding `np.random` before `generate_parray` was reverted — the double
`seed(42)` made x and px use identical random draws (spurious x–x′
correlation, ε→0).  The pre-existing RNG structure is preserved.

## Files changed
- `shared/beamline_config.yaml` (new lattice schema; removed solenoid/rf sections)
- `shared/params.py` (lattice helpers, multi-instance)
- `validation/backend.py` (lattice-only parsing, multi-instance regions/elements/kicks)
- `validation/test_solenoid.py`, `validation/test_rf.py` (lattice reads)
- `AG/run_shared.py` (lattice-only regions + multi-RF thin lens)
- `GPT模拟/ued_beamline_v2.py` (flat_elem views)
