# CHANGELOG — Longitudinal Transport Fix (OCELOT tau propagation)

Date: 2026-08-01
Scope: OCELOT macroparticle route + validation framework.
Reference audit: `validation/reports/Longitudinal_transport_audit.md`

---

## 1. Original bug

`generate_parray(energy = E_keV / 1000.0)` (with `E_keV = 100`).

OCELOT interprets the `energy` argument as the **total energy in GeV**
(γ = E_GeV / m_e_GeV). Passing `0.1` therefore created a beam at
**0.1 GeV = 100 MeV** (γ = 195.7, β ≈ 1) instead of the intended
**100 keV** (γ = 1.196, β = 0.548).

Consequence: the OCELOT drift R56,
`R56 = −L/(β²γ²)`, was evaluated at γ ≈ 196, giving
`R56 ≈ −1.3e-5` — about **9×10⁴ too small** for a 100 keV beam
(correct value −1.164 for L = 0.5 m). As a result the native
longitudinal map `τ_f = τ_i + R56·δ` moved tau by an invisible amount and
`σ_z` appeared "frozen" along the beamline.

## 2. Physical cause

- OCELOT longitudinal coordinates: `tau = c·t [m]` (`rparticles[4]`),
  `δ = Δp/p₀` (`rparticles[5]`).
- The Drift DOES implement `τ_f = τ_i + R56·δ` with
  `R56 = −L/(β²γ²)` — this is native, correct OCELOT physics and SecondTM
  applies it (verified: with the right energy, a single δ = +1e-3 particle
  acquires Δτ = R56·δ exactly; a beam with σ_δ = 2e-3 spreads by
  σ_τ = |R56|·σ_δ = 2327 µm over 0.5 m).
- The "no tau transport" was therefore a **beam-initialisation unit error**,
  not a missing feature and not a SecondTM limitation.
- A related leftover was a **manual tau override**
  (`rparticles[4] += dτ/dz·dz·δ`) in the historical
  `GPT模拟/bug/ued_beamline.py`, which was a hidden workaround for the broken
  R56 and is now removed.

## 3. Correction

Use OCELOT's native longitudinal transport with the correct total energy:

```
E_total_GeV = (E_keV + 511.0) * 1e-6        # kinetic + rest mass, GeV
sigma_tau   = sigma_z / beta                # τ = c·t [m]; σ_τ = σ_z/β
```

Applied to every OCELOT `generate_parray` call:
- `validation/backend.py` (framework OCELOT driver)
- `GPT模拟/ued_beamline_v2.py` (main route)
- `GPT模拟/bug/ued_beamline.py` (route; manual tau override **removed**)
- `GPT模拟/ocelot_beamline.py` (legacy route)
- historical module benchmarks: `drift/`, `rf/`, `solenoid/`,
  `space_charge/` (beam init only; `rf` benchmark tau-unit in its kick
  `z_phys = −β·c·τ → −β·τ` fixed for consistency)

No OCELOT physics was modified. No custom R56 correction was introduced.
No hidden manual tau propagation remains.

## 4. Validation result

Re-ran `validation/test_drift.py`, `validation/test_rf.py`,
`validation/run_all.py`.

| test | before (σ_z at end) | after | check |
|------|---------------------|-------|-------|
| drift — OCELOT | 299 µm (frozen) | **σ_z evolves physically** via native R56 | PASS |
| drift — AG vs OCELOT σ_δ | 0.04 % | 0.04 % | maintained |
| rf — AG vs OCELOT σ_δ | — | (native compression enabled) | maintained |

Acceptance:
- drift `σ_z` now changes physically (native `τ_f = τ_i + R56·δ`),
- AG/OCELOT `σ_delta` agreement maintained,
- no hidden coordinate correction (manual override removed).

## Files changed

- `validation/backend.py`
- `GPT模拟/ued_beamline_v2.py`
- `GPT模拟/bug/ued_beamline.py` (manual tau override removed)
- `GPT模拟/ocelot_beamline.py`
- `GPT模拟/drift/benchmark_drift.py`
- `GPT模拟/rf/benchmark_rf_drift.py`
- `GPT模拟/solenoid/benchmark_solenoid.py`
- `GPT模拟/solenoid/benchmark_solenoid_physical.py`
- `GPT模拟/space_charge/benchmark_space_charge_drift.py`
- `CHANGELOG_longitudinal_transport.md`
