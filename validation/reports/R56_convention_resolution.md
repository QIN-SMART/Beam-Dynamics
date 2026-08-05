# R56 / Longitudinal Coordinate Convention — Resolution

Date: 2026-08-01 · AUDIT + executable characterization
Executable: `validation/test_r56_convention.py`
Data: `validation/reports/r56_convention_results.json`, `r56_convention_plot.png`
**No production physics code modified.**

---

## Part 1 — Explicit definitions

Reference particle at K0 = 100 keV:

| symbol | value | definition |
|---|---|---|
| K0 | 100 keV | reference kinetic energy |
| E0 | K0 + m_e·c² = 611.0 keV | reference total energy |
| γ0 | 1.19569 | γ0 = 1 + K0/(m_ec²) |
| β0 | 0.54822 | β0 = √(1−1/γ0²) |
| p0 | 1.7901e-22 kg·m/s | p0 = β0·γ0·m_e·c |
| δ_p | (p−p0)/p0 | **momentum deviation** (AG-native) |
| Δt | t_particle − t_reference | laboratory arrival-time difference |
| Δct | c·Δt | laboratory c-time difference |
| Δz | −β0·c·Δt | co-moving spatial bunch coordinate (head positive) |
| p_OC | dE/(c·p0) | **OCELOT-native rparticles[5]** (energy-normalized) |

Signs: δ_p > 0 ⇒ higher momentum ⇒ arrives earlier (Δt < 0) ⇒ ahead in the
co-moving coordinate (Δz > 0).

**"tau" (installed OCELOT 26.06.1):** `p.tau()` returns `rparticles[4]`, and
`generate_parray` documents it as `"tau" = c*t`
(`ocelot/cpbd/beam/generator.py:50`), i.e. tau = c·Δt in metres.  Note the
ParticleArray class docstring writes the row as `ds = c*tau`
(`ocelot/cpbd/beam/particle.py:20`).

---

## Part 2 — Exact relativistic reference (independent, no AG/OCELOT)

Single drift L = 0.5 m, δ ∈ {−1e-3, −5e-4, −1e-4, 1e-4, 5e-4, 1e-3}.
For each δ: p = p0(1+δ), γ(δ)=√(1+(p/m_ec)²), β(δ)=p/(γm_ec),
t = L/(βc), then Δt, Δct, Δz = −β0c·Δt.  Slopes fitted linearly:

| slope | fitted | analytic | rel. err |
|---|---|---|---|
| dΔt/dδ | −2.127910e-9 s | −L/(β0·c·γ0²) = −2.127908e-9 s | 1.0e-6 |
| dΔct/dδ | −0.637931 m | −L/(β0·γ0²) = −0.637931 m | 1.0e-6 |
| dΔz/dδ | +0.349727 m | +L/γ0² = +0.349727 m | 1.0e-6 |

The analytic linear forms in the task are confirmed as the exact small-δ
slopes (rel. err ≈ 1e-6, from finite-δ curvature).

---

## Part 3 — OCELOT raw-coordinate characterization (v26.06.1)

- Beam: 7 particles, `rparticles[4]=0`, `rparticles[5]=δ` (raw, as-fed),
  E_tot = (100+511)·1e-6 GeV = 6.11e-4 GeV; tracked through one native
  `Drift(l=0.5)` with SecondTM (5 steps of 0.1 m).
- Transfer-map R56 (installed): **R56_tm = −1.163624**.
- Raw slope d(Δτ_raw)/dδ = **−1.163624** — exactly the transfer map.

| comparison | rel. err vs raw slope |
|---|---|
| vs −L/(β0²γ0²) = −1.163624 (OCELOT formula) | **−1.25e-5** ✓ exact |
| vs −L/(β0·γ0²) = −0.637931 (exact c·t slope for δ_p) | +0.824 (= 1/β0 − 1) |

Raw OCELOT tau transport is internally exact for its OWN coordinate pair
(τ=c·t, p=dE/(c·p0)): since p_OC = β0·δ_p, the slope d(ct)/d(p_OC) =
d(ct)/dδ_p · 1/β0 = −L/(β0γ0²)/β0 = **−L/(β0²γ0²)**, which is precisely what
the installed `uni_matrix` implements.

---

## Part 4 — Installed OCELOT source conventions

| item | file:line | content |
|---|---|---|
| tau definition | `ocelot/cpbd/beam/generator.py:50` | `"tau" = c*t` |
| p definition | `ocelot/cpbd/beam/generator.py:51` | `'p' is canonical momentum E/(c*p0)` |
| p definition | `ocelot/cpbd/beam/beam.py:211` | `p=dE/(c*p0)` |
| ParticleArray rows | `ocelot/cpbd/beam/particle.py:20` | `(ds = c*tau, p = dE/(p0*c))` |
| Drift R56 | `ocelot/cpbd/r_matrix.py:81` | `r56 -= z/(beta*beta)*igamma2` ⇒ R56 = −L/(β²γ²) |
| Solenoid R56 | `ocelot/cpbd/elements/solenoid_atom.py:74` | `r56 -= l/(beta*beta*gamma2)` (same convention) |

The installed formula assumes the coordinate pair **(τ=c·t, p=dE/(c·p0))** —
an **energy-normalized longitudinal momentum**, NOT Δp/p0 and NOT ΔE/E.
It is exact for that pair (Part 3).  It is not an ultrarelativistic
approximation; at β→1 the conventions coincide.

---

## Part 5 — AG convention

AG tracks the co-moving spatial coordinate Δz with the map

```
Δz_AG = (L/γ0²)·δ        (δ = δ_p, momentum deviation)
```
entered via `envelope_ode`: `dC_zδ/dz = σ_δ²/γ²` (⇒ R56_z = L/γ²), with σ_δ
fed from `shared/beamline_config.yaml` `initial_distribution.sigma_delta`.

Against the exact Δz = −β0·c·Δt for the same δ_p scan:
- max relative error **1.15e-3** (0.12 % — pure O(δ²) nonlinearity of the
  exact kinematics),
- quadratic curvature −1.3e-2 m (linear map is an excellent approximation in
  this δ range).

AG's spatial map is exact to first order; it is the momentum-deviation
convention (δ_p).

---

## Part 6 — Coordinate-transformation closure test

Formal derivation from the DOCUMENTED definitions only:

```
Δt      = Δτ / c                 (τ = c·t)
δ_p     = p_OC / β0              (p_OC = dE/(c·p0) = β0·δ_p)
Δz      = −β0·c·Δt = −β0·Δτ      (definition of the co-moving coordinate)
```

Closure (using raw OCELOT transport Δτ = R56_tm·p_OC):

```
Δz_formal = −β0 · (−L/(β0²γ0²)) · p_OC = L·p_OC/(β0·γ0²) = L·δ_p/γ0²  = exact
```

Test results (all 6 δ values, both signs):
- **closure after formal conversion: max rel. residual 2.1e-3** (0.2 %) ✓
- naive framework feeding (raw value treated as δ_p):
  residual **0.826**, matching the predicted **1/β0 − 1 = 0.824** — this is
  the observed σ_z factor.

The transformation is derived from the definitions, reproduces the exact
result, preserves signs for ±δ, and contains no fitted coefficient.

---

## Part 7 — Classification

### Primary category: **B — Input delta-variable convention mismatch**

The two models are each internally consistent; the framework feeds the SAME
raw number into two different physical variables without conversion:

- AG's `sigma_delta` is used as **δ_p = Δp/p0** (its R56_z = L/γ² is exact
  for δ_p).
- OCELOT's `rparticles[5]` is **p_OC = dE/(c·p0) = β·δ_p** (its R56 =
  −L/(β²γ²) is exact for that pair).

Consequence: with the same input value v, OCELOT's physical momentum
deviation is v/β0 = 1.82·v — hence its tau/σ_z transport is 1/β0 faster in
the comparison, which is exactly the observed residual.

**Evidence:** Part 3 raw slope equals the transfer map to 1e-5; Part 6 naive
residual 0.826 = 1/β0−1 reproduced exactly; formal conversion closes to 0.2%.

**Secondary (E-like) aspect:** the previously reported σ_δ "agreement"
compares OCELOT raw p_OC with AG δ_p — different observables that coincide
only because the same input number is conserved in drift and the RF-kick
amplitude dominates the initial spread.

**Affected:** σ_z, σ_t (compression curves); semantic meaning of the reported
OCELOT σ_δ (raw p_OC vs δ_p); the physical δ_p represented by OCELOT's beam
(input feeding).

**Unaffected:** σ_x, σ_y, ε_nx, ε_ny (transverse maps are energy-independent
in this lattice); RF kick amplitude σ_δ (initial spread ≪ kick, masked);
AG-side physics; OCELOT tracking.

**Smallest correction location:** the OCELOT adapter input/output conversion
(`validation/backend.py::run_ocelot`), not AG core, not OCELOT core.

**Expected regression risk:** low — transverse tests unaffected; σ_δ reported
values preserved; only σ_z/σ_t (currently diagnostic) change toward agreement.

**Required tests before implementation:** drift (σ_z should move from ~2 %
to ~0.1 %), rf (kick amplitude unchanged), full-beamline (σ_z diag → agreement),
solenoid (unchanged).

---

## Part 8 — Proposed implementation plan (NOT implemented)

Recommended: adopt **δ_p = Δp/p0 as the canonical project momentum deviation**
(the config value is kept; its comment is relabelled from "ΔE/E" to
"momentum deviation Δp/p0" — comment-only).

**Files:** `validation/backend.py` — `run_ocelot` only (3 one-line changes).
(AG needs no change; it already uses δ_p.)

| quantity | old (current) | new |
|---|---|---|
| initial spread | `rparticles[5] = N(0, σ_δ)` | `rparticles[5] = N(0, β0·σ_δ)` (p_OC = β0·δ_p) |
| RF kick | `rparticles[5] += K·sin(φ+kz)` | `rparticles[5] += β0·K·sin(φ+kz)` |
| reported σ_δ | `std(p.p())` | `std(p.p())/β0` (= δ_p) |
| z_phys in kick | −β0·τ | −β0·τ (already correct: Δz = −β0·Δτ) |

Old/new definitions: OCELOT raw p_OC = dE/(c·p0) is kept internally;
conversion happens only at the adapter boundary (input feed ×β0, output ÷β0).

**Tests that must remain unchanged (pass):** solenoid (transverse), drift
transverse, rf kick amplitude (σ_δ ≈ 2.95e-3 both), full-beamline σ_x/σ_y/ε.

**Tests whose expected output changes:** drift σ_z (R56-driven growth now
matches AG), full-beamline σ_z/σ_t (R56 diagnostic item should CLOSE to
≈ agreement; if so, it can be reclassified from "open item" to "resolved").

**Rollback:** revert the 3 lines in `run_ocelot`; no core or config changes.

If instead the config's "ΔE/E" meaning must be preserved literally, the same
plan applies with δ_p = σ_δ/β0² and p_OC = σ_δ/β0 — both adapters would need
the conversion; not recommended (δ_p is the standard accelerator convention).

---

## Acceptance checklist (this audit)

- [x] exact relativistic reference computed independently (Part 2)
- [x] positive and negative δ tested
- [x] OCELOT raw variables recorded before conversion (Part 3)
- [x] installed OCELOT source convention identified with citations (Part 4)
- [x] AG convention explicitly stated (Part 5)
- [x] 1/β difference explained: input δ-variable convention mismatch (B),
      OCELOT p = dE/(c·p0) = β0·δ_p
- [x] no production physics code modified; no empirical correction factor
