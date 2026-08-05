# RF Model Physics Audit — AG vs OCELOT

Date: 2026-08-01 · Audit only, **no code modified**.
Scope: compare the AG RF chirp model (`AG/external_forces.py` `rf_chirp_rate`,
`make_rf_chirp_func`; `AG/beam_dynamics_6d.py` `envelope_ode`) with the OCELOT
RF kick model (`GPT模拟/ued_beamline_v2.py` `apply_rf_kick`; shared config).

Shared numbers used throughout (100 keV):
γ = 1.1957, β = 0.5482, β² = 0.3004, γ² = 1.4297, γ³ = 1.7096,
E₀ = γmc² = 611.0 keV, V = 30 kV, f = 2.856 GHz, k = 2πf/c = 59.86 m⁻¹,
φ = π, L_cav = 0.022 m, E_rf = V/L = 1.364 MV/m, m_ec² = 8.187e-14 J,
e·V = 4.806e-15 J.

---

## 1. Definition of δ (momentum deviation)

**OCELOT** — `δ ≡ (p − p₀)/p₀` (momentum deviation).
Confirmed by OCELOT source ("delta = (p-p0)/p0") and by the kick formula, which
carries a β² factor (the ΔE→Δp/p conversion, §4):

```
d_delta = (V_RF / (β²·E_total)) · sin(φ + k·z_phys)     # E_total = γmc²
```

**AG** — the shared config labels `sigma_delta` as "relative energy spread
ΔE/E", but AG's RF chirp coefficient contains a β² factor (see §5) which only
makes sense for a *momentum* deviation. So AG's δ is *nominally* ΔE/E but is
*dynamically* treated as a momentum deviation — an internal inconsistency.

**Audit statement:** the physically correct definition for a tracking code that
evolves `δ` through R56/drift and RF is the **momentum deviation δ = Δp/p₀**.
OCELOT is consistent; AG is ambiguous.

---

## 2. Chirp h — definition, equation, unit

Chirp = correlated energy/momentum slope across the bunch:
h = ∂δ/∂z, **unit m⁻¹**.

| | OCELOT | AG |
|---|---|---|
| quantity | momentum chirp ∂δ/∂z | envelope-focusing coefficient |
| formula | eV·k·cosφ / (β²·E₀) | −e·E_rf·k / (β²·γ³·m_ec²) |
| unit | **m⁻¹** | **m⁻²** |
| value | **−9.78 m⁻¹** | **−310.9 m⁻²** |

AG's coefficient is *not* ∂δ/∂z; it is used as `dC_zd/ds = h·σ_z²` in the 6D
envelope ODE (hence the extra 1/m from σ_z² and the m⁻² unit). This is a
different physical quantity, not merely a different number.

---

## 3. RF energy modulation ΔE(z)

Standing-wave cavity, phase φ w.r.t. the zero-crossing, bunch coordinate z
(head positive, z = 0 at the reference particle):

```
ΔE(z) = e·V·sin(φ + k·z)
dΔE/dz = e·V·k·cos(φ + k·z)
```

At the shared setting φ = π (max-chirp operating point):
```
ΔE(z) = −e·V·sin(k·z) ≈ −e·V·k·z        (linear ramp across the bunch)
dΔE/dz = −e·V·k = −1.80 MeV/m             (energy gradient across the bunch)
```
Note AG's `rf_chirp_rate` contains **no φ dependence at all** — it hard-codes
the max-chirp case (cosφ = −1). OCELOT retains sin(φ + kz) explicitly.

---

## 4. Conversion ΔE → Δp/p

From the relativistic dispersion relation E² = (pc)² + (m_ec²)²:

```
E dE = p c² dp  ⇒  dE = (pc²/E) dp = β·c·dp        (p = γmβc, E = γmc²)
⇒  Δp = ΔE/(βc)
⇒  δ = Δp/p₀ = ΔE / (βc · γmβc) = ΔE / (β²·γ·m_ec²) = ΔE / (β²·E₀)
```

So the RF kick in momentum terms is:

```
δ(z) = e·V·sin(φ + k·z) / (β²·γ·m_ec²)              (dimensionless)
```

This is exactly OCELOT's `V/(β²·E_total)` with `E_total = γmc²`.

---

## 5. Relativistic factors

In the correct momentum-deviation chirp:

```
h = ∂δ/∂z = e·V·k·cosφ / (β²·γ·m_ec²)
```
- **β²** — from the ΔE→Δp/p conversion (§4), `dE = βc·dp`.
- **γ¹** — from the total energy `E₀ = γmc²` in the same conversion.
- **γ²**, **γ³** — do NOT appear in the momentum chirp.

Where γ² legitimately appears in the *same* framework: the drift R56 for the
path-length coordinate, `R56 = L/γ²`, which AG uses for the z–δ correlation
growth (`dC_zd/dz = σ_δ²/γ²`). That γ² belongs to *time-of-flight*, not to the
RF energy kick.

AG's formula has **γ³ = γ¹ · γ²**. The γ¹ is the (correct) energy
normalization; the extra **γ² is spurious** — it is the R56/time-of-flight
factor applied to the wrong place, or an artifact of mixing a time-domain /
energy-deviation derivation into a momentum-deviation code.

| factor | OCELOT kick | AG chirp | source |
|---|---|---|---|
| β² | ✓ | ✓ | ΔE→Δp/p |
| γ¹ | ✓ | ✓ | E₀ = γmc² |
| γ² | — | ✓ (extra) | spurious (time-of-flight factor misplaced) |
| γ³ | — | ✓ | γ¹ × spurious γ² |

---

## 6. Why AG gives −310.9 m⁻² and OCELOT gives −9.78 m⁻¹

Two distinct effects combine:

**(a) Different quantity and units (1/L).**
OCELOT applies the whole kick at the cavity entrance (thin lens) → h is the
chirp across the bunch, m⁻¹. AG spreads the chirp over the cavity as a rate
per unit beamline length (`dC_zd/ds = h·σ_z²`, then `h·L` recovers the thin
lens) → the coefficient carries 1/L_cav.

**(b) Extra 1/γ² from the γ³ denominator.**
AG `h = −e·E_rf·k/(β²γ³mc²)`; OCELOT `h = eV·k·cosφ/(β²γmc²)` with E_rf = V/L.

Exact relation (verified numerically to 4 digits):

```
h_AG = h_OCELOT / (γ² · L_cav) = −9.78 / (1.4297 × 0.022) = −310.9 m⁻²
```
- factor 1/γ² = 0.70  ← γ³ vs γ¹ (§5)
- factor 1/L = 45.5   ← distributed-rate vs thin-kick convention (§6a)
- product 31.8, unit m⁻¹ → m⁻².

**(c) Additional ODE mis-integration (why σ_δ is ~10×, not ~32×).**
AG's `envelope_ode` keeps the thin-lens "sudden-kick" term `h²σ_z²` in
`dσ_δ/dz = h·(C_zd + h·σ_z²)/σ_δ` while simultaneously treating `h` as a
per-unit-length rate. The correct distributed-lens coupling is
`dσ_δ²/dz = 2h·C_zd`; the `h²σ_z²` term double-counts and inflates σ_δ.
Net AG σ_δ after the cavity = 28.1e-3 vs OCELOT 2.9e-3 (**≈9.6×**).

---

## 7. Minimal correction (proposal only — not applied)

Two one-line changes in AG:

**(7a) `external_forces.rf_chirp_rate` — fix the formula**
```
h = −e·E_rf·k·cosφ / (β²·γ·m_e·c²)          # γ³ → γ, insert cosφ
```
→ −444.5 m⁻² for this config.

**(7b) `beam_dynamics_6d.envelope_ode` — remove the spurious h² term**
```
dsigma_delta_dz = h_rate * beam.C_zd / sigma_delta    # was: h·(C_zd + h·σ_z²)/σ_δ
```

**Validation** (audit-only re-implementation of the longitudinal subsystem,
/tmp): the corrected pair reproduces the thin-lens reference to **0.04%**
(σ_δ = 4.936e-3 vs 4.938e-3); the current AG produces 33.2e-3 in the same
subsystem, matching the full-model ~28e-3.

**Residual σ_δ difference (4.9e-3 vs OCELOT 2.9e-3) is NOT the chirp** — it is
the σ_z at the cavity entrance: AG's σ_z has grown to ~505 µm because the shared
`epsilon_nz_mm_mrad: 0.2` is 10× the physically-consistent value
(βγ·σ_z·σ_δ = 0.02 µm). OCELOT's σ_z stays 298 µm (no tau transport, known GPT
limitation). With a consistent ε_nz, the corrected AG would match OCELOT.

**Checklist of proposed changes (for a future, separate task):**
1. `rf_chirp_rate`: γ³ → γ, multiply by cos(φ).  [external_forces.py]
2. `envelope_ode`: drop the `h²σ_z²` term in `dσ_δ/dz`.  [beam_dynamics_6d.py]
3. (separate issue) reconcile `epsilon_nz_mm_mrad` with σ_z·σ_δ (0.02, not 0.2).
