# RF Model Standardization — Decision & Physics Notes

Date: 2026-08-01 · **Physics notes only — no code modified.**
Applies to: AG (continuous cavity model) vs OCELOT (thin RF kick).
Goal: make AG and OCELOT physically comparable in the shared framework.

---

## 1. Should the project adopt thin-lens RF as the common model?

**Decision: YES — thin-lens (impulse) RF is adopted as the canonical shared
model.**

Rationale:
- **Standard physics.** UED RF bunchers are conventionally treated as thin
  longitudinal lenses: the cavity imprints a momentum modulation `δ(z)` across
  the bunch in a single impulse; the downstream drift compresses it. OCELOT
  already implements this; AG has an (existing) thin-lens path
  (`apply_rf_thin_lens`) alongside its continuous chirp path.
- **Validity here.** σ_z = 300 µm ≪ λ_RF/4 = 26 mm (bunch samples a small,
  nearly-linear phase region); the cavity acts as a lens, not an accelerator
  (φ = π ⇒ reference-particle net gain ≈ 0).
- **Comparability.** The two codes can only agree if they use the SAME RF
  definition. The continuous model is currently implemented with an incorrect
  coefficient (γ³) and a mis-integrated ODE term, and adds transverse-RF and
  acceleration effects OCELOT does not have — making direct comparison
  meaningless.
- Thin-lens caveat (documented, not blocking): cavity transit angle
  θ = ω·L/(βc) ≈ 2.4 rad ⇒ a transit-time factor TT ≈ sin(θ/2)/(θ/2) ≈ 0.78
  applies to the *effective* voltage. Both codes must use the same convention;
  TT is an optional *shared* refinement, not a per-code difference.

---

## 2. Standardized definitions (thin-lens)

**δ (momentum deviation)**
```
δ ≡ (p − p₀)/p₀ = Δp/p₀
```
Conversion from energy gain (from E² = (pc)² + (m_ec²)² ⇒ dE = βc·dp):
```
δ = ΔE/(β²·γ·m_e·c²) = ΔE/(β²·E₀),   E₀ = γm_ec²
```

**Chirp h**  (unit m⁻¹)
```
h = ∂δ/∂z = e·V·k·cosφ / (β²·γ·m_e·c²) = eV·k·cosφ/(β²E₀)
shared config:  h = 30000·59.86·(−1)/(0.3004·611000) = −9.78 m⁻¹
```

**RF kick equation**  (thin lens at z_rf)
```
δ_out(z) = δ_in(z) + K·sin(φ + k·z),     K = eV/(β²E₀) = 0.1637
```
Linear (buncher) form near bunch centre:
```
δ_out ≈ δ_in + h·z
```
Envelope-map form (used by both AG's `apply_rf_thin_lens` and the OCELOT
moment equations):
```
σ_z  → σ_z                 (unchanged by a thin lens)
C_zδ → C_zδ + h·σ_z²
σ_δ² → σ_δ² + 2h·C_zδ + h²·σ_z²
```

---

## 3. AG terms to remove / change

The AG RF region currently injects **three** effects. For thin-lens
standardization:

| # | AG term | location | action |
|---|---------|----------|--------|
| 1 | Continuous chirp coefficient `h = −e·E_rf·k/(β²γ³·m_ec²)` (m⁻², γ³, no φ) | `external_forces.rf_chirp_rate` / `build_rf_chirp_func` | **Remove** from the standardized run; use the existing `beam_dynamics_6d.apply_rf_thin_lens(beam, H)` with `H = eV·k·cosφ/(β²E₀)` passed from the framework |
| 2 | ODE chirp integration `dσ_δ/dz = h(C_zδ + h·σ_z²)/σ_δ` (spurious h²σ_z² term) | `beam_dynamics_6d.envelope_ode` | **Remove** (replaced by the impulse map) |
| 3 | Transverse RF force `Fe_x,y = −η(β·k·E_rf)/(2γβ²)·σ_x,y` (dE_cdt = k·E_rf) | `external_forces.rf_cavity_force` | **Remove** (zero `dE_cdt`; OCELOT thin kick has no transverse effect) |
| 4 | Acceleration `dγ/dz = η·E_rf` (phase-independent, +30 keV) | `external_forces.build_gamma_prime_func` / `rf_acceleration_gradient` | **Remove** (at φ = π the net gain ≈ 0; OCELOT leaves energy unchanged) |
| 5 | (config) `epsilon_nz_mm_mrad: 0.2` — 10× the consistent value | `shared/beamline_config.yaml` | **Change to 0.02** = βγ·σ_z·σ_δ (else σ_z inflates to 505 µm at the cavity and σ_δ scales with it) |

Note: AG already ships a correct thin-lens envelope map (`apply_rf_thin_lens`),
so the standardization reuses existing code rather than adding new physics.

---

## 4. Migration plan

Framework-adapter driven — **no modification of AG core or OCELOT tracking**:

1. **Shared RF module** (framework): single source of `δ`, `h`, `K`, and the
   thin-lens envelope map; `H = eV·k·cosφ/(β²E₀)` computed from the shared
   config.
2. **OCELOT backend**: already thin-kick; assert it uses the shared `K`, `φ`,
   `k` (no change needed).
3. **AG backend adapter** (`validation/backend.run_ag`, `AG/run_shared.py`):
   - propagate 0 → z_rf with the RF-region forces zeroed (`dE_cdt = 0`, no
     `gamma_prime` from the rf region);
   - apply `apply_rf_thin_lens(beam, H)` at z_rf (existing function, H from
     the shared module);
   - propagate z_rf → end.
4. **Fix ε_nz** in the shared config (0.02) so σ_z at the cavity matches
   (303 µm), removing the σ_δ scale mismatch.
5. **Re-run** `validation/test_rf.py`; acceptance: σ_δ AG vs OCELOT < 2%.
6. **Record** in `validation/CHECKPOINTS.md`.

Verification already performed (audit, /tmp, read-only):
- AG `apply_rf_thin_lens(H=−9.78)` gives σ_δ = H·σ_z exactly (4.938 vs 4.943 e-3).
- With consistent ε_nz (σ_z@cavity = 303 µm): AG σ_δ = 2.953e-3 vs OCELOT
  2.921e-3 → **1.1 % agreement**.

---

## 5. Notes / open items

- **σ_z at sample remains model-limited**: with the thin-lens chirp,
  AG's σ_z at z = 0.777 → 477 µm (compression/overshoot, R56·h = −2.4 ⇒
  the bunch focuses then re-expands). OCELOT's σ_z stays 298 µm because
  `SecondTM` drift does not transport tau (known GPT limitation) — separate
  from the RF model, tracked in CHECKPOINTS.
- **Transit-time factor** TT ≈ 0.78 is a shared refinement option for both
  codes; do not apply to one code only.
- **Phase dependence**: the shared kick keeps `sin(φ + k·z)`; AG's old
  continuous coefficient had no φ at all. The thin-lens form restores the
  correct phase behaviour for free.

---

## Physics note (2026-08-01): shared RF physics switches

To let the two backends be compared with the RF transverse effect either on or
off, a shared `physics_switches` block was added to `shared/beamline_config.yaml`
(read identically by both adapters):

- `rf_longitudinal_kick` (default true): standardized thin-lens longitudinal
  kick `δ += K·sin(φ+kz)`, `H = eV·k·cosφ/(β²E₀)`. This is the compression
  physics and is always the same in both models.
- `rf_transverse_kick` (default false): Panofsky-Wenzel transverse RF kick
  `x' += K_trans·x`, `K_trans = −e·k·V/(2γβ·m_ec²)` (≈ −2.68 m⁻¹ here).
  When OFF, no RF transverse focusing/kick is applied by either model.

Default OFF keeps the baseline comparison clean (transverse-only RF force is a
second-order effect for the buncher at φ≈π); switching ON uses the SAME
equation in both models — no change to the formula, only to which models apply
it. The switch does not touch R56 or any longitudinal coordinate convention.
