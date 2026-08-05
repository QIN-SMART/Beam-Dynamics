# Longitudinal Transport Audit — AG vs OCELOT (z–δ)

Date: 2026-08-01 · Audit only, **no physics code modified**.
Goal: understand and standardize the z–δ (longitudinal) transport between the
AG envelope model and the OCELOT macroparticle model.

---

## 1. OCELOT longitudinal coordinates

`ParticleArray.rparticles` is 6×N:

| row | symbol | meaning |
|-----|--------|---------|
| 0,1 | x, px | transverse (canonical) |
| 2,3 | y, py | transverse (canonical) |
| 4 | **tau** | longitudinal coordinate, **τ = c·t [m]** (arrival-time × c) |
| 5 | **p / δ** | momentum deviation **δ = Δp/p₀** |

- `p.tau()` returns `rparticles[4]` in **metres** (τ = c·t), confirmed against
  OCELOT source ("tau = c*t") and by `sigma_tau` → physical σ_z = β·σ_tau.
- `p.p()` / `rparticles[5]` is the **momentum deviation δ = (p−p₀)/p₀**
  (OCELOT docs: "delta = (p-p0)/p0", not ΔE/E).

Reference energy convention (critical): OCELOT's beam energy is **total energy
in GeV** (γ = E_GeV/m_e_GeV). `generate_parray(energy=0.1)` ⇒ γ = 195.7,
β ≈ 1, i.e. **100 MeV**, not 100 keV.

---

## 2. Does the Drift apply  z_f = z_i + R56·δ ?

**Yes — natively, in tau coordinates.** The Drift (via `uni_matrix`, k1=hx=0)
has the (τ,δ) block

```
R = | 1   R56 |        R56 = −L/(β²·γ²)   < 0
    | 0   1   |
```

and SecondTM applies it: `τ_f = τ_i + R56·δ`.

Verified numerically at the correct 100 keV energy (γ=1.1957, β=0.5482,
L=0.5 m): R56 = −1.164; a single particle with δ = +1e-3 acquires
Δτ = −1163.6 µm = R56·δ exactly.

Sign: δ > 0 (faster particle) ⇒ τ decreases ⇒ earlier arrival ⇒ travels ahead
in space. Physically consistent.

---

## 3. Why "SecondTM does not change tau"

**It is NOT a SecondTM limitation.** SecondTM applies R56 correctly. The reason
tau appeared frozen in our runs is a **beam-energy unit bug** in the GPT route:

```
generate_parray(..., energy = E_keV/1000.0)      # E_keV=100 → 0.1
```
OCELOT reads 0.1 **GeV total energy** → γ = 195.7, β ≈ 1, so

```
R56 = −L/(β²γ²) ≈ −0.5/195.7² = −1.31e-5      # ~9×10⁴ too small
```
⇒ Δτ = R56·δ is numerically invisible (0.026 µm over 0.5 m in the test), so
σ_z stays "frozen".

With the **correct total energy** `E = (E_keV + 511)·1e-6 GeV = 6.11e-4 GeV`
(γ = 1.1957, β = 0.5482): R56 = −1.164 and a beam with σ_δ = 2e-3 acquires
σ_tau = 2327 µm over 0.5 m = |R56|·σ_δ exactly; σ_z(z) evolves and compresses
natively. **The longitudinal transport was never missing — it was using the
wrong energy.**

---

## 4. AG drift matrix comparison

AG uses the physical longitudinal offset z (head positive, in m) and
δ = Δp/p₀. Its drift map is

```
R_AG = | 1  R56_AG |        R56_AG = +L/γ²   > 0
       | 0  1      |
```
(dC_zδ/dz = σ_δ²/γ² in `envelope_ode` ⇒ R56 = L/γ²).

| quantity | AG (z) | OCELOT (τ=c·t) |
|----------|--------|----------------|
| coordinate | spatial offset z, head + | τ = c·t (arrival time) |
| R56 | +L/γ² | −L/(β²γ²) |
| fast particle (δ>0) | Δz > 0 (ahead) | Δτ < 0 (earlier) |

Coordinate relation:  z = β·c·t = β·τ, head-positive z ⇔ smaller τ (z =
−β(τ−τ_ref)). Mapping OCELOT to the spatial frame:

```
Δz = −β·Δτ = −β·R56_τ·δ = +L·δ/(β·γ²)      (OCELOT)
Δz = +L·δ/γ²                               (AG)
```

**Consistent in sign** (fast particles travel ahead in both), but the mapped
magnitudes differ by **1/β = 1.82** because OCELOT's R56_τ carries β², whereas
the self-consistent arrival-time R56 is −L/(βγ²). Flag for validation:
either OCELOT's `uni_matrix` R56 has a factor-1/β convention quirk, or a
different definition of τ is intended; the two codes are not 1:1 on the
absolute longitudinal kick yet.

---

## 5. Minimal solution (no hidden manual correction)

**Adopt Option A: use native OCELOT longitudinal transport.**

Single change in the beam-generation call of the GPT route
(`GPT模拟/ued_beamline_v2.py`, `validation/backend.run_ocelot`, and the
benchmark scripts that call `generate_parray`):

```
energy = (E_keV + 511.0) * 1e-6        # total energy in GeV   (was E_keV/1000)
sigma_tau = sigma_z / beta             # τ = c·t [m]           (already correct)
```

This activates the native `τ_f = τ_i + R56·δ` map:
- σ_z(z) evolves correctly in every drift (verified: R56 = −1.164 at 100 keV);
- the RF thin-lens chirp is then transported downstream by OCELOT itself,
  giving real bunch compression;
- **no manual `rparticles[4] += ...` override** anywhere (removes the hidden
  correction in `GPT模拟/bug/ued_beamline.py`).

Option B (fallback, documented): a **longitudinal map adapter** — apply the
AG-style drift map `z_f = z_i + R56·δ`, `R56 = +L/γ²` on the envelope/beam
after each drift step. This is explicit and documented (not hidden), but it
duplicates OCELOT's native physics and re-introduces the τ/z convention
subtleties. Only choose B if Option A proves insufficient.

**Recommended follow-up (separate task):** settle the factor-1/β between
OCELOT's `R56 = −L/(β²γ²)` and the exact arrival-time value −L/(βγ²), and
standardize the τ↔z mapping so AG and OCELOT agree on the absolute σ_z
evolution (the shared-config transverse/energy numbers already agree to
<1 %).

---

## References
- OCELOT `r_matrix.uni_matrix`: R56 = −z/(β²γ²) for a drift.
- OCELOT `SecondTM` / `SecondOrderMult.tmat_multip`: applies R (and T) to all 6
  coordinates including τ.
- OCELOT `generate_parray` / `m_e_GeV`: energy is total energy in GeV.
- AG `beam_dynamics_6d.envelope_ode`: `dC_zδ/dz = σ_δ²/γ²` ⇒ R56_z = L/γ².
- AG `beam_dynamics_6d.compute_R56(L,γ,β)`: R56 = L/(γ²β²) (drift, spatial).
