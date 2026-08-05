# BASELINE PHYSICS STATUS — v0.10 no-SC validated baseline

Freeze date: 2026-08-05 · config SHA: dd8ada3d4cb2

## Physics switches (shared config)
| switch | value |
|---|---|
| space_charge | **OFF** (not yet formally validated) |
| RF longitudinal kick | **ON** (standardized thin lens, H = eV·k·cosφ/(β²E₀)) |
| RF transverse kick | **OFF** (default; formula exists, independently unvalidated) |
| Solenoid reduced coupling | **OFF** (round-beam correction, adapter-level) |
| R56 convention | **RESOLVED** (classification B; OCELOT longitudinal adapter) |
| OCELOT native R56 | **unmodified** |
| shared longitudinal variable | δ_p = Δp/p0 |
| OCELOT native 6th coordinate | p_oc = ΔE/(c·p0) |
| conversion | p_oc = β0·δ_p (input), δ_p = p_oc/β0 (reported) |

## Beam parameters
| parameter | value |
|---|---|
| beam energy K0 | 100 keV |
| γ0 / β0 | 1.1957 / 0.5482 |
| initial σ_x = σ_y | 85 µm |
| initial σ_z | 300 µm |
| initial σ_δ_p | 1e-4 (momentum deviation) |
| ε_nx = ε_ny | 0.08 mm·mrad |
| ε_nz | 0.02 mm·mrad (consistent = βγ·σ_z·σ_δ) |
| charge | 100 fC |
| particles | 50000 |

## Beamline (lattice, 7 elements)
cathode(0) → drift(0.1) → solenoid1(0.06, B=0.05 T) → drift(0.24) →
rf1(0.022, 2.856 GHz, 30 kV, φ=π) → drift(0.355) → sample(0.777 m)

## RF parameters
V = 30 kV, f = 2.856 GHz, φ = π, L = 22 mm, k = 59.86 m⁻¹,
H = −9.78 m⁻¹, K_trans = −2.68 m⁻¹ (switch OFF).

## Solenoid parameters
B = 0.05 T, L = 60 mm, k_s = 22.38 m⁻¹ (k_s² = 500.65 m⁻²).

## Validation status (all four acceptance tests PASS)
- drift: PASS (σ_x 0.2 %, σ_δ_p 0.04 %, σ_z analytic 0.3–1.8 %)
- solenoid: PASS (coupling OFF vs OCELOT 0.40 %)
- rf: PASS (σ_δ_p 0.3–1 %, kick semantics 2.7e-14, routing {0,0,1,1})
- full beamline: PASS — 7 quantitative metrics (σ_x 0.26 %, σ_y 0.39 %,
  σ_δ_p 0.98 %, ε_nx 0.14 %, ε_ny 0.09 %, σ_z 1.17 %, σ_t 1.17 %;
  waist Δz = 0.6 mm)
- R56 characterization: unchanged (raw slope = transfer map = −1.163624;
  closure residual 0.2 %)
