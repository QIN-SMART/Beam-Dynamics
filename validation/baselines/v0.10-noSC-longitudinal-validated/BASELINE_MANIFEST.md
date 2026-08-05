# BASELINE MANIFEST — v0.10 no-SC longitudinal-validated

Freeze date: 2026-08-05

| item | value |
|---|---|
| Tag | v0.10-noSC-longitudinal-validated |
| Code commit | 9657b525e2e078f18fbbc0fcbe26ef44f3a9aa51 |
| Tag commit | (see `git rev-parse v0.10-noSC-longitudinal-validated`) |
| Branch | main |
| config SHA | dd8ada3d4cb2 |
| OCELOT | 26.06.1 |

## Files in this commit set
57 tracked files (first commit of the repository): AG/ (core + run_shared),
GPT模拟/ (route + benchmarks), shared/ (params/config/output), validation/
(backend, tests, reports), scripts/ (phone sync), root changelogs/docs,
.gitignore.  Excluded via .gitignore: __pycache__, .DS_Store, .vscode,
Okemos1.22b_imagSystem/, AG-prime/, notes-codes/, codes-result/,
code-update-log/, log/, 文献/ (PDFs), *.png (regenerable, except
validation/baselines/**), *_results.json and shared/results/*.json
(regenerable).

## Acceptance results (this freeze)
| test | result | key numbers |
|---|---|---|
| drift | PASS | σ_x 0.2%, σ_δ_p 0.04%, σ_z vs analytic 0.33% (OCELOT) |
| solenoid | PASS | coupling OFF vs OCELOT σ_x 0.40% |
| rf | PASS | σ_δ_p ~0.3-1%, kick residual 2.7e-14, routing {0,0,1,1} |
| full_beamline | PASS | σ_x 0.26%, σ_y 0.39%, σ_δ_p 0.98%, ε_nx 0.14%, ε_ny 0.09%, σ_z 1.17%, σ_t 1.17% |
| r56_convention | unchanged | raw slope = R56_tm = −1.163624, closure 0.2% |

## Key sample-plane values (full beamline, z = 777 mm)
| quantity | AG | OCELOT |
|---|---|---|
| σ_x [µm] | 1984.2 | 1989.4 |
| σ_y [µm] | 1984.2 | 1991.9 |
| σ_δ_p [e-3] | 2.953 | 2.924 |
| ε_nx [mm·mrad] | 0.0800 | 0.0800 |
| σ_z [µm] | 477.0 | 471.5 |
| σ_t [fs] | 2902.3 | 2868.8 |

Compression waist: AG z = 546 mm, σ_z = 14.8 µm; OCELOT z = 547 mm,
σ_z = 10.3 µm (Δz = 0.6 mm).

## Physics switches (frozen)
SC: OFF · RF longitudinal: ON · RF transverse: OFF · Solenoid coupling:
OFF (round beam) · R56: RESOLVED · OCELOT native R56: unmodified ·
δ_p = Δp/p0 · p_oc = ΔE/(c·p0) · conversion p_oc = β0·δ_p.

## Resolved items
- Solenoid reduced-order coupling spurious for round beams (exact 4×4 ⇒ σ_xy ≡ 0).
- OCELOT beam-energy unit (total GeV vs kinetic).
- Manual tau override removed (native OCELOT transport).
- RF thin-lens standardization (longitudinal + transverse, switch-gated).
- ε_nz made consistent (0.02 = βγ·σ_z·σ_δ).
- R56 input δ-variable convention mismatch (classification B) — adapter
  p_oc = β0·δ_p.
- Geometry single-source lattice refactor (multi-instance).

## Open items
- Space charge not yet validated (both models have implementations).
- RF transverse formula not independently validated; default OFF.
- AG emittance model cannot describe nonlinear growth — limitation.
- OCELOT full per-process RNG reproducibility not handled (first
  generate_parray per process unseeded).
- GPT main route (ued_beamline_v2.py) still builds its step lattices with
  hardcoded geometry rather than reading lattice.elements for steps 1–3
  (registered, not part of the framework tests).

## Rollback
```
git switch --detach v0.10-noSC-longitudinal-validated
# return to development:
git switch main
```
