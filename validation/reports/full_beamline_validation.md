# Full-Beamline Validation — AG vs OCELOT

Date: 2026-08-01 · shared config · RF transverse kick OFF (`physics_switches`)

## Sample-plane metrics (z = 777 mm)

| quantity | AG | OCELOT | dev [%] | verdict |
|---|---|---|---|---|
| sigma_x_um | 1984.191 | 1989.402 | 0.26 | PASS |
| sigma_y_um | 1984.191 | 1991.873 | 0.39 | PASS |
| sigma_delta_e3 | 2.953 | 2.924 | 0.98 | PASS |
| eps_nx_mm_mrad | 0.080 | 0.080 | 0.14 | PASS |
| eps_ny_mm_mrad | 0.080 | 0.080 | 0.09 | PASS |
| sigma_z_um | 477.001 | 471.497 | 1.17 | PASS |
| time_res_fs | 2902.307 | 2868.815 | 1.17 | PASS |

## Passed transverse / RF-kick metrics

- **sigma_x_um**: dev = 0.26% (threshold 5.0%) — PASS
- **sigma_y_um**: dev = 0.39% (threshold 5.0%) — PASS
- **sigma_delta_e3**: dev = 0.98% (threshold 2.0%) — PASS
- **eps_nx_mm_mrad**: dev = 0.14% (threshold 5.0%) — PASS
- **eps_ny_mm_mrad**: dev = 0.09% (threshold 5.0%) — PASS
- **sigma_z_um**: dev = 1.17% (threshold 5.0%) — PASS
- **time_res_fs**: dev = 1.17% (threshold 5.0%) — PASS

## Longitudinal (R56 adapter applied)

- σ_z at sample: AG 477.0 µm vs OCELOT 471.5 µm (dev 1.17%, threshold 5%).
- σ_t at sample: AG 2902.3 fs vs OCELOT 2868.8 fs (dev 1.17%).
- compression waist: AG z=546 mm vs OCELOT z=547 mm (Δz=0.6 mm).
- The input δ-variable convention mismatch (classification B) was resolved by the OCELOT longitudinal adapter (p_oc = β0·δ_p input, δ_p = p_oc/β0 output); see `CHANGELOG_r56_variable_adapter.md` and `R56_convention_resolution.md`.

## Integration status

- Both backends ran the complete lattice (cathode → sample, 7 elements) with identical shared parameters (config_sha dd8ada3d4cb2).
- Integration failures: **NONE**.
- Figure: `full_AG_vs_OCELOT.png` · results: `full_results.json`
