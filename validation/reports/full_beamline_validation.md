# Full-Beamline Validation — AG vs OCELOT

Date: 2026-08-01 · shared config · RF transverse kick OFF (`physics_switches`)

## Sample-plane metrics (z = 777 mm)

| quantity | AG | OCELOT | dev [%] | verdict |
|---|---|---|---|---|
| sigma_x_um | 1984.191 | 1996.205 | 0.60 | PASS |
| sigma_y_um | 1984.191 | 1991.415 | 0.36 | PASS |
| sigma_delta_e3 | 2.953 | 2.939 | 0.46 | PASS |
| eps_nx_mm_mrad | 0.080 | 0.080 | 0.53 | PASS |
| eps_ny_mm_mrad | 0.080 | 0.080 | 0.26 | PASS |
| sigma_z_um | 477.001 | 474.022 | 0.63 | PASS |
| time_res_fs | 2902.307 | 2884.180 | 0.63 | PASS |

## Passed transverse / RF-kick metrics

- **sigma_x_um**: dev = 0.60% (threshold 5.0%) — PASS
- **sigma_y_um**: dev = 0.36% (threshold 5.0%) — PASS
- **sigma_delta_e3**: dev = 0.46% (threshold 2.0%) — PASS
- **eps_nx_mm_mrad**: dev = 0.53% (threshold 5.0%) — PASS
- **eps_ny_mm_mrad**: dev = 0.26% (threshold 5.0%) — PASS
- **sigma_z_um**: dev = 0.63% (threshold 5.0%) — PASS
- **time_res_fs**: dev = 0.63% (threshold 5.0%) — PASS

## Longitudinal (R56 adapter applied)

- σ_z at sample: AG 477.0 µm vs OCELOT 474.0 µm (dev 0.63%, threshold 5%).
- σ_t at sample: AG 2902.3 fs vs OCELOT 2884.2 fs (dev 0.63%).
- compression waist: AG z=546 mm vs OCELOT z=546 mm (Δz=0.4 mm).
- The input δ-variable convention mismatch (classification B) was resolved by the OCELOT longitudinal adapter (p_oc = β0·δ_p input, δ_p = p_oc/β0 output); see `CHANGELOG_r56_variable_adapter.md` and `R56_convention_resolution.md`.

## Integration status

- Both backends ran the complete lattice (cathode → sample, 7 elements) with identical shared parameters (config_sha 23e97edb733b).
- Integration failures: **NONE**.
- Figure: `full_AG_vs_OCELOT.png` · results: `full_results.json`
