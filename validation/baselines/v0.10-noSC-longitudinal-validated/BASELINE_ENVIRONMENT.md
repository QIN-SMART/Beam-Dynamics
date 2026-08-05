# BASELINE ENVIRONMENT — v0.10 no-SC validated baseline

Freeze date: 2026-08-05

| item | value |
|---|---|
| OS | macOS 13.7.1 (x86_64) |
| Python | 3.12.2 (conda-forge) |
| NumPy | 1.26.4 |
| SciPy | 1.13.1 |
| Matplotlib | 3.9.2 |
| PyYAML | 6.0.1 |
| OCELOT | 26.06.1 (installed at /opt/anaconda3/lib/python3.12/site-packages/ocelot) |
| Git commit | 9657b525e2e078f18fbbc0fcbe26ef44f3a9aa51 (code baseline; tag commit adds the baseline folder, see BASELINE_MANIFEST.md) |
| Git branch | main |
| shared config SHA | dd8ada3d4cb2 |
| n_particles | 50000 |

Random-number policy: unchanged from the validated state — `generate_parray`
uses the unseeded global `np.random` for the first per-process beam (x/y/tau;
pre-existing), then `np.random.seed(42)` before overwriting px/py/δ_p
(deterministic).  Full per-process reproducibility of the FIRST beam remains
an open item.

Notes:
- The interpreter is /opt/anaconda3/bin/python3; the system python3
  (/Library/Frameworks/.../3.6) must NOT be used (missing yaml/numpy/ocelot).
- OCELOT source is not part of this repository; the exact version/commit is
  pinned above for reproducibility.
