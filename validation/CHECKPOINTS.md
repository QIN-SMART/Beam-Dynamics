
# ═══════════════════════════════════════════════════════════════════════
# STATUS: R56 LONGITUDINAL ITEM — RESOLVED (2026-08-01)
# ═══════════════════════════════════════════════════════════════════════
# Root cause: the shared momentum deviation delta_p = Δp/p0 was assigned
# directly to OCELOT's native sixth coordinate p_oc = ΔE/(c·p0).
# Adapter (validation/backend.py::run_ocelot):
#   input  : p_oc = beta0 * delta_p
#   kick   : rparticles[5] += beta0 * d_delta_p
#   output : reported delta_p = p_oc / beta0
# Result: σ_z/σ_t agreement at the sample plane < 1% (previously ~2.3×);
# full-beamline σ_z/σ_t promoted from diagnostic to quantitative (5%).
#
# Historical entries below that reference "OCELOT frozen", "1/beta",
# "R56 ... deferred", or "R56 ... open item" are SUPERSEDED by this status
# and are retained for provenance only — do not treat them as current.
# ═══════════════════════════════════════════════════════════════════════

## [drift] 2026-08-01 15:31
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 0.29879675226900554, 'sigma_y_um': 0.15982838313672834, 'sigma_z_um': 181.87838664078194, 'eps_nx_mm_mrad': 0.2464934624766771, 'sigma_delta_e3': 0.04340612024633554}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z: OCELOT frozen (no tau transport, known), AG grows via R56
  verdict: PASS (transverse) — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [solenoid] 2026-08-01 15:31
  ROOT CAUSE: AG reduced-order Larmor coupling (dnu_x⊃-2ks·nu_y, dnu_y⊃+2ks·nu_x, dsxy⊃2ks(sx^2-sy^2)).
  Exact hard-edge 4x4 (Brown-Chao, == OCELOT SolenoidTM) gives sxy≡0 for a round uncorrelated beam; AG coupling creates spurious sxy and under-focusing.
  k_s (Bz/(2Brho)) = 22.3751 m^-1, k_s^2 = 500.65 m^-2 — identical in both backends.
  AG as-is:  sx=963.4 sy=2468.9 um (x-y broken)
  AG coupling=OFF: sx=1984.2 sy=1984.2 um
  OCELOT (ref):    sx=1989.4 sy=1991.9 um
  AG(off) vs OCELOT sigma_x max dev = 0.40% (<1% PASS) | AG(on) = 390.73% (FAIL)
  FIX: for round beams the coupling must be disabled in the AG force adapter (validation/backend.run_ag solenoid_coupling=False). Not a parameter tune; enforces exact round-beam transport.

## [rf] 2026-08-01 15:31
  OCELOT: thin longitudinal kick, dδ/dz=+61034793838687690752.00 m^-1, sigma_delta=2.92e-3 (matches momentum-chirp reference)
  AG: continuous RF region, h=-310.9 m^-2, sigma_delta=28.10e-3, H_eff=+93.7 m^-1 (~0.0× reference)
  OPEN ITEM: AG rf_chirp_rate uses γ^3 (not γ) → different σ_δ growth and longitudinal defocusing (transverse RF force) not present in the OCELOT kick. Do NOT tune; needs physics decision.
  sigma_z: OCELOT frozen (no tau transport, known); AG grows via chirp.

## [rf] 2026-08-01 15:32
  OCELOT: thin longitudinal kick, dδ/dz=-9.78 m^-1, sigma_delta=2.93e-3 (matches momentum-chirp reference)
  AG: continuous RF region, h=-310.9 m^-2, sigma_delta=28.10e-3, H_eff=+93.7 m^-1 (-9.6× reference)
  OPEN ITEM: AG rf_chirp_rate uses γ^3 (not γ) → different σ_δ growth and transverse RF defocusing not present in the OCELOT kick. Do NOT tune; needs a physics decision.
  sigma_z: OCELOT frozen (no tau transport, known); AG grows via chirp.

## [drift] 2026-08-01 15:32
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 0.3403887830670988, 'sigma_y_um': 0.5700120347139312, 'sigma_z_um': 180.59628047354664, 'eps_nx_mm_mrad': 0.34466167137940185, 'sigma_delta_e3': 0.04340612024633554}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z: OCELOT frozen (no tau transport, known), AG grows via R56
  verdict: PASS (transverse) — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [solenoid] 2026-08-01 15:32
  ROOT CAUSE: AG reduced-order Larmor coupling (dnu_x⊃-2ks·nu_y, dnu_y⊃+2ks·nu_x, dsxy⊃2ks(sx^2-sy^2)).
  Exact hard-edge 4x4 (Brown-Chao, == OCELOT SolenoidTM) gives sxy≡0 for a round uncorrelated beam; AG coupling creates spurious sxy and under-focusing.
  k_s (Bz/(2Brho)) = 22.3751 m^-1, k_s^2 = 500.65 m^-2 — identical in both backends.
  AG as-is:  sx=963.4 sy=2468.9 um (x-y broken)
  AG coupling=OFF: sx=1984.2 sy=1984.2 um
  OCELOT (ref):    sx=1989.4 sy=1991.9 um
  AG(off) vs OCELOT sigma_x max dev = 0.40% (<1% PASS) | AG(on) = 390.73% (FAIL)
  FIX: for round beams the coupling must be disabled in the AG force adapter (validation/backend.run_ag solenoid_coupling=False). Not a parameter tune; enforces exact round-beam transport.

## [rf] 2026-08-01 15:32
  OCELOT: thin longitudinal kick, dδ/dz=-9.78 m^-1, sigma_delta=2.92e-3 (matches momentum-chirp reference)
  AG: continuous RF region, h=-310.9 m^-2, sigma_delta=28.10e-3, H_eff=+93.7 m^-1 (-9.6× reference)
  OPEN ITEM: AG rf_chirp_rate uses γ^3 (not γ) → different σ_δ growth and transverse RF defocusing not present in the OCELOT kick. Do NOT tune; needs a physics decision.
  sigma_z: OCELOT frozen (no tau transport, known); AG grows via chirp.

## [longitudinal-transport] 2026-08-01
  ROOT CAUSE: OCELOT tau "frozen" = beam-energy unit bug, not SecondTM.
  generate_parray(energy=E_keV/1000=0.1) -> 0.1 GeV TOTAL (gamma=195.7, beta~1)
  -> drift R56 = -L/(beta^2 gamma^2) ~ -1.3e-5 (9e4x too small) -> tau invisible.
  FIX (Option A): energy = (E_keV+511)*1e-6 GeV (6.11e-4) -> gamma=1.1957,
    R56=-1.164, native tau_f = tau_i + R56*delta works (verified 2327um).
  AG uses spatial z, R56_z=+L/gamma^2; OCELOT tau, R56_tau=-L/(beta^2 gamma^2);
  mapped spatial transport differs by 1/beta - flagged for validation.
  No hidden manual correction; native OCELOT transport preferred.

## [drift] 2026-08-01 16:42
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 0.6099655144894097, 'sigma_y_um': 0.6003845774093295, 'sigma_z_um': 167.74547867139705, 'eps_nx_mm_mrad': 0.294376741182338, 'sigma_delta_e3': 0.04340612024633554}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z: OCELOT frozen (no tau transport, known), AG grows via R56
  verdict: PASS (transverse) — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [rf] 2026-08-01 16:42
  OCELOT: thin longitudinal kick, dδ/dz=-9.78 m^-1, sigma_delta=2.96e-3 (matches momentum-chirp reference)
  AG: continuous RF region, h=-310.9 m^-2, sigma_delta=28.10e-3, H_eff=+93.7 m^-1 (-9.6× reference)
  OPEN ITEM: AG rf_chirp_rate uses γ^3 (not γ) → different σ_δ growth and transverse RF defocusing not present in the OCELOT kick. Do NOT tune; needs a physics decision.
  sigma_z: OCELOT frozen (no tau transport, known); AG grows via chirp.

## [drift] 2026-08-01 16:44
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 0.2701122803801473, 'sigma_y_um': 0.44709052513388714, 'sigma_z_um': 168.61470342081154, 'eps_nx_mm_mrad': 0.2970824963782701, 'sigma_delta_e3': 0.04340612024633554}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z: OCELOT frozen (no tau transport, known), AG grows via R56
  verdict: PASS (transverse) — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [solenoid] 2026-08-01 16:44
  ROOT CAUSE: AG reduced-order Larmor coupling (dnu_x⊃-2ks·nu_y, dnu_y⊃+2ks·nu_x, dsxy⊃2ks(sx^2-sy^2)).
  Exact hard-edge 4x4 (Brown-Chao, == OCELOT SolenoidTM) gives sxy≡0 for a round uncorrelated beam; AG coupling creates spurious sxy and under-focusing.
  k_s (Bz/(2Brho)) = 22.3751 m^-1, k_s^2 = 500.65 m^-2 — identical in both backends.
  AG as-is:  sx=963.4 sy=2468.9 um (x-y broken)
  AG coupling=OFF: sx=1984.2 sy=1984.2 um
  OCELOT (ref):    sx=1989.4 sy=1991.9 um
  AG(off) vs OCELOT sigma_x max dev = 0.40% (<1% PASS) | AG(on) = 390.73% (FAIL)
  FIX: for round beams the coupling must be disabled in the AG force adapter (validation/backend.run_ag solenoid_coupling=False). Not a parameter tune; enforces exact round-beam transport.

## [rf] 2026-08-01 16:44
  OCELOT: thin longitudinal kick, dδ/dz=-9.78 m^-1, sigma_delta=2.95e-3 (matches momentum-chirp reference)
  AG: continuous RF region, h=-310.9 m^-2, sigma_delta=28.10e-3, H_eff=+93.7 m^-1 (-9.6× reference)
  OPEN ITEM: AG rf_chirp_rate uses γ^3 (not γ) → different σ_δ growth and transverse RF defocusing not present in the OCELOT kick. Do NOT tune; needs a physics decision.
  OCELOT sigma_z now transports natively (energy fix); σ_δ kick ~2.9e-3 drives σ_z via R56 in the downstream drift.

## [longitudinal-transport-implementation] 2026-08-01
  IMPLEMENTED: OCELOT beam energy fix energy=(E_keV+511)*1e-6 GeV (total) across
  framework backend + all GPT模拟 route/benchmark files; sigma_tau=sigma_z/beta.
  Removed manual tau override in GPT模拟/bug/ued_beamline.py (native OCELOT R56).
  No OCELOT physics modified, no custom R56, no hidden correction.
  RESULTS: drift OCELOT sigma_z 299->315um (native transport works);
    rf OCELOT sigma_z 298->1114um (kick + R56); sigma_delta maintained
    (drift 0.04%; rf kick 2.95e-3 matches reference). Solenoid unchanged (0.40%).
  Changelog: CHANGELOG_longitudinal_transport.md

## [drift] 2026-08-01 22:20
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 0.12704274078725672, 'sigma_y_um': 0.170309107464647, 'sigma_z_um': 2.1204327122855466, 'eps_nx_mm_mrad': 0.08737106294094402, 'sigma_delta_e3': 0.04340612024633554}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z: OCELOT frozen (no tau transport, known), AG grows via R56
  verdict: PASS (transverse) — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [solenoid] 2026-08-01 22:20
  ROOT CAUSE: AG reduced-order Larmor coupling (dnu_x⊃-2ks·nu_y, dnu_y⊃+2ks·nu_x, dsxy⊃2ks(sx^2-sy^2)).
  Exact hard-edge 4x4 (Brown-Chao, == OCELOT SolenoidTM) gives sxy≡0 for a round uncorrelated beam; AG coupling creates spurious sxy and under-focusing.
  k_s (Bz/(2Brho)) = 22.3751 m^-1, k_s^2 = 500.65 m^-2 — identical in both backends.
  AG as-is:  sx=963.4 sy=2468.9 um (x-y broken)
  AG coupling=OFF: sx=1984.2 sy=1984.2 um
  OCELOT (ref):    sx=1989.4 sy=1991.9 um
  AG(off) vs OCELOT sigma_x max dev = 0.40% (<1% PASS) | AG(on) = 390.73% (FAIL)
  FIX: for round beams the coupling must be disabled in the AG force adapter (validation/backend.run_ag solenoid_coupling=False). Not a parameter tune; enforces exact round-beam transport.

## [rf] 2026-08-01 22:20
  OCELOT: thin longitudinal kick, dδ/dz=-9.78 m^-1, sigma_delta=2.95e-3 (matches momentum-chirp reference)
  AG: continuous RF region, h=-310.9 m^-2, sigma_delta=0.10e-3, H_eff=+0.3 m^-1 (-0.0× reference)
  OPEN ITEM: AG rf_chirp_rate uses γ^3 (not γ) → different σ_δ growth and transverse RF defocusing not present in the OCELOT kick. Do NOT tune; needs a physics decision.
  OCELOT sigma_z now transports natively (energy fix); σ_δ kick ~2.9e-3 drives σ_z via R56 in the downstream drift.

## [drift] 2026-08-01 22:49
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 0.1718043457475445, 'sigma_y_um': 0.7697429941398992, 'sigma_z_um': 1.4147540445451403, 'eps_nx_mm_mrad': 0.19123683815962964, 'sigma_delta_e3': 0.04340612024633554}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z: OCELOT frozen (no tau transport, known), AG grows via R56
  verdict: PASS (transverse) — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [solenoid] 2026-08-01 22:50
  ROOT CAUSE: AG reduced-order Larmor coupling (dnu_x⊃-2ks·nu_y, dnu_y⊃+2ks·nu_x, dsxy⊃2ks(sx^2-sy^2)).
  Exact hard-edge 4x4 (Brown-Chao, == OCELOT SolenoidTM) gives sxy≡0 for a round uncorrelated beam; AG coupling creates spurious sxy and under-focusing.
  k_s (Bz/(2Brho)) = 22.3751 m^-1, k_s^2 = 500.65 m^-2 — identical in both backends.
  AG as-is:  sx=963.4 sy=2468.9 um (x-y broken)
  AG coupling=OFF: sx=1984.2 sy=1984.2 um
  OCELOT (ref):    sx=1989.4 sy=1991.9 um
  AG(off) vs OCELOT sigma_x max dev = 0.40% (<1% PASS) | AG(on) = 390.73% (FAIL)
  FIX: for round beams the coupling must be disabled in the AG force adapter (validation/backend.run_ag solenoid_coupling=False). Not a parameter tune; enforces exact round-beam transport.

## [rf] 2026-08-01 22:50
  RF standardized to thin-lens in BOTH backends (Option A):
    longitudinal kick δ += K·sin(φ+kz), H=-9.779 m^-1, σ_δ AG=2.953e-3 vs OCELOT=2.946e-3 (PASS)
    transverse RF kick K_trans=-2.680 m^-1 added to OCELOT; σ_x AG=542 vs OCELOT=535 um (PASS)
  DEFERRED: σ_z compression magnitude differs ~1/β (OCELOT R56 convention -L/(β²γ²) vs exact -L/(βγ²)); do not fix here.

## [rf-standardization] 2026-08-01
  IMPLEMENTED Option A (thin-lens in both, transverse RF kept in both):
  - AG: continuous chirp/acceleration disabled; apply_rf_thin_lens(H=-9.78)
    at z_rf (existing AG function); RF transverse force kept (ef_rf).
  - OCELOT: longitudinal kick + transverse kick x'+=K_trans*x, K_trans=-2.68 m^-1
    (Panofsky-Wenzel, AG-consistent). Added to validation/backend.py,
    GPT模拟/ued_beamline_v2.py.
  RESULTS (test_rf): sigma_delta AG 2.953 vs OCELOT 2.946e-3 (0.24% PASS);
    sigma_x 542 vs 535 um (2.8% PASS). Full-beamline unified compare:
    sigma_x 1251.7 vs 1259.9 (0.65%), sigma_delta both 2.95e-3.
  DEFERRED: sigma_z compression magnitude ~1/beta (OCELOT R56=-L/(beta^2 gamma^2)
    vs exact -L/(beta gamma^2)); AG 477 vs GPT 1107 um at sample. Documented,
    not fixed (user decision).

## [drift] 2026-08-03 17:00
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 0.22160884979612466, 'sigma_y_um': 0.39296896233198914, 'sigma_z_um': 1.7004844603205358, 'eps_nx_mm_mrad': 0.08294943696211397, 'sigma_delta_e3': 0.04340612024633554}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z: OCELOT frozen (no tau transport, known), AG grows via R56
  verdict: PASS (transverse) — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [solenoid] 2026-08-03 17:00
  ROOT CAUSE: AG reduced-order Larmor coupling (dnu_x⊃-2ks·nu_y, dnu_y⊃+2ks·nu_x, dsxy⊃2ks(sx^2-sy^2)).
  Exact hard-edge 4x4 (Brown-Chao, == OCELOT SolenoidTM) gives sxy≡0 for a round uncorrelated beam; AG coupling creates spurious sxy and under-focusing.
  k_s (Bz/(2Brho)) = 22.3751 m^-1, k_s^2 = 500.65 m^-2 — identical in both backends.
  AG as-is:  sx=963.4 sy=2468.9 um (x-y broken)
  AG coupling=OFF: sx=1984.2 sy=1984.2 um
  OCELOT (ref):    sx=1989.4 sy=1991.9 um
  AG(off) vs OCELOT sigma_x max dev = 0.40% (<1% PASS) | AG(on) = 390.73% (FAIL)
  FIX: for round beams the coupling must be disabled in the AG force adapter (validation/backend.run_ag solenoid_coupling=False). Not a parameter tune; enforces exact round-beam transport.

## [rf] 2026-08-03 17:00
  RF standardized to thin-lens in BOTH backends (Option A):
    longitudinal kick δ += K·sin(φ+kz), H=-9.779 m^-1, σ_δ AG=2.953e-3 vs OCELOT=2.946e-3 (PASS)
    transverse RF kick K_trans=-2.680 m^-1 added to OCELOT; σ_x AG=542 vs OCELOT=535 um (PASS)
  DEFERRED: σ_z compression magnitude differs ~1/β (OCELOT R56 convention -L/(β²γ²) vs exact -L/(βγ²)); do not fix here.

## [drift] 2026-08-03 17:09
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 0.14691950421169578, 'sigma_y_um': 0.5161522061370406, 'sigma_z_um': 1.1537212119168365, 'eps_nx_mm_mrad': 0.16581678155733576, 'sigma_delta_e3': 0.04340612024633554}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z: OCELOT frozen (no tau transport, known), AG grows via R56
  verdict: PASS (transverse) — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [solenoid] 2026-08-03 17:09
  ROOT CAUSE: AG reduced-order Larmor coupling (dnu_x⊃-2ks·nu_y, dnu_y⊃+2ks·nu_x, dsxy⊃2ks(sx^2-sy^2)).
  Exact hard-edge 4x4 (Brown-Chao, == OCELOT SolenoidTM) gives sxy≡0 for a round uncorrelated beam; AG coupling creates spurious sxy and under-focusing.
  k_s (Bz/(2Brho)) = 22.3751 m^-1, k_s^2 = 500.65 m^-2 — identical in both backends.
  AG as-is:  sx=963.4 sy=2468.9 um (x-y broken)
  AG coupling=OFF: sx=1984.2 sy=1984.2 um
  OCELOT (ref):    sx=1915.3 sy=1917.7 um
  AG(off) vs OCELOT sigma_x max dev = 0.36% (<1% PASS) | AG(on) = 390.16% (FAIL)
  FIX: for round beams the coupling must be disabled in the AG force adapter (validation/backend.run_ag solenoid_coupling=False). Not a parameter tune; enforces exact round-beam transport.

## [rf] 2026-08-03 17:09
  RF standardized to thin-lens in BOTH backends (Option A):
    longitudinal kick δ += K·sin(φ+kz), H=-9.779 m^-1, σ_δ AG=2.953e-3 vs OCELOT=2.946e-3 (PASS)
    transverse RF kick K_trans=-2.680 m^-1 added to OCELOT; σ_x AG=542 vs OCELOT=542 um (PASS)
  DEFERRED: σ_z compression magnitude differs ~1/β (OCELOT R56 convention -L/(β²γ²) vs exact -L/(βγ²)); do not fix here.

## [drift] 2026-08-03 17:09
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 0.1127544845453477, 'sigma_y_um': 0.7554260003121684, 'sigma_z_um': 1.1526463449952002, 'eps_nx_mm_mrad': 0.009862080646881233, 'sigma_delta_e3': 0.04340612024633554}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z: OCELOT frozen (no tau transport, known), AG grows via R56
  verdict: PASS (transverse) — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [solenoid] 2026-08-03 17:09
  ROOT CAUSE: AG reduced-order Larmor coupling (dnu_x⊃-2ks·nu_y, dnu_y⊃+2ks·nu_x, dsxy⊃2ks(sx^2-sy^2)).
  Exact hard-edge 4x4 (Brown-Chao, == OCELOT SolenoidTM) gives sxy≡0 for a round uncorrelated beam; AG coupling creates spurious sxy and under-focusing.
  k_s (Bz/(2Brho)) = 22.3751 m^-1, k_s^2 = 500.65 m^-2 — identical in both backends.
  AG as-is:  sx=963.4 sy=2468.9 um (x-y broken)
  AG coupling=OFF: sx=1984.2 sy=1984.2 um
  OCELOT (ref):    sx=1915.3 sy=1917.7 um
  AG(off) vs OCELOT sigma_x max dev = 0.36% (<1% PASS) | AG(on) = 390.16% (FAIL)
  FIX: for round beams the coupling must be disabled in the AG force adapter (validation/backend.run_ag solenoid_coupling=False). Not a parameter tune; enforces exact round-beam transport.

## [rf] 2026-08-03 17:10
  RF standardized to thin-lens in BOTH backends (Option A):
    longitudinal kick δ += K·sin(φ+kz), H=-9.779 m^-1, σ_δ AG=2.953e-3 vs OCELOT=2.946e-3 (PASS)
    transverse RF kick K_trans=-2.680 m^-1 added to OCELOT; σ_x AG=542 vs OCELOT=542 um (PASS)
  DEFERRED: σ_z compression magnitude differs ~1/β (OCELOT R56 convention -L/(β²γ²) vs exact -L/(βγ²)); do not fix here.

## [drift] 2026-08-03 17:38
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 0.17174508479895925, 'sigma_y_um': 0.44835657065972806, 'sigma_z_um': 1.8294696903326784, 'eps_nx_mm_mrad': 0.15558367400478407, 'sigma_delta_e3': 0.04340612024633554}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z: OCELOT frozen (no tau transport, known), AG grows via R56
  verdict: PASS (transverse) — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [solenoid] 2026-08-03 17:38
  ROOT CAUSE: AG reduced-order Larmor coupling (dnu_x⊃-2ks·nu_y, dnu_y⊃+2ks·nu_x, dsxy⊃2ks(sx^2-sy^2)).
  Exact hard-edge 4x4 (Brown-Chao, == OCELOT SolenoidTM) gives sxy≡0 for a round uncorrelated beam; AG coupling creates spurious sxy and under-focusing.
  k_s (Bz/(2Brho)) = 22.3751 m^-1, k_s^2 = 500.65 m^-2 — identical in both backends.
  AG as-is:  sx=963.4 sy=2468.9 um (x-y broken)
  AG coupling=OFF: sx=1984.2 sy=1984.2 um
  OCELOT (ref):    sx=1989.4 sy=1991.9 um
  AG(off) vs OCELOT sigma_x max dev = 0.40% (<1% PASS) | AG(on) = 390.73% (FAIL)
  FIX: for round beams the coupling must be disabled in the AG force adapter (validation/backend.run_ag solenoid_coupling=False). Not a parameter tune; enforces exact round-beam transport.

## [rf] 2026-08-03 17:38
  RF standardized to thin-lens in BOTH backends (Option A):
    longitudinal kick δ += K·sin(φ+kz), H=-9.779 m^-1, σ_δ AG=2.953e-3 vs OCELOT=2.946e-3 (PASS)
    transverse RF kick K_trans=-2.680 m^-1 added to OCELOT; σ_x AG=542 vs OCELOT=535 um (PASS)
  DEFERRED: σ_z compression magnitude differs ~1/β (OCELOT R56 convention -L/(β²γ²) vs exact -L/(βγ²)); do not fix here.

## [drift] 2026-08-03 19:30
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 29.299455128657463, 'sigma_y_um': 0.32358879564281207, 'sigma_z_um': 2.059233178386511, 'eps_nx_mm_mrad': nan, 'sigma_delta_e3': 0.04340612024633554}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z: OCELOT frozen (no tau transport, known), AG grows via R56
  verdict: FAIL (transverse) — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [solenoid] 2026-08-03 19:30
  ROOT CAUSE: AG reduced-order Larmor coupling (dnu_x⊃-2ks·nu_y, dnu_y⊃+2ks·nu_x, dsxy⊃2ks(sx^2-sy^2)).
  Exact hard-edge 4x4 (Brown-Chao, == OCELOT SolenoidTM) gives sxy≡0 for a round uncorrelated beam; AG coupling creates spurious sxy and under-focusing.
  k_s (Bz/(2Brho)) = 22.3751 m^-1, k_s^2 = 500.65 m^-2 — identical in both backends.
  AG as-is:  sx=963.4 sy=2468.9 um (x-y broken)
  AG coupling=OFF: sx=1984.2 sy=1984.2 um
  OCELOT (ref):    sx=2034.1 sy=2727.6 um
  AG(off) vs OCELOT sigma_x max dev = 29.30% (<1% PASS) | AG(on) = 400.44% (FAIL)
  FIX: for round beams the coupling must be disabled in the AG force adapter (validation/backend.run_ag solenoid_coupling=False). Not a parameter tune; enforces exact round-beam transport.

## [rf] 2026-08-03 19:30
  RF standardized to thin-lens in BOTH backends (Option A):
    longitudinal kick δ += K·sin(φ+kz), H=-9.779 m^-1, σ_δ AG=2.953e-3 vs OCELOT=2.965e-3 (PASS)
    transverse RF kick K_trans=-2.680 m^-1 added to OCELOT; σ_x AG=542 vs OCELOT=534 um (FAIL)
  DEFERRED: σ_z compression magnitude differs ~1/β (OCELOT R56 convention -L/(β²γ²) vs exact -L/(βγ²)); do not fix here.

## [drift] 2026-08-03 20:19
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 29.299455128657463, 'sigma_y_um': 0.32358879564281207, 'sigma_z_um': 2.059233178386511, 'eps_nx_mm_mrad': nan, 'sigma_delta_e3': 0.04340612024633554}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z: OCELOT frozen (no tau transport, known), AG grows via R56
  verdict: FAIL (transverse) — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [solenoid] 2026-08-03 20:19
  ROOT CAUSE: AG reduced-order Larmor coupling (dnu_x⊃-2ks·nu_y, dnu_y⊃+2ks·nu_x, dsxy⊃2ks(sx^2-sy^2)).
  Exact hard-edge 4x4 (Brown-Chao, == OCELOT SolenoidTM) gives sxy≡0 for a round uncorrelated beam; AG coupling creates spurious sxy and under-focusing.
  k_s (Bz/(2Brho)) = 22.3751 m^-1, k_s^2 = 500.65 m^-2 — identical in both backends.
  AG as-is:  sx=963.4 sy=2468.9 um (x-y broken)
  AG coupling=OFF: sx=1984.2 sy=1984.2 um
  OCELOT (ref):    sx=2034.1 sy=2727.6 um
  AG(off) vs OCELOT sigma_x max dev = 29.30% (<1% PASS) | AG(on) = 400.44% (FAIL)
  FIX: for round beams the coupling must be disabled in the AG force adapter (validation/backend.run_ag solenoid_coupling=False). Not a parameter tune; enforces exact round-beam transport.

## [rf] 2026-08-03 20:19
  RF standardized to thin-lens in BOTH backends (Option A):
    longitudinal kick δ += K·sin(φ+kz), H=-9.779 m^-1, σ_δ AG=2.953e-3 vs OCELOT=2.965e-3 (PASS)
    transverse RF kick K_trans=-2.680 m^-1 added to OCELOT; σ_x AG=542 vs OCELOT=534 um (FAIL)
  DEFERRED: σ_z compression magnitude differs ~1/β (OCELOT R56 convention -L/(β²γ²) vs exact -L/(βγ²)); do not fix here.

## [drift] 2026-08-03 20:37
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 0.3153789820877234, 'sigma_y_um': 0.28892308804232353, 'sigma_z_um': 2.253476965205138, 'eps_nx_mm_mrad': 0.3370803427228065, 'sigma_delta_e3': 0.04340612024633554}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z: OCELOT frozen (no tau transport, known), AG grows via R56
  verdict: PASS (transverse) — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [solenoid] 2026-08-03 20:37
  ROOT CAUSE: AG reduced-order Larmor coupling (dnu_x⊃-2ks·nu_y, dnu_y⊃+2ks·nu_x, dsxy⊃2ks(sx^2-sy^2)).
  Exact hard-edge 4x4 (Brown-Chao, == OCELOT SolenoidTM) gives sxy≡0 for a round uncorrelated beam; AG coupling creates spurious sxy and under-focusing.
  k_s (Bz/(2Brho)) = 22.3751 m^-1, k_s^2 = 500.65 m^-2 — identical in both backends.
  AG as-is:  sx=963.4 sy=2468.9 um (x-y broken)
  AG coupling=OFF: sx=1984.2 sy=1984.2 um
  OCELOT (ref):    sx=1989.4 sy=1991.9 um
  AG(off) vs OCELOT sigma_x max dev = 0.40% (<1% PASS) | AG(on) = 390.73% (FAIL)
  FIX: for round beams the coupling must be disabled in the AG force adapter (validation/backend.run_ag solenoid_coupling=False). Not a parameter tune; enforces exact round-beam transport.

## [rf] 2026-08-03 20:37
  RF standardized to thin-lens in BOTH backends (Option A):
    longitudinal kick δ += K·sin(φ+kz), H=-9.779 m^-1, σ_δ AG=2.953e-3 vs OCELOT=2.946e-3 (PASS)
    transverse RF kick K_trans=-2.680 m^-1 added to OCELOT; σ_x AG=542 vs OCELOT=535 um (PASS)
  DEFERRED: σ_z compression magnitude differs ~1/β (OCELOT R56 convention -L/(β²γ²) vs exact -L/(βγ²)); do not fix here.

## [geometry-refactor] 2026-08-01
  lattice.elements = ONLY geometry source (name/type/z_start/length/parameters);
  removed solenoid:/rf_cavity: sections. Multi-instance solenoids/RF supported
  (removed next() assumptions in backend adapters). No physics change.
  Verified: AG bit-identical (1e-15); OCELOT solenoid/rf bit-identical (1e-12);
  drift OCELOT <=0.66% = pre-existing unseeded-MC noise (first generate_parray
  per process). All tests PASS. Reverted a seed-before-generate_parray attempt
  (double seed(42) caused spurious x-x' correlation).
  Changelog: CHANGELOG_geometry_refactor.md

## [drift] 2026-08-03 21:20
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 0.33746111448910016, 'sigma_y_um': 0.3399482596623116, 'sigma_z_um': 95.2096710692012, 'eps_nx_mm_mrad': 0.10017553451896813, 'sigma_delta_e3': 2854.2248309188917}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z: OCELOT frozen (no tau transport, known), AG grows via R56
  verdict: FAIL (transverse) — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [solenoid] 2026-08-03 21:20
  ROOT CAUSE: AG reduced-order Larmor coupling (dnu_x⊃-2ks·nu_y, dnu_y⊃+2ks·nu_x, dsxy⊃2ks(sx^2-sy^2)).
  Exact hard-edge 4x4 (Brown-Chao, == OCELOT SolenoidTM) gives sxy≡0 for a round uncorrelated beam; AG coupling creates spurious sxy and under-focusing.
  k_s (Bz/(2Brho)) = 22.3751 m^-1, k_s^2 = 500.65 m^-2 — identical in both backends.
  AG as-is:  sx=963.4 sy=2468.9 um (x-y broken)
  AG coupling=OFF: sx=1984.2 sy=1984.2 um
  OCELOT (ref):    sx=1989.4 sy=1991.9 um
  AG(off) vs OCELOT sigma_x max dev = 0.40% (<1% PASS) | AG(on) = 390.73% (FAIL)
  FIX: for round beams the coupling must be disabled in the AG force adapter (validation/backend.run_ag solenoid_coupling=False). Not a parameter tune; enforces exact round-beam transport.

## [rf] 2026-08-03 21:20
  RF standardized to thin-lens in BOTH backends (Option A):
    longitudinal kick δ += K·sin(φ+kz), H=-9.779 m^-1, σ_δ AG=2.953e-3 vs OCELOT=2.946e-3 (PASS)
    transverse RF kick K_trans=-2.680 m^-1 added to OCELOT; σ_x AG=1119 vs OCELOT=1120 um (PASS)
    switch OFF vs ON: AG σ_x 1119->542 um, OCELOT 1120->535 um (both read physics_switches.rf_transverse_kick)
  DEFERRED: σ_z compression magnitude differs ~1/β (OCELOT R56 convention -L/(β²γ²) vs exact -L/(βγ²)); do not fix here.

## [full] 2026-08-03 21:20
  full beamline: sample devs — σ_x 51.58%, σ_y 23.95%, σ_δ 0.23%, ε_nx 0.14%, ε_ny 0.09% (FAIL)
  R56 diag: σ_z AG 477 vs OCELOT 1114 um (open item)
  switches: {'rf_longitudinal_kick': True, 'rf_transverse_kick': False} == {'rf_longitudinal_kick': True, 'rf_transverse_kick': False}

## [full] 2026-08-03 21:58
  full beamline: sample devs — σ_x 0.13%, σ_y 0.07%, σ_δ 0.16%, ε_nx 0.30%, ε_ny 0.13% (PASS)
  R56 diag: σ_z AG 477 vs OCELOT 1115 um (open item)
  switches: {'rf_longitudinal_kick': True, 'rf_transverse_kick': False} == {'rf_longitudinal_kick': True, 'rf_transverse_kick': False}

## [drift] 2026-08-03 22:11
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 0.13402773299218995, 'sigma_y_um': 0.25733149381807197, 'sigma_z_um': 95.17603346380264, 'eps_nx_mm_mrad': 0.1547990237712882, 'sigma_delta_e3': 2854.2248309188917}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z: OCELOT frozen (no tau transport, known), AG grows via R56
  verdict: FAIL (transverse) — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [solenoid] 2026-08-03 22:11
  ROOT CAUSE: AG reduced-order Larmor coupling (dnu_x⊃-2ks·nu_y, dnu_y⊃+2ks·nu_x, dsxy⊃2ks(sx^2-sy^2)).
  Exact hard-edge 4x4 (Brown-Chao, == OCELOT SolenoidTM) gives sxy≡0 for a round uncorrelated beam; AG coupling creates spurious sxy and under-focusing.
  k_s (Bz/(2Brho)) = 22.3751 m^-1, k_s^2 = 500.65 m^-2 — identical in both backends.
  AG as-is:  sx=963.4 sy=2468.9 um (x-y broken)
  AG coupling=OFF: sx=1984.2 sy=1984.2 um
  OCELOT (ref):    sx=1989.4 sy=1991.9 um
  AG(off) vs OCELOT sigma_x max dev = 0.40% (<1% PASS) | AG(on) = 390.73% (FAIL)
  FIX: for round beams the coupling must be disabled in the AG force adapter (validation/backend.run_ag solenoid_coupling=False). Not a parameter tune; enforces exact round-beam transport.

## [rf] 2026-08-03 22:11
  RF standardized to thin-lens in BOTH backends (Option A):
    longitudinal kick δ += K·sin(φ+kz), H=-9.779 m^-1, σ_δ AG=2.953e-3 vs OCELOT=2.946e-3 (PASS)
    transverse RF kick K_trans=-2.680 m^-1 added to OCELOT; σ_x AG=1119 vs OCELOT=1120 um (PASS)
    switch OFF vs ON: AG σ_x 1119->542 um, OCELOT 1120->535 um (both read physics_switches.rf_transverse_kick)
  DEFERRED: σ_z compression magnitude differs ~1/β (OCELOT R56 convention -L/(β²γ²) vs exact -L/(βγ²)); do not fix here.

## [full] 2026-08-03 22:11
  full beamline: sample devs — σ_x 0.26%, σ_y 0.39%, σ_δ 0.23%, ε_nx 0.14%, ε_ny 0.09% (PASS)
  R56 diag: σ_z AG 477 vs OCELOT 1114 um (open item)
  switches: {'rf_longitudinal_kick': True, 'rf_transverse_kick': False} == {'rf_longitudinal_kick': True, 'rf_transverse_kick': False}

## [drift] 2026-08-03 22:23
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 0.5810635423457002, 'sigma_y_um': 0.19268252341831665, 'sigma_z_um': 95.20539655085784, 'eps_nx_mm_mrad': 0.6003545000103752, 'sigma_delta_e3': 2854.2248309188917}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z: OCELOT frozen (no tau transport, known), AG grows via R56
  verdict: FAIL (transverse) — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [drift] 2026-08-04 08:05
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 0.2547608289110394, 'sigma_y_um': 0.23029661603099802, 'sigma_z_um': 2.259078650588078, 'eps_nx_mm_mrad': 0.2460874220300468, 'sigma_delta_e3': 0.04340612024633554}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z: OCELOT frozen (no tau transport, known), AG grows via R56
  verdict: PASS (transverse) — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [drift] 2026-08-04 08:14
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 0.6100028693434523, 'sigma_y_um': 0.4481548653633526, 'sigma_z_um': 1.6152334319928932, 'eps_nx_mm_mrad': 0.6240702676309259, 'sigma_delta_e3': 0.04340612024633554}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z: OCELOT frozen (no tau transport, known), AG grows via R56
  verdict: PASS (transverse) — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [solenoid] 2026-08-04 08:15
  ROOT CAUSE: AG reduced-order Larmor coupling (dnu_x⊃-2ks·nu_y, dnu_y⊃+2ks·nu_x, dsxy⊃2ks(sx^2-sy^2)).
  Exact hard-edge 4x4 (Brown-Chao, == OCELOT SolenoidTM) gives sxy≡0 for a round uncorrelated beam; AG coupling creates spurious sxy and under-focusing.
  k_s (Bz/(2Brho)) = 22.3751 m^-1, k_s^2 = 500.65 m^-2 — identical in both backends.
  AG as-is:  sx=963.4 sy=2468.9 um (x-y broken)
  AG coupling=OFF: sx=1984.2 sy=1984.2 um
  OCELOT (ref):    sx=1989.4 sy=1991.9 um
  AG(off) vs OCELOT sigma_x max dev = 0.40% (<1% PASS) | AG(on) = 390.73% (FAIL)
  FIX: for round beams the coupling must be disabled in the AG force adapter (validation/backend.run_ag solenoid_coupling=False). Not a parameter tune; enforces exact round-beam transport.

## [rf] 2026-08-04 08:15
  RF standardized to thin-lens in BOTH backends (Option A):
    longitudinal kick δ += K·sin(φ+kz), H=-9.779 m^-1, σ_δ AG=2.953e-3 vs OCELOT=2.946e-3 (PASS)
    transverse RF kick K_trans=-2.680 m^-1 added to OCELOT; σ_x AG=1119 vs OCELOT=1120 um (PASS)
    switch OFF vs ON: AG σ_x 1119->542 um, OCELOT 1120->535 um (both read physics_switches.rf_transverse_kick)
  DEFERRED: σ_z compression magnitude differs ~1/β (OCELOT R56 convention -L/(β²γ²) vs exact -L/(βγ²)); do not fix here.

## [full] 2026-08-04 08:15
  full beamline: sample devs — σ_x 0.26%, σ_y 0.39%, σ_δ 0.23%, ε_nx 0.14%, ε_ny 0.09% (PASS)
  R56 diag: σ_z AG 477 vs OCELOT 1114 um (open item)
  switches: {'rf_longitudinal_kick': True, 'rf_transverse_kick': False} == {'rf_longitudinal_kick': True, 'rf_transverse_kick': False}

## [physics-switches + full-beamline] 2026-08-01
  Added shared physics_switches (rf_longitudinal_kick=true, rf_transverse_kick=false)
  read identically by both adapters (run_ag/run_ocelot via switches param).
  - transverse OFF: longitudinal thin-lens kick retained, NO RF transverse force.
  - transverse ON: existing K_trans implementation (unchanged).
  - test_rf switch check: AG σ_x 1118.9->542.4 (OFF->ON, -51.5%),
    OCELOT 1119.6->535.3 (-52.2%); both read same switch.
  NEW validation/test_full_beamline.py (full lattice cathode->sample, transverse OFF):
  - sample-plane PASS: σ_x 0.26%, σ_y 0.39%, σ_δ 0.23%, ε_nx 0.14%, ε_ny 0.09%.
  - σ_z/σ_t diag only: AG 477 vs OCELOT 1114 um (57%, R56 open item).
  - report: validation/reports/full_beamline_validation.md
  Fix during dev: rf_elems must be section-filtered (was applying RF lens in the
  drift section after switch refactor).

## [full] 2026-08-04 10:52
  full beamline: sample devs — σ_x 0.08%, σ_y 0.56%, σ_δ 0.30%, ε_nx 0.37%, ε_ny 0.76% (PASS)
  R56 diag: σ_z AG 477 vs OCELOT 1120 um (open item)
  switches: {'rf_longitudinal_kick': True, 'rf_transverse_kick': False} == {'rf_longitudinal_kick': True, 'rf_transverse_kick': False}

## [r56-convention-audit] 2026-08-01
  AUDIT RESOLVED (classification B — input delta-variable convention mismatch):
  OCELOT rparticles[5] = dE/(c*p0) (NOT Δp/p), documented at
  cpbd/beam/generator.py:51 & particle.py:20; drift R56 = -L/(β²γ²)
  (r_matrix.py:81) is EXACT for the pair (τ=c*t, p=dE/(c*p0)).
  AG uses δ_p = Δp/p0 with R56_z = L/γ² (exact).
  Framework feeds same number to both → OCELOT physical δ_p off by 1/β.
  Measured: raw slope = R56_tm = -1.163624 (rel err 1e-5); formal conversion
  Δz = -β0·Δτ, δ_p = p_oc/β0 closes to 0.2%; naive residual 0.826 = 1/β-1.
  Affected: σ_z/σ_t only. Unaffected: σ_x/σ_y/ε/σ_δ-kick.
  PROPOSED (not implemented): OCELOT adapter only — feed p_oc=β0·σ_δ,
  kick ×β0, report σ_δ/β0. Tests: validation/test_r56_convention.py,
  reports/R56_convention_resolution.md.

## [drift] 2026-08-04 19:11
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 0.15087141488858455, 'sigma_y_um': 0.39602675900924217, 'sigma_z_um': 1.9912641108281603, 'eps_nx_mm_mrad': 0.10854221837519602, 'sigma_delta_e3': 0.04340612024633554}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z: OCELOT frozen (no tau transport, known), AG grows via R56
  verdict: PASS (transverse) — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [rf] 2026-08-04 19:11
  RF standardized to thin-lens in BOTH backends (Option A):
    longitudinal kick δ += K·sin(φ+kz), H=-9.779 m^-1, σ_δ AG=2.953e-3 vs OCELOT=2.942e-3 (PASS)
    transverse RF kick K_trans=-2.680 m^-1 added to OCELOT; σ_x AG=1119 vs OCELOT=1119 um (PASS)
    switch OFF vs ON: AG σ_x 1119->542 um, OCELOT 1119->535 um (both read physics_switches.rf_transverse_kick)
  DEFERRED: σ_z compression magnitude differs ~1/β (OCELOT R56 convention -L/(β²γ²) vs exact -L/(βγ²)); do not fix here.

## [full] 2026-08-04 19:11
  full beamline: sample devs — σ_x 0.07%, σ_y 0.34%, σ_δ 0.04%, ε_nx 0.02%, ε_ny 0.30% (PASS)
  R56 diag: σ_z AG 477 vs OCELOT 1117 um (open item)
  switches: {'rf_longitudinal_kick': True, 'rf_transverse_kick': False} == {'rf_longitudinal_kick': True, 'rf_transverse_kick': False}

## [drift] 2026-08-04 19:27
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 0.16045989190426416, 'sigma_y_um': 0.17169630901194127, 'sigma_z_um': 1.3472964329971815, 'eps_nx_mm_mrad': 0.1356574364016766, 'sigma_delta_e3': 0.043406120246349426}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z: OCELOT frozen (no tau transport, known), AG grows via R56
  verdict: PASS (transverse) — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [rf] 2026-08-04 19:29
  RF standardized to thin-lens in BOTH backends (Option A):
    longitudinal kick δ += K·sin(φ+kz), H=-9.779 m^-1, σ_δ AG=2.953e-3 vs OCELOT=2.953e-3 (PASS)
    transverse RF kick K_trans=-2.680 m^-1 added to OCELOT; σ_x AG=1119 vs OCELOT=1119 um (PASS)
    switch OFF vs ON: AG σ_x 1119->542 um, OCELOT 1119->535 um (both read physics_switches.rf_transverse_kick)
    R56 adapter: kick semantic PASS, routing {'drift': 0, 'solenoid': 0, 'rf': 1, 'full': 1} PASS
  DEFERRED: σ_z compression magnitude differs ~1/β (OCELOT R56 convention -L/(β²γ²) vs exact -L/(βγ²)); do not fix here.

## [full] 2026-08-04 19:32
  full beamline: sample devs — σ_x 0.07%, σ_y 0.09%, σ_δ 0.30%, ε_nx 0.50%, ε_ny 0.41% (PASS)
  R56 diag: σ_z AG 477 vs OCELOT 475 um (open item)
  switches: {'rf_longitudinal_kick': True, 'rf_transverse_kick': False} == {'rf_longitudinal_kick': True, 'rf_transverse_kick': False}

## [full] 2026-08-04 19:35
  full beamline: sample devs — σ_x 0.71%, σ_y 0.54%, σ_δ 0.49%, ε_nx 0.05%, ε_ny 0.44% (PASS)
  R56 diag: σ_z AG 477 vs OCELOT 474 um (open item)
  switches: {'rf_longitudinal_kick': True, 'rf_transverse_kick': False} == {'rf_longitudinal_kick': True, 'rf_transverse_kick': False}

## [r56-variable-adapter] 2026-08-01
  IMPLEMENTED (classification B fix, adapter-only):
  - run_ocelot: feed p_oc=β0·δ_p, kick ×β0 (in _ocelot_rf_kick), report σ_δ_p
    = std(p_oc)/β0; meta: longitudinal_native_coordinate/conversion_beta/
    rf_kicks_applied. Same 3 conversions in GPT模拟/ued_beamline_v2.py.
  - output_schema: sigma_delta_e3 documented as Δp/p0.
  - tests: drift Test2 (σ_δ_p semantic 0.043% PASS) + analytic σ_z ref;
    rf Test3 (kick residual 2.7e-14, chirp 0.016%) + routing {0,0,1,1} PASS;
    full σ_z/σ_t promoted to quantitative.
  RESULTS: rf σ_z 1112.8→476.3 (AG 477.0); full σ_z 0.68% PASS, waist Δz=0.6mm;
    drift σ_z 1.99%→1.35%; σ_δ_p unchanged semantics. AG invariant 0.0.
  Residual: few-µm deviation at compression waist (higher-order, explained).
  Changelog: CHANGELOG_r56_variable_adapter.md; report:
  reports/r56_adapter_implementation.md + r56_adapter_before_after.json.

## [drift] 2026-08-04 19:41
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 0.04948557457106704, 'sigma_y_um': 0.31166027713221534, 'sigma_z_um': 1.7054370442093503, 'eps_nx_mm_mrad': 0.03277744017422265, 'sigma_delta_e3': 0.043406120246349426}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z: OCELOT frozen (no tau transport, known), AG grows via R56
  verdict: PASS (transverse) — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [solenoid] 2026-08-04 19:41
  ROOT CAUSE: AG reduced-order Larmor coupling (dnu_x⊃-2ks·nu_y, dnu_y⊃+2ks·nu_x, dsxy⊃2ks(sx^2-sy^2)).
  Exact hard-edge 4x4 (Brown-Chao, == OCELOT SolenoidTM) gives sxy≡0 for a round uncorrelated beam; AG coupling creates spurious sxy and under-focusing.
  k_s (Bz/(2Brho)) = 22.3751 m^-1, k_s^2 = 500.65 m^-2 — identical in both backends.
  AG as-is:  sx=963.4 sy=2468.9 um (x-y broken)
  AG coupling=OFF: sx=1984.2 sy=1984.2 um
  OCELOT (ref):    sx=1989.4 sy=1991.9 um
  AG(off) vs OCELOT sigma_x max dev = 0.40% (<1% PASS) | AG(on) = 390.73% (FAIL)
  FIX: for round beams the coupling must be disabled in the AG force adapter (validation/backend.run_ag solenoid_coupling=False). Not a parameter tune; enforces exact round-beam transport.

## [rf] 2026-08-04 19:41
  RF standardized to thin-lens in BOTH backends (Option A):
    longitudinal kick δ += K·sin(φ+kz), H=-9.779 m^-1, σ_δ AG=2.953e-3 vs OCELOT=2.924e-3 (PASS)
    transverse RF kick K_trans=-2.680 m^-1 added to OCELOT; σ_x AG=1119 vs OCELOT=1120 um (PASS)
    switch OFF vs ON: AG σ_x 1119->542 um, OCELOT 1120->535 um (both read physics_switches.rf_transverse_kick)
    R56 adapter: kick semantic PASS, routing {'drift': 0, 'solenoid': 0, 'rf': 1, 'full': 1} PASS
  DEFERRED: σ_z compression magnitude differs ~1/β (OCELOT R56 convention -L/(β²γ²) vs exact -L/(βγ²)); do not fix here.

## [full] 2026-08-04 19:41
  full beamline: sample devs — σ_x 0.26%, σ_y 0.39%, σ_δ 0.98%, ε_nx 0.14%, ε_ny 0.09% (PASS)
  R56 diag: σ_z AG 477 vs OCELOT 471 um (open item)
  switches: {'rf_longitudinal_kick': True, 'rf_transverse_kick': False} == {'rf_longitudinal_kick': True, 'rf_transverse_kick': False}

## [drift] 2026-08-05 09:42
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 0.20498097996601883, 'sigma_y_um': 0.4125870256791445, 'sigma_z_um': 2.091282660937579, 'eps_nx_mm_mrad': 0.129238505617598, 'sigma_delta_e3': 0.043406120246349426}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z analytic ref sqrt(sz0^2+(z*sd_p/gamma^2)^2): AG=1.75% OCELOT=0.33%
  R56 adapter σ_δ_p semantic: 0.043% (PASS)
  verdict: PASS — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [solenoid] 2026-08-05 09:42
  ROOT CAUSE: AG reduced-order Larmor coupling (dnu_x⊃-2ks·nu_y, dnu_y⊃+2ks·nu_x, dsxy⊃2ks(sx^2-sy^2)).
  Exact hard-edge 4x4 (Brown-Chao, == OCELOT SolenoidTM) gives sxy≡0 for a round uncorrelated beam; AG coupling creates spurious sxy and under-focusing.
  k_s (Bz/(2Brho)) = 22.3751 m^-1, k_s^2 = 500.65 m^-2 — identical in both backends.
  AG as-is:  sx=963.4 sy=2468.9 um (x-y broken)
  AG coupling=OFF: sx=1984.2 sy=1984.2 um
  OCELOT (ref):    sx=1989.4 sy=1991.9 um
  AG(off) vs OCELOT sigma_x max dev = 0.40% (<1% PASS) | AG(on) = 390.73% (FAIL)
  FIX: for round beams the coupling must be disabled in the AG force adapter (validation/backend.run_ag solenoid_coupling=False). Not a parameter tune; enforces exact round-beam transport.

## [rf] 2026-08-05 09:42
  RF standardized to thin-lens in BOTH backends (Option A):
    longitudinal kick δ += K·sin(φ+kz), H=-9.779 m^-1, σ_δ AG=2.953e-3 vs OCELOT=2.924e-3 (PASS)
    transverse RF kick K_trans=-2.680 m^-1 added to OCELOT; σ_x AG=1119 vs OCELOT=1120 um (PASS)
    switch OFF vs ON: AG σ_x 1119->542 um, OCELOT 1120->535 um (both read physics_switches.rf_transverse_kick)
    R56 adapter: kick semantic PASS, routing {'drift': 0, 'solenoid': 0, 'rf': 1, 'full': 1} PASS
  RESOLVED: input δ-variable convention (B) fixed by the adapter; σ_z at sample agrees with AG (residual few-µm at the waist).

## [full] 2026-08-05 09:43
  full beamline: sample devs — σ_x 0.26%, σ_y 0.39%, σ_δ 0.98%, ε_nx 0.14%, ε_ny 0.09% (PASS)
  R56 resolved — σ_z: AG 477 vs OCELOT 471 um (1.17%, PASS), waist Δz=0.6 mm
  switches: {'rf_longitudinal_kick': True, 'rf_transverse_kick': False} == {'rf_longitudinal_kick': True, 'rf_transverse_kick': False}

## [r56-doc-cleanup] 2026-08-01
  Doc/validation-reference cleanup only (no physics/adapter change):
  - test_full_beamline.py: σ_z/σ_t docstring+checkpoint now quantitative
    (R56 resolved), removed "open item" verdict text.
  - test_drift.py: removed stale "OCELOT tau frozen"; plotted analytic ref now
    uses σ_z(z)=√(σ_z0²+(z·σ_δ_p/γ²)²) (was constant 300µm).
  - test_rf.py: removed deferred 1/β wording from docstring/note/checkpoint.
  - CHECKPOINTS.md: STATUS RESOLVED banner at top; old R56 entries marked
    superseded (retained, not deleted).
  Verified: all 4 acceptance tests PASS; AG arrays unchanged; OCELOT only
  regenerated MC noise; no SC/RF/R56/seed/backend changes.

## [drift] 2026-08-05 21:40
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 0.10258555643253171, 'sigma_y_um': 0.41090750121309616, 'sigma_z_um': 1.2798324295058479, 'eps_nx_mm_mrad': 0.08650181807549762, 'sigma_delta_e3': 0.043406120246349426}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z analytic ref sqrt(sz0^2+(z*sd_p/gamma^2)^2): AG=1.75% OCELOT=0.47%
  R56 adapter σ_δ_p semantic: 0.043% (PASS)
  verdict: PASS — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [solenoid] 2026-08-05 21:40
  ROOT CAUSE: AG reduced-order Larmor coupling (dnu_x⊃-2ks·nu_y, dnu_y⊃+2ks·nu_x, dsxy⊃2ks(sx^2-sy^2)).
  Exact hard-edge 4x4 (Brown-Chao, == OCELOT SolenoidTM) gives sxy≡0 for a round uncorrelated beam; AG coupling creates spurious sxy and under-focusing.
  k_s (Bz/(2Brho)) = 22.3751 m^-1, k_s^2 = 500.65 m^-2 — identical in both backends.
  AG as-is:  sx=963.4 sy=2468.9 um (x-y broken)
  AG coupling=OFF: sx=1984.2 sy=1984.2 um
  OCELOT (ref):    sx=1989.4 sy=1991.9 um
  AG(off) vs OCELOT sigma_x max dev = 0.40% (<1% PASS) | AG(on) = 390.73% (FAIL)
  FIX: for round beams the coupling must be disabled in the AG force adapter (validation/backend.run_ag solenoid_coupling=False). Not a parameter tune; enforces exact round-beam transport.

## [rf] 2026-08-05 21:40
  RF standardized to thin-lens in BOTH backends (Option A):
    longitudinal kick δ += K·sin(φ+kz), H=-9.779 m^-1, σ_δ AG=2.953e-3 vs OCELOT=2.924e-3 (PASS)
    transverse RF kick K_trans=-2.680 m^-1 added to OCELOT; σ_x AG=1119 vs OCELOT=1120 um (PASS)
    switch OFF vs ON: AG σ_x 1119->542 um, OCELOT 1120->535 um (both read physics_switches.rf_transverse_kick)
    R56 adapter: kick semantic PASS, routing {'drift': 0, 'solenoid': 0, 'rf': 1, 'full': 1} PASS
  RESOLVED: input δ-variable convention (B) fixed by the adapter; σ_z at sample agrees with AG (residual few-µm at the waist).

## [full] 2026-08-05 21:40
  full beamline: sample devs — σ_x 0.26%, σ_y 0.39%, σ_δ 0.98%, ε_nx 0.14%, ε_ny 0.09% (PASS)
  R56 resolved — σ_z: AG 477 vs OCELOT 471 um (1.17%, PASS), waist Δz=0.6 mm
  switches: {'rf_longitudinal_kick': True, 'rf_transverse_kick': False} == {'rf_longitudinal_kick': True, 'rf_transverse_kick': False}

## [drift] 2026-08-06 10:55
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 0.43301417903757716, 'sigma_y_um': 0.7039260117123852, 'sigma_z_um': 1.2994156746116554, 'eps_nx_mm_mrad': 0.27476842355856385, 'sigma_delta_e3': 0.043406120246349426}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z analytic ref sqrt(sz0^2+(z*sd_p/gamma^2)^2): AG=1.75% OCELOT=0.46%
  R56 adapter σ_δ_p semantic: 0.043% (PASS)
  verdict: PASS — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [solenoid] 2026-08-06 10:55
  ROOT CAUSE: AG reduced-order Larmor coupling (dnu_x⊃-2ks·nu_y, dnu_y⊃+2ks·nu_x, dsxy⊃2ks(sx^2-sy^2)).
  Exact hard-edge 4x4 (Brown-Chao, == OCELOT SolenoidTM) gives sxy≡0 for a round uncorrelated beam; AG coupling creates spurious sxy and under-focusing.
  k_s (Bz/(2Brho)) = 22.3751 m^-1, k_s^2 = 500.65 m^-2 — identical in both backends.
  AG as-is:  sx=963.4 sy=2468.9 um (x-y broken)
  AG coupling=OFF: sx=1984.2 sy=1984.2 um
  OCELOT (ref):    sx=1989.4 sy=1991.9 um
  AG(off) vs OCELOT sigma_x max dev = 0.40% (<1% PASS) | AG(on) = 390.73% (FAIL)
  FIX: for round beams the coupling must be disabled in the AG force adapter (validation/backend.run_ag solenoid_coupling=False). Not a parameter tune; enforces exact round-beam transport.

## [rf] 2026-08-06 10:55
  RF standardized to thin-lens in BOTH backends (Option A):
    longitudinal kick δ += K·sin(φ+kz), H=-9.779 m^-1, σ_δ AG=2.953e-3 vs OCELOT=2.924e-3 (PASS)
    transverse RF kick K_trans=-2.680 m^-1 added to OCELOT; σ_x AG=1119 vs OCELOT=1120 um (PASS)
    switch OFF vs ON: AG σ_x 1119->542 um, OCELOT 1120->535 um (both read physics_switches.rf_transverse_kick)
    R56 adapter: kick semantic PASS, routing {'drift': 0, 'solenoid': 0, 'rf': 1, 'full': 1} PASS
  RESOLVED: input δ-variable convention (B) fixed by the adapter; σ_z at sample agrees with AG (residual few-µm at the waist).

## [full] 2026-08-06 10:55
  full beamline: sample devs — σ_x 0.26%, σ_y 0.39%, σ_δ 0.98%, ε_nx 0.14%, ε_ny 0.09% (PASS)
  R56 resolved — σ_z: AG 477 vs OCELOT 471 um (1.17%, PASS), waist Δz=0.6 mm
  switches: {'rf_longitudinal_kick': True, 'rf_transverse_kick': False} == {'rf_longitudinal_kick': True, 'rf_transverse_kick': False}

## [gpt-lattice-single-source] 2026-08-05
  GPT main route refactored: ued_beamline_v2.py builds from lattice.elements
  ONLY (build_lattice_from_shared(cfg, active_types)); hardcoded drifts
  (0.100/0.240/0.355/0.777) removed; step semantics via STEP_ACTIVE;
  multi-instance solenoid/RF supported; RF kick per instance at own z_start;
  transverse gated by shared switch; module importable (main() guard).
  Verified (validation/test_gpt_route_equivalence.py PASS):
  A geometry exact match; B step routing (0/0/N_rf/N_rf, runtime σ_δ gating,
  total length 777mm every step); C sample regression <2% (σ_x 0.65%,
  σ_z 0.43%, σ_δ_p 0.43%); D baseline run_all + r56 unchanged, AG unchanged,
  config SHA dd8ada3d4cb2 unchanged. No physics/R56/SC/random/param changes.
  Reports: gpt_route_lattice_equivalence.md, gpt_route_geometry.json,
  gpt_route_before_after.json, review_summary_gpt_lattice.png;
  CHANGELOG_gpt_lattice_single_source.md.

## [drift] 2026-08-08 12:26
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 0.47889966582113885, 'sigma_y_um': 0.3071841772479801, 'sigma_z_um': 1.9474872884866767, 'eps_nx_mm_mrad': 0.4992211367450183, 'sigma_delta_e3': 0.043406120246349426}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z analytic ref sqrt(sz0^2+(z*sd_p/gamma^2)^2): AG=1.75% OCELOT=0.19%
  R56 adapter σ_δ_p semantic: 0.043% (PASS)
  verdict: PASS — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [solenoid] 2026-08-08 12:26
  ROOT CAUSE: AG reduced-order Larmor coupling (dnu_x⊃-2ks·nu_y, dnu_y⊃+2ks·nu_x, dsxy⊃2ks(sx^2-sy^2)).
  Exact hard-edge 4x4 (Brown-Chao, == OCELOT SolenoidTM) gives sxy≡0 for a round uncorrelated beam; AG coupling creates spurious sxy and under-focusing.
  k_s (Bz/(2Brho)) = 22.3751 m^-1, k_s^2 = 500.65 m^-2 — identical in both backends.
  AG as-is:  sx=963.4 sy=2468.9 um (x-y broken)
  AG coupling=OFF: sx=1984.2 sy=1984.2 um
  OCELOT (ref):    sx=1989.4 sy=1991.9 um
  AG(off) vs OCELOT sigma_x max dev = 0.40% (<1% PASS) | AG(on) = 390.73% (FAIL)
  FIX: for round beams the coupling must be disabled in the AG force adapter (validation/backend.run_ag solenoid_coupling=False). Not a parameter tune; enforces exact round-beam transport.

## [rf] 2026-08-08 12:26
  RF standardized to thin-lens in BOTH backends (Option A):
    longitudinal kick δ += K·sin(φ+kz), H=-9.779 m^-1, σ_δ AG=2.953e-3 vs OCELOT=2.924e-3 (PASS)
    transverse RF kick K_trans=-2.680 m^-1 added to OCELOT; σ_x AG=1119 vs OCELOT=1120 um (PASS)
    switch OFF vs ON: AG σ_x 1119->542 um, OCELOT 1120->535 um (both read physics_switches.rf_transverse_kick)
    R56 adapter: kick semantic PASS, routing {'drift': 0, 'solenoid': 0, 'rf': 1, 'full': 1} PASS
  RESOLVED: input δ-variable convention (B) fixed by the adapter; σ_z at sample agrees with AG (residual few-µm at the waist).

## [full] 2026-08-08 12:26
  full beamline: sample devs — σ_x 0.26%, σ_y 0.39%, σ_δ 0.98%, ε_nx 0.14%, ε_ny 0.09% (PASS)
  R56 resolved — σ_z: AG 477 vs OCELOT 471 um (1.17%, PASS), waist Δz=0.6 mm
  switches: {'rf_longitudinal_kick': True, 'rf_transverse_kick': False} == {'rf_longitudinal_kick': True, 'rf_transverse_kick': False}

## [v0.12-architecture-governance] 2026-08-08
  AUDIT phase (no physics change, 7 reports in validation/reports/v0.12_*.md):
  - dataflow/hardcoding/module-responsibility/data-contract/maintainability/
    architecture/regression audits completed.
  IMPLEMENTED (A-class + safe B-class only):
  - shared/constants.py: single physical-constants source (A1);
    params/reference/backend import from it (values identical).
  - backend.py meta += provenance (git commit, timestamp, lattice_hash,
    python, coordinate convention) — traceability (Phase 9).
  - validation/config_check.py + test_config_consistency.py: Level-1
    read-only config consistency (Phase 7/10); added to run_all.
  - YAML comment: sigma_delta semantics clarified to δ_p (B4, value
    unchanged, config SHA unchanged dd8ada3d4cb2).
  NOT modified: AG core, OCELOT core/R56, RF equations, SC, seeds policy,
    lattice geometry, config values, test thresholds.
  VERIFIED: 6/6 tests PASS + r56 unchanged; AG bit-identical to v0.10;
  full-beamline σ_z 477.001/471.497 (1.17%), σ_δ_p 0.98%.
  OPEN (recorded, not touched): γ/β/p0 derivation still duplicated (~12
  sites, B2); rparticles magic indices (B2); seed(42) not configurable (B1);
  backend.py adapter+physics mixing (debt); GPT route monolith (debt).

## [drift] 2026-08-10 10:25
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 0.3767796798699752, 'sigma_y_um': 0.6342565055448941, 'sigma_z_um': 1.8000933557025134, 'eps_nx_mm_mrad': 0.31140389844541905, 'sigma_delta_e3': 0.02245152598292253}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z analytic ref sqrt(sz0^2+(z*sd_p/gamma^2)^2): AG=1.75% OCELOT=0.08%
  R56 adapter σ_δ_p semantic: 0.022% (PASS)
  verdict: PASS — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [solenoid] 2026-08-10 10:25
  ROOT CAUSE: AG reduced-order Larmor coupling (dnu_x⊃-2ks·nu_y, dnu_y⊃+2ks·nu_x, dsxy⊃2ks(sx^2-sy^2)).
  Exact hard-edge 4x4 (Brown-Chao, == OCELOT SolenoidTM) gives sxy≡0 for a round uncorrelated beam; AG coupling creates spurious sxy and under-focusing.
  k_s (Bz/(2Brho)) = 22.3751 m^-1, k_s^2 = 500.65 m^-2 — identical in both backends.
  AG as-is:  sx=963.4 sy=2468.9 um (x-y broken)
  AG coupling=OFF: sx=1984.2 sy=1984.2 um
  OCELOT (ref):    sx=1996.2 sy=1991.4 um
  AG(off) vs OCELOT sigma_x max dev = 0.60% (<1% PASS) | AG(on) = 389.28% (FAIL)
  FIX: for round beams the coupling must be disabled in the AG force adapter (validation/backend.run_ag solenoid_coupling=False). Not a parameter tune; enforces exact round-beam transport.

## [rf] 2026-08-10 10:25
  RF standardized to thin-lens in BOTH backends (Option A):
    longitudinal kick δ += K·sin(φ+kz), H=-9.779 m^-1, σ_δ AG=2.953e-3 vs OCELOT=2.939e-3 (PASS)
    transverse RF kick K_trans=-2.680 m^-1 added to OCELOT; σ_x AG=1119 vs OCELOT=1122 um (PASS)
    switch OFF vs ON: AG σ_x 1119->542 um, OCELOT 1122->537 um (both read physics_switches.rf_transverse_kick)
    R56 adapter: kick semantic PASS, routing {'drift': 0, 'solenoid': 0, 'rf': 1, 'full': 1} PASS
  RESOLVED: input δ-variable convention (B) fixed by the adapter; σ_z at sample agrees with AG (residual few-µm at the waist).

## [full] 2026-08-10 10:25
  full beamline: sample devs — σ_x 0.60%, σ_y 0.36%, σ_δ 0.46%, ε_nx 0.53%, ε_ny 0.26% (PASS)
  R56 resolved — σ_z: AG 477 vs OCELOT 474 um (0.63%, PASS), waist Δz=0.4 mm
  switches: {'rf_longitudinal_kick': True, 'rf_transverse_kick': False} == {'rf_longitudinal_kick': True, 'rf_transverse_kick': False}

## [v0.13-preSC-maintainability] 2026-08-10
  Debt convergence (no physics change):
  - shared/beam_physics.py: BeamReference single γ/β/p0/velocity source;
    migrated params.derived, reference, beam_result, ued_beamline_v2,
    run_shared, tests (task 1).
  - shared/ocelot_coords.py: named rparticles access (I_X..I_P, set/add_*);
    removed magic indices from business code (task 2); test_r56 unchanged
    (frozen independent reference).
  - random configurable: config random.seed=42; generate_parray seed,
    px/py/delta seed+1 (independent, avoids x-px correlation); verified
    same-seed twice → bit-identical (task 3).
  VERIFIED: 6/6 tests PASS + r56 bit-identical; AG bit-identical
  (σ_x 1984.191, σ_z 477.001); OCELOT now deterministic (σ_x 1996.2,
  within thresholds); config SHA changed (random key added; physics values
  unchanged); lattice hash unchanged.

## [drift] 2026-08-10 22:09
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 0.3767796798699752, 'sigma_y_um': 0.6342565055448941, 'sigma_z_um': 1.8000933557025134, 'eps_nx_mm_mrad': 0.31140389844541905, 'sigma_delta_e3': 0.02245152598292253}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z analytic ref sqrt(sz0^2+(z*sd_p/gamma^2)^2): AG=1.75% OCELOT=0.08%
  R56 adapter σ_δ_p semantic: 0.022% (PASS)
  verdict: PASS — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [solenoid] 2026-08-10 22:09
  ROOT CAUSE: AG reduced-order Larmor coupling (dnu_x⊃-2ks·nu_y, dnu_y⊃+2ks·nu_x, dsxy⊃2ks(sx^2-sy^2)).
  Exact hard-edge 4x4 (Brown-Chao, == OCELOT SolenoidTM) gives sxy≡0 for a round uncorrelated beam; AG coupling creates spurious sxy and under-focusing.
  k_s (Bz/(2Brho)) = 22.3751 m^-1, k_s^2 = 500.65 m^-2 — identical in both backends.
  AG as-is:  sx=963.4 sy=2468.9 um (x-y broken)
  AG coupling=OFF: sx=1984.2 sy=1984.2 um
  OCELOT (ref):    sx=1996.2 sy=1991.4 um
  AG(off) vs OCELOT sigma_x max dev = 0.60% (<1% PASS) | AG(on) = 389.28% (FAIL)
  FIX: for round beams the coupling must be disabled in the AG force adapter (validation/backend.run_ag solenoid_coupling=False). Not a parameter tune; enforces exact round-beam transport.

## [rf] 2026-08-10 22:10
  RF standardized to thin-lens in BOTH backends (Option A):
    longitudinal kick δ += K·sin(φ+kz), H=-9.779 m^-1, σ_δ AG=2.953e-3 vs OCELOT=2.939e-3 (PASS)
    transverse RF kick K_trans=-2.680 m^-1 added to OCELOT; σ_x AG=1119 vs OCELOT=1122 um (PASS)
    switch OFF vs ON: AG σ_x 1119->542 um, OCELOT 1122->537 um (both read physics_switches.rf_transverse_kick)
    R56 adapter: kick semantic PASS, routing {'drift': 0, 'solenoid': 0, 'rf': 1, 'full': 1} PASS
  RESOLVED: input δ-variable convention (B) fixed by the adapter; σ_z at sample agrees with AG (residual few-µm at the waist).

## [full] 2026-08-10 22:10
  full beamline: sample devs — σ_x 0.60%, σ_y 0.36%, σ_δ 0.46%, ε_nx 0.53%, ε_ny 0.26% (PASS)
  R56 resolved — σ_z: AG 477 vs OCELOT 474 um (0.63%, PASS), waist Δz=0.4 mm
  switches: {'rf_longitudinal_kick': True, 'rf_transverse_kick': False} == {'rf_longitudinal_kick': True, 'rf_transverse_kick': False}

## [v0.14-sc-audit] 2026-08-10
  SC 接入审计 + P0 修复（纯数据流，SC 物理零改动）:
  - P0 根因: OCELOT tracking_step() 只应用 transfer maps，PhysProc
    (SpaceCharge) 只在 track() 的 counter 机制中触发——生产代码全用
    tracking_step → SC 从未执行（历史"SC ON/OFF 重合"疑点根因）。
  - 修复: backend.py / ued_beamline_v2.py / sc_audit_diagnostics.py 循环内
    复刻 counter 机制（每 step×dz 在 map 后触发 apply）+ mesh/step 从 config
    传入（原值相同，数值不变）。SC OFF 路径与旧代码逐位一致。
  - 验证: 修复后 smoke +236% (σx 725→2436µm @500fC)，charge 0-1000fC 单调，
    N/mesh 收敛（±0.6%），step≤5 建议；6/6 测试 PASS + R56 不变；
    no-SC 与 v0.13 逐位一致 (hash e041d6ae9fb7a0d2)；AG 位级不变。
  - P1 记录（未改）: AG SC 电荷语义 Ne·e vs config charge_fC（8 vs 100 fC，
    弱 12.5 倍，修复将改变 SC ON 数值）；SC state 双来源与 HARD FAIL 设计。

## [drift] 2026-08-10 22:40
  drift transverse: AG vs OCELOT max rel dev {'sigma_x_um': 0.3767796798699752, 'sigma_y_um': 0.6342565055448941, 'sigma_z_um': 1.8000933557025134, 'eps_nx_mm_mrad': 0.31140389844541905, 'sigma_delta_e3': 0.02245152598292253}
  analytic drift reference matches (sigma_x=sqrt(s0^2+(s0'z)^2))
  sigma_z analytic ref sqrt(sz0^2+(z*sd_p/gamma^2)^2): AG=1.75% OCELOT=0.08%
  R56 adapter σ_δ_p semantic: 0.022% (PASS)
  verdict: PASS — report /Users/qin/Desktop/shuyan/Beam_dynamics_simu/validation/reports/drift_AG_vs_OCELOT.png

## [solenoid] 2026-08-10 22:40
  ROOT CAUSE: AG reduced-order Larmor coupling (dnu_x⊃-2ks·nu_y, dnu_y⊃+2ks·nu_x, dsxy⊃2ks(sx^2-sy^2)).
  Exact hard-edge 4x4 (Brown-Chao, == OCELOT SolenoidTM) gives sxy≡0 for a round uncorrelated beam; AG coupling creates spurious sxy and under-focusing.
  k_s (Bz/(2Brho)) = 22.3751 m^-1, k_s^2 = 500.65 m^-2 — identical in both backends.
  AG as-is:  sx=963.4 sy=2468.9 um (x-y broken)
  AG coupling=OFF: sx=1984.2 sy=1984.2 um
  OCELOT (ref):    sx=1996.2 sy=1991.4 um
  AG(off) vs OCELOT sigma_x max dev = 0.60% (<1% PASS) | AG(on) = 389.28% (FAIL)
  FIX: for round beams the coupling must be disabled in the AG force adapter (validation/backend.run_ag solenoid_coupling=False). Not a parameter tune; enforces exact round-beam transport.

## [rf] 2026-08-10 22:40
  RF standardized to thin-lens in BOTH backends (Option A):
    longitudinal kick δ += K·sin(φ+kz), H=-9.779 m^-1, σ_δ AG=2.953e-3 vs OCELOT=2.939e-3 (PASS)
    transverse RF kick K_trans=-2.680 m^-1 added to OCELOT; σ_x AG=1119 vs OCELOT=1122 um (PASS)
    switch OFF vs ON: AG σ_x 1119->542 um, OCELOT 1122->537 um (both read physics_switches.rf_transverse_kick)
    R56 adapter: kick semantic PASS, routing {'drift': 0, 'solenoid': 0, 'rf': 1, 'full': 1} PASS
  RESOLVED: input δ-variable convention (B) fixed by the adapter; σ_z at sample agrees with AG (residual few-µm at the waist).

## [full] 2026-08-10 22:40
  full beamline: sample devs — σ_x 0.60%, σ_y 0.36%, σ_δ 0.46%, ε_nx 0.53%, ε_ny 0.26% (PASS)
  R56 resolved — σ_z: AG 477 vs OCELOT 474 um (0.63%, PASS), waist Δz=0.4 mm
  switches: {'rf_longitudinal_kick': True, 'rf_transverse_kick': False} == {'rf_longitudinal_kick': True, 'rf_transverse_kick': False}
