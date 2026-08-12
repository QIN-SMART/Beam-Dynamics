"""
Backend drivers for the validation framework.

Both backends are driven through their EXISTING public APIs only:
  - AG : AG/beam_dynamics_6d.propagate + AG/external_forces
         (make_beam_100keV, ExtFieldRegion, propagate, build_all_external)
  - OCELOT : ocelot Drift/Solenoid/Cavity elements + generate_parray + Navigator
             (same pattern as GPT模拟/ued_beamline_v2.py / solenoid benchmark)

No physics kernel is re-implemented here.

Sections:
  'drift'    drift only            (cathode → sample, no elements)
  'solenoid' drift + solenoid + drift
  'rf'       drift + rf + drift
"""

import os
import sys

import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_THIS_DIR)
for p in (_REPO, _THIS_DIR, os.path.join(_REPO, "AG")):
    if p not in sys.path:
        sys.path.insert(0, p)

from shared.params import load_config, parse, derived, _lattice_elements, config_sha  # noqa: E402
from shared.constants import C_SI, M_E_SI, E_SI, MEC2_KEV  # noqa: E402
from shared.ocelot_coords import (add_p_oc, add_px, add_py,  # noqa: E402
                                  set_px, set_py, set_p_oc)
from beam_result import BeamResult  # noqa: E402


def _provenance(cfg):
    """Provenance metadata for a simulation run (v0.12 traceability).

    Recorded in every BeamResult.meta so any figure/result answers
    'which version, which parameters, which backend produced this?'.
    """
    import subprocess
    from datetime import datetime
    import json as _json
    try:
        commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                                capture_output=True, text=True, timeout=5)
        git_commit = commit.stdout.strip() or "unknown"
    except Exception:
        git_commit = "unknown"
    try:
        lattice_hash = _json.dumps(_lattice_elements(cfg),
                                   sort_keys=True, default=str)
        lattice_hash = hashlib_hex(lattice_hash)
    except Exception:
        lattice_hash = "unknown"
    return {
        "git_commit": git_commit,
        "config_sha": config_sha(cfg),
        "lattice_hash": lattice_hash,
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "python": ".".join(__import__("sys").version.split()[:1]),
        "random_seed": int(cfg.get("random", {}).get("seed", 0)),
        "coordinate_convention": "delta_p = dp/p0; OCELOT native p_oc = dE/(c*p0)",
    }


def hashlib_hex(s):
    import hashlib
    return hashlib.sha1(s.encode("utf-8")).hexdigest()[:12]


# ════════════════════════════════════════════════════════════════════════
# Standardized RF (thin-lens) model — shared by both backends.
#   δ = Δp/p₀,  kick: δ_out = δ_in + K·sin(φ+kz),  K = eV/(β²E₀)
#   chirp H = ∂δ/∂z = eV·k·cosφ/(β²E₀)   [m⁻¹]
#   transverse RF kick K_trans (Panofsky-Wenzel, AG-consistent)  [m⁻¹]
# ════════════════════════════════════════════════════════════════════════

def _rf_constants(cfg, d, rf_elem=None):
    """Return (H, K_trans, k_rf, E_rf) for ONE rf cavity element.

    rf_elem : lattice element dict of type 'rf_cavity' (default: first one).
    """
    from shared.params import elem_params, first_of_type
    rf = elem_params(rf_elem if rf_elem is not None
                     else first_of_type(cfg, "rf_cavity"))
    k_rf = 2.0 * np.pi * rf["frequency_GHz"] * 1e9 / C_SI
    V = rf["voltage_kV"] * 1e3
    E_rf = V / rf["length_m"]
    gamma, beta = d["gamma"], d["beta"]
    ETA = E_SI / (M_E_SI * C_SI**2)
    H = E_SI * V * k_rf * np.cos(rf["phase_rad"]) / (beta**2 * gamma * M_E_SI * C_SI**2)
    K_trans = -ETA * k_rf * V / (2.0 * gamma * beta)
    return H, K_trans, k_rf, E_rf


# ════════════════════════════════════════════════════════════════════════
# shared section geometry derived from the shared config (lattice only)
# ════════════════════════════════════════════════════════════════════════

def _section_regions(cfg, section):
    """z-spans (z0, z1, etype) for a section, derived from lattice.elements.

    Multi-instance safe: 'solenoid' activates ALL solenoids, 'rf' activates
    ALL rf cavities, 'full' activates every element with its real type,
    'drift' has no active elements.  Drift spans fill the gaps between active
    elements (cathode/sample markers are zero-length).
    """
    from shared.params import _lattice_elements, z_sample
    spans = []
    z_prev = 0.0
    for e in _lattice_elements(cfg):
        z0, z1, L = e["z_start"], e["z_start"] + e["length"], e["length"]
        if L <= 0:
            continue                       # cathode / sample markers
        etype = e["type"]
        if section == "solenoid":
            active = "solenoid" if etype == "solenoid" else "drift"
        elif section == "rf":
            active = "rf" if etype == "rf_cavity" else "drift"
        elif section == "full":
            active = "rf" if etype == "rf_cavity" else etype
        else:
            active = "drift"
        if z0 > z_prev + 1e-12:
            spans.append((z_prev, z0, "drift"))
        spans.append((z0, z1, active))
        z_prev = z1
    if z_sample(cfg) > z_prev + 1e-12:
        spans.append((z_prev, z_sample(cfg), "drift"))
    return spans


# ════════════════════════════════════════════════════════════════════════
# AG backend
# ════════════════════════════════════════════════════════════════════════

def run_ag(cfg, section, n_points=2000, sc_enabled=None,
           solenoid_coupling=True, switches=None):
    """Drive the AG 6D envelope ODE through the section.

    RF is handled by the STANDARDIZED THIN-LENS model (see _rf_constants):
    the continuous cavity chirp is replaced by apply_rf_thin_lens(H) at the
    cavity entrance; the RF TRANSVERSE force follows the shared
    physics_switches (default OFF).
    """
    from beam_dynamics_6d import (make_beam_100keV, ExtFieldRegion,
                                  get_alpha_interpolators, propagate,
                                  apply_rf_thin_lens, Beam6D)
    from external_forces import build_all_external, build_external_force_func

    P = parse(cfg)
    d = derived(cfg)
    ib = cfg["initial_distribution"]
    sc_enabled = cfg["space_charge"]["enabled"] if sc_enabled is None else sc_enabled
    sw = dict(P.switches.as_dict())
    sw.update(switches or {})

    get_alpha_interpolators()

    # v0.14.1 task 2 (AG charge semantics): the AG SC force is ∝ Ne·e
    # (beam_dynamics_6d.py:379 fb = η·Ne·e/(8π√π·ε₀)).  Ne must be the
    # PHYSICAL bunch charge Q/e — beam.n_particles is only the OCELOT
    # macroparticle numerical resolution.  (Previously Ne=n_particles=5e4
    # → 8 fC equivalent charge instead of the configured 100 fC: 12.5× too
    # weak.  SC OFF is unaffected: Ne is forced to 0.0 below.)
    ne_phys = abs(P.beam.Q_C) / E_SI

    beam0 = make_beam_100keV(
        Ne=ne_phys,
        beamK_eV=P.beam.energy_keV * 1e3,
        sigma_x0_um=ib["spot_rms_um"], sigma_y0_um=ib["spot_rms_um"],
        sigma_z0_um=ib["bunch_length_um"],
        sigma_delta=ib["sigma_delta"],
        eps_nx_um=ib["epsilon_n_mm_mrad"], eps_ny_um=ib["epsilon_n_mm_mrad"],
        eps_nz_um=ib["epsilon_nz_mm_mrad"],
    )

    regions = []
    active_types = {"drift": {"drift"},
                    "solenoid": {"solenoid", "drift"},
                    "rf": {"rf_cavity", "drift"},
                    "full": {"solenoid", "rf_cavity", "drift"}}[section]
    for e in _lattice_elements(cfg):
        etype, z0, z1, L = e["type"], e["z_start"], e["z_start"] + e["length"], e["length"]
        if L <= 0:
            continue                       # cathode / sample markers
        if etype not in active_types:
            continue                       # non-active → field-free drift
        if etype == "solenoid":
            regions.append(ExtFieldRegion(z0, z1, "solenoid",
                                          {"Bz": e["parameters"]["B_field_T"],
                                           "dBz_dz": 0.0}))
        elif etype == "rf_cavity":
            p = e["parameters"]
            k_rf = 2.0 * np.pi * p["frequency_GHz"] * 1e9 / C_SI
            E_rf = p["voltage_kV"] * 1e3 / L
            regions.append(ExtFieldRegion(z0, z1, "rf",
                                          {"E_rf": E_rf, "k_rf": k_rf,
                                           "dE_dz": 0.0, "dE_cdt": k_rf * E_rf}))

    # core external fields (solenoid only); RF kept ONLY for its transverse
    # force when the shared switch is ON — the continuous RF chirp/acceleration
    # is replaced by the thin lens (also switch-gated).
    regions_core = [r for r in regions if r.ftype != "rf"]
    rf_regions = ([r for r in regions if r.ftype == "rf"]
                  if sw["rf_transverse_kick"] else [])
    ef_core, gp_core, ch_core = build_all_external(regions_core)

    ef = ef_core
    if rf_regions:
        ef_rf = build_external_force_func(rf_regions)
        ef = lambda z, beam, _core=ef_core, _rf=ef_rf: tuple(
            a + b for a, b in zip(_core(z, beam), _rf(z, beam)))

    if not solenoid_coupling and any(r.ftype == "solenoid" for r in regions_core):
        base = ef
        ef = (lambda z, beam, _base=base:
              _base(z, beam)[:3] + (0.0, 0.0, 0.0))

    z_max = _section_regions(cfg, section)[-1][1]
    beam_k = beam0.copy()
    if not sc_enabled:
        beam_k.Ne = 0.0
    sc_model = "ellipsoid" if sc_enabled else "gaussian"

    rf_elems = ([e for e in _lattice_elements(cfg)
                 if e["type"] == "rf_cavity" and e["length"] > 0
                 and "rf_cavity" in active_types]
                if sw["rf_longitudinal_kick"] else [])
    if rf_elems:
        # thin-lens RF (multi-instance): propagate to each cavity entrance,
        # apply H·z chirp, continue to the next one.
        z_parts, st_parts = [], []
        z_cur = 0.0
        beam_cur = beam_k
        n_per = max(200, n_points // (len(rf_elems) + 1))
        for e in rf_elems:
            z_rf = e["z_start"]
            if z_rf > z_cur:
                z1, st1 = propagate(beam_cur, (z_cur, z_rf), n_points=n_per,
                                    external_force_func=ef, gamma_prime_func=gp_core,
                                    rf_chirp_func=ch_core, sc_model=sc_model)
                z_parts.append(z1); st_parts.append(st1)
                beam_cur = Beam6D.from_state(st1[-1, :11], st1[-1, 11], beam_k.Ne,
                                             beam_k.eps_nx, beam_k.eps_ny, beam_k.eps_nz)
            H, _, _, _ = _rf_constants(cfg, d, e)
            beam_cur = apply_rf_thin_lens(beam_cur, H)
            z_cur = z_rf
        if z_cur < z_max:
            z2, st2 = propagate(beam_cur, (z_cur, z_max), n_points=n_per,
                                external_force_func=ef, gamma_prime_func=gp_core,
                                rf_chirp_func=ch_core, sc_model=sc_model)
            z_parts.append(z2); st_parts.append(st2)
        z_arr = np.concatenate([zz if i == 0 else zz[1:]
                                for i, zz in enumerate(z_parts)])
        st = np.concatenate([ss if i == 0 else ss[1:]
                             for i, ss in enumerate(st_parts)])
    else:
        z_arr, st = propagate(beam_k, (0.0, z_max), n_points=n_points,
                              external_force_func=ef, gamma_prime_func=gp_core,
                              rf_chirp_func=ch_core, sc_model=sc_model)

    from beamline_sim import compute_emittance
    em = compute_emittance(st, beam0.eps_nx, beam0.eps_ny, beam0.eps_nz)

    result = BeamResult(
        route="AG",
        z_mm=z_arr * 1e3,
        sigma_x_um=st[:, 0] * 1e6,
        sigma_y_um=st[:, 1] * 1e6,
        sigma_z_um=st[:, 2] * 1e6,
        eps_nx_mm_mrad=em["eps_n_x"] * 1e6,
        eps_ny_mm_mrad=em["eps_n_y"] * 1e6,
        energy_keV=np.full_like(z_arr, P.beam.energy_keV),
        sigma_delta_e3=st[:, 4] * 1e3,
        meta={"section": section, "solenoid_coupling": solenoid_coupling,
              "sc_enabled": sc_enabled, "switches": sw, "config_sha": config_sha(cfg),
              "provenance": _provenance(cfg),
              # v0.14.1 task 3 AG state metadata (formal contract):
              "sc_requested": bool(sc_enabled),
              "sc_effective": bool(sc_enabled),   # AG: requested==effective
                                                   # (adapter injects Ne=Q/e)
              "physical_charge_C": float(P.beam.Q_C),
              "physical_electron_number": float(ne_phys),
              "ag_ne_phys": float(ne_phys),       # alias (task 2 compat, NOT
                                                   # the formal contract)
              "rf": "thin-lens (H) + transverse kick" if sw["rf_transverse_kick"]
              else "thin-lens (H) only"},
    )
    return result


# ════════════════════════════════════════════════════════════════════════
# OCELOT adapter helpers
# ════════════════════════════════════════════════════════════════════════

def _ocelot_rf_kick(p, rf_elem, cfg, d, sw):
    """Standardized thin-lens RF kick applied to an OCELOT ParticleArray.

    Longitudinal: the shared formula gives a momentum-relative kick
      d_delta_p = Δp/p0 = K·sin(φ + k·z_phys)
    which is stored in the OCELOT native coordinate
      p_oc = ΔE/(c·p0) = β0·δ_p
    hence the conversion  d_p_oc = β0 · d_delta_p.  (R56 audit: classification B.)
    Transverse (switch-gated, unchanged equation): x' += K_trans·x.
    """
    pv = rf_elem["parameters"]
    H_rf, K_trans, k_rf, E_rf = _rf_constants(cfg, d, rf_elem)
    E_tot = d["gamma"] * MEC2_KEV * 1e3          # eV (from derived, single source)
    tau = p.tau()
    z_phys = -d["beta"] * tau
    d_delta_p = (pv["voltage_kV"] * 1e3 / (d["beta"]**2 * E_tot)) \
        * np.sin(pv["phase_rad"] + k_rf * z_phys)
    add_p_oc(p, d["beta"] * d_delta_p)             # p_oc += β0·δ_p (adapter boundary)
    if sw["rf_transverse_kick"]:
        x = p.x(); y = p.y()
        add_px(p, K_trans * x)
        add_py(p, K_trans * y)


# ════════════════════════════════════════════════════════════════════════
# OCELOT backend
# ════════════════════════════════════════════════════════════════════════

def run_ocelot(cfg, section, n_particles=None, dz=0.001, sc_enabled=None,
               switches=None, sc_step=None, sc_mesh=None):
    """Drive OCELOT macroparticle tracking through the section.

    Uses ocelot's Drift / Solenoid / Cavity elements and Navigator —
    identical element physics to GPT模拟/ued_beamline_v2.py.
    """
    from ocelot.cpbd.elements import Drift, Solenoid
    from ocelot.cpbd.magnetic_lattice import MagneticLattice
    from ocelot.cpbd.beam import generate_parray
    from ocelot.cpbd.navi import Navigator
    from ocelot.cpbd.track import tracking_step

    P = parse(cfg)
    d = derived(cfg)
    ib = cfg["initial_distribution"]

    n_particles = n_particles or P.beam.n_particles
    sc_enabled = cfg["space_charge"]["enabled"] if sc_enabled is None else sc_enabled
    sw = dict(P.switches.as_dict())
    sw.update(switches or {})

    def emit(x, xp):
        return np.sqrt(np.mean(x**2) * np.mean(xp**2) - np.mean(x * xp)**2)

    # ── lattice for the section (built from lattice.elements only) ──
    # Non-active elements keep their length as a plain drift so the total
    # beamline length is preserved.
    elems = []
    active_types = {"drift": {"drift"},
                    "solenoid": {"solenoid", "drift"},
                    "rf": {"rf_cavity", "drift"},
                    "full": {"solenoid", "rf_cavity", "drift"}}[section]
    for e in _lattice_elements(cfg):
        etype, z0, z1, L = e["type"], e["z_start"], e["z_start"] + e["length"], e["length"]
        if L <= 0:
            # v0.14.1 task 1: with SC ON the zero-length cathode/sample
            # markers are KEPT in the runtime sequence so the PhysProc can be
            # attached as cathode → sample (SC coverage ends at the sample
            # position, NOT at the start of the last non-zero element).
            # SC OFF keeps the previous sequence (markers skipped) — the
            # tracking_step loop below stays bitwise identical to v0.13.
            if sc_enabled and etype in ("cathode", "sample"):
                elems.append(Drift(l=0.0, eid=f"{e['name']}_MARKER"))
            continue
        if etype == "solenoid" and "solenoid" in active_types:
            k_sol = 1.602176634e-19 * e["parameters"]["B_field_T"] / (2.0 * d["p_SI"])
            elems.append(Solenoid(l=L, k=k_sol, eid=e["name"]))
        elif etype == "rf_cavity" and "rf_cavity" in active_types:
            elems.append(Drift(l=L, eid=f"{e['name']}_BODY"))  # analytic kick below
        else:                              # drift, or non-active element
            elems.append(Drift(l=L, eid=e["name"]))
    lat = MagneticLattice(elems, method={"global": "SecondTM"})
    lat.update_transfer_maps()
    total_length = sum(e.l for e in elems)

    # ── beam ──
    # Random policy (v0.13): configured seed for the ENTIRE beam, so the same
    # seed reproduces the same beam bit-for-bit.  seed → generate_parray
    # (x/y/tau), seed+1 → px/py/δ_p (independent draws; avoids the x–px
    # correlation trap of reusing the identical sequence).
    rng_seed = int(cfg["random"]["seed"])
    np.random.seed(rng_seed)
    p = generate_parray(
        sigma_x=ib["spot_rms_um"] * 1e-6, sigma_y=ib["spot_rms_um"] * 1e-6,
        sigma_tau=ib["bunch_length_um"] * 1e-6 / d["beta"],   # τ=c·t [m]
        energy=(P.beam.energy_keV + 511.0) * 1e-6,  # TOTAL energy in GeV (E_kin+mc²)
        charge=P.beam.Q_C,
        nparticles=n_particles,
    )
    np.random.seed(rng_seed + 1)
    N = p.rparticles.shape[1]
    set_px(p, np.random.normal(0.0, d["sigma_xp"], N))
    set_py(p, np.random.normal(0.0, d["sigma_yp"], N))
    # OCELOT p() = ΔE/(c·p0), NOT Δp/p0.  Convert the shared momentum
    # deviation to the native coordinate:  p_oc = β0·δ_p  (R56 audit: B).
    delta_p = np.random.normal(0.0, ib["sigma_delta"], N)
    set_p_oc(p, d["beta"] * delta_p)

    # ── SC runtime state contract (v0.14.1 task 3) ──
    # sc_requested → available → configured → attached → apply_count → effective.
    # HARD FAIL on every broken link when requested: NO silent fallback.
    from shared.sc_state import SCState  # noqa: E402
    sc_state = SCState(requested=bool(sc_enabled))
    sc_proc = None
    if sc_enabled:
        # HARD FAIL on import: `except ImportError: pass` is FORBIDDEN.
        from ocelot.cpbd.sc import SpaceCharge   # ImportError propagates
        sc_state.available = True
        sc_cfg = cfg["space_charge"]
        try:
            # v0.14 P1 fix: mesh/step must come from shared config, not
            # defaults; sc_step/sc_mesh are optional overrides (scheduler
            # tests only)
            sc_proc = SpaceCharge(step=sc_step or sc_cfg.get("step", 1),
                                  nmesh_xyz=list(sc_mesh or sc_cfg.get("mesh", [63, 63, 63])))
        except Exception as e:
            sc_state.fail(f"SpaceCharge construction failed: {e!r}")
        sc_state.configured = True

    navi = Navigator(lat, unit_step=dz)
    if sc_proc is not None:
        try:
            navi.add_physics_proc(sc_proc, lat.sequence[0], lat.sequence[-1])
        except Exception as e:
            sc_state.fail(f"add_physics_proc failed: {e!r}")
        sc_state.attached = True
        sc_state.coverage_start_m = float(sc_proc.s_start)
        sc_state.coverage_stop_m = float(sc_proc.s_stop)

    # RF kicks at each cavity entrance — standardized thin-lens model
    # (see _ocelot_rf_kick; longitudinal in δ_p converted to p_oc = β0·δ_p,
    #  transverse part follows the shared physics_switches).
    rf_elems = ([e for e in _lattice_elements(cfg)
                 if e["type"] == "rf_cavity" and e["length"] > 0
                 and "rf_cavity" in active_types]
                if sw["rf_longitudinal_kick"] else [])
    rf_kicks_applied = 0
    rf_done = set()
    sc_apply_count = 0
    if sc_proc is not None:
        # ── SC ON: OCELOT NATIVE scheduler (v0.14.1 task 1) ──
        # Navigator.get_next_step() is the exact mechanism of ocelot.track():
        # PhysProc (SpaceCharge) is triggered by the Navigator process counter,
        # and its coverage [s_start, s_stop) is defined by the add_physics_proc
        # anchors (cathode → sample).  The former manual counter clone is
        # retired (scheduler characterization T5–T7 showed it is not generally
        # equivalent: lost tail / out-of-coverage applies).
        z_list, sx_list, sy_list, sz_list = [], [], [], []
        enx_list, eny_list, sd_list = [], [], []
        sc_events = []                       # (z, zstep) per SC apply
        for t_maps, dz_step, proc_list, phys_steps in navi.get_next_step():
            # RF kick at the step START (same semantics as the SC OFF loop:
            # kick is applied to the particle state at the cavity entrance).
            z_start = navi.z0 - dz_step
            for rf_elem in rf_elems:
                z_rf = rf_elem["z_start"]
                if z_rf not in rf_done and z_start >= z_rf - 1e-12:
                    _ocelot_rf_kick(p, rf_elem, cfg, d, sw)
                    rf_done.add(z_rf)
                    rf_kicks_applied += 1
            for tm in t_maps:
                tm.apply(p)
            for proc, zstep in zip(proc_list, phys_steps):
                proc.z0 = navi.z0
                proc.apply(p, zstep)
                if proc is sc_proc:
                    sc_apply_count += 1
                    sc_events.append((float(navi.z0), float(zstep)))
            z = navi.z0
            x = p.x(); xp = p.px(); y = p.y(); yp = p.py()
            z_list.append(z)
            sx_list.append(np.std(x)); sy_list.append(np.std(y))
            sz_list.append(np.std(p.tau()) * d["beta"])
            enx_list.append(emit(x, xp) * d["beta_gamma"])
            eny_list.append(emit(y, yp) * d["beta_gamma"])
            sd_list.append(np.std(p.p()) / d["beta"])
        z_arr = np.array(z_list)
        sx = np.array(sx_list); sy = np.array(sy_list); sz = np.array(sz_list)
        enx = np.array(enx_list); eny = np.array(eny_list); sd = np.array(sd_list)
        # HARD FAIL (post-tracking): apply_count>0 and coverage == cathode→sample
        sc_state.apply_count = sc_apply_count
        sample_elems = [e for e in _lattice_elements(cfg) if e["type"] == "sample"]
        exp_stop = float(sample_elems[0]["z_start"]) if sample_elems else 0.0
        sc_state.verify_final(0.0, exp_stop)
        meta = {"section": section, "sc_enabled": sc_enabled, "switches": sw,
                "config_sha": config_sha(cfg), "provenance": _provenance(cfg),
                "longitudinal_native_coordinate": "p_oc = dE/(c*p0)",
                "reported_delta": "delta_p = dp/p0",
                "conversion_beta": float(d["beta"]),
                "rf_kicks_applied": rf_kicks_applied,
                "sc_scheduler": "ocelot_native",
                "sc_events": sc_events,
                "rf": "thin-lens (K·sin) + transverse kick K_trans"
                if sw["rf_transverse_kick"] else "thin-lens (K·sin) only"}
        meta.update(sc_state.to_meta())
    else:
        # ── SC OFF: UNCHANGED loop (tracking_step) — bitwise identical to
        # v0.13 (no PhysProc attached, sc_proc is None).
        n_steps = int(total_length / dz) if total_length > 0 else 1
        z_arr = np.zeros(n_steps); sx = np.zeros(n_steps); sy = np.zeros(n_steps)
        sz = np.zeros(n_steps); enx = np.zeros(n_steps); eny = np.zeros(n_steps)
        sd = np.zeros(n_steps)
        for i in range(n_steps):
            z_before = navi.z0
            for rf_elem in rf_elems:
                z_rf = rf_elem["z_start"]
                if z_rf not in rf_done and z_before >= z_rf - 1e-12:
                    _ocelot_rf_kick(p, rf_elem, cfg, d, sw)
                    rf_done.add(z_rf)
                    rf_kicks_applied += 1
            tracking_step(lat, p, dz, navi)
            z = navi.z0
            x = p.x(); xp = p.px(); y = p.y(); yp = p.py()
            z_arr[i] = z
            sx[i] = np.std(x); sy[i] = np.std(y)
            sz[i] = np.std(p.tau()) * d["beta"]
            enx[i] = emit(x, xp) * d["beta_gamma"]
            eny[i] = emit(y, yp) * d["beta_gamma"]
            sd[i] = np.std(p.p()) / d["beta"]   # report δ_p = p_oc/β0 (Δp/p0)
        meta = {"section": section, "sc_enabled": sc_enabled, "switches": sw,
                "config_sha": config_sha(cfg), "provenance": _provenance(cfg),
                "longitudinal_native_coordinate": "p_oc = dE/(c*p0)",
                "reported_delta": "delta_p = dp/p0",
                "conversion_beta": float(d["beta"]),
                "rf_kicks_applied": rf_kicks_applied,
                "sc_requested": False,
                "sc_effective": False,
                "rf": "thin-lens (K·sin) + transverse kick K_trans"
                if sw["rf_transverse_kick"] else "thin-lens (K·sin) only"}

    return BeamResult(
        route="OCELOT",
        z_mm=z_arr * 1e3,
        sigma_x_um=sx * 1e6,
        sigma_y_um=sy * 1e6,
        sigma_z_um=sz * 1e6,
        eps_nx_mm_mrad=enx * 1e6,
        eps_ny_mm_mrad=eny * 1e6,
        energy_keV=np.full_like(z_arr, P.beam.energy_keV),
        sigma_delta_e3=sd * 1e3,
        meta=meta,
    )
