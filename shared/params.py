"""
Shared parameter loader — single source of truth: shared/beamline_config.yaml.

Both independent routes read all beamline / initial-beam parameters through
this module (pure-python + yaml, no numpy/ocelot/matplotlib dependency).

Usage:
    from shared.params import load_config, parse, config_sha, derived

    cfg  = load_config()          # raw dict (same keys as before)
    d    = derived(cfg)           # relativistic + derived quantities (SI)
    P    = parse(cfg)             # typed dataclass bundle
"""

import os
import hashlib

import yaml

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CONFIG = os.path.join(_THIS_DIR, "beamline_config.yaml")

# ── Physical constants (SI) ────────────────────────────────────────────────
C_SI       = 2.99792458e8
M_E_SI     = 9.10938356e-31
E_SI       = 1.602176634e-19
EPSILON_0  = 8.854187817e-12
MEC2_KEV   = 511.0              # m_e c² in keV


def config_path(path: str = None) -> str:
    return path or DEFAULT_CONFIG


def load_config(path: str = None) -> dict:
    """Load the shared YAML config (raw dict, keys identical to the old
    GPT模拟/beamline_config.yaml)."""
    with open(config_path(path), "r") as f:
        return yaml.safe_load(f)


def config_sha(cfg: dict) -> str:
    """Short hash of the config, recorded in result files to guarantee that
    both routes were run on the SAME parameters."""
    return hashlib.sha1(str(cfg).encode("utf-8")).hexdigest()[:12]


# ── Derived (SI) quantities from the shared config ────────────────────────

def derived(cfg: dict) -> dict:
    """Relativistic kinematics + beam-derived quantities, all SI."""
    b  = cfg["beam"]
    ib = cfg["initial_distribution"]

    E_keV    = b["energy_keV"]
    gamma    = 1.0 + E_keV / MEC2_KEV
    beta     = (1.0 - 1.0 / gamma**2) ** 0.5
    beta_gamma = beta * gamma
    p_SI     = gamma * M_E_SI * beta * C_SI
    v_e      = beta * C_SI

    eps_n     = ib["epsilon_n_mm_mrad"] * 1e-6       # [m·rad]
    spot      = ib["spot_rms_um"] * 1e-6             # [m]
    eps_geom  = eps_n / beta_gamma                   # geometric [m·rad]

    return {
        "E_keV": E_keV,
        "gamma": gamma,
        "beta": beta,
        "beta_gamma": beta_gamma,
        "p_SI": p_SI,
        "v_e": v_e,
        "eps_n_m": eps_n,
        "eps_geom_m": eps_geom,
        "sigma_xp": eps_geom / spot,                 # [rad]
        "sigma_yp": eps_geom / spot,                 # [rad]
    }


# ── Lattice (single source of truth for geometry + element parameters) ────

def _lattice_elements(cfg: dict) -> list:
    """Return list of element dicts {name, type, z_start, length, parameters}."""
    return list(cfg["lattice"]["elements"])


def elements_of_type(cfg: dict, etype: str) -> list:
    """All lattice elements of a given type (multi-instance safe)."""
    return [e for e in _lattice_elements(cfg) if e["type"] == etype]


def first_of_type(cfg: dict, etype: str) -> dict:
    """First element of a type, or None."""
    els = elements_of_type(cfg, etype)
    return els[0] if els else None


def elem_geom(e: dict):
    """Return (z_start, z_end, length) of an element."""
    return e["z_start"], e["z_start"] + e["length"], e["length"]


def elem_params(e: dict) -> dict:
    """Flat parameter dict of an element (parameters + z_start/length)."""
    p = dict(e.get("parameters") or {})
    p["z_start_m"] = e["z_start"]
    p["length_m"] = e["length"]
    return p


def flat_elem(cfg: dict, etype: str) -> dict:
    """Flat dict of the FIRST element of a type (backward-compatible view:
    keys B_field_T / frequency_GHz / voltage_kV / phase_rad / z_start_m /
    length_m as in the pre-refactor config sections).  None if absent."""
    e = first_of_type(cfg, etype)
    return elem_params(e) if e is not None else None


def z_sample(cfg: dict) -> float:
    """Sample-plane longitudinal position [m] (last lattice element end)."""
    return max(e["z_start"] + e["length"] for e in _lattice_elements(cfg))


def lattice_active(cfg: dict) -> list:
    """Active (non-drift) elements with resolved (start, end) spans [m].
    Returns list of (z_start, z_end, etype)."""
    out = []
    for e in _lattice_elements(cfg):
        if e["type"] in ("solenoid", "rf_cavity"):
            z0, z1, _ = elem_geom(e)
            out.append((z0, z1, e["type"]))
    return out


# ── Typed parameter bundle ────────────────────────────────────────────────

class BeamParams:
    """Initial beam / bunch parameters (shared keys, converted to SI where
    the raw value is not already SI)."""
    def __init__(self, cfg):
        b  = cfg["beam"]
        ib = cfg["initial_distribution"]
        self.energy_keV         = b["energy_keV"]
        self.charge_fC          = b["charge_fC"]
        self.n_particles        = b["n_particles"]
        self.Q_C                = b["charge_fC"] * 1e-15
        self.spot_rms_um        = ib["spot_rms_um"]
        self.bunch_length_um    = ib["bunch_length_um"]
        self.epsilon_n_mm_mrad  = ib["epsilon_n_mm_mrad"]
        self.epsilon_nz_mm_mrad = ib["epsilon_nz_mm_mrad"]
        self.sigma_delta        = ib["sigma_delta"]
        # SI
        self.spot_rms           = ib["spot_rms_um"] * 1e-6
        self.sig_z0             = ib["bunch_length_um"] * 1e-6
        self.eps_n              = ib["epsilon_n_mm_mrad"] * 1e-6
        self.eps_nz             = ib["epsilon_nz_mm_mrad"] * 1e-6


class SolenoidParams:
    """First solenoid element (geometry + parameters from lattice)."""

    def __init__(self, cfg, elem=None):
        s = elem_params(elem if elem is not None else first_of_type(cfg, "solenoid"))
        self.B_field_T = s["B_field_T"]
        self.length_m  = s["length_m"]
        self.z_start_m = s["z_start_m"]


class RFParams:
    """First RF cavity element (geometry + parameters from lattice)."""

    def __init__(self, cfg, elem=None):
        r = elem_params(elem if elem is not None else first_of_type(cfg, "rf_cavity"))
        self.frequency_GHz = r["frequency_GHz"]
        self.voltage_kV    = r["voltage_kV"]
        self.phase_rad     = r["phase_rad"]
        self.length_m      = r["length_m"]
        self.z_start_m     = r["z_start_m"]
        # SI / derived
        self.f_RF          = r["frequency_GHz"] * 1e9
        self.V_RF          = r["voltage_kV"] * 1e3
        self.k_rf          = 2.0 * 3.141592653589793 * self.f_RF / C_SI
        self.E_rf          = self.V_RF / self.length_m   # peak field [V/m]


class SCParams:
    def __init__(self, cfg):
        s = cfg["space_charge"]
        self.enabled = bool(s["enabled"])
        self.mesh    = list(s["mesh"])
        self.step    = s["step"]


class PhysicsSwitchParams:
    """Shared physics switches (read identically by both backends)."""

    def __init__(self, cfg):
        s = cfg.get("physics_switches", {})
        self.rf_longitudinal_kick = bool(s.get("rf_longitudinal_kick", True))
        self.rf_transverse_kick   = bool(s.get("rf_transverse_kick", False))

    def as_dict(self):
        return {"rf_longitudinal_kick": self.rf_longitudinal_kick,
                "rf_transverse_kick": self.rf_transverse_kick}


class OutputParams:
    def __init__(self, cfg):
        o = cfg["output"]
        self.z_diagnostics_mm = list(o["z_diagnostics_mm"])
        self.step_size_m      = o["step_size_m"]


class Config:
    """All shared parameters grouped for convenient access."""
    def __init__(self, cfg):
        self.raw      = cfg
        self.beam     = BeamParams(cfg)
        self.sc       = SCParams(cfg)
        self.switches = PhysicsSwitchParams(cfg)
        self.output   = OutputParams(cfg)
        self.elements = _lattice_elements(cfg)
        # first-instance convenience (multi-instance via elements_of_type)
        self.solenoid = SolenoidParams(cfg) if first_of_type(cfg, "solenoid") else None
        self.rf       = RFParams(cfg) if first_of_type(cfg, "rf_cavity") else None


def parse(cfg: dict = None) -> Config:
    return Config(cfg if cfg is not None else load_config())
