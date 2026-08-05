"""
Unified output schema — identical for both routes.

Each route writes one JSON file (shared/results/<route>_results.json) with the
same keys and units so shared/compare.py can produce overlay plots and a
comparison table.

File structure:
{
  "route": "AG" | "GPT",
  "config_sha": str,          # params.config_sha() of the config used
  "sc_enabled": bool,
  "units": {...},             # human-readable unit record
  "meta": {...},              # optional route-specific notes
  "probes": [ {...}, ... ],   # one record per z_diagnostics_mm
  "history": { "z_mm": [...], "sigma_x_um": [...], ... }   # full evolution
}

Probe / history record keys (units fixed):
  z_mm             [mm]
  sigma_x_um       [μm]
  sigma_y_um       [μm]
  sigma_z_um       [μm]
  sigma_delta_e3   [×10⁻³]  (σ_δ of the MOMENTUM deviation δ_p = Δp/p₀;
                             OCELOT raw p = ΔE/(c·p₀) is divided by β₀ before
                             reporting — see R56_convention_resolution.md)
  eps_nx_mm_mrad   [mm·mrad]  normalized emittance
  eps_ny_mm_mrad   [mm·mrad]  normalized emittance
"""

import os
import json

_RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")

UNITS = {
    "z_mm": "mm",
    "sigma_x_um": "um",
    "sigma_y_um": "um",
    "sigma_z_um": "um",
    "sigma_delta_e3": "e-3",
    "eps_nx_mm_mrad": "mm.mrad",
    "eps_ny_mm_mrad": "mm.mrad",
}


def make_probe(z_mm, sigma_x_um, sigma_y_um, sigma_z_um,
               sigma_delta_e3, eps_nx_mm_mrad, eps_ny_mm_mrad, **extra):
    """Build a standard probe record dict (extra fields allowed)."""
    rec = {
        "z_mm": float(z_mm),
        "sigma_x_um": float(sigma_x_um),
        "sigma_y_um": float(sigma_y_um),
        "sigma_z_um": float(sigma_z_um),
        "sigma_delta_e3": float(sigma_delta_e3),
        "eps_nx_mm_mrad": float(eps_nx_mm_mrad),
        "eps_ny_mm_mrad": float(eps_ny_mm_mrad),
    }
    rec.update(extra)
    return rec


def write_results(route, probes, history, config_sha, sc_enabled,
                  out_dir=None, meta=None):
    """Write <route>_results.json to shared/results/ (default). Returns path."""
    out_dir = out_dir or _RESULTS_DIR
    os.makedirs(out_dir, exist_ok=True)

    payload = {
        "route": route,
        "config_sha": config_sha,
        "sc_enabled": bool(sc_enabled),
        "units": UNITS,
        "probes": probes,
        "history": history,
    }
    if meta:
        payload["meta"] = meta

    path = os.path.join(out_dir, f"{route}_results.json")
    with open(path, "w") as f:
        json.dump(payload, f, indent=2)
    return path


def load_results(path):
    """Load a result file written by write_results()."""
    with open(path, "r") as f:
        return json.load(f)


def list_results(out_dir=None):
    """Return sorted list of *_results.json paths in the results dir."""
    out_dir = out_dir or _RESULTS_DIR
    if not os.path.isdir(out_dir):
        return []
    return sorted(
        os.path.join(out_dir, p) for p in os.listdir(out_dir)
        if p.endswith("_results.json")
    )
