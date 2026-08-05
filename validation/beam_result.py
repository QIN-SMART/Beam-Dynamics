"""
Unified result container shared by every validation test.

Both backends (AG, OCELOT) are reduced to this identical structure so the
framework can compare them field-by-field and plot AG vs OCELOT overlays.

Fields (arrays are sampled on the same z grid):
  z_mm             [mm]
  sigma_x_um       [μm]  projected RMS
  sigma_y_um       [μm]
  sigma_z_um       [μm]  longitudinal RMS (physical, z = β·c·t)
  eps_nx_mm_mrad   [mm·mrad]  normalized transverse emittance (x)
  eps_ny_mm_mrad   [mm·mrad]  normalized transverse emittance (y)
  energy_keV       [keV]  reference-particle kinetic energy
  sigma_delta_e3   [×10⁻³]
  chirp_m          [m⁻¹]  dδ/dz of the correlated energy spread
  time_res_fs      [fs]   RMS time spread σ_t = σ_z / (βc)
"""

import os
import json

import numpy as np

C_SI = 2.99792458e8


class BeamResult:
    def __init__(self, route, z_mm, sigma_x_um, sigma_y_um, sigma_z_um,
                 eps_nx_mm_mrad, eps_ny_mm_mrad,
                 energy_keV, sigma_delta_e3, chirp_m=None,
                 time_res_fs=None, meta=None):
        self.route = route
        self.z_mm = np.asarray(z_mm, dtype=float)
        self.sigma_x_um = np.asarray(sigma_x_um, dtype=float)
        self.sigma_y_um = np.asarray(sigma_y_um, dtype=float)
        self.sigma_z_um = np.asarray(sigma_z_um, dtype=float)
        self.eps_nx_mm_mrad = np.asarray(eps_nx_mm_mrad, dtype=float)
        self.eps_ny_mm_mrad = np.asarray(eps_ny_mm_mrad, dtype=float)
        self.energy_keV = np.asarray(energy_keV, dtype=float)
        self.sigma_delta_e3 = np.asarray(sigma_delta_e3, dtype=float)
        self.chirp_m = chirp_m                      # scalar [m⁻¹]
        self.time_res_fs = (time_res_fs if time_res_fs is not None
                            else self._time_from_sigma_z())
        self.meta = meta or {}

    def _beta(self):
        # beta from energy: gamma = 1 + E[keV]/511, beta = sqrt(1-1/gamma^2)
        gamma = 1.0 + self.energy_keV / 511.0
        return np.sqrt(np.maximum(1.0 - 1.0 / gamma**2, 0.0))

    def _time_from_sigma_z(self):
        beta = np.mean(self._beta())
        if beta < 1e-6:
            return np.full_like(self.sigma_z_um, np.nan)
        return self.sigma_z_um * 1e-6 / (beta * C_SI) * 1e15   # fs

    # ── comparison helpers ────────────────────────────────────────────────
    def at(self, field, z_mm):
        """Linear interpolation of a field onto query z positions."""
        return np.interp(z_mm, self.z_mm, getattr(self, field))

    def metrics_vs(self, other, fields=None):
        """Max-relative deviation vs another BeamResult on the common z range."""
        zlo = max(self.z_mm[0], other.z_mm[0])
        zhi = min(self.z_mm[-1], other.z_mm[-1])
        if zhi <= zlo:
            return {}
        zq = np.linspace(zlo, zhi, 400)
        fields = fields or ["sigma_x_um", "sigma_y_um", "sigma_z_um",
                            "eps_nx_mm_mrad", "sigma_delta_e3"]
        out = {}
        for f in fields:
            a = np.interp(zq, self.z_mm, getattr(self, f))
            b = np.interp(zq, other.z_mm, getattr(other, f))
            scale = np.maximum(np.abs(b), 1e-12)
            out[f] = float(np.max(np.abs(a - b) / scale) * 100.0)
        return out

    def summary(self):
        i = len(self.z_mm) - 1
        return (f"{self.route}: z_end={self.z_mm[-1]:.0f}mm  "
                f"σ_x={self.sigma_x_um[-1]:.1f}  σ_y={self.sigma_y_um[-1]:.1f}  "
                f"σ_z={self.sigma_z_um[-1]:.1f} μm  "
                f"ε_nx={self.eps_nx_mm_mrad[-1]:.4f} mm·mrad  "
                f"σ_δ={self.sigma_delta_e3[-1]:.3f}e-3  "
                f"Δt={self.time_res_fs[-1]:.1f} fs")

    # ── IO ────────────────────────────────────────────────────────────────
    def to_dict(self):
        return {
            "route": self.route,
            "z_mm": self.z_mm.tolist(),
            "sigma_x_um": self.sigma_x_um.tolist(),
            "sigma_y_um": self.sigma_y_um.tolist(),
            "sigma_z_um": self.sigma_z_um.tolist(),
            "eps_nx_mm_mrad": self.eps_nx_mm_mrad.tolist(),
            "eps_ny_mm_mrad": self.eps_ny_mm_mrad.tolist(),
            "energy_keV": self.energy_keV.tolist(),
            "sigma_delta_e3": self.sigma_delta_e3.tolist(),
            "chirp_m": self.chirp_m,
            "time_res_fs": self.time_res_fs.tolist(),
            "meta": self.meta,
        }

    def save(self, path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)
        return path
