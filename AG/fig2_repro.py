#!/usr/bin/env python3
"""
Reproduce Kelisani 2023, Fig. 2:
  (a) bunch size, (b) beam divergence along a drift section
  due to only space-charge effect.
  Solid lines = beam-envelope equations.

Parameters (confirmed from paper text, p.7):
  5 MeV beam kinetic energy, 1 nC bunch charge,
  σ_x = √2 mm, σ_y = 2√2 mm, σ_z = 0.030 mm,
  ε_nx = 0.05 μm, ε_ny = 0.10 μm, ε_nz = 2.93 μm,
  energy spread 1%.

Paper layout (Fig 4, which uses same parameters as Fig 2):
  (a) transvers σ_x,y (mm) [left axis] + longitudinal σ_z (mm) [right axis]
  (b) transverse Λ_x,y (mrad) [left] + longitudinal Λ_z (mrad) [right]
  X axis: Z (cm)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from beam_dynamics_6d import (Beam6D, propagate, get_alpha_interpolators,
                               beamK_to_gamma)

# ═════════════════════════════════════════════════════════════════════════════
# Paper parameters  (Kelisani 2023, Fig. 2 / Fig. 4)
# ═════════════════════════════════════════════════════════════════════════════
BEAMK_EV    = 5.0e6                                    # 5 MeV
NE          = 1e-9 / 1.60217662e-19                    # 1 nC → ~6.24e9 e⁻
SIGMA_X0    = np.sqrt(2) * 1e-3                        # √2 mm
SIGMA_Y0    = 2.0 * np.sqrt(2) * 1e-3                  # 2√2 mm
SIGMA_Z0    = 30e-6                                    # 30 μm
EPS_NX      = 0.05e-6                                  # 0.05 μm
EPS_NY      = 0.10e-6                                  # 0.10 μm
EPS_NZ      = 2.93e-6                                  # 2.93 μm
SIGMA_DELTA = 0.01                                     # 1% energy spread
Z_DRIFT     = 0.10                                     # drift length

# ═════════════════════════════════════════════════════════════════════════════
# Build initial beam
# ═════════════════════════════════════════════════════════════════════════════
gamma0 = beamK_to_gamma(BEAMK_EV)
beta0  = np.sqrt(1.0 - 1.0 / gamma0**2)
p0     = gamma0 * beta0

beam0 = Beam6D(
    sigma_x=SIGMA_X0, sigma_y=SIGMA_Y0, sigma_z=SIGMA_Z0,
    sigma_xy=0.0, sigma_delta=SIGMA_DELTA, C_zd=0.0,
    nu_x=0.0, nu_y=0.0, nu_z=0.0,
    nu_xy=0.0, nu_delta=0.0,
    gamma=gamma0, Ne=NE,
    eps_nx=EPS_NX, eps_ny=EPS_NY, eps_nz=EPS_NZ,
)
print(beam0.summary("Initial beam (Kelisani 2023, Fig.2)"))
print(f"  χ_τ ≈ {NE * 2.818e-15 * SIGMA_X0 / (EPS_NX * EPS_NZ):.1e} (space-charge parameter)")

# ═════════════════════════════════════════════════════════════════════════════
# Propagate: pure drift with space charge  (Gaussian 9-coeff model, fig.2)
# ═════════════════════════════════════════════════════════════════════════════
_ = get_alpha_interpolators()

z_arr, st = propagate(beam0, (0.0, Z_DRIFT), n_points=400,
                      sc_model='gaussian', rtol=1e-10, atol=1e-12)

# Also propagate with ellipsoid model for Fig. 4 reference
z_arr4, st4 = propagate(beam0, (0.0, Z_DRIFT), n_points=400,
                        sc_model='ellipsoid', rtol=1e-10, atol=1e-12)

# Extract state arrays
def extract_state(s, eps_nx, eps_ny, eps_nz, p0_const):
    sx, sy, sz = s[:, 0], s[:, 1], s[:, 2]
    nx, ny, nz = s[:, 6], s[:, 7], s[:, 8]
    egx, egy, egz = eps_nx / p0_const, eps_ny / p0_const, eps_nz / p0_const
    lx = np.sqrt(egx**2 / np.maximum(sx, 1e-20)**2 + nx**2)
    ly = np.sqrt(egy**2 / np.maximum(sy, 1e-20)**2 + ny**2)
    lz = np.sqrt(egz**2 / np.maximum(sz, 1e-20)**2 + nz**2)
    return sx, sy, sz, lx, ly, lz

sigma_x, sigma_y, sigma_z, Lambda_x, Lambda_y, Lambda_z = \
    extract_state(st, EPS_NX, EPS_NY, EPS_NZ, p0)

sigma_x4, sigma_y4, sigma_z4, Lambda_x4, Lambda_y4, Lambda_z4 = \
    extract_state(st4, EPS_NX, EPS_NY, EPS_NZ, p0)

# ═════════════════════════════════════════════════════════════════════════════
# Plot — matching the paper's multi-axis layout
# ═════════════════════════════════════════════════════════════════════════════
z_cm = z_arr * 1e2

fig, ((ax_sz_xy, ax_sz_z), (ax_div_xy, ax_div_z)) = plt.subplots(
    2, 2, figsize=(13, 9))

# ── (a1) Transverse bunch size  σ_x, σ_y ──
# Dual y-axes: left for σ_x, right for σ_y
color_x, color_y = '#1f77b4', '#d62728'
ax_sz_xy.plot(z_cm, sigma_x * 1e3, '-', color=color_x, lw=2.0, label=r'$\sigma_x$ (Gauss)')
ax_sz_xy.plot(z_cm, sigma_x4 * 1e3, ':', color=color_x, lw=1.2, alpha=0.5, label=r'$\sigma_x$ (Ellip)')
ax_sz_xy.set_xlabel('Z (cm)', fontsize=12)
ax_sz_xy.set_ylabel(r'$\sigma_x$  (mm)', fontsize=12, color=color_x)
ax_sz_xy.tick_params(axis='y', labelcolor=color_x)
ax_sz_xy.grid(True, alpha=0.3)

ax_sz_xy2 = ax_sz_xy.twinx()
ax_sz_xy2.plot(z_cm, sigma_y * 1e3, '--', color=color_y, lw=2.0, label=r'$\sigma_y$ (Gauss)')
ax_sz_xy2.plot(z_cm, sigma_y4 * 1e3, ':', color=color_y, lw=1.2, alpha=0.5, label=r'$\sigma_y$ (Ellip)')
ax_sz_xy2.set_ylabel(r'$\sigma_y$  (mm)', fontsize=12, color=color_y)
ax_sz_xy2.tick_params(axis='y', labelcolor=color_y)

lines1, labels1 = ax_sz_xy.get_legend_handles_labels()
lines2, labels2 = ax_sz_xy2.get_legend_handles_labels()
ax_sz_xy.legend(lines1 + lines2, labels1 + labels2, fontsize=7, loc='upper left', ncol=2)
ax_sz_xy.set_title('(a) Transverse bunch size', fontsize=13, fontweight='bold')
ax_sz_xy.set_xlim(0, Z_DRIFT * 1e2)

# ── (a2) Longitudinal bunch size  σ_z ──
ax_sz_z.plot(z_cm, sigma_z * 1e3, 'g-', lw=1.8, label=r'$\sigma_z$ (Gauss)')
ax_sz_z.plot(z_cm, sigma_z4 * 1e3, 'g:', lw=1.2, alpha=0.6, label=r'$\sigma_z$ (Ellip)')
ax_sz_z.set_ylabel(r'$\sigma_z$  (mm)', fontsize=12)
ax_sz_z.legend(fontsize=9, loc='upper left')
ax_sz_z.grid(True, alpha=0.3)
ax_sz_z.set_title('(a) Longitudinal bunch size', fontsize=13, fontweight='bold')

# ── (b1) Transverse divergence  Λ_x, Λ_y ──
ax_div_xy.plot(z_cm, Lambda_x * 1e3, 'b-', lw=1.8, label=r'$\Lambda_x$ (Gauss)')
ax_div_xy.plot(z_cm, Lambda_y * 1e3, 'r--', lw=1.8, label=r'$\Lambda_y$ (Gauss)')
ax_div_xy.plot(z_cm, Lambda_x4 * 1e3, 'b:', lw=1.2, alpha=0.6, label=r'$\Lambda_x$ (Ellip)')
ax_div_xy.plot(z_cm, Lambda_y4 * 1e3, 'r:', lw=1.2, alpha=0.6, label=r'$\Lambda_y$ (Ellip)')
ax_div_xy.set_xlabel('Z (cm)', fontsize=12)
ax_div_xy.set_ylabel(r'$\Lambda_{x,y}$  (mrad)', fontsize=12)
ax_div_xy.legend(fontsize=8, loc='upper left', ncol=2)
ax_div_xy.grid(True, alpha=0.3)
ax_div_xy.set_title('(b) Transverse divergence', fontsize=13, fontweight='bold')

# ── (b2) Longitudinal divergence  Λ_z ──
ax_div_z.plot(z_cm, Lambda_z * 1e3, 'm-', lw=1.8, label=r'$\Lambda_z$ (Gauss)')
ax_div_z.plot(z_cm, Lambda_z4 * 1e3, 'm:', lw=1.2, alpha=0.6, label=r'$\Lambda_z$ (Ellip)')
ax_div_z.set_xlabel('Z (cm)', fontsize=12)
ax_div_z.set_ylabel(r'$\Lambda_z$  (mrad)', fontsize=12)
ax_div_z.legend(fontsize=9, loc='upper left')
ax_div_z.grid(True, alpha=0.3)
ax_div_z.set_title('(b) Longitudinal divergence', fontsize=13, fontweight='bold')

fig.suptitle(
    'Kelisani 2023 Fig. 2 — Space-charge driven evolution in 10 cm drift\n'
    r'5 MeV, 1 nC,  $\sigma_x{=}\sqrt{2}$ mm, $\sigma_y{=}2\sqrt{2}$ mm, '
    r'$\sigma_z{=}30\ \mu$m, $\varepsilon_{nx}{=}0.05\ \mu$m, '
    r'$\varepsilon_{ny}{=}0.10\ \mu$m, $\varepsilon_{nz}{=}2.93\ \mu$m',
    fontsize=11, fontweight='bold')
plt.tight_layout()
plt.savefig("kelisani2023_fig2_repro.png", dpi=150)
print("\nSaved: kelisani2023_fig2_repro.png")

# Print final values
b_end = Beam6D.from_state(st[-1, :11], st[-1, 11], NE, EPS_NX, EPS_NY, EPS_NZ)
print(b_end.summary("z=10 cm (Gaussian SC model)"))
b_end4 = Beam6D.from_state(st4[-1, :11], st4[-1, 11], NE, EPS_NX, EPS_NY, EPS_NZ)
print(b_end4.summary("z=10 cm (Ellipsoid SC model)"))
