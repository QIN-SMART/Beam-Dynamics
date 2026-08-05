#!/usr/bin/env python3
"""
External Field Forces for 6D Beam Envelope Equations
=====================================================
Based on Kelisani 2023, Phys. Rev. Applied 19, 054011, Section IV.

Provides F-type force functions for the 6D envelope ODE:
  σ_u'' + (γ γ')/(γ² β²) σ_u' = F_u^e + F_u^s + F_u^ε

Implemented force types
-----------------------
  Drift        — free space, zero forces
  Solenoid     — Larmor focusing (k_s = eB/2p) + x-y coupling
  RF Cavity    — longitudinal acceleration γ'(z) + RF focusing δ' = −k_RF z

All implementations use linear / hard-edge approximations.
No complex field distributions are assumed.
"""

import numpy as np
from typing import Tuple, Optional, Callable, List
from dataclasses import dataclass

# ── Physical constants (SI) ──
# Duplicated here so the module is self-contained;
# values match beam_dynamics_6d.py exactly.
M_E = 9.10938356e-31             # electron mass  [kg]
M_E_EV = 5.109989461e5           # electron rest energy  [eV]
C_LIGHT = 2.99792458e8           # speed of light  [m/s]
E_CHARGE = 1.60217662e-19        # elementary charge  [C]
EPSILON_0 = 8.854187817e-12      # vacuum permittivity  [F/m]
R_E = 2.8179403262e-15           # classical electron radius  [m]

# η = e / (m_e c²)  — charge-to-mass-energy ratio  [1 / V·m]
ETA = E_CHARGE / (M_E * C_LIGHT**2)

# ── Import Beam6D from the main module ──
from beam_dynamics_6d import Beam6D, ExtFieldRegion, gamma_to_beta


# ═════════════════════════════════════════════════════════════════════════════
# 1.  Drift  —  Fe_x = Fe_y = Fe_z = 0
# ═════════════════════════════════════════════════════════════════════════════

def drift_force(beam: Beam6D
                ) -> Tuple[float, float, float, float, float, float]:
    """
    Drift space: no external forces, no coupling.

    Returns (Fe_x, Fe_y, Fe_z, dσ_xy/dz, dν_x_sol, dν_y_sol).
    All zero for a field-free drift.
    """
    return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0


# ═════════════════════════════════════════════════════════════════════════════
# 2.  Solenoid  —  Larmor focusing  k_s = e B / (2 p)  +  x-y coupling
# ═════════════════════════════════════════════════════════════════════════════

def solenoid_force(beam: Beam6D, Bz: float, dBz_dz: float = 0.0
                   ) -> Tuple[float, float, float, float, float, float]:
    """
    Solenoid force with Larmor focusing and x-y coupling.

    Based on Kelisani 2023 Eqs. (31)-(33), Section IV.A.

    **Larmor wavenumber**
        k_s = e·Bz / (2 p) = η·c·Bz / (2 p₀)          p₀ = γ β

    **Envelope equation**   (linear / hard-edge approximation)
        σ_x'' + (γ γ')/(γ² β²) σ_x' + k_s² σ_x = F^s_x + F^ε_x
        σ_y'' + (γ γ')/(γ² β²) σ_y' + k_s² σ_y = F^s_y + F^ε_y

        For this module the focusing term is returned as external force:
        F^sol_x = −k_s² σ_x,   F^sol_y = −k_s² σ_y

    **x-y coupling**  (Larmor rotation)
        d⟨xy⟩/dz  = 2 k_s (σ_x² − σ_y²)
        dν_x/dz  ⊃  −2 k_s ν_y
        dν_y/dz  ⊃  +2 k_s ν_x

    Parameters
    ----------
    beam : Beam6D
        Current beam state.
    Bz : float [T]
        On-axis longitudinal magnetic field  μ.
    dBz_dz : float [T/m]
        Axial derivative ∂μ/∂z.  Use 0 for hard-edge (linear) approximation.

    Returns
    -------
    (Fe_x, Fe_y, Fe_z, dσ_xy/dz, dν_x_sol, dν_y_sol)
        Fe_x, Fe_y   — solenoid focusing forces
        Fe_z         — always 0 (solenoid does not focus longitudinally)
        dσ_xy/dz     — x-y correlation growth from Larmor rotation
        dν_x_sol     — Larmor contribution to  dν_x/dz  (couples ν_y)
        dν_y_sol     — Larmor contribution to  dν_y/dz  (couples ν_x)
    """
    p0 = beam.p0
    if p0 < 1e-15:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

    # Larmor wavenumber  k_s = η c Bz / (2 p₀)
    k_s = ETA * C_LIGHT * Bz / (2.0 * p0)

    # Focusing strength includes fringe contribution  (k_s² + k_fringe²)
    k_s_fringe = ETA * C_LIGHT * dBz_dz * beam.sigma_z / (2.0 * p0)
    k_s2 = k_s**2 + k_s_fringe**2

    Fe_x = -k_s2 * beam.sigma_x
    Fe_y = -k_s2 * beam.sigma_y
    Fe_z = 0.0

    # x-y correlation growth:  d⟨xy⟩/dz = 2 k_s (σ_x² − σ_y²)
    d_sigma_xy = 2.0 * k_s * (beam.sigma_x**2 - beam.sigma_y**2)

    # Larmor coupling to ν-derivatives
    dnx_sol = -2.0 * k_s * beam.nu_y
    dny_sol = +2.0 * k_s * beam.nu_x

    return Fe_x, Fe_y, Fe_z, d_sigma_xy, dnx_sol, dny_sol


def solenoid_larmor_wavenumber(beam: Beam6D, Bz: float) -> float:
    """Return the Larmor wavenumber  k_s = e Bz / (2 p)  [1/m]."""
    p0 = beam.p0
    if p0 < 1e-15:
        return 0.0
    return ETA * C_LIGHT * Bz / (2.0 * p0)


# ═════════════════════════════════════════════════════════════════════════════
# 3.  RF Cavity  —  longitudinal acceleration + RF focusing
# ═════════════════════════════════════════════════════════════════════════════

def rf_cavity_force(beam: Beam6D, E_rf: float, dE_dz: float = 0.0,
                    dE_cdt: float = 0.0
                    ) -> Tuple[float, float, float]:
    """
    RF cavity F-type forces (transverse + longitudinal).

    Based on Kelisani 2023 Eqs. (53)-(55), Section IV.D.

    **Transverse** — RF defocusing from combined E + v×B
        F^rf_{x,y} = −η/(2γβ²) · (∂E/∂z + β·∂E/c∂t) · σ_{x,y}

        For a standing-wave cavity:  ∂E/c∂t ≈ k_rf·E_rf
        For a traveling-wave cavity: ∂E/c∂t ≈ 0

    **Longitudinal** — focusing from axial field gradient
        F^rf_z = +η/(γ p₀²) · (∂E/∂z) · σ_z

    Parameters
    ----------
    beam : Beam6D
    E_rf : float [V/m]
        On-axis peak RF electric field.
    dE_dz : float [V/m²]
        Axial derivative ∂E_rf/∂z.  ≈ 0 for TW, ≠ 0 for SW.
    dE_cdt : float [V/m²]
        Time derivative (1/c)·∂E_rf/∂t.
        For standing wave on crest: dE_cdt ≈ k_rf·E_rf.

    Returns (Fe_x, Fe_y, Fe_z).
    """
    p0 = beam.p0
    gamma = beam.gamma
    beta = beam.beta

    gamma_beta2 = gamma * beta**2
    if abs(gamma_beta2) < 1e-15:
        return 0.0, 0.0, 0.0

    Erf_trans = dE_dz + beta * dE_cdt
    factor_trans = -ETA * Erf_trans / (2.0 * gamma_beta2)
    Fe_x = factor_trans * beam.sigma_x
    Fe_y = factor_trans * beam.sigma_y

    Fe_z = +ETA * dE_dz / (gamma * p0**2) * beam.sigma_z

    return Fe_x, Fe_y, Fe_z


def rf_acceleration_gradient(E_rf: float) -> float:
    """
    Energy gain rate from RF accelerating field.

        dγ/dz = η · E_rf = e · E_rf / (m_e · c²)

    Parameters
    ----------
    E_rf : float [V/m]
        On-axis accelerating field (signed for electron beam).

    Returns dγ/dz [1/m].
    """
    return ETA * E_rf


def rf_chirp_rate(beam: Beam6D, E_rf: float, k_rf: float) -> float:
    """
    RF longitudinal focusing chirp rate  h = dδ/dz.

    For an RF cavity the correlated energy modulation is:

        δ'  ≈  h · z                                  (within cavity)

        h = − e · E_rf · k_rf / (β² γ³ · m_e · c²)

    where  k_rf = 2π · f_rf / c  is the RF wavenumber.

    **Sign convention**
        h > 0  →  head (z > 0) gains more energy  →  positive chirp
                  →  bunch compresses in subsequent drift.
        h < 0  →  tail gains energy  →  decompression.

    Parameters
    ----------
    beam : Beam6D
    E_rf : float [V/m]
        Peak on-axis RF field.
    k_rf : float [1/m]
        RF wavenumber  2π·f_rf / c.

    Returns h [1/m].
    """
    beta2 = max(beam.beta**2, 1e-10)
    gamma3 = beam.gamma**3
    return -E_CHARGE * E_rf * k_rf / (beta2 * gamma3 * M_E * C_LIGHT**2)


def rf_compute_chirp_coefficient(beam: Beam6D, E_rf: float, L_cav: float,
                                 f_rf: float = 2.856e9) -> float:
    """
    Total RF chirp coefficient H over a cavity of length L_cav.

        H = − e · E_rf · k_rf · L_cav / (β² γ³ · m_e · c²)

    After passing through the cavity (thin-lens approximation):
        δ_out = δ_in + H · z

    Parameters
    ----------
    beam : Beam6D
    E_rf : float [V/m]
    L_cav : float [m]
    f_rf : float [Hz]

    Returns H [1/m].
    """
    k_rf = 2.0 * np.pi * f_rf / C_LIGHT
    return rf_chirp_rate(beam, E_rf, k_rf) * L_cav


# ═════════════════════════════════════════════════════════════════════════════
# Composite function builders for beamline simulation
# ═════════════════════════════════════════════════════════════════════════════

def build_external_force_func(
    regions: List[ExtFieldRegion],
) -> Callable[[float, Beam6D], Tuple[float, float, float, float, float, float]]:
    """
    Build a composite external-force function from a list of field regions.

    Each region is active on  z ∈ [z_start, z_end]  and specifies its
    field type ('solenoid', 'rf', 'drift') and parameters.

    Parameters
    ----------
    regions : list of ExtFieldRegion

    Returns
    -------
    force_func(z, beam) → (Fe_x, Fe_y, Fe_z, dσ_xy/dz, dν_x_sol, dν_y_sol)
        All six contributions, summed over active regions.
    """
    def force_func(z: float, beam: Beam6D
                   ) -> Tuple[float, float, float, float, float, float]:
        fx, fy, fz = 0.0, 0.0, 0.0
        dsxy, dnx, dny = 0.0, 0.0, 0.0
        for reg in regions:
            if reg.z_start <= z <= reg.z_end:
                p = reg.params
                if reg.ftype == 'solenoid':
                    fx1, fy1, fz1, dsxy1, dnx1, dny1 = solenoid_force(
                        beam, p['Bz'], p.get('dBz_dz', 0.0))
                    fx += fx1
                    fy += fy1
                    fz += fz1
                    dsxy += dsxy1
                    dnx += dnx1
                    dny += dny1
                elif reg.ftype == 'rf':
                    k_rf_default = 2.0 * np.pi * 2.856e9 / C_LIGHT
                    fx1, fy1, fz1 = rf_cavity_force(
                        beam, p['E_rf'],
                        p.get('dE_dz', 0.0),
                        p.get('dE_cdt', 0.0))
                    fx += fx1
                    fy += fy1
                    fz += fz1
                elif reg.ftype == 'drift':
                    pass    # zero contributions
        return fx, fy, fz, dsxy, dnx, dny

    return force_func


def build_gamma_prime_func(
    regions: List[ExtFieldRegion],
) -> Callable[[float, Beam6D], float]:
    """
    Build a function that returns dγ/dz at position z.

    Energy gain from RF cavities:  dγ/dz = η · E_rf.

    Parameters
    ----------
    regions : list of ExtFieldRegion

    Returns
    -------
    gamma_prime_func(z, beam) → dγ/dz [1/m]
    """
    def gamma_prime_func(z: float, beam: Beam6D) -> float:
        dgamma = 0.0
        for reg in regions:
            if reg.z_start <= z <= reg.z_end:
                if reg.ftype == 'rf':
                    dgamma += rf_acceleration_gradient(
                        reg.params.get('E_rf', 0.0))
        return dgamma

    return gamma_prime_func


def build_rf_chirp_func(
    regions: List[ExtFieldRegion],
) -> Callable[[float, Beam6D], float]:
    """
    Build a function that returns the RF chirp rate h = dδ/dz.

    Parameters
    ----------
    regions : list of ExtFieldRegion

    Returns
    -------
    chirp_func(z, beam) → h [1/m]
    """
    def chirp_func(z: float, beam: Beam6D) -> float:
        h_tot = 0.0
        for reg in regions:
            if reg.z_start <= z <= reg.z_end and reg.ftype == 'rf':
                E_rf = reg.params.get('E_rf', 0.0)
                k_rf = reg.params.get(
                    'k_rf', 2.0 * np.pi * 2.856e9 / C_LIGHT)
                h_tot += rf_chirp_rate(beam, E_rf, k_rf)
        return h_tot

    return chirp_func


def build_all_external(
    regions: List[ExtFieldRegion],
) -> Tuple[Callable, Callable, Callable]:
    """
    Convenience: build all three composite functions at once.

    Returns (external_force_func, gamma_prime_func, rf_chirp_func).
    Ready to pass to  beam_dynamics_6d.propagate().
    """
    return (
        build_external_force_func(regions),
        build_gamma_prime_func(regions),
        build_rf_chirp_func(regions),
    )


# ═════════════════════════════════════════════════════════════════════════════
# Demo / smoke test
# ═════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    from beam_dynamics_6d import (Beam6D, make_beam_100keV, propagate,
                                   beamK_to_gamma, get_alpha_interpolators)

    print("=" * 60)
    print("  External Forces Demo")
    print("=" * 60)

    _ = get_alpha_interpolators()

    # ── Beam at 100 keV ──
    beam0 = make_beam_100keV(Ne=1e5, beamK_eV=100_000,
                             sigma_z0_um=200, sigma_delta=1e-3)
    print(beam0.summary("Initial beam (z=0)"))

    # ── 1. Drift-only propagation ──
    print("\n── 1. Drift propagation ──")
    z_drift, st_drift = propagate(beam0, (0, 0.5), n_points=100,
                                  sc_model='ellipsoid')
    b1 = Beam6D.from_state(st_drift[-1, :11], st_drift[-1, 11],
                            beam0.Ne, beam0.eps_nx, beam0.eps_ny, beam0.eps_nz)
    print(f"  σ_x={b1.sigma_x*1e6:.2f} μm  σ_z={b1.sigma_z*1e6:.2f} μm")

    # ── 2. Solenoid beamline ──
    print("\n── 2. Solenoid (Bz=0.2 T, z=0.1-0.4 m) ──")
    regions_sol = [
        ExtFieldRegion(0.0, 0.1, 'drift', {}),
        ExtFieldRegion(0.1, 0.4, 'solenoid', {'Bz': 0.2, 'dBz_dz': 0.0}),
        ExtFieldRegion(0.4, 0.5, 'drift', {}),
    ]
    ef_sol, gp_sol, _ = build_all_external(regions_sol)
    z_sol, st_sol = propagate(beam0, (0, 0.5), n_points=200,
                              external_force_func=ef_sol,
                              sc_model='ellipsoid')
    b2 = Beam6D.from_state(st_sol[-1, :11], st_sol[-1, 11],
                            beam0.Ne, beam0.eps_nx, beam0.eps_ny, beam0.eps_nz)
    print(f"  σ_x={b2.sigma_x*1e6:.2f} μm  σ_y={b2.sigma_y*1e6:.2f} μm  "
          f"σ_z={b2.sigma_z*1e6:.2f} μm")
    if abs(b2.sigma_xy) > 1e-20:
        print(f"  σ_xy={b2.sigma_xy*1e12:.2f} pm²  (x-y coupling from Larmor)")

    # ── 3. RF cavity beamline ──
    print("\n── 3. RF cavity (E_rf=10 MV/m, z=0.1-0.3 m) ──")
    k_sband = 2.0 * np.pi * 2.856e9 / C_LIGHT
    regions_rf = [
        ExtFieldRegion(0.0, 0.1, 'drift', {}),
        ExtFieldRegion(0.1, 0.3, 'rf',
                       {'E_rf': 10e6, 'k_rf': k_sband,
                        'dE_dz': 0.0, 'dE_cdt': k_sband * 10e6}),
        ExtFieldRegion(0.3, 0.5, 'drift', {}),
    ]
    ef_rf, gp_rf, ch_rf = build_all_external(regions_rf)
    z_rf, st_rf = propagate(beam0, (0, 0.5), n_points=200,
                            external_force_func=ef_rf,
                            gamma_prime_func=gp_rf,
                            rf_chirp_func=ch_rf,
                            sc_model='ellipsoid')
    b3 = Beam6D.from_state(st_rf[-1, :11], st_rf[-1, 11],
                            beam0.Ne, beam0.eps_nx, beam0.eps_ny, beam0.eps_nz)
    print(f"  γ_end={b3.gamma:.4f}  (Δγ≈{b3.gamma-beam0.gamma:.4f})")
    print(f"  σ_x={b3.sigma_x*1e6:.2f} μm  σ_z={b3.sigma_z*1e6:.2f} μm  "
          f"σ_δ={b3.sigma_delta*1e3:.2f}‰")
    print(f"  C_zd={b3.C_zd*1e6:.2f} μm  Δt={b3.time_spread_ps:.1f} ps")

    # ── Plot ──
    fig, axes = plt.subplots(2, 3, figsize=(16, 8))

    # Row 0: transverse sizes
    labels = ['Drift', 'Solenoid (0.2 T)', 'RF (10 MV/m)']
    datasets = [(z_drift, st_drift), (z_sol, st_sol), (z_rf, st_rf)]
    for ax, (label, (zz, ss)) in zip(axes[0], zip(labels, datasets)):
        ax.plot(zz * 1e3, ss[:, 0] * 1e6, 'b-', lw=1, label=r'$\sigma_x$')
        ax.plot(zz * 1e3, ss[:, 1] * 1e6, 'r--', lw=1, label=r'$\sigma_y$')
        ax.set_title(label, fontsize=11)
        ax.set_ylabel(r'$\sigma_{x,y}$ (μm)')
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

    # Row 1: longitudinal
    for ax, (label, (zz, ss)) in zip(axes[1], zip(labels, datasets)):
        ax.plot(zz * 1e3, ss[:, 2] * 1e6, 'g-', lw=1.5)
        ax.set_xlabel('z (mm)')
        ax.set_ylabel(r'$\sigma_z$ (μm)')
        ax.grid(True, alpha=0.3)

    fig.suptitle("External Forces — Drift / Solenoid / RF Cavity", fontsize=13,
                 fontweight='bold')
    plt.tight_layout()
    plt.savefig("external_forces_demo.png", dpi=150)
    print("\nSaved: external_forces_demo.png")
    print("\n=== Demo complete ===")
