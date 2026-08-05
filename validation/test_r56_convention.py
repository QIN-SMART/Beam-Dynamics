#!/usr/bin/env python3
"""
R56 / longitudinal-coordinate convention — AUDIT + CHARACTERIZATION TEST.

Parts:
  1. Explicit definitions (see report)
  2. Exact relativistic reference (independent, no AG/OCELOT)
  3. OCELOT raw-coordinate characterization (installed OCELOT)
  4. Installed OCELOT source-convention inspection
  5. AG convention characterization
  6. Coordinate-transformation closure test

Outputs (validation/reports/):
  R56_convention_resolution.md  (report)
  r56_convention_results.json   (all measured numbers)
  r56_convention_plot.png       (6 panels)

No production physics code is modified.
"""

import json
import os
import sys

import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPORTS = os.path.join(_THIS_DIR, "reports")
for p in (_THIS_DIR, os.path.dirname(_THIS_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# ──────────────────────────────────────────────────────────────────────────
# Part 1 — definitions (constants)
# ──────────────────────────────────────────────────────────────────────────
K0_keV = 100.0
M_E = 9.10938356e-31
C = 2.99792458e8
L = 0.5
D_ARR = np.array([-1e-3, -5e-4, -1e-4, 1e-4, 5e-4, 1e-3])

# reference particle
E0_J = K0_keV * 1e3 * 1.602176634e-19
GAMMA0 = 1.0 + E0_J / (M_E * C**2)
BETA0 = np.sqrt(1.0 - 1.0 / GAMMA0**2)
P0 = GAMMA0 * BETA0 * M_E * C


def kin_from_p(p):
    gamma = np.sqrt(1.0 + (p / (M_E * C))**2)
    beta = p / (gamma * M_E * C)
    return gamma, beta


def main():
    out = {"part1": {
        "K0_keV": K0_keV, "E0_J": E0_J, "gamma0": GAMMA0, "beta0": BETA0,
        "p0_kg_m_s": P0, "deltas": D_ARR.tolist(),
        "definitions": {
            "delta_p": "(p-p0)/p0 momentum deviation",
            "delta_E": "dE/E0 energy deviation",
            "p_ocelot": "dE/(c*p0) — OCELOT native rparticles[5]",
            "Delta_t": "t_particle - t_reference (lab)",
            "Delta_ct": "c*Delta_t",
            "Delta_z": "-beta0*c*Delta_t (co-moving spatial, head positive)",
            "signs": "delta_p>0 -> arrives earlier (Delta_t<0), ahead (Delta_z>0)",
        }}}

    # ──────────────────────────────────────────────────────────────────────
    # Part 2 — exact relativistic reference
    # ──────────────────────────────────────────────────────────────────────
    dt_ex = np.zeros_like(D_ARR)
    dct_ex = np.zeros_like(D_ARR)
    dz_ex = np.zeros_like(D_ARR)
    for i, d in enumerate(D_ARR):
        p = P0 * (1.0 + d)
        g, b = kin_from_p(p)
        t0 = L / (BETA0 * C)
        t = L / (b * C)
        dt = t - t0
        dt_ex[i] = dt
        dct_ex[i] = C * dt
        dz_ex[i] = -BETA0 * C * dt

    def slope(x, y):
        return float(np.polyfit(x, y, 1)[0])

    s_dt, s_dct, s_dz = slope(D_ARR, dt_ex), slope(D_ARR, dct_ex), slope(D_ARR, dz_ex)
    R56_t = -L / (BETA0 * C * GAMMA0**2)
    R56_ct = -L / (BETA0 * GAMMA0**2)
    R56_z = L / GAMMA0**2

    def cmp(fit, ref):
        return {"fit": fit, "ref": ref,
                "abs_err": float(fit - ref), "rel_err": float((fit - ref) / ref)}

    part2 = {
        "R56_t_exact_fit": cmp(s_dt, R56_t),
        "R56_ct_exact_fit": cmp(s_dct, R56_ct),
        "R56_z_exact_fit": cmp(s_dz, R56_z),
    }
    out["part2"] = part2

    # ──────────────────────────────────────────────────────────────────────
    # Part 3 — OCELOT raw characterization (installed OCELOT 26.06.1)
    # ──────────────────────────────────────────────────────────────────────
    import ocelot
    from ocelot.cpbd.elements import Drift
    from ocelot.cpbd.magnetic_lattice import MagneticLattice
    from ocelot.cpbd.beam import generate_parray
    from ocelot.cpbd.navi import Navigator
    from ocelot.cpbd.track import tracking_step

    E_tot_GeV = (K0_keV + 511.0) * 1e-6
    lat = MagneticLattice([Drift(l=L)], method={"global": "SecondTM"})
    lat.update_transfer_maps()
    R = lat.sequence[0].first_order_tms[0].get_params(E_tot_GeV).R
    R56_ocelot = float(R[4, 5])

    pa = generate_parray(sigma_x=1e-6, sigma_y=1e-6, sigma_tau=1e-9,
                         energy=E_tot_GeV, charge=1e-15, nparticles=len(D_ARR))
    pa.rparticles[:] = 0.0
    pa.rparticles[4, :] = 0.0                 # raw tau = 0
    pa.rparticles[5, :] = D_ARR               # raw p = dE/(c·p0) convention
    pa.E = E_tot_GeV

    tau0 = pa.tau().copy()
    p0_arr = pa.p().copy()
    navi = Navigator(lat)
    for _ in range(5):
        tracking_step(lat, pa, L / 5.0, navi)
    d_tau_raw = pa.tau() - tau0
    slope_oc_raw = slope(D_ARR, d_tau_raw)

    part3 = {
        "ocelot_version": getattr(ocelot, "__version__", "unknown"),
        "ocelot_file": ocelot.__file__,
        "E_tot_GeV": E_tot_GeV, "beta0": BETA0, "gamma0": GAMMA0,
        "R56_transfer_map": R56_ocelot,
        "raw_tau_before": tau0.tolist(), "raw_p_before": p0_arr.tolist(),
        "raw_tau_after": pa.tau().tolist(), "raw_p_after": pa.p().tolist(),
        "d_tau_raw": d_tau_raw.tolist(),
        "slope_d_tau_raw": slope_oc_raw,
        "vs_ct_exact": cmp(slope_oc_raw, R56_ct),          # -L/(βγ²)
        "vs_ct_beta2": cmp(slope_oc_raw, -L / (BETA0**2 * GAMMA0**2)),
        "vs_transfer_map": cmp(slope_oc_raw, R56_ocelot),
    }
    out["part3"] = part3

    # ──────────────────────────────────────────────────────────────────────
    # Part 4 — installed OCELOT source-convention inspection
    # ──────────────────────────────────────────────────────────────────────
    ocelot_base = os.path.dirname(ocelot.__file__)
    src = {}
    def grab(path, keys):
        full = os.path.join(ocelot_base, path)
        lines = open(full).read().splitlines()
        hits = []
        for i, ln in enumerate(lines, 1):
            if any(k in ln for k in keys):
                hits.append((i, ln.strip()))
        return {"file": full, "hits": hits}
    src["beam_docstring"] = grab("cpbd/beam/beam.py",
                                 ["ds = c*tau", "dE/(p0*c)", "c*tau",
                                  "E/(c*p0)", "sigma_tau: std", "canonical"])
    src["r_matrix_r56"] = grab("cpbd/r_matrix.py",
                               ["igamma2", "r56 -= z", "r56 = hx"])
    src["solenoid_r56"] = grab("cpbd/elements/solenoid_atom.py", ["r56"])
    out["part4_source"] = src

    # ──────────────────────────────────────────────────────────────────────
    # Part 5 — AG convention characterization
    # ──────────────────────────────────────────────────────────────────────
    dz_ag = (L / GAMMA0**2) * D_ARR          # AG map: Δz = (L/γ²)·δ
    rel_err_ag = np.abs(dz_ag - dz_ex) / np.maximum(np.abs(dz_ex), 1e-30)
    # linearity: fit exact Δz vs δ with a quadratic, report curvature
    c2 = np.polyfit(D_ARR, dz_ex, 2)[0]
    part5 = {
        "map": "Delta_z_AG = (L/gamma0^2)*delta  [spatial, delta=delta_p]",
        "delta_z_ag": dz_ag.tolist(),
        "delta_z_exact": dz_ex.tolist(),
        "max_rel_err": float(np.max(rel_err_ag)),
        "quadratic_curvature": float(c2),
        "where_used": ("AG envelope_ode: dC_zd/dz = sigma_delta^2/gamma^2 "
                       "(R56_z = L/gamma^2); entered via make_beam_100keV "
                       "sigma_delta and the drift dynamics"),
    }
    out["part5"] = part5

    # ──────────────────────────────────────────────────────────────────────
    # Part 6 — coordinate-transformation closure test
    # ──────────────────────────────────────────────────────────────────────
    # Formal derivation (from the DOCUMENTED definitions):
    #   tau = c*t  (OCELOT)   →  Delta_t = Delta_tau / c
    #   p_ocelot = dE/(c·p0)  →  delta_p = p_ocelot / beta0
    #   Delta_z = -beta0*c*Delta_t = -beta0*Delta_tau
    # Closure:  Delta_z_formal = -beta0 * R56_ocelot * p_ocelot
    #         = -beta0 * (-L/(beta0^2*gamma0^2)) * p_ocelot
    #         = L * p_ocelot / (beta0*gamma0^2)  =  L*delta_p/gamma0^2  = exact
    delta_p_from_oc = D_ARR / BETA0                    # p_oc → δ_p
    dz_formal = -BETA0 * (R56_ocelot * D_ARR)          # from raw OCELOT transport
    # exact Δz at the same PHYSICAL δ_p
    dz_ex_p = np.zeros_like(D_ARR)
    for i, dp in enumerate(delta_p_from_oc):
        p = P0 * (1.0 + dp)
        g, b = kin_from_p(p)
        dz_ex_p[i] = -BETA0 * C * (L / (b * C) - L / (BETA0 * C))
    residual_formal = np.abs(dz_formal - dz_ex_p) / np.maximum(np.abs(dz_ex_p), 1e-30)
    # naive framework feeding (p_oc raw value used as delta_p): residual = 1/β − 1
    dz_naive = -BETA0 * (R56_ocelot * D_ARR)           # same as formal, but δ_p = raw
    dz_ex_naive = np.zeros_like(D_ARR)
    for i, dp in enumerate(D_ARR):                     # δ_p wrongly = raw value
        p = P0 * (1.0 + dp)
        g, b = kin_from_p(p)
        dz_ex_naive[i] = -BETA0 * C * (L / (b * C) - L / (BETA0 * C))
    residual_naive = np.abs(dz_naive - dz_ex_naive) / np.maximum(np.abs(dz_ex_naive), 1e-30)

    part6 = {
        "derived_transforms": {
            "delta_t": "Delta_tau / c",
            "delta_p": "p_ocelot / beta0",
            "delta_z": "-beta0 * Delta_tau",
        },
        "closure_after_formal_conversion_max_rel_residual":
            float(np.max(residual_formal)),
        "closure_naive_feeding_max_rel_residual": float(np.max(residual_naive)),
        "expected_naive_residual_1_over_beta_minus_1": 1.0 / BETA0 - 1.0,
        "residual_formal_per_delta": residual_formal.tolist(),
    }
    out["part6"] = part6

    # ──────────────────────────────────────────────────────────────────────
    # plot — 6 panels
    # ──────────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(3, 2, figsize=(11, 11))
    ax = axes[0, 0]
    ax.plot(D_ARR * 1e3, dt_ex * 1e12, "o-")
    ax.set_xlabel(r"$\delta$ [$10^{-3}$]"); ax.set_ylabel(r"$\Delta t$ [ps]")
    ax.set_title("1. Exact lab arrival-time difference")
    ax.grid(alpha=0.3)

    ax = axes[0, 1]
    ax.plot(D_ARR * 1e3, dct_ex * 1e6, "o-")
    ax.set_xlabel(r"$\delta$ [$10^{-3}$]"); ax.set_ylabel(r"$\Delta(ct)$ [$\mu$m]")
    ax.set_title("2. Exact c-time difference")
    ax.grid(alpha=0.3)

    ax = axes[1, 0]
    ax.plot(D_ARR * 1e3, dz_ex * 1e6, "o-", label="exact")
    ax.plot(D_ARR * 1e3, dz_ag * 1e6, "s--", label="AG map (L/γ²)")
    ax.set_xlabel(r"$\delta_p$ [$10^{-3}$]"); ax.set_ylabel(r"$\Delta z$ [$\mu$m]")
    ax.set_title("3. Exact co-moving Δz vs AG map")
    ax.legend(); ax.grid(alpha=0.3)

    ax = axes[1, 1]
    ax.plot(D_ARR * 1e3, d_tau_raw * 1e6, "o-", label="raw OCELOT")
    fit = np.poly1d(np.polyfit(D_ARR, d_tau_raw, 1))
    ax.plot(D_ARR * 1e3, fit(D_ARR) * 1e6, "--", label=f"fit slope {slope_oc_raw:.3f} m")
    ax.set_xlabel(r"$p_{\mathrm{ocelot}}$ [$10^{-3}$]")
    ax.set_ylabel(r"$\Delta\tau_{\mathrm{raw}}$ [$\mu$m]")
    ax.set_title("4. OCELOT raw tau change")
    ax.legend(); ax.grid(alpha=0.3)

    ax = axes[2, 0]
    ax.plot(D_ARR * 1e3, dz_ex * 1e6, "o-", label="exact Δz")
    ax.plot(D_ARR * 1e3, dz_ex_p * 1e6, "s--", label="exact Δz at δ_p=p_oc/β")
    ax.set_xlabel(r"$p_{\mathrm{ocelot}}$ [$10^{-3}$]"); ax.set_ylabel(r"$\Delta z$ [$\mu$m]")
    ax.set_title("5. AG Δz vs exact (δ_p convention)")
    ax.legend(); ax.grid(alpha=0.3)

    ax = axes[2, 1]
    ax.semilogy(D_ARR * 1e3, np.maximum(residual_formal, 1e-12), "o-",
                label="formal conversion residual")
    ax.semilogy(D_ARR * 1e3, np.maximum(residual_naive, 1e-12), "s--",
                label="naive feeding residual (1/β−1)")
    ax.set_xlabel(r"$p_{\mathrm{ocelot}}$ [$10^{-3}$]")
    ax.set_ylabel("relative residual")
    ax.set_title("6. Residuals after derived coordinate conversion")
    ax.legend(); ax.grid(alpha=0.3)

    fig.suptitle(f"R56 coordinate convention — K0={K0_keV:.0f} keV, L={L} m",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    os.makedirs(_REPORTS, exist_ok=True)
    plot_path = os.path.join(_REPORTS, "r56_convention_plot.png")
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)

    # ──────────────────────────────────────────────────────────────────────
    # save results
    # ──────────────────────────────────────────────────────────────────────
    with open(os.path.join(_REPORTS, "r56_convention_results.json"), "w") as f:
        json.dump(out, f, indent=2)

    # console summary
    print("Part 2 (exact reference, independent):")
    for k, v in part2.items():
        print(f"  {k}: fit={v['fit']:.6e} ref={v['ref']:.6e} "
              f"rel_err={v['rel_err']:.2e}")
    print(f"Part 3 (OCELOT raw): version={part3['ocelot_version']}  "
          f"R56_tm={R56_ocelot:.6f}  slope_raw={slope_oc_raw:.6f}")
    print(f"  raw slope vs ct_exact(-L/βγ²)={part3['vs_ct_exact']['rel_err']:.2e}  "
          f"vs -L/(β²γ²)={part3['vs_ct_beta2']['rel_err']:.2e}")
    print(f"Part 5 (AG): map Δz=(L/γ²)δ, max rel err vs exact = "
          f"{part5['max_rel_err']:.2e}")
    print(f"Part 6 (closure): formal-conversion residual = "
          f"{part6['closure_after_formal_conversion_max_rel_residual']:.2e}  |  "
          f"naive-feeding residual = "
          f"{part6['closure_naive_feeding_max_rel_residual']:.3f} "
          f"(expected 1/β−1 = {part6['expected_naive_residual_1_over_beta_minus_1']:.3f})")
    print(f"plot  -> {plot_path}")
    print(f"json  -> {os.path.join(_REPORTS, 'r56_convention_results.json')}")
    return out


if __name__ == "__main__":
    main()
