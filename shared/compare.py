#!/usr/bin/env python3
"""
Cross-route comparison — AG (6D envelope ODE) vs GPT模拟 (OCELOT tracking).

Reads shared/results/AG_results.json + shared/results/GPT_results.json
(any two <route>_results.json files) and produces:
  1. terminal comparison table at the shared probe positions
  2. overlay figures  σ_x(z), σ_z(z), ε_nx(z), σ_δ(z)

Usage:
  python3 shared/compare.py
  python3 shared/compare.py --out shared/results/compare_overlay.png
  python3 shared/compare.py --results shared/results
"""

import os
import sys
import argparse

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from output_schema import list_results, load_results  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

PROBE_FIELDS = [
    ("sigma_x_um", "σ_x (μm)"),
    ("sigma_y_um", "σ_y (μm)"),
    ("sigma_z_um", "σ_z (μm)"),
    ("sigma_delta_e3", "σ_δ (×10⁻³)"),
    ("eps_nx_mm_mrad", "ε_nx (mm·mrad)"),
    ("eps_ny_mm_mrad", "ε_ny (mm·mrad)"),
]

HIST_FIELDS = [
    ("sigma_x_um", "σ_x (μm)", "#1f77b4"),
    ("sigma_z_um", "σ_z (μm)", "#d62728"),
    ("eps_nx_mm_mrad", "ε_nx (mm·mrad)", "#2ca02c"),
    ("sigma_delta_e3", "σ_δ (×10⁻³)", "#9467bd"),
]


def _routes(results_dir):
    return {load_results(p)["route"]: load_results(p)
            for p in list_results(results_dir)}


def _interp_history(hist, field, z_query_mm):
    """Linear interpolation of a history field at query z positions."""
    z = np.asarray(hist["z_mm"])
    v = np.asarray(hist[field])
    return np.interp(z_query_mm, z, v)


def print_probe_table(routes):
    route_names = list(routes.keys())
    zs = sorted({round(p["z_mm"]) for r in routes.values() for p in r["probes"]})
    hdr = "z(mm)"
    for name in route_names:
        hdr += f"   σ_x({name})   ε_nx({name})"
    print("Probe comparison (AG = envelope, GPT = tracking):")
    print(f"  {'z(mm)':>6s}" + "".join(f"   {'σ_x(um)':>9s}   {'εnx(mm.mrad)':>12s}" for _ in route_names))
    for z in zs:
        row = f"  {z:6.0f}"
        for r in routes.values():
            p = next((p for p in r["probes"] if abs(round(p["z_mm"]) - z) < 1e-9), None)
            if p is None:
                row += f"   {'—':>9s}   {'—':>12s}"
            else:
                row += f"   {p['sigma_x_um']:9.1f}   {p['eps_nx_mm_mrad']:12.4f}"
        print(row)


def print_full_history_stats(routes):
    """Per-route final-plane summary from history."""
    for name, r in routes.items():
        h = r["history"]
        print(f"  {name:>4s}: z_final={h['z_mm'][-1]:6.0f} mm  "
              f"σ_x={h['sigma_x_um'][-1]:8.1f} μm  σ_z={h['sigma_z_um'][-1]:8.1f} μm  "
              f"σ_δ={h['sigma_delta_e3'][-1]:6.2f}e-3  "
              f"ε_nx={h['eps_nx_mm_mrad'][-1]:7.4f} mm·mrad")


def plot_overlay(routes, out_png):
    n = len(HIST_FIELDS)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    axes = axes.ravel()
    for ax, (field, label, color) in zip(axes, HIST_FIELDS):
        for name, r in routes.items():
            h = r["history"]
            ax.plot(h["z_mm"], h[field], label=f"{name}", lw=1.6)
        ax.set_xlabel("z (mm)")
        ax.set_ylabel(label)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    fig.suptitle("UED Beamline — AG vs GPT  (shared config)",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    if out_png:
        os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)
        plt.savefig(out_png, dpi=150)
        print(f"  saved -> {out_png}")
    return fig


def main():
    ap = argparse.ArgumentParser(description="Compare AG vs GPT results")
    ap.add_argument("--results", default=None, help="results dir")
    ap.add_argument("--out", default=None,
                    help="overlay png path (default: <results>/compare_overlay.png)")
    args = ap.parse_args()

    results_dir = args.results
    routes = _routes(results_dir)

    if len(routes) < 2:
        print(f"  Need at least 2 result files in {results_dir or 'shared/results'}")
        print(f"  found: {list(routes.keys())}")
        sys.exit(1)

    shas = {r["config_sha"] for r in routes.values()}
    scs = {r["sc_enabled"] for r in routes.values()}
    print("routes:", {k: v for k, v in routes.items()})
    print(f"  config_sha  = {shas}   {'(MATCH ✓)' if len(shas)==1 else '(MISMATCH!)'}")
    print(f"  sc_enabled  = {scs}")

    print_full_history_stats(routes)
    print()
    print_probe_table(routes)

    out_png = args.out or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "results", "compare_overlay.png")
    plot_overlay(routes, out_png)


if __name__ == "__main__":
    main()
