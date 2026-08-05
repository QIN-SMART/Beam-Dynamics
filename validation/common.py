"""
Shared helpers for validation tests: overlay plots, comparison metrics,
checkpoint logging.
"""

import os
import json

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPORTS_DIR = os.path.join(_THIS_DIR, "reports")
CHECKPOINTS = os.path.join(_THIS_DIR, "CHECKPOINTS.md")

PANELS = [
    ("sigma_x_um", r"$\sigma_x$ [$\mu$m]"),
    ("sigma_y_um", r"$\sigma_y$ [$\mu$m]"),
    ("sigma_z_um", r"$\sigma_z$ [$\mu$m]"),
    ("eps_nx_mm_mrad", r"$\varepsilon_{nx}$ [mm$\cdot$mrad]"),
    ("sigma_delta_e3", r"$\sigma_\delta$ [$10^{-3}$]"),
]


def plot_compare(section, results, out_png, title=None):
    """Overlay plot of one field per panel for several BeamResults.

    results : list of (label, BeamResult)
    """
    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)
    n = len(PANELS)
    fig, axes = plt.subplots(n, 1, figsize=(10, 2.4 * n), sharex=True)
    for ax, (field, label) in zip(axes, PANELS):
        for name, r in results:
            ax.plot(r.z_mm, getattr(r, field), lw=1.6, label=name)
        ax.set_ylabel(label)
        ax.legend(fontsize=8, loc="upper right")
        ax.grid(alpha=0.3)
    axes[-1].set_xlabel("z [mm]")
    fig.suptitle(title or f"{section} — AG vs OCELOT", fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    return out_png


def compare_metrics(a, b, fields=None):
    """Max-relative deviation (%) of BeamResult a vs b on common z range."""
    fields = fields or [f for f, _ in PANELS]
    zlo = max(a.z_mm[0], b.z_mm[0])
    zhi = min(a.z_mm[-1], b.z_mm[-1])
    if zhi <= zlo:
        return {}
    zq = np.linspace(zlo, zhi, 500)
    out = {}
    for f in fields:
        av = np.interp(zq, a.z_mm, getattr(a, f))
        bv = np.interp(zq, b.z_mm, getattr(b, f))
        out[f] = float(np.max(np.abs(av - bv) / np.maximum(np.abs(bv), 1e-12)) * 100)
    return out


def save_results(section, results, json_path=None):
    """Save one JSON per result under validation/reports/."""
    json_path = json_path or os.path.join(REPORTS_DIR, f"{section}_results.json")
    os.makedirs(os.path.dirname(json_path), exist_ok=True)
    payload = {r.route: r.to_dict() for r in results}
    with open(json_path, "w") as f:
        json.dump(payload, f, indent=2)
    return json_path


def log_checkpoint(section, lines):
    """Append a dated checkpoint entry to validation/CHECKPOINTS.md."""
    os.makedirs(os.path.dirname(CHECKPOINTS), exist_ok=True)
    header = f"\n## [{section}] {__import__('datetime').datetime.now():%Y-%m-%d %H:%M}\n"
    with open(CHECKPOINTS, "a") as f:
        f.write(header)
        for ln in lines:
            f.write("  " + str(ln) + "\n")
    return CHECKPOINTS


def print_summary(section, results, metrics=None):
    print(f"\n== {section} ==")
    for name, r in results:
        print(f"  {name:<8s} {r.summary()}")
    if metrics:
        print("  AG vs OCELOT max rel dev (%):",
              {k: round(v, 2) for k, v in metrics.items()})
