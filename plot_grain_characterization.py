#!/usr/bin/env python3
"""
Plot one grain realization's surface energetics and site-type mix.

Example:
  python plot_grain_characterization.py \
    --config config_astro_full_paperfit.yaml
"""

from __future__ import annotations

import argparse
import os
from typing import Any, Dict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from kmc_simulation import KineticMonteCarlo
from paper_plot_style import DOUBLE_COLUMN_WIDTH, add_panel_labels, apply_publication_style, save_figure, style_axes, style_legend


NON_SIM_KEYS = {
    "mode",
    "output_filename",
    "raw_runs_output",
    "mrn_output_filename",
    "parameter_sweeps",
    "explicit_conditions",
    "save_raw_runs",
    "ensemble_runs",
    "min_ensemble_runs",
    "max_time_s",
    "max_steps",
    "use_mrn",
    "aggregate_across_sizes",
    "mrn_min_um",
    "mrn_max_um",
    "mrn_bins",
    "burnin_arrivals",
    "measure_arrivals",
}

SITE_LABELS = {
    1: "Physisorption",
    2: "Chemisorption",
    3: "Defect",
}


def _load_sim_params(config_path: str, seed: int | None) -> Dict[str, Any]:
    with open(config_path, "r") as f:
        config = yaml.safe_load(f) or {}

    params: Dict[str, Any] = {}
    for key, value in config.items():
        if key in NON_SIM_KEYS:
            continue
        if isinstance(value, list):
            continue
        if isinstance(value, dict) and key != "temp_ramp":
            continue
        params[key] = value

    if seed is not None:
        params["rng_seed"] = int(seed)
    elif params.get("rng_seed", None) is None:
        params["rng_seed"] = 0
    return params


def main() -> None:
    p = argparse.ArgumentParser(description="Plot a grain realization's binding/diffusion-energy distributions.")
    p.add_argument("--config", default="config_astro_full_paperfit.yaml", help="YAML config used to build the grain")
    p.add_argument("--seed", type=int, default=None, help="Optional RNG seed override for the grain realization")
    p.add_argument("--include-bulk", action="store_true", help="Use all non-void lattice cells instead of only the top layer")
    p.add_argument("--out-dir", default="results/plots/grain_characterization", help="Output directory")
    args = p.parse_args()

    sim_params = _load_sim_params(args.config, args.seed)
    kmc = KineticMonteCarlo(sim_params)

    if args.include_bulk:
        mask = kmc.lattice != None
        layer_tag = "bulk"
    else:
        mask = kmc.lattice[0] != None
        layer_tag = "surface"

    if not np.any(mask):
        raise SystemExit("No accessible lattice cells found for the requested grain realization.")

    if args.include_bulk:
        bind = np.asarray(kmc.E_bind_eV_map[mask], dtype=float)
        diff = np.asarray(kmc.E_diff_eV_map[mask], dtype=float)
        site_types = np.asarray(kmc.site_types[mask], dtype=int)
    else:
        bind = np.asarray(kmc.E_bind_eV_map[0][mask], dtype=float)
        diff = np.asarray(kmc.E_diff_eV_map[0][mask], dtype=float)
        site_types = np.asarray(kmc.site_types[0][mask], dtype=int)

    valid_ratio = bind > 0.0
    ratio = diff[valid_ratio] / bind[valid_ratio]
    ratio = ratio[np.isfinite(ratio)]

    summary_rows = []
    total_sites = int(site_types.size)
    for site_code, label in SITE_LABELS.items():
        site_mask = site_types == int(site_code)
        count = int(np.sum(site_mask))
        frac = float(count / total_sites) if total_sites > 0 else 0.0
        bind_mean = float(np.mean(bind[site_mask])) if count > 0 else 0.0
        diff_mean = float(np.mean(diff[site_mask])) if count > 0 else 0.0
        ratio_mean = float(np.mean((diff[site_mask] / bind[site_mask])[bind[site_mask] > 0])) if count > 0 and np.any(bind[site_mask] > 0) else 0.0
        summary_rows.append(
            {
                "site_type": label,
                "count": count,
                "fraction": frac,
                "mean_E_bind_eV": bind_mean,
                "mean_E_diff_eV": diff_mean,
                "mean_Ediff_over_Ebind": ratio_mean,
            }
        )

    os.makedirs(args.out_dir, exist_ok=True)

    apply_publication_style()
    fig, axes = plt.subplots(2, 2, figsize=(DOUBLE_COLUMN_WIDTH, 4.9))

    axes[0, 0].hist(bind, bins=40, color="#1d3557", alpha=0.84)
    axes[0, 0].set_xlabel(r"$E_{\rm bind}$ (eV)")
    axes[0, 0].set_ylabel("Count")
    style_axes(axes[0, 0])

    axes[0, 1].hist(diff, bins=40, color="#c05640", alpha=0.84)
    axes[0, 1].set_xlabel(r"$E_{\rm hop}$ (eV)")
    axes[0, 1].set_ylabel("Count")
    style_axes(axes[0, 1])

    if ratio.size >= 2 and (float(np.nanmax(ratio)) - float(np.nanmin(ratio))) > 1e-12:
        axes[1, 0].hist(ratio, bins=min(30, max(8, int(np.sqrt(ratio.size)))), color="#2a9d8f", alpha=0.84)
        axes[1, 0].axvline(float(np.mean(ratio)), color="black", linestyle="--", linewidth=1.5, label=f"mean={np.mean(ratio):.3f}")
        style_legend(axes[1, 0], loc="upper right")
    elif ratio.size >= 2:
        axes[1, 0].axvline(float(ratio[0]), color="#2a9d8f", linewidth=2.0, label=f"all values={float(ratio[0]):.3f}")
        style_legend(axes[1, 0], loc="upper right")
        axes[1, 0].text(0.5, 0.55, "Uniform across accessible sites", ha="center", va="center", transform=axes[1, 0].transAxes, fontsize=7.1, color="#4a4a4a")
    elif ratio.size == 1:
        axes[1, 0].axvline(float(ratio[0]), color="#2a9d8f", linewidth=2.0, label=f"value={float(ratio[0]):.3f}")
        style_legend(axes[1, 0], loc="upper right")
        axes[1, 0].text(0.5, 0.55, "Single finite ratio value", ha="center", va="center", transform=axes[1, 0].transAxes, fontsize=7.1, color="#4a4a4a")
    else:
        axes[1, 0].text(0.5, 0.55, r"No finite $E_{\rm hop}/E_{\rm bind}$ values", ha="center", va="center", transform=axes[1, 0].transAxes, fontsize=7.1, color="#4a4a4a")
    axes[1, 0].set_xlabel(r"$E_{\rm hop}/E_{\rm bind}$")
    axes[1, 0].set_ylabel("Count")
    style_axes(axes[1, 0])

    site_counts = [int(np.sum(site_types == code)) for code in SITE_LABELS]
    site_fracs = [count / total_sites if total_sites else 0.0 for count in site_counts]
    site_names = [SITE_LABELS[code] for code in SITE_LABELS]
    bars = axes[1, 1].bar(site_names, site_fracs, color=["#1d3557", "#c05640", "#6d597a"])
    axes[1, 1].set_ylabel("Fraction of accessible sites")
    axes[1, 1].tick_params(axis="x", rotation=15)
    style_axes(axes[1, 1])
    for bar, frac in zip(bars, site_fracs):
        axes[1, 1].text(bar.get_x() + bar.get_width() / 2.0, frac + 0.015, f"{100.0 * frac:.1f}%", ha="center", va="bottom", fontsize=7.2)
    axes[1, 1].set_ylim(0.0, max(site_fracs + [0.45]) + 0.08)

    add_panel_labels(axes.flatten(), labels="ABCD")
    grain_radius_um = float(sim_params.get("grain_radius_um", 0.0) or 0.0)
    fig.text(
        0.985,
        0.99,
        f"{layer_tag} lattice, r={grain_radius_um:g} μm",
        ha="right",
        va="top",
        fontsize=6.9,
        color="#4a4a4a",
    )
    fig.subplots_adjust(left=0.08, right=0.985, bottom=0.1, top=0.95, hspace=0.28, wspace=0.22)

    stem = f"grain_characterization_{layer_tag}"
    out_png = os.path.join(args.out_dir, stem)
    save_figure(fig, out_png)
    plt.close(fig)

    summary = pd.DataFrame(summary_rows)
    summary.loc[len(summary)] = {
        "site_type": "All",
        "count": total_sites,
        "fraction": 1.0,
        "mean_E_bind_eV": float(np.mean(bind)),
        "mean_E_diff_eV": float(np.mean(diff)),
        "mean_Ediff_over_Ebind": float(np.mean(ratio)),
    }
    out_csv = os.path.join(args.out_dir, f"{stem}.csv")
    summary.to_csv(out_csv, index=False)

    print(f"Wrote {out_png}.png/.pdf")
    print(f"Wrote {out_csv}")


if __name__ == "__main__":
    main()
