#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import yaml
from matplotlib.transforms import blended_transform_factory

from kmc_simulation import KineticMonteCarlo
from paper_plot_style import SINGLE_COLUMN_WIDTH, apply_publication_style, save_figure, style_legend, style_axes


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


def _load_sim_params(config_path: str, seed: int | None) -> dict[str, Any]:
    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f) or {}
    params: dict[str, Any] = {}
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
    elif params.get("rng_seed") is None:
        params["rng_seed"] = 0
    return params


def main() -> None:
    p = argparse.ArgumentParser(description="Plot the manuscript binding-energy landscape figure.")
    p.add_argument("--config", default="config_astro_full_paperfit.yaml")
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--out-dir", default="results/plots/manuscript")
    args = p.parse_args()

    params = _load_sim_params(args.config, args.seed)
    kmc = KineticMonteCarlo(params)
    mask = kmc.lattice[0] != None
    site_types = np.asarray(kmc.site_types[0][mask], dtype=int)
    ebind = np.asarray(kmc.E_bind_eV_map[0][mask], dtype=float)

    phys_mask = site_types == 1
    chem_mask = site_types == 2
    defect_mask = site_types == 3

    bins = np.logspace(-3, 0.5, 45)
    apply_publication_style()
    fig, ax = plt.subplots(figsize=(SINGLE_COLUMN_WIDTH, 2.55))
    ax.hist(ebind[phys_mask], bins=bins, color="#0072B2", alpha=0.72, label="Physisorption")
    ax.hist(ebind[chem_mask], bins=bins, color="#D55E00", alpha=0.72, label="Chemisorption")
    if np.any(defect_mask):
        ax.hist(ebind[defect_mask], bins=bins, histtype="step", color="#6d597a", linewidth=1.2, label="Defect sites")

    trans = blended_transform_factory(ax.transData, ax.transAxes)
    for temp, label in [(20, "20 K"), (100, "100 K"), (250, "250 K")]:
        energy = 8.617e-5 * temp
        ax.axvline(energy, color="#7f7f7f", ls="--", lw=0.7, alpha=0.8)
        ax.text(energy * 1.06, 0.92, label, rotation=90, fontsize=6.3, color="#6a6a6a", va="top", ha="left", transform=trans)

    style_axes(ax)
    ax.set_xscale("log")
    ax.set_xlim(1e-3, 3.2)
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel(r"Binding energy $E_{\mathrm{bind}}$ (eV)")
    ax.set_ylabel("Number of sites")
    style_legend(ax, loc="upper right", fontsize=7.0)

    os.makedirs(args.out_dir, exist_ok=True)
    out_base = os.path.join(args.out_dir, "fig02_binding_energy_distribution")
    save_figure(fig, out_base)
    plt.close(fig)
    print(f"Wrote {out_base}.png/.pdf")
