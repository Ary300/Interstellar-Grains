#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import yaml
from matplotlib.patches import Patch

from kmc_simulation import KineticMonteCarlo
from paper_plot_style import SINGLE_COLUMN_WIDTH, apply_publication_style, save_figure


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


def _surface_shell_mask(lattice: np.ndarray) -> np.ndarray:
    occupied = lattice != None
    mask = np.zeros_like(occupied, dtype=bool)
    depth, rows, cols = occupied.shape
    neighbor_offsets = [
        (-1, 0, 0),
        (1, 0, 0),
        (0, -1, 0),
        (0, 1, 0),
        (0, 0, -1),
        (0, 0, 1),
    ]
    for d, r, c in np.argwhere(occupied):
        for dd, dr, dc in neighbor_offsets:
            nd, nr, nc = d + dd, r + dr, c + dc
            if nd < 0 or nr < 0 or nc < 0 or nd >= depth or nr >= rows or nc >= cols or not occupied[nd, nr, nc]:
                mask[d, r, c] = True
                break
    return mask


def main() -> None:
    p = argparse.ArgumentParser(description="Render the outer grain lattice shell for manuscript Figure 1.")
    p.add_argument("--config", default="config_astro_full_paperfit.yaml")
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--out-dir", default="results/plots/manuscript")
    args = p.parse_args()

    params = _load_sim_params(args.config, args.seed)
    kmc = KineticMonteCarlo(params)
    lattice = np.asarray(kmc.lattice, dtype=object)
    site_types = np.asarray(kmc.site_types, dtype=int)
    shell = _surface_shell_mask(lattice)

    spacing_a = float(np.sqrt(float(params.get("site_area_angstroms_sq", 25.0))))
    zz, yy, xx = np.where(shell)
    xs = xx.astype(float) * spacing_a
    ys = yy.astype(float) * spacing_a
    zs = zz.astype(float) * spacing_a

    colors = []
    sizes = []
    for d, r, c in zip(zz, yy, xx):
        st = int(site_types[d, r, c])
        if st == 2:
            colors.append("#0072B2")
            sizes.append(18.0)
        elif st == 3:
            colors.append("#D55E00")
            sizes.append(18.0)
        else:
            colors.append("#CCCCCC")
            sizes.append(11.0)

    apply_publication_style()
    fig = plt.figure(figsize=(SINGLE_COLUMN_WIDTH, 3.42))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_proj_type("persp", focal_length=1.55)

    for color, size in sorted(set(zip(colors, sizes)), key=lambda x: x[1]):
        sel = [(c == color and s == size) for c, s in zip(colors, sizes)]
        alpha = 0.92 if color != "#CCCCCC" else 0.35
        ax.scatter(xs[sel], ys[sel], zs[sel], c=color, alpha=alpha, s=size, marker="s", depthshade=False, linewidths=0)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.grid(False)
    ax.view_init(elev=29, azim=-50)
    ax.set_xlim(xs.min() - spacing_a, xs.max() + spacing_a)
    ax.set_ylim(ys.min() - spacing_a, ys.max() + spacing_a)
    ax.set_zlim(zs.min() - 0.5 * spacing_a, zs.max() + 1.5 * spacing_a)
    ax.set_box_aspect((1.0, 1.0, 0.34))
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor((1.0, 1.0, 1.0, 0.0))
        axis.pane.set_edgecolor((0.85, 0.85, 0.85, 0.25))
    ax.text2D(0.27, 0.08, "x (A)", transform=ax.transAxes, fontsize=7.4, rotation=-17, color="#222222")
    ax.text2D(0.75, 0.09, "y (A)", transform=ax.transAxes, fontsize=7.4, rotation=16, color="#222222")
    ax.text2D(0.05, 0.56, "z (A)", transform=ax.transAxes, fontsize=7.7, rotation=90, color="#222222")

    legend_handles = [
        Patch(color="#CCCCCC", label="Regular site"),
        Patch(color="#D55E00", label="Surface defect"),
        Patch(color="#0072B2", label="Chemisorption site"),
    ]
    ax.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(0.03, 0.97),
        fontsize=6.5,
        frameon=False,
        borderpad=0.2,
        handlelength=0.9,
    )
    ax.set_position([0.02, 0.03, 0.96, 0.92])
    fig.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)

    os.makedirs(args.out_dir, exist_ok=True)
    out_base = os.path.join(args.out_dir, "fig01_grain_structure")
    save_figure(fig, out_base)
    plt.close(fig)
    print(f"Wrote {out_base}.png/.pdf")


if __name__ == "__main__":
    main()
