#!/usr/bin/env python3
"""
Plot cumulative ensemble convergence for one raw-run condition.

Example:
  python plot_ensemble_convergence.py \
    --raw results/astro_transition_deep_raw.csv \
    --temp 100 --nH 1000
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from paper_plot_style import SINGLE_COLUMN_WIDTH, apply_publication_style, save_figure, style_axes, style_legend


def _sem(values: np.ndarray) -> float:
    if values.size <= 1:
        return 0.0
    return float(np.std(values, ddof=1) / np.sqrt(values.size))


def _pick_points(n_runs: int, requested: list[int] | None) -> list[int]:
    if requested:
        return sorted({n for n in requested if 1 <= n <= n_runs})
    default = [10, 20, 50, 100, 200, 500, 1000]
    picked = [n for n in default if n <= n_runs]
    if n_runs not in picked:
        picked.append(n_runs)
    return picked


def main() -> None:
    p = argparse.ArgumentParser(description="Plot running mean ± SEM versus ensemble size.")
    p.add_argument("--raw", default="results/astro_transition_deep_raw.csv", help="Raw runs CSV")
    p.add_argument("--temp", type=float, required=True, help="Surface temperature to select")
    p.add_argument("--nH", type=float, required=True, help="Gas density to select")
    p.add_argument("--metric", default="epsilon", help="Per-run metric column to track")
    p.add_argument("--run-points", nargs="*", type=int, default=None, help="Specific ensemble sizes to sample")
    p.add_argument("--out-dir", default="results/plots/astro_transition_deep", help="Output plot directory")
    p.add_argument("--out-table", default=None, help="Optional output CSV for sampled convergence values")
    args = p.parse_args()

    df = pd.read_csv(args.raw)
    need = {"surface_temperature_k", "h_gas_density_cm3", args.metric}
    missing = sorted(need - set(df.columns))
    if missing:
        raise SystemExit(f"Missing columns in {args.raw}: {missing}")

    block = df[
        (pd.to_numeric(df["surface_temperature_k"], errors="coerce") == float(args.temp))
        & (pd.to_numeric(df["h_gas_density_cm3"], errors="coerce") == float(args.nH))
    ].copy()
    if block.empty:
        raise SystemExit(f"No rows found for T={args.temp} K, nH={args.nH} cm^-3 in {args.raw}")

    values = pd.to_numeric(block[args.metric], errors="coerce").dropna().to_numpy(dtype=float)
    if values.size == 0:
        raise SystemExit(f"No numeric values found in column '{args.metric}' for selected condition")

    cumulative_mean = np.cumsum(values) / np.arange(1, values.size + 1, dtype=float)
    cumulative_sem = np.array([_sem(values[:i]) for i in range(1, values.size + 1)], dtype=float)
    run_points = _pick_points(values.size, args.run_points)

    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)
    if args.out_table:
        os.makedirs(os.path.dirname(args.out_table) or ".", exist_ok=True)

    apply_publication_style()
    fig, ax = plt.subplots(figsize=(SINGLE_COLUMN_WIDTH, 2.65))
    x = np.arange(1, values.size + 1, dtype=int)
    ax.plot(x, cumulative_mean, color="#1d3557", linewidth=2.3, label="Running mean")
    ax.fill_between(
        x,
        cumulative_mean - cumulative_sem,
        cumulative_mean + cumulative_sem,
        color="#1d3557",
        alpha=0.18,
        label="±1 SEM",
    )
    ax.scatter(run_points, cumulative_mean[np.array(run_points) - 1], color="#c05640", zorder=3, s=24, label="Highlighted N")
    ax.axvline(20, color="#666666", lw=0.8, ls="--", alpha=0.8)
    style_axes(ax)
    ax.set_xlabel("Number of ensemble realizations")
    ax.set_ylabel(args.metric)
    ax.set_xscale("log")
    style_legend(ax, loc="best")
    fig.subplots_adjust(left=0.2, right=0.98, bottom=0.18, top=0.98)

    metric_slug = args.metric.replace("/", "_")
    stem = f"ensemble_convergence_{metric_slug}_T{args.temp:g}_nH{args.nH:g}".replace(".", "p")
    out_png = os.path.join(out_dir, stem)
    save_figure(fig, out_png)
    plt.close(fig)
    print(f"Wrote {out_png}.png/.pdf")

    sampled = pd.DataFrame(
        {
            "n_runs": run_points,
            f"{args.metric}_mean": [float(cumulative_mean[n - 1]) for n in run_points],
            f"{args.metric}_sem": [float(cumulative_sem[n - 1]) for n in run_points],
        }
    )

    out_table = args.out_table or os.path.join(out_dir, f"{stem}.csv")
    sampled.to_csv(out_table, index=False)
    print(f"Wrote {out_table}")


if __name__ == "__main__":
    main()
