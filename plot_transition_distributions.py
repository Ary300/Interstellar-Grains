#!/usr/bin/env python3
"""
Plot per-run distributions for a deep transition-region campaign.

Intended input is the *raw runs* CSV from run_sweep.py (save_raw_runs: true),
e.g.:

  python run_sweep.py config_astro_transition_deep.yaml

Then:

  python plot_transition_distributions.py --raw results/astro_transition_deep_raw.csv

Outputs:
  - results/plots/astro_transition_deep/epsilon_boxplot.png
  - results/plots/astro_transition_deep/h2_release_rate_boxplot.png
  - results/tables/astro_transition_deep_distribution_summary.csv
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from paper_plot_style import SINGLE_COLUMN_WIDTH, apply_publication_style, save_figure, style_axes, style_legend


def _ensure_out_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def _ci95_halfwidth(x: pd.Series) -> float:
    x = pd.to_numeric(x, errors="coerce").dropna()
    if len(x) <= 1:
        return float("nan")
    return float(1.96 * x.std(ddof=1) / np.sqrt(len(x)))


def main() -> None:
    p = argparse.ArgumentParser(description="Plot raw-run distributions for the 80–120 K transition campaign.")
    p.add_argument("--raw", default="results/astro_transition_deep_raw.csv", help="Raw runs CSV from run_sweep.py")
    p.add_argument("--out-plots", default="results/plots/astro_transition_deep", help="Output directory for plots")
    p.add_argument(
        "--out-table",
        default="results/tables/astro_transition_deep_distribution_summary.csv",
        help="Output CSV for grouped distribution summary",
    )
    args = p.parse_args()

    df = pd.read_csv(args.raw)
    needed = {"surface_temperature_k", "h_gas_density_cm3", "epsilon", "h2_release_rate_cm2_s"}
    missing = sorted(needed - set(df.columns))
    if missing:
        raise SystemExit(f"Missing required columns in {args.raw}: {missing}")

    _ensure_out_dir(args.out_plots)
    _ensure_out_dir(os.path.dirname(args.out_table) or ".")

    df["surface_temperature_k"] = pd.to_numeric(df["surface_temperature_k"], errors="coerce")
    df["h_gas_density_cm3"] = pd.to_numeric(df["h_gas_density_cm3"], errors="coerce")
    df["epsilon"] = pd.to_numeric(df["epsilon"], errors="coerce")
    df["h2_release_rate_cm2_s"] = pd.to_numeric(df["h2_release_rate_cm2_s"], errors="coerce")

    # Grouped summary table.
    gcols = ["surface_temperature_k", "h_gas_density_cm3"]
    summary = (
        df.groupby(gcols, dropna=False)
        .agg(
            n_runs=("epsilon", "count"),
            epsilon_mean=("epsilon", "mean"),
            epsilon_std=("epsilon", "std"),
            epsilon_ci95=("epsilon", _ci95_halfwidth),
            h2_release_rate_cm2_s_mean=("h2_release_rate_cm2_s", "mean"),
            h2_release_rate_cm2_s_std=("h2_release_rate_cm2_s", "std"),
            h2_release_rate_cm2_s_ci95=("h2_release_rate_cm2_s", _ci95_halfwidth),
        )
        .reset_index()
        .sort_values(gcols)
    )
    summary.to_csv(args.out_table, index=False)
    print(f"Wrote {args.out_table}")

    apply_publication_style()

    df_plot = df.dropna(subset=["surface_temperature_k", "h_gas_density_cm3"])
    df_plot["surface_temperature_k"] = df_plot["surface_temperature_k"].astype(float)
    df_plot["h_gas_density_cm3"] = df_plot["h_gas_density_cm3"].astype(float)
    df_plot["h_gas_density_cm3"] = df_plot["h_gas_density_cm3"].astype(int)

    # ε distribution
    fig, ax = plt.subplots(figsize=(SINGLE_COLUMN_WIDTH, 2.8))
    sns.boxplot(
        data=df_plot,
        x="surface_temperature_k",
        y="epsilon",
        hue="h_gas_density_cm3",
        palette=["#1d3557", "#c05640"],
        showfliers=False,
        ax=ax,
    )
    style_axes(ax)
    ax.set_xlabel("Surface temperature (K)")
    ax.set_ylabel(r"Per-run efficiency $\epsilon$")
    style_legend(ax, title=r"$n_{\rm H}$ (cm$^{-3}$)", loc="upper right")
    fig.subplots_adjust(left=0.2, right=0.98, bottom=0.16, top=0.98)
    out_eps = os.path.join(args.out_plots, "epsilon_boxplot")
    save_figure(fig, out_eps)
    plt.close(fig)
    print(f"Wrote {out_eps}.png/.pdf")

    # H2 release rate distribution
    fig, ax = plt.subplots(figsize=(SINGLE_COLUMN_WIDTH, 2.8))
    sns.boxplot(
        data=df_plot,
        x="surface_temperature_k",
        y="h2_release_rate_cm2_s",
        hue="h_gas_density_cm3",
        palette=["#1d3557", "#c05640"],
        showfliers=False,
        ax=ax,
    )
    style_axes(ax)
    ax.set_xlabel("Surface temperature (K)")
    ax.set_ylabel(r"Per-run H$_2$ release rate (cm$^{-2}$ s$^{-1}$)")
    style_legend(ax, title=r"$n_{\rm H}$ (cm$^{-3}$)", loc="upper right")
    fig.subplots_adjust(left=0.2, right=0.98, bottom=0.16, top=0.98)
    out_rate = os.path.join(args.out_plots, "h2_release_rate_boxplot")
    save_figure(fig, out_rate)
    plt.close(fig)
    print(f"Wrote {out_rate}.png/.pdf")


if __name__ == "__main__":
    main()
