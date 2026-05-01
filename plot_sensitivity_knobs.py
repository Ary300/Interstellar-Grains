#!/usr/bin/env python3
"""
Plot a simple parameter-sensitivity envelope for ISM sweeps.

Designed for outputs from config_astro_sensitivity_knobs.yaml (or similar), where
chemisorption_fraction and er_reaction_probability are swept at fixed (nH, UV).

Outputs:
  - results/tables/astro_sensitivity_knobs_summary.csv
  - results/plots/astro_sensitivity_knobs/k_eff_envelope.png
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from paper_plot_style import SINGLE_COLUMN_WIDTH, apply_publication_style, save_figure, style_axes, style_legend


def _ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def _pick_col(df: pd.DataFrame, base: str) -> str:
    if base in df.columns:
        return base
    m = f"{base}_mean"
    if m in df.columns:
        return m
    raise KeyError(f"Missing column: {base} or {base}_mean")


def _has_positive(*arrays: np.ndarray) -> bool:
    for arr in arrays:
        vals = np.asarray(arr, dtype=float)
        if np.any(np.isfinite(vals) & (vals > 0.0)):
            return True
    return False


def main() -> None:
    p = argparse.ArgumentParser(description="Plot chemisorption/ER knob sensitivity envelope.")
    p.add_argument("--input", default="results/astro_sensitivity_knobs.csv", help="Aggregated sweep CSV")
    p.add_argument("--sigma-H", dest="sigma_H", type=float, default=1e-21, help="Cross-section per H [cm^2/H] used for k_eff")
    p.add_argument("--out-plots", default="results/plots/astro_sensitivity_knobs", help="Output directory for plots")
    p.add_argument("--out-table", default="results/tables/astro_sensitivity_knobs_summary.csv", help="Output CSV summary")
    p.add_argument("--baseline-fchem", type=float, default=0.4, help="Baseline chemisorption_fraction to highlight")
    p.add_argument("--baseline-per", type=float, default=0.9, help="Baseline er_reaction_probability to highlight")
    args = p.parse_args()

    df = pd.read_csv(args.input)

    need = {"surface_temperature_k", "h_gas_density_cm3", "uv_flux_factor", "chemisorption_fraction", "er_reaction_probability"}
    missing = sorted(need - set(df.columns))
    if missing:
        raise SystemExit(f"Missing required columns in {args.input}: {missing}")

    rate_col = _pick_col(df, "h2_release_rate_cm2_s")

    df["surface_temperature_k"] = pd.to_numeric(df["surface_temperature_k"], errors="coerce")
    df["h_gas_density_cm3"] = pd.to_numeric(df["h_gas_density_cm3"], errors="coerce")
    df["chemisorption_fraction"] = pd.to_numeric(df["chemisorption_fraction"], errors="coerce")
    df["er_reaction_probability"] = pd.to_numeric(df["er_reaction_probability"], errors="coerce")
    df[rate_col] = pd.to_numeric(df[rate_col], errors="coerce")
    df = df.dropna(subset=["surface_temperature_k", "h_gas_density_cm3", "chemisorption_fraction", "er_reaction_probability", rate_col])

    if df[["h_gas_density_cm3", "uv_flux_factor"]].drop_duplicates().shape[0] != 1:
        raise SystemExit("Expected a single (nH, uv) slice in the input. Filter your CSV or rerun with fixed nH/uv.")

    df["k_eff_cm3_s"] = 4.0 * float(args.sigma_H) * df[rate_col].astype(float) / df["h_gas_density_cm3"].astype(float)

    g = df.groupby("surface_temperature_k", dropna=False)
    summary = (
        g["k_eff_cm3_s"]
        .agg(k_eff_min="min", k_eff_median="median", k_eff_max="max", n="count")
        .reset_index()
        .sort_values("surface_temperature_k")
    )

    # Baseline line (if present)
    base = df[
        np.isclose(df["chemisorption_fraction"], float(args.baseline_fchem))
        & np.isclose(df["er_reaction_probability"], float(args.baseline_per))
    ][["surface_temperature_k", "k_eff_cm3_s"]].sort_values("surface_temperature_k")

    _ensure_dir(os.path.dirname(args.out_table) or ".")
    summary.to_csv(args.out_table, index=False)
    print(f"Wrote {args.out_table}")

    _ensure_dir(args.out_plots)
    out_png = os.path.join(args.out_plots, "k_eff_envelope.png")

    apply_publication_style()

    x = summary["surface_temperature_k"].astype(float).values
    y_min = summary["k_eff_min"].astype(float).values
    y_med = summary["k_eff_median"].astype(float).values
    y_max = summary["k_eff_max"].astype(float).values

    fig, ax = plt.subplots(figsize=(SINGLE_COLUMN_WIDTH, 2.85))
    ax.fill_between(x, y_min, y_max, color="#d9c27a", alpha=0.28, label="Sweep envelope")
    ax.plot(x, y_med, marker="o", ms=5, lw=2.2, color="#1d3557", label="Median response")
    if not base.empty:
        ax.plot(base["surface_temperature_k"], base["k_eff_cm3_s"], marker="s", ms=5, lw=2.0, linestyle="--", color="#c05640", label="Paperfit baseline")
    style_axes(ax)
    if _has_positive(y_min, y_med, y_max, base["k_eff_cm3_s"].to_numpy(dtype=float) if not base.empty else np.array([])):
        floor = min(v for v in np.concatenate([y_min, y_med, y_max, base["k_eff_cm3_s"].to_numpy(dtype=float) if not base.empty else np.array([])]) if np.isfinite(v) and v > 0.0)
        ax.set_yscale("symlog", linthresh=max(floor / 5.0, 1e-30))
    ax.set_xlabel("Surface temperature (K)")
    ax.set_ylabel(r"Effective $k_{\rm eff}$ (cm$^3$ s$^{-1}$)")
    style_legend(ax, loc="best")
    fig.subplots_adjust(left=0.2, right=0.98, bottom=0.18, top=0.98)
    save_figure(fig, os.path.splitext(out_png)[0])
    plt.close(fig)
    print(f"Wrote {os.path.splitext(out_png)[0]}.png/.pdf")


if __name__ == "__main__":
    main()
