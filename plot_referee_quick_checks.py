#!/usr/bin/env python3
"""
Summarize two small referee-facing checks:

1) Grain-size dependence at fixed microphysics
2) Constant vs temperature-dependent sticking model
3) Porosity consistency (0.0 vs 0.2)

Each check can be run independently by providing only the relevant input.
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


def _pick_metric(df: pd.DataFrame, base: str) -> str:
    if base in df.columns:
        return base
    mean_col = f"{base}_mean"
    if mean_col in df.columns:
        return mean_col
    raise KeyError(f"Missing column: {base} or {base}_mean")


def _plot_grain(grain_input: str, out_dir: str) -> None:
    df = pd.read_csv(grain_input)
    eps_col = _pick_metric(df, "epsilon")
    rate_col = _pick_metric(df, "h2_release_rate_cm2_s")

    need = {"surface_temperature_k", "grain_radius_um", eps_col, rate_col}
    missing = sorted(need - set(df.columns))
    if missing:
        raise SystemExit(f"Missing columns in {grain_input}: {missing}")

    df = df.copy()
    df["surface_temperature_k"] = pd.to_numeric(df["surface_temperature_k"], errors="coerce")
    df["grain_radius_um"] = pd.to_numeric(df["grain_radius_um"], errors="coerce")
    df[eps_col] = pd.to_numeric(df[eps_col], errors="coerce")
    df[rate_col] = pd.to_numeric(df[rate_col], errors="coerce")
    df = df.dropna(subset=["surface_temperature_k", "grain_radius_um", eps_col, rate_col])

    baseline_radius = 0.005
    base = (
        df[np.isclose(df["grain_radius_um"], baseline_radius)]
        .set_index("surface_temperature_k")[[eps_col, rate_col]]
        .rename(columns={eps_col: "epsilon_baseline", rate_col: "rate_baseline"})
    )
    merged = df.join(base, on="surface_temperature_k", how="left")
    merged["epsilon_ratio_to_0p005um"] = merged[eps_col] / merged["epsilon_baseline"]
    merged["rate_ratio_to_0p005um"] = merged[rate_col] / merged["rate_baseline"]
    summary_path = os.path.join(out_dir, "grain_size_summary.csv")
    merged.sort_values(["surface_temperature_k", "grain_radius_um"]).to_csv(summary_path, index=False)

    apply_publication_style()

    for metric_col, ylabel, out_name in [
        (eps_col, "Recombination efficiency ε", "grain_size_epsilon.png"),
        (rate_col, r"H$_2$ release rate (cm$^{-2}$ s$^{-1}$)", "grain_size_rate.png"),
    ]:
        fig, ax = plt.subplots(figsize=(SINGLE_COLUMN_WIDTH, 2.55))
        for temp, group in merged.groupby("surface_temperature_k", dropna=False):
            group = group.sort_values("grain_radius_um")
            ax.plot(group["grain_radius_um"], group[metric_col], marker="o", ms=5, lw=2.0, label=f"{temp:.0f} K")
        ax.set_xscale("log")
        style_axes(ax)
        ax.set_xlabel("Grain radius (μm)")
        ax.set_ylabel(ylabel)
        style_legend(ax, title="Surface T", loc="best")
        fig.subplots_adjust(left=0.2, right=0.98, bottom=0.18, top=0.98)
        save_figure(fig, os.path.join(out_dir, os.path.splitext(out_name)[0]))
        plt.close(fig)

    print(f"Wrote {summary_path}")


def _plot_sticking(sticking_input: str, out_dir: str) -> None:
    df = pd.read_csv(sticking_input)
    eps_col = _pick_metric(df, "epsilon")
    rate_col = _pick_metric(df, "h2_release_rate_cm2_s")

    need = {"surface_temperature_k", "gas_temperature_k", "sticking_temp_model", eps_col, rate_col}
    missing = sorted(need - set(df.columns))
    if missing:
        raise SystemExit(f"Missing columns in {sticking_input}: {missing}")

    df = df.copy()
    for c in ["surface_temperature_k", "gas_temperature_k", eps_col, rate_col]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["surface_temperature_k", "gas_temperature_k", eps_col, rate_col])

    pivot_eps = df.pivot_table(
        index=["surface_temperature_k", "gas_temperature_k"],
        columns="sticking_temp_model",
        values=eps_col,
        aggfunc="first",
    )
    pivot_rate = df.pivot_table(
        index=["surface_temperature_k", "gas_temperature_k"],
        columns="sticking_temp_model",
        values=rate_col,
        aggfunc="first",
    )

    if "constant" not in pivot_eps.columns:
        raise SystemExit("Expected sticking_temp_model='constant' in sticking check CSV")

    alt_eps_cols = [c for c in pivot_eps.columns if c != "constant"]
    alt_rate_cols = [c for c in pivot_rate.columns if c != "constant"]
    if not alt_eps_cols or not alt_rate_cols:
        raise SystemExit("Expected at least one non-constant sticking model in sticking check CSV")

    alt_eps = alt_eps_cols[0]
    alt_rate = alt_rate_cols[0]

    summary = pivot_eps[["constant", alt_eps]].rename(columns={alt_eps: "epsilon_alt"})
    summary["epsilon_frac_change"] = summary["epsilon_alt"] / summary["constant"] - 1.0
    summary = summary.join(
        pivot_rate[["constant", alt_rate]].rename(columns={"constant": "rate_constant", alt_rate: "rate_alt"})
    )
    summary["rate_frac_change"] = summary["rate_alt"] / summary["rate_constant"] - 1.0
    summary = summary.reset_index()
    summary_path = os.path.join(out_dir, "sticking_model_summary.csv")
    summary.to_csv(summary_path, index=False)

    xlabels = [f"Ts={ts:.0f}, Tg={tg:.0f}" for ts, tg in zip(summary["surface_temperature_k"], summary["gas_temperature_k"])]
    x = np.arange(len(summary))
    width = 0.35

    apply_publication_style()
    fig, ax = plt.subplots(figsize=(SINGLE_COLUMN_WIDTH, 2.75))
    ax.bar(x - width / 2, summary["epsilon_frac_change"], width=width, color="#1d3557", label=r"$\epsilon$")
    ax.bar(x + width / 2, summary["rate_frac_change"], width=width, color="#c05640", label="Rate")
    ax.axhline(0.0, color="k", linewidth=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(xlabels, rotation=20, ha="right")
    style_axes(ax)
    ax.set_ylabel("Fractional change vs constant sticking")
    style_legend(ax, loc="best")
    fig.subplots_adjust(left=0.2, right=0.98, bottom=0.3, top=0.98)
    save_figure(fig, os.path.join(out_dir, "sticking_model_fractional_change"))
    plt.close(fig)

    print(f"Wrote {summary_path}")


def _plot_porosity(porosity_input: str, out_dir: str) -> None:
    df = pd.read_csv(porosity_input)
    eps_col = _pick_metric(df, "epsilon")
    rate_col = _pick_metric(df, "h2_release_rate_cm2_s")

    need = {"surface_temperature_k", "h_gas_density_cm3", "porosity_fraction", eps_col, rate_col}
    missing = sorted(need - set(df.columns))
    if missing:
        raise SystemExit(f"Missing columns in {porosity_input}: {missing}")

    df = df.copy()
    for c in ["surface_temperature_k", "h_gas_density_cm3", "porosity_fraction", eps_col, rate_col]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["surface_temperature_k", "h_gas_density_cm3", "porosity_fraction", eps_col, rate_col])

    pivot_eps = df.pivot_table(
        index=["surface_temperature_k", "h_gas_density_cm3"],
        columns="porosity_fraction",
        values=eps_col,
        aggfunc="first",
    )
    pivot_rate = df.pivot_table(
        index=["surface_temperature_k", "h_gas_density_cm3"],
        columns="porosity_fraction",
        values=rate_col,
        aggfunc="first",
    )

    if 0.0 not in pivot_eps.columns or 0.2 not in pivot_eps.columns:
        raise SystemExit("Expected porosity_fraction values 0.0 and 0.2 in porosity check CSV")

    summary = pivot_eps[[0.0, 0.2]].rename(columns={0.0: "epsilon_porosity_0p0", 0.2: "epsilon_porosity_0p2"})
    summary["epsilon_frac_change"] = summary["epsilon_porosity_0p2"] / summary["epsilon_porosity_0p0"] - 1.0
    summary = summary.join(
        pivot_rate[[0.0, 0.2]].rename(columns={0.0: "rate_porosity_0p0", 0.2: "rate_porosity_0p2"})
    )
    summary["rate_frac_change"] = summary["rate_porosity_0p2"] / summary["rate_porosity_0p0"] - 1.0
    summary = summary.reset_index()
    summary_path = os.path.join(out_dir, "porosity_summary.csv")
    summary.to_csv(summary_path, index=False)
    print(f"Wrote {summary_path}")

    apply_publication_style()
    fig, ax = plt.subplots(figsize=(SINGLE_COLUMN_WIDTH, 2.65))
    labels = [f"T={t:.0f}, n={n:.0f}" for t, n in zip(summary["surface_temperature_k"], summary["h_gas_density_cm3"])]
    x = np.arange(len(summary))
    width = 0.35
    ax.bar(x - width / 2, summary["epsilon_frac_change"], width=width, color="#1d3557", label=r"$\epsilon$")
    ax.bar(x + width / 2, summary["rate_frac_change"], width=width, color="#c05640", label="Rate")
    ax.axhline(0.0, color="k", linewidth=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=18, ha="right")
    style_axes(ax)
    ax.set_ylabel("Fractional change from porosity 0.0 to 0.2")
    style_legend(ax, loc="best")
    fig.subplots_adjust(left=0.2, right=0.98, bottom=0.26, top=0.98)
    save_figure(fig, os.path.join(out_dir, "porosity_fractional_change"))
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser(description="Plot quick referee-facing sensitivity checks.")
    p.add_argument("--grain-input", default=None, help="Aggregated CSV from config_astro_grain_size_check.yaml")
    p.add_argument("--sticking-input", default=None, help="Aggregated CSV from config_astro_sticking_model_check.yaml")
    p.add_argument("--porosity-input", default=None, help="Aggregated CSV from config_astro_porosity_check.yaml")
    p.add_argument("--out-dir", default="results/plots/referee_quick_checks", help="Output directory")
    args = p.parse_args()

    _ensure_dir(args.out_dir)
    if not args.grain_input and not args.sticking_input and not args.porosity_input:
        raise SystemExit("Provide at least one of --grain-input, --sticking-input, or --porosity-input")

    if args.grain_input:
        _plot_grain(args.grain_input, args.out_dir)
    if args.sticking_input:
        _plot_sticking(args.sticking_input, args.out_dir)
    if args.porosity_input:
        _plot_porosity(args.porosity_input, args.out_dir)


if __name__ == "__main__":
    main()
