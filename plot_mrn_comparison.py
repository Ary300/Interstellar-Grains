#!/usr/bin/env python3
"""
Compare single-grain and MRN-integrated paperfit outputs.
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from paper_plot_style import SINGLE_COLUMN_WIDTH, add_panel_labels, add_regime_band, apply_publication_style, save_figure, style_axes


def _pick_col(df: pd.DataFrame, base: str) -> str:
    if base in df.columns:
        return base
    mean_col = f"{base}_mean"
    if mean_col in df.columns:
        return mean_col
    raise KeyError(f"Missing column: {base} or {mean_col}")


def main() -> None:
    p = argparse.ArgumentParser(description="Compare single-grain and MRN-integrated ISM outputs.")
    p.add_argument("--single-input", default="results/jhub_full_merged.csv")
    p.add_argument("--mrn-input", default="results/astro_mrn_integration.csv")
    p.add_argument("--nH", type=float, default=1000.0, help="Density slice to compare")
    p.add_argument("--uv", type=float, default=0.0, help="UV slice to compare")
    p.add_argument("--plot-min-temp", type=float, default=20.0, help="Minimum temperature to display in the figure")
    p.add_argument("--out-table", default="results/tables/astro_mrn_comparison.csv")
    p.add_argument("--out-dir", default="results/plots/astro_mrn_comparison")
    args = p.parse_args()

    single = pd.read_csv(args.single_input)
    mrn = pd.read_csv(args.mrn_input)
    eps_single = _pick_col(single, "epsilon")
    rate_single = _pick_col(single, "h2_release_rate_cm2_s")
    eps_mrn = _pick_col(mrn, "epsilon")
    rate_mrn = _pick_col(mrn, "h2_release_rate_cm2_s")

    def _slice(df: pd.DataFrame, eps_col: str, rate_col: str) -> pd.DataFrame:
        out = df.copy()
        for c in ["surface_temperature_k", "h_gas_density_cm3", "uv_flux_factor", eps_col, rate_col]:
            out[c] = pd.to_numeric(out[c], errors="coerce")
        out = out.dropna(subset=["surface_temperature_k", "h_gas_density_cm3", "uv_flux_factor", eps_col, rate_col])
        out = out[np.isclose(out["h_gas_density_cm3"], float(args.nH)) & np.isclose(out["uv_flux_factor"], float(args.uv))]
        return out[["surface_temperature_k", eps_col, rate_col]].rename(columns={eps_col: "epsilon", rate_col: "h2_release_rate_cm2_s"})

    single_slice = _slice(single, eps_single, rate_single).rename(columns={"epsilon": "epsilon_single", "h2_release_rate_cm2_s": "rate_single"})
    mrn_slice = _slice(mrn, eps_mrn, rate_mrn).rename(columns={"epsilon": "epsilon_mrn", "h2_release_rate_cm2_s": "rate_mrn"})
    merged = pd.merge(single_slice, mrn_slice, on="surface_temperature_k", how="inner").sort_values("surface_temperature_k")
    if merged.empty:
        raise SystemExit("No overlapping single-grain / MRN rows for the requested nH and UV slice.")

    merged["epsilon_ratio_mrn_over_single"] = merged["epsilon_mrn"] / merged["epsilon_single"]
    merged["rate_ratio_mrn_over_single"] = merged["rate_mrn"] / merged["rate_single"]

    os.makedirs(os.path.dirname(args.out_table) or ".", exist_ok=True)
    merged.to_csv(args.out_table, index=False)
    print(f"Wrote {args.out_table}")

    os.makedirs(args.out_dir, exist_ok=True)
    plot_df = merged[merged["surface_temperature_k"] >= float(args.plot_min_temp)].copy()
    if plot_df.empty:
        raise SystemExit("No rows remain after applying --plot-min-temp to the MRN comparison.")

    apply_publication_style()

    fig, axes = plt.subplots(2, 1, figsize=(SINGLE_COLUMN_WIDTH, 4.1), sharex=True)
    axes[0].plot(plot_df["surface_temperature_k"], plot_df["epsilon_ratio_mrn_over_single"], marker="o", ms=5, lw=2.2, color="#1d3557")
    axes[0].axhline(1.0, color="k", linewidth=1.0, alpha=0.5)
    axes[0].set_ylabel(r"$\epsilon_{\rm MRN}/\epsilon_{\rm single}$")
    style_axes(axes[0])
    axes[0].grid(False)
    add_regime_band(axes[0], 150.0, 250.0, "Warm regime", color="#f8e1d8", alpha=0.08)

    axes[1].plot(plot_df["surface_temperature_k"], plot_df["rate_ratio_mrn_over_single"], marker="o", ms=5, lw=2.2, color="#c05640")
    axes[1].axhline(1.0, color="k", linewidth=1.0, alpha=0.5)
    axes[1].set_xlabel("Surface temperature (K)")
    axes[1].set_ylabel(r"$R_{\rm MRN}/R_{\rm single}$")
    style_axes(axes[1])
    axes[1].grid(False)
    add_regime_band(axes[1], 150.0, 250.0, "Warm regime", color="#f8e1d8", alpha=0.08)
    add_panel_labels(axes, labels="AB")
    fig.text(
        0.985,
        0.99,
        fr"$n_{{\rm H}}={args.nH:g}$ cm$^{{-3}}$, $G_0={args.uv:g}$",
        ha="right",
        va="top",
        fontsize=6.9,
        color="#4a4a4a",
    )
    fig.subplots_adjust(left=0.18, right=0.98, bottom=0.12, top=0.96, hspace=0.12)
    out_png = os.path.join(args.out_dir, "mrn_correction_factors")
    save_figure(fig, out_png)
    plt.close(fig)
    print(f"Wrote {out_png}.png/.pdf")


if __name__ == "__main__":
    main()
