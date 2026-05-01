#!/usr/bin/env python3
"""
Create a polished multi-panel summary figure from the current ISM baseline run.
"""

from __future__ import annotations

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from paper_plot_style import (
    DOUBLE_COLUMN_WIDTH,
    add_panel_labels,
    add_regime_band,
    apply_publication_style,
    density_style,
    save_figure,
    style_axes,
    style_legend,
)


def _pick_col(df: pd.DataFrame, base: str) -> str:
    if base in df.columns:
        return base
    mean_col = f"{base}_mean"
    if mean_col in df.columns:
        return mean_col
    raise KeyError(f"Missing column: {base} or {mean_col}")


def _ci_col(df: pd.DataFrame, base: str) -> str | None:
    col = f"{base}_ci95"
    return col if col in df.columns else None


def main() -> None:
    p = argparse.ArgumentParser(description="Create polished paper figures from the baseline ISM outputs.")
    p.add_argument("--input", default="results/jhub_full_merged.csv")
    p.add_argument("--ct-input", default="results/tables/ct02_comparison_jhub_full.csv")
    p.add_argument("--uv", type=float, default=0.0, help="UV slice to plot")
    p.add_argument("--mechanism-nh", type=float, default=1000.0, help="Density slice for mechanism decomposition")
    p.add_argument("--out", default="results/plots/paper/ism_main_results")
    args = p.parse_args()

    apply_publication_style()

    df = pd.read_csv(args.input)
    eps_col = _pick_col(df, "epsilon")
    rate_col = _pick_col(df, "h2_release_rate_cm2_s")
    eps_ci = _ci_col(df, "epsilon")
    rate_ci = _ci_col(df, "h2_release_rate_cm2_s")

    for col in [
        "surface_temperature_k",
        "h_gas_density_cm3",
        "uv_flux_factor",
        eps_col,
        rate_col,
        "h2_formed_LH_mean",
        "h2_formed_ER_mean",
    ]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df[np.isclose(df["uv_flux_factor"], float(args.uv))].copy()
    densities = sorted(df["h_gas_density_cm3"].dropna().unique())

    ct = pd.read_csv(args.ct_input)
    for col in ["surface_temperature_k", "h_gas_density_cm3", "uv_flux_factor", "ratio_kmc_over_ct02"]:
        ct[col] = pd.to_numeric(ct[col], errors="coerce")
    ct = ct[np.isclose(ct["uv_flux_factor"], float(args.uv))].copy()

    fig, axes = plt.subplots(2, 2, figsize=(DOUBLE_COLUMN_WIDTH, 5.05))
    ax_eps, ax_rate, ax_mech, ax_ct = axes.flatten()

    for n_h in densities:
        sub = df[np.isclose(df["h_gas_density_cm3"], float(n_h))].sort_values("surface_temperature_k")
        style = density_style(float(n_h))
        color = style["color"]
        marker = style["marker"]
        linestyle = style["linestyle"]
        ax_eps.plot(
            sub["surface_temperature_k"],
            sub[eps_col],
            color=color,
            lw=2.2,
            marker=marker,
            ms=4.8,
            linestyle=linestyle,
            label=fr"$n_{{\rm H}}={int(n_h):,}$",
        )
        if eps_ci:
            ci = pd.to_numeric(sub[eps_ci], errors="coerce").fillna(0.0)
            ax_eps.fill_between(sub["surface_temperature_k"], sub[eps_col] - ci, sub[eps_col] + ci, color=color, alpha=0.13)

        ax_rate.plot(
            sub["surface_temperature_k"],
            sub[rate_col],
            color=color,
            lw=2.2,
            marker=marker,
            ms=4.8,
            linestyle=linestyle,
            label=fr"$n_{{\rm H}}={int(n_h):,}$",
        )
        if rate_ci:
            ci = pd.to_numeric(sub[rate_ci], errors="coerce").fillna(0.0)
            ax_rate.fill_between(sub["surface_temperature_k"], sub[rate_col] - ci, sub[rate_col] + ci, color=color, alpha=0.13)

    mech = df[np.isclose(df["h_gas_density_cm3"], float(args.mechanism_nh))].sort_values("surface_temperature_k").copy()
    total = mech["h2_formed_LH_mean"].fillna(0.0) + mech["h2_formed_ER_mean"].fillna(0.0)
    total = total.replace(0.0, np.nan)
    mech["lh_frac"] = mech["h2_formed_LH_mean"] / total
    mech["er_frac"] = mech["h2_formed_ER_mean"] / total
    ax_mech.plot(mech["surface_temperature_k"], mech["lh_frac"], color="#2a9d8f", lw=2.4, marker="o", ms=4.5, label="LH fraction")
    ax_mech.plot(mech["surface_temperature_k"], mech["er_frac"], color="#c05640", lw=2.4, marker="s", ms=4.5, label="ER fraction")
    ax_mech.fill_between(mech["surface_temperature_k"], 0.0, mech["lh_frac"], color="#2a9d8f", alpha=0.10)
    ax_mech.fill_between(mech["surface_temperature_k"], mech["lh_frac"], mech["lh_frac"] + mech["er_frac"], color="#c05640", alpha=0.08)

    for n_h in densities:
        sub = ct[np.isclose(ct["h_gas_density_cm3"], float(n_h))].sort_values("surface_temperature_k")
        if sub.empty:
            continue
        style = density_style(float(n_h))
        ax_ct.plot(
            sub["surface_temperature_k"],
            sub["ratio_kmc_over_ct02"],
            color=style["color"],
            lw=2.1,
            marker=style["marker"],
            ms=4.5,
            linestyle=style["linestyle"],
            label=fr"$n_{{\rm H}}={int(n_h):,}$",
        )
    ax_ct.axhline(1.0, color="#555555", lw=1.0, linestyle="--", alpha=0.8)

    for ax in [ax_eps, ax_rate, ax_mech, ax_ct]:
        style_axes(ax)
        ax.set_xlim(8, 252)
        ax.grid(False)

    add_regime_band(ax_eps, 20, 100, "LH window", color="#fff4bf", alpha=0.07)
    add_regime_band(ax_eps, 150, 250, "ER plateau", color="#f8e1d8", alpha=0.07)
    add_regime_band(ax_rate, 150, 250, "ER plateau", color="#f8e1d8", alpha=0.07)
    add_regime_band(ax_mech, 100, 120, "Transition", color="#b8d8ba", alpha=0.08)

    ax_eps.set_xlabel("Surface temperature (K)")
    ax_eps.set_ylabel(r"Formation efficiency $\epsilon$")
    ax_eps.set_ylim(0.0, 0.33)
    style_legend(ax_eps, loc="upper right")
    ax_eps.axhline(0.19, color="#444444", lw=0.8, ls="--", alpha=0.8)

    ax_rate.set_xlabel("Surface temperature (K)")
    ax_rate.set_ylabel(r"H$_2$ release rate (cm$^{-2}$ s$^{-1}$)")
    ax_rate.set_yscale("log")

    ax_mech.set_xlabel("Surface temperature (K)")
    ax_mech.set_ylabel(r"H$_2$ fraction")
    ax_mech.set_ylim(0.0, 1.05)
    style_legend(ax_mech, loc="center right")
    crossover = mech.iloc[np.abs(mech["lh_frac"] - mech["er_frac"]).to_numpy().argsort()[:1]]
    if not crossover.empty:
        t_cross = float(crossover["surface_temperature_k"].iloc[0])
        ax_mech.annotate(
            fr"{t_cross:.0f} K",
            xy=(t_cross, 0.5),
            xytext=(t_cross + 17.0, 0.23),
            arrowprops={"arrowstyle": "->", "lw": 0.9, "color": "#444444"},
            fontsize=7.2,
            color="#444444",
        )

    ax_ct.set_xlabel("Surface temperature (K)")
    ax_ct.set_ylabel("Rate ratio (KMC / CT10)")
    ax_ct.set_yscale("log")

    add_panel_labels([ax_eps, ax_rate, ax_mech, ax_ct], labels="ABCD")
    fig.subplots_adjust(left=0.09, right=0.985, bottom=0.1, top=0.98, wspace=0.22, hspace=0.3)
    save_figure(fig, args.out)
    plt.close(fig)
    print(f"Wrote {args.out}.png/.pdf")


if __name__ == "__main__":
    main()
