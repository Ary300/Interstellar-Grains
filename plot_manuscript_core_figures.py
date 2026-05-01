#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from paper_plot_style import (
    DOUBLE_COLUMN_WIDTH,
    OKABE_ITO,
    SINGLE_COLUMN_WIDTH,
    add_panel_labels,
    add_regime_band,
    apply_publication_style,
    density_style,
    interpolate_phase_map,
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


def _maybe_ci(df: pd.DataFrame, base: str) -> str | None:
    name = f"{base}_ci95"
    return name if name in df.columns else None


def _read_digitized(path: str, method: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df[df["method"] == method].copy()
    for col in ["T_K", "epsilon", "epsilon_err_low", "epsilon_err_high"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df.dropna(subset=["T_K", "epsilon"]).sort_values("T_K")


def _axis_edges(values: np.ndarray, *, log: bool = False) -> np.ndarray:
    vals = np.asarray(sorted(values), dtype=float)
    if vals.size == 1:
        delta = 0.5 * vals[0] if vals[0] else 0.5
        return np.array([vals[0] - delta, vals[0] + delta], dtype=float)
    if log:
        lv = np.log10(vals)
        mids = 0.5 * (lv[:-1] + lv[1:])
        first = lv[0] - (mids[0] - lv[0])
        last = lv[-1] + (lv[-1] - mids[-1])
        return 10 ** np.concatenate([[first], mids, [last]])
    mids = 0.5 * (vals[:-1] + vals[1:])
    first = vals[0] - (mids[0] - vals[0])
    last = vals[-1] + (vals[-1] - mids[-1])
    return np.concatenate([[first], mids, [last]])


def _load_baseline(path: str, uv: float) -> tuple[pd.DataFrame, str, str, str | None, str | None]:
    df = pd.read_csv(path)
    eps_col = _pick_col(df, "epsilon")
    rate_col = _pick_col(df, "h2_release_rate_cm2_s")
    eps_ci = _maybe_ci(df, "epsilon")
    rate_ci = _maybe_ci(df, "h2_release_rate_cm2_s")
    numeric_cols = [
        "surface_temperature_k",
        "h_gas_density_cm3",
        "uv_flux_factor",
        eps_col,
        rate_col,
        "h2_formed_LH_mean",
        "h2_formed_ER_mean",
        "final_h_atoms_on_surface_mean",
    ]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    if eps_ci:
        df[eps_ci] = pd.to_numeric(df[eps_ci], errors="coerce")
    if rate_ci:
        df[rate_ci] = pd.to_numeric(df[rate_ci], errors="coerce")
    df = df[np.isclose(df["uv_flux_factor"], float(uv))].copy()
    return df, eps_col, rate_col, eps_ci, rate_ci


def make_fig04_validation(iso_path: str, digitized_path: str, out_dir: Path) -> None:
    iso = pd.read_csv(iso_path)
    for col in ["surface_temperature_k", "epsilon_mean", "epsilon_ci95"]:
        iso[col] = pd.to_numeric(iso[col], errors="coerce")
    iso = iso.dropna(subset=["surface_temperature_k", "epsilon_mean"]).sort_values("surface_temperature_k")
    expt = _read_digitized(digitized_path, "isothermal")

    model_at_expt = np.interp(expt["T_K"], iso["surface_temperature_k"], iso["epsilon_mean"])
    sigma = 0.5 * (expt["epsilon_err_low"].to_numpy(dtype=float) + expt["epsilon_err_high"].to_numpy(dtype=float))
    sigma[sigma <= 0.0] = np.nanmedian(sigma[sigma > 0.0])
    residual = model_at_expt - expt["epsilon"].to_numpy(dtype=float)
    warm_mask = expt["T_K"].between(100.0, 250.0).to_numpy(dtype=bool)
    dof = max(int(np.sum(warm_mask)) - 1, 1)
    chi2_red = float(np.nansum((residual[warm_mask] / sigma[warm_mask]) ** 2) / dof)

    fig, (ax1, ax2) = plt.subplots(
        1,
        2,
        figsize=(DOUBLE_COLUMN_WIDTH, 3.15),
        gridspec_kw={"width_ratios": [2.2, 1.0], "wspace": 0.24},
    )

    ax1.plot(iso["surface_temperature_k"], iso["epsilon_mean"], color=OKABE_ITO["blue"], lw=1.9, label="This work (KMC)")
    if "epsilon_ci95" in iso.columns:
        ax1.fill_between(
            iso["surface_temperature_k"],
            iso["epsilon_mean"] - iso["epsilon_ci95"],
            iso["epsilon_mean"] + iso["epsilon_ci95"],
            color=OKABE_ITO["blue"],
            alpha=0.22,
            lw=0,
        )
    ax1.errorbar(
        expt["T_K"],
        expt["epsilon"],
        yerr=[expt["epsilon_err_low"], expt["epsilon_err_high"]],
        fmt="o",
        mfc="white",
        mec="black",
        mew=1.0,
        ms=5.0,
        capsize=2,
        elinewidth=0.8,
        color="black",
        label="Grieco et al. (2023)",
        zorder=4,
    )
    style_axes(ax1)
    ax1.set_xlim(0, 270)
    ax1.set_ylim(0.0, 0.35)
    ax1.set_xlabel("Grain temperature (K)")
    ax1.set_ylabel(r"Formation efficiency $\epsilon$")
    ax1.text(
        0.97,
        0.08,
        fr"$\chi^2_{{\rm red}}={chi2_red:.3f}$" + "\n(100-250 K)",
        transform=ax1.transAxes,
        ha="right",
        va="bottom",
        fontsize=7.4,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 1.5},
    )
    style_legend(ax1, loc="upper right")

    ax2.axhline(0.0, color="black", lw=0.8, zorder=0)
    colors = np.where(residual >= 0.0, OKABE_ITO["blue"], OKABE_ITO["vermillion"])
    for temp, resid, err, color in zip(expt["T_K"], residual, sigma, colors):
        ax2.vlines(float(temp), 0.0, float(resid), color=color, lw=1.1)
        ax2.errorbar(
            [float(temp)],
            [float(resid)],
            yerr=[[float(err)], [float(err)]],
            fmt="o",
            mfc="white",
            mec="black",
            mew=1.0,
            ms=5.0,
            capsize=2,
            elinewidth=0.8,
            color="black",
        )
    style_axes(ax2)
    ax2.set_xlim(0, 270)
    ax2.set_ylim(-0.08, 0.18)
    ax2.set_xlabel("Grain temperature (K)")
    ax2.set_ylabel(r"$\epsilon_{\rm model}-\epsilon_{\rm expt}$")
    add_panel_labels([ax1, ax2], labels="AB")
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.16, top=0.98, wspace=0.26)
    save_figure(fig, str(out_dir / "fig04_grieco_validation"))
    plt.close(fig)


def make_fig07_efficiency(df: pd.DataFrame, eps_col: str, eps_ci: str | None, out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(DOUBLE_COLUMN_WIDTH, 3.45))
    add_regime_band(ax, 5, 20, "Diffusion-limited", color="#d9d9d9", alpha=0.09, y=0.972)
    add_regime_band(ax, 20, 100, "Peak LH window", color="#fff4bf", alpha=0.15, y=0.972)
    add_regime_band(ax, 120, 270, "ER plateau", color="#fde1de", alpha=0.13, y=0.972)
    for n_h in sorted(df["h_gas_density_cm3"].dropna().unique()):
        sub = df[np.isclose(df["h_gas_density_cm3"], float(n_h))].sort_values("surface_temperature_k")
        style = density_style(float(n_h))
        ax.plot(
            sub["surface_temperature_k"],
            sub[eps_col],
            color=style["color"],
            lw=1.6,
            marker=style["marker"],
            ms=4.2,
            label=fr"$n_{{\rm H}}=10^{{{int(np.log10(n_h))}}}\,\rm cm^{{-3}}$",
        )
        if eps_ci:
            ci = sub[eps_ci].fillna(0.0)
            ax.fill_between(sub["surface_temperature_k"], sub[eps_col] - ci, sub[eps_col] + ci, color=style["color"], alpha=0.14, lw=0)
    ax.hlines(0.18, xmin=120.0, xmax=270.0, colors="black", linestyles="--", lw=0.8)
    ax.text(266.0, 0.183, r"analytic limit", ha="right", va="bottom", fontsize=7.0, color="#444444")
    style_axes(ax)
    ax.set_xlim(5.0, 270.0)
    ax.set_ylim(0.0, 0.35)
    ax.set_xlabel("Grain surface temperature (K)")
    ax.set_ylabel(r"H$_2$ formation efficiency $\epsilon$")
    style_legend(ax, loc="lower left", bbox_to_anchor=(0.015, 0.02), ncol=2)
    save_figure(fig, str(out_dir / "fig07_efficiency_density"))
    plt.close(fig)


def make_fig08_mechanisms(df: pd.DataFrame, out_dir: Path) -> None:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(DOUBLE_COLUMN_WIDTH, 3.0), gridspec_kw={"wspace": 0.28})
    left = df[np.isclose(df["h_gas_density_cm3"], 100.0)].sort_values("surface_temperature_k").copy()
    total = left["h2_formed_LH_mean"] + left["h2_formed_ER_mean"]
    left["lh_frac"] = left["h2_formed_LH_mean"] / total.replace(0.0, np.nan)
    left["er_frac"] = left["h2_formed_ER_mean"] / total.replace(0.0, np.nan)
    ax1.stackplot(
        left["surface_temperature_k"],
        left["lh_frac"],
        left["er_frac"],
        colors=[OKABE_ITO["blue"], OKABE_ITO["vermillion"]],
        alpha=0.82,
    )
    ax1.axvline(110.0, color="#666666", lw=0.8, ls="--")
    ax1.text(42.0, 0.74, "LH", color="white", fontsize=7.6, fontweight="bold")
    ax1.text(175.0, 0.76, "ER", color="white", fontsize=7.6, fontweight="bold")
    ax1.text(112.0, 0.08, "crossover", rotation=90, va="bottom", ha="left", fontsize=6.4, color="#666666")
    style_axes(ax1)
    ax1.set_xlim(10.0, 250.0)
    ax1.set_ylim(0.0, 1.0)
    ax1.set_xlabel("Grain surface temperature (K)")
    ax1.set_ylabel("Fraction of formed H$_2$")

    for n_h in sorted(df["h_gas_density_cm3"].dropna().unique()):
        sub = df[np.isclose(df["h_gas_density_cm3"], float(n_h))].sort_values("surface_temperature_k").copy()
        total = sub["h2_formed_LH_mean"] + sub["h2_formed_ER_mean"]
        er_frac = 100.0 * sub["h2_formed_ER_mean"] / total.replace(0.0, np.nan)
        mask = sub["surface_temperature_k"].between(80.0, 150.0)
        style = density_style(float(n_h))
        ax2.plot(
            sub.loc[mask, "surface_temperature_k"],
            er_frac.loc[mask],
            color=style["color"],
            lw=1.5,
            marker=style["marker"],
            ms=4.2,
            label=fr"$10^{{{int(np.log10(n_h))}}}$",
        )
    style_axes(ax2)
    ax2.set_xlim(80.0, 150.0)
    ax2.set_ylim(0.0, 100.0)
    ax2.set_xlabel("Grain surface temperature (K)")
    ax2.set_ylabel("ER fraction (%)")
    style_legend(ax2, title=r"$n_{\rm H}$", loc="lower right")
    add_panel_labels([ax1, ax2], labels="AB")
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.18, top=0.98, wspace=0.28)
    save_figure(fig, str(out_dir / "fig08_mechanism_decomposition"))
    plt.close(fig)


def make_fig09_surface_inventory(df: pd.DataFrame, out_dir: Path) -> None:
    sub = df[np.isclose(df["h_gas_density_cm3"], 100.0)].sort_values("surface_temperature_k")
    ycol = "final_h_atoms_on_surface_mean"
    fig, ax = plt.subplots(figsize=(SINGLE_COLUMN_WIDTH, 2.7))
    add_regime_band(ax, 10, 80, "LH-depleted", color="#d9d9d9", alpha=0.10, y=0.08)
    add_regime_band(ax, 80, 140, "Transition", color="#fff4bf", alpha=0.18, y=0.08)
    add_regime_band(ax, 140, 250, "ER reservoir", color="#fde1de", alpha=0.16, y=0.08)
    ax.plot(sub["surface_temperature_k"], sub[ycol], color=OKABE_ITO["blue"], lw=1.8, marker="o", ms=4.4)
    sat = float(sub.loc[sub["surface_temperature_k"] >= 150.0, ycol].iloc[0])
    ax.axhline(sat, color="#666666", lw=0.8, ls=":")
    ax.text(247.0, sat * 1.03, "chemisorption saturation", ha="right", va="bottom", fontsize=6.8, color="#555555")
    style_axes(ax)
    ax.set_xlim(10.0, 250.0)
    ax.set_yscale("log")
    ax.set_ylim(1.0, 500.0)
    ax.set_xlabel("Grain surface temperature (K)")
    ax.set_ylabel(r"$\langle N_{\rm H}\rangle$ on grain")
    save_figure(fig, str(out_dir / "fig09_surface_inventory"))
    plt.close(fig)


def make_fig10_transition_zoom(df: pd.DataFrame, eps_col: str, eps_ci: str | None, out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(SINGLE_COLUMN_WIDTH, 2.8))
    for n_h in sorted(df["h_gas_density_cm3"].dropna().unique()):
        sub = df[np.isclose(df["h_gas_density_cm3"], float(n_h))].sort_values("surface_temperature_k")
        mask = sub["surface_temperature_k"].between(80.0, 140.0)
        style = density_style(float(n_h))
        ax.plot(sub.loc[mask, "surface_temperature_k"], sub.loc[mask, eps_col], color=style["color"], lw=1.5, marker=style["marker"], ms=4.2, label=fr"$10^{{{int(np.log10(n_h))}}}$")
        if eps_ci:
            ci = sub.loc[mask, eps_ci].fillna(0.0)
            ax.fill_between(sub.loc[mask, "surface_temperature_k"], sub.loc[mask, eps_col] - ci, sub.loc[mask, eps_col] + ci, color=style["color"], alpha=0.12, lw=0)
    ax.axvline(100.0, color="#666666", lw=0.8, ls="--")
    low = float(df[(np.isclose(df["h_gas_density_cm3"], 10.0)) & (np.isclose(df["surface_temperature_k"], 100.0))][eps_col].iloc[0])
    high = float(df[(np.isclose(df["h_gas_density_cm3"], 10000.0)) & (np.isclose(df["surface_temperature_k"], 100.0))][eps_col].iloc[0])
    ax.annotate("", xy=(100.0, high), xytext=(100.0, low), arrowprops={"arrowstyle": "<->", "lw": 0.9, "color": "#444444"})
    ax.text(102.0, 0.5 * (low + high), f"{100.0*(high/low-1.0):.0f}%", fontsize=7.0, va="center", color="#444444")
    ax.text(83.0, 0.284, "stochastic enhancement", fontsize=7.0, style="italic", color="#555555")
    style_axes(ax)
    ax.set_xlim(80.0, 140.0)
    ax.set_ylim(0.17, 0.30)
    ax.set_xlabel("Grain surface temperature (K)")
    ax.set_ylabel(r"Formation efficiency $\epsilon$")
    style_legend(ax, title=r"$n_{\rm H}$", loc="upper right")
    save_figure(fig, str(out_dir / "fig10_transition_zoom"))
    plt.close(fig)


def make_fig11_phase_map(df: pd.DataFrame, eps_col: str, out_dir: Path) -> None:
    pivot = df.pivot_table(index="h_gas_density_cm3", columns="surface_temperature_k", values=eps_col, aggfunc="first").sort_index().sort_index(axis=1)
    dens = pivot.index.to_numpy(dtype=float)
    temps = pivot.columns.to_numpy(dtype=float)
    z = pivot.to_numpy(dtype=float)
    temp_dense, dens_dense, z_dense = interpolate_phase_map(temps, dens, z, nx=260, ny=220, log_y=True)

    fig, ax = plt.subplots(figsize=(SINGLE_COLUMN_WIDTH, 3.05))
    mesh = ax.pcolormesh(temp_dense, dens_dense, z_dense, cmap="viridis", vmin=0.0, vmax=0.32, shading="nearest")
    tt, nn = np.meshgrid(temp_dense, dens_dense)
    contours = ax.contour(tt, nn, z_dense, levels=[0.19, 0.25], colors="white", linewidths=0.85, linestyles="--")
    ax.clabel(
        contours,
        fmt=lambda v: f"{v:.2f}",
        inline=True,
        fontsize=6.0,
        colors="white",
        manual=[(84.0, 16.0), (150.0, 3500.0)],
    )
    style_axes(ax)
    ax.set_yscale("log")
    ax.set_xlim(10.0, 250.0)
    ax.set_yticks([10, 100, 1000, 10000])
    ax.set_xlabel("Grain surface temperature (K)")
    ax.set_ylabel(r"$n_{\rm H}$ (cm$^{-3}$)")
    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label(r"$\epsilon$")
    cbar.ax.tick_params(direction="in", length=3)
    save_figure(fig, str(out_dir / "fig11_phase_map"))
    plt.close(fig)


def make_fig12_rates(df: pd.DataFrame, rate_col: str, rate_ci: str | None, out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(SINGLE_COLUMN_WIDTH, 2.65))
    for n_h in sorted(df["h_gas_density_cm3"].dropna().unique()):
        sub = df[np.isclose(df["h_gas_density_cm3"], float(n_h))].sort_values("surface_temperature_k")
        style = density_style(float(n_h))
        ax.plot(sub["surface_temperature_k"], sub[rate_col], color=style["color"], lw=1.6, marker=style["marker"], ms=4.2, label=fr"$10^{{{int(np.log10(n_h))}}}$")
        if rate_ci:
            ci = sub[rate_ci].fillna(0.0)
            ax.fill_between(sub["surface_temperature_k"], np.maximum(sub[rate_col] - ci, 1e-30), sub[rate_col] + ci, color=style["color"], alpha=0.14, lw=0)
    style_axes(ax)
    ax.set_yscale("log")
    ax.set_xlim(10.0, 250.0)
    ax.set_xlabel("Grain surface temperature (K)")
    ax.set_ylabel(r"$R({\rm H}_2)$ (cm$^{-2}$ s$^{-1}$)")
    style_legend(ax, title=r"$n_{\rm H}$", loc="upper right")
    save_figure(fig, str(out_dir / "fig12_release_rate"))
    plt.close(fig)


def make_fig13_ct10(df: pd.DataFrame, rate_col: str, rate_ci: str | None, ct_path: str, out_dir: Path) -> None:
    ct = pd.read_csv(ct_path)
    for col in ["surface_temperature_k", "h_gas_density_cm3", "uv_flux_factor", "ct02_h2_release_rate_cm2_s", "ratio_kmc_over_ct02"]:
        ct[col] = pd.to_numeric(ct[col], errors="coerce")
    ct = ct[np.isclose(ct["uv_flux_factor"], 0.0) & np.isclose(ct["h_gas_density_cm3"], 1000.0)].sort_values("surface_temperature_k")
    sub = df[np.isclose(df["h_gas_density_cm3"], 1000.0)].sort_values("surface_temperature_k")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(DOUBLE_COLUMN_WIDTH, 3.0), gridspec_kw={"wspace": 0.28})
    ax1.plot(sub["surface_temperature_k"], sub[rate_col], color=OKABE_ITO["blue"], lw=1.8, marker="o", ms=4.4, label="KMC")
    if rate_ci:
        ci = sub[rate_ci].fillna(0.0)
        ax1.fill_between(sub["surface_temperature_k"], np.maximum(sub[rate_col] - ci, 1e-30), sub[rate_col] + ci, color=OKABE_ITO["blue"], alpha=0.18, lw=0)
    ax1.plot(ct["surface_temperature_k"], ct["ct02_h2_release_rate_cm2_s"], color=OKABE_ITO["orange"], lw=1.6, ls="--", label="CT10")
    y_100 = float(sub.loc[np.isclose(sub["surface_temperature_k"], 100.0), rate_col].iloc[0])
    ax1.annotate(
        "back-diffusion\nsuppression",
        xy=(100.0, y_100),
        xytext=(145.0, 4.8e6),
        fontsize=6.9,
        color="#555555",
        ha="left",
        va="center",
        arrowprops={"arrowstyle": "->", "lw": 0.8, "color": "#666666"},
    )
    style_axes(ax1)
    ax1.set_yscale("log")
    ax1.set_xlim(10.0, 250.0)
    ax1.set_xlabel("Grain surface temperature (K)")
    ax1.set_ylabel(r"$R({\rm H}_2)$ (cm$^{-2}$ s$^{-1}$)")
    style_legend(ax1, loc="lower right")

    ax2.axhspan(0.5, 0.9, color="#e9e9e9", alpha=0.6, zorder=0)
    ax2.axhline(1.0, color="#666666", lw=0.8, ls="--")
    ax2.plot(ct["surface_temperature_k"], ct["ratio_kmc_over_ct02"], color="black", lw=1.6)
    style_axes(ax2)
    ax2.set_xlim(10.0, 250.0)
    ax2.set_ylim(0.0, 2.0)
    ax2.set_xlabel("Grain surface temperature (K)")
    ax2.set_ylabel(r"$R_{\rm KMC}/R_{\rm CT10}$")
    add_panel_labels([ax1, ax2], labels="AB")
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.18, top=0.98, wspace=0.28)
    save_figure(fig, str(out_dir / "fig13_ct10_comparison"))
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser(description="Generate publication-style core manuscript figures.")
    p.add_argument("--baseline", default="results/jhub_full_merged.csv")
    p.add_argument("--ct", default="results/tables/ct02_comparison_jhub_full.csv")
    p.add_argument("--iso-validation", default="results/grieco_validation_paper_iso_paperfit.csv")
    p.add_argument("--digitized", default="grieco_fig2_digitized.csv")
    p.add_argument("--out-dir", default="results/plots/manuscript")
    p.add_argument("--uv", type=float, default=0.0)
    args = p.parse_args()

    apply_publication_style()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df, eps_col, rate_col, eps_ci, rate_ci = _load_baseline(args.baseline, args.uv)

    make_fig04_validation(args.iso_validation, args.digitized, out_dir)
    make_fig07_efficiency(df, eps_col, eps_ci, out_dir)
    make_fig08_mechanisms(df, out_dir)
    make_fig09_surface_inventory(df, out_dir)
    make_fig10_transition_zoom(df, eps_col, eps_ci, out_dir)
    make_fig11_phase_map(df, eps_col, out_dir)
    make_fig12_rates(df, rate_col, rate_ci, out_dir)
    make_fig13_ct10(df, rate_col, rate_ci, args.ct, out_dir)


if __name__ == "__main__":
    main()
