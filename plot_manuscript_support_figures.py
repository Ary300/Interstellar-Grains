#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch
from matplotlib.colors import TwoSlopeNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from paper_plot_style import (
    DOUBLE_COLUMN_WIDTH,
    OKABE_ITO,
    SINGLE_COLUMN_WIDTH,
    add_panel_labels,
    apply_publication_style,
    interpolate_phase_map,
    save_figure,
    style_axes,
    style_legend,
)


def _pick_col(df: pd.DataFrame, base: str) -> str:
    if base in df.columns:
        return base
    m = f"{base}_mean"
    if m in df.columns:
        return m
    raise KeyError(f"Missing column: {base} or {m}")


def _maybe_ci(df: pd.DataFrame, base: str) -> str | None:
    name = f"{base}_ci95"
    return name if name in df.columns else None


def make_fig03_schematic(out_dir: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE_COLUMN_WIDTH, 2.4), gridspec_kw={"wspace": 0.18})
    titles = ["Langmuir-Hinshelwood", "Eley-Rideal", "UV (future work)"]
    for i, (ax, title) in enumerate(zip(axes, titles)):
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        ax.text(0.02, 0.96, f"({chr(97+i)})", fontsize=10, fontweight="bold", va="top", ha="left")
        ax.text(0.5, 0.93, title, fontsize=8.2, ha="center", va="top", color="#222222")
        if i == 2:
            ax.add_patch(
                FancyBboxPatch(
                    (0.04, 0.08),
                    0.92,
                    0.78,
                    boxstyle="round,pad=0.015,rounding_size=0.02",
                    ec="#aaaaaa",
                    fc="#fafafa",
                    lw=0.9,
                    linestyle=(0, (3, 2)),
                    alpha=0.8,
                )
            )
        ax.plot([0.08, 0.92], [0.28, 0.28], color="#7a7a7a", lw=1.4, solid_capstyle="round")

    # (a) LH cartoon
    ax = axes[0]
    for x in [0.24, 0.48, 0.72]:
        ax.plot([x - 0.05, x, x + 0.05], [0.28, 0.22, 0.28], color="#9a9a9a", lw=0.9)
    h1 = Circle((0.28, 0.24), 0.04, facecolor=OKABE_ITO["blue"], edgecolor="white", lw=0.8)
    h2 = Circle((0.68, 0.24), 0.04, facecolor=OKABE_ITO["blue"], edgecolor="white", lw=0.8)
    ax.add_patch(h1)
    ax.add_patch(h2)
    ax.text(0.28, 0.24, "H", ha="center", va="center", fontsize=7.0, color="white")
    ax.text(0.68, 0.24, "H", ha="center", va="center", fontsize=7.0, color="white")
    ax.add_patch(FancyArrowPatch((0.34, 0.36), (0.46, 0.30), arrowstyle="->", mutation_scale=10, lw=0.9, color="#444444"))
    ax.add_patch(FancyArrowPatch((0.62, 0.36), (0.54, 0.30), arrowstyle="->", mutation_scale=10, lw=0.9, color="#444444"))
    ax.add_patch(FancyArrowPatch((0.50, 0.33), (0.50, 0.60), arrowstyle="->", mutation_scale=11, lw=1.0, color="#444444"))
    ax.plot([0.46, 0.54], [0.58, 0.58], color="#444444", lw=0.8)
    ax.add_patch(Circle((0.44, 0.58), 0.03, facecolor=OKABE_ITO["blue"], edgecolor="white", lw=0.7))
    ax.add_patch(Circle((0.56, 0.58), 0.03, facecolor=OKABE_ITO["blue"], edgecolor="white", lw=0.7))
    ax.text(0.50, 0.68, r"H$_2$", ha="center", va="center", fontsize=7.4, color="#333333")
    ax.text(0.50, 0.11, "surface diffusion", ha="center", va="center", fontsize=6.5, color="#666666", style="italic")

    # (b) ER cartoon
    ax = axes[1]
    ax.plot([0.42, 0.50, 0.58], [0.28, 0.16, 0.28], color="#9a9a9a", lw=1.0)
    ax.add_patch(Circle((0.50, 0.20), 0.04, facecolor=OKABE_ITO["vermillion"], edgecolor="white", lw=0.8))
    ax.text(0.50, 0.20, "H", ha="center", va="center", fontsize=7.0, color="white")
    ax.add_patch(Circle((0.50, 0.72), 0.04, facecolor=OKABE_ITO["blue"], edgecolor="white", lw=0.8))
    ax.text(0.50, 0.72, "H", ha="center", va="center", fontsize=7.0, color="white")
    ax.add_patch(FancyArrowPatch((0.50, 0.66), (0.50, 0.33), arrowstyle="->", mutation_scale=10, lw=0.9, color="#444444"))
    ax.add_patch(FancyArrowPatch((0.56, 0.32), (0.74, 0.56), arrowstyle="->", mutation_scale=11, lw=1.0, color="#444444"))
    ax.add_patch(Circle((0.76, 0.58), 0.03, facecolor=OKABE_ITO["blue"], edgecolor="white", lw=0.7))
    ax.add_patch(Circle((0.86, 0.58), 0.03, facecolor=OKABE_ITO["vermillion"], edgecolor="white", lw=0.7))
    ax.plot([0.79, 0.83], [0.58, 0.58], color="#444444", lw=0.8)
    ax.text(0.81, 0.67, r"H$_2$", ha="center", va="center", fontsize=7.4, color="#333333")
    ax.text(0.50, 0.11, "direct gas-surface reaction", ha="center", va="center", fontsize=6.5, color="#666666", style="italic")

    # (c) UV future-work panel
    ax = axes[2]
    for theta in np.linspace(0, 2 * np.pi, 6, endpoint=False):
        ax.add_patch(Circle((0.50 + 0.12 * np.cos(theta), 0.48 + 0.12 * np.sin(theta)), 0.028, facecolor="#999999", edgecolor="white", lw=0.6, alpha=0.65))
    ax.add_patch(Circle((0.50, 0.48), 0.028, facecolor="#999999", edgecolor="white", lw=0.6, alpha=0.65))
    ax.add_patch(Circle((0.66, 0.58), 0.03, facecolor=OKABE_ITO["blue"], edgecolor="white", lw=0.6, alpha=0.75))
    ax.add_patch(Circle((0.36, 0.63), 0.03, facecolor=OKABE_ITO["blue"], edgecolor="white", lw=0.6, alpha=0.75))
    ax.text(0.18, 0.76, r"$h\nu$", fontsize=8.0, color="#777777")
    ax.plot([0.22, 0.27, 0.30, 0.35, 0.38], [0.72, 0.75, 0.69, 0.73, 0.68], color="#777777", lw=1.0, ls="--")
    ax.add_patch(FancyArrowPatch((0.63, 0.55), (0.82, 0.72), arrowstyle="->", mutation_scale=10, lw=0.9, color="#777777"))
    ax.text(0.84, 0.74, r"H$_2$", ha="center", va="center", fontsize=7.3, color="#666666")
    for artist in ax.get_children():
        if hasattr(artist, "set_alpha"):
            try:
                artist.set_alpha(0.82 if getattr(artist, "get_alpha", lambda: None)() is None else artist.get_alpha())
            except Exception:
                pass

    save_figure(fig, str(out_dir / "fig03_mechanism_schematic"))
    plt.close(fig)


def make_fig05_ablation(full_iso: str, full_ded: str, phys_only: str, digitized: str, out_dir: Path) -> None:
    full = pd.read_csv(full_iso)
    ded = pd.read_csv(full_ded)
    phys = pd.read_csv(phys_only)
    expt = pd.read_csv(digitized)
    for df, cols in (
        (full, ["surface_temperature_k", "epsilon_mean", "epsilon_ci95"]),
        (ded, ["temperature_k", "epsilon_released_total_mean", "epsilon_released_total_ci95"]),
        (phys, ["surface_temperature_k", "epsilon_released_total_mean", "epsilon_released_total_ci95"]),
        (expt, ["T_K", "epsilon", "epsilon_err_low", "epsilon_err_high"]),
    ):
        for col in cols:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    fig, ax = plt.subplots(figsize=(SINGLE_COLUMN_WIDTH, 2.65))
    ax.plot(ded["temperature_k"], ded["epsilon_released_total_mean"], color=OKABE_ITO["blue"], lw=1.8)
    ax.fill_between(ded["temperature_k"], ded["epsilon_released_total_mean"] - ded["epsilon_released_total_ci95"], ded["epsilon_released_total_mean"] + ded["epsilon_released_total_ci95"], color=OKABE_ITO["blue"], alpha=0.18, lw=0)
    ax.plot(full["surface_temperature_k"], full["epsilon_mean"], color=OKABE_ITO["blue"], lw=1.8, label="Full model")
    ax.fill_between(full["surface_temperature_k"], full["epsilon_mean"] - full["epsilon_ci95"], full["epsilon_mean"] + full["epsilon_ci95"], color=OKABE_ITO["blue"], alpha=0.18, lw=0)
    ax.plot(phys["surface_temperature_k"], phys["epsilon_released_total_mean"], color=OKABE_ITO["vermillion"], lw=1.7, ls="--", label="Physisorption only")
    ax.fill_between(phys["surface_temperature_k"], phys["epsilon_released_total_mean"] - phys["epsilon_released_total_ci95"], phys["epsilon_released_total_mean"] + phys["epsilon_released_total_ci95"], color=OKABE_ITO["vermillion"], alpha=0.12, lw=0)
    ax.errorbar(
        expt["T_K"],
        expt["epsilon"],
        yerr=[expt["epsilon_err_low"], expt["epsilon_err_high"]],
        fmt="o",
        mfc="white",
        mec="black",
        mew=1.0,
        ms=4.6,
        capsize=2,
        elinewidth=0.8,
        color="black",
        label="Grieco et al. (2023)",
    )
    ax.annotate(
        "chemisorption stabilizes\nwarm formation",
        xy=(20, float(phys.loc[np.isclose(phys["surface_temperature_k"], 20.0), "epsilon_released_total_mean"].iloc[0])),
        xytext=(60, 0.08),
        fontsize=6.8,
        color="#555555",
        style="italic",
        arrowprops={"arrowstyle": "->", "lw": 0.8, "color": "#666666"},
    )
    style_axes(ax)
    ax.set_xlim(0.0, 270.0)
    ax.set_ylim(0.0, 0.50)
    ax.set_xlabel("Grain temperature (K)")
    ax.set_ylabel(r"Formation efficiency $\epsilon$")
    style_legend(ax, loc="upper right")
    save_figure(fig, str(out_dir / "fig05_physisorption_ablation"))
    plt.close(fig)


def make_fig06_tau(table_path: str, digitized: str, out_dir: Path) -> None:
    df = pd.read_csv(table_path)
    expt = pd.read_csv(digitized)
    expt_20 = expt[np.isclose(expt["T_K"], 20.0)]["epsilon"].iloc[0]
    expt_200 = expt[np.isclose(expt["T_K"], 200.0)]["epsilon"].iloc[0]
    expt_20_err = expt[np.isclose(expt["T_K"], 20.0)]["epsilon_err_high"].iloc[0]
    expt_200_err = expt[np.isclose(expt["T_K"], 200.0)]["epsilon_err_high"].iloc[0]

    fig, ax = plt.subplots(figsize=(SINGLE_COLUMN_WIDTH, 2.4))
    ax.axhspan(expt_20 - expt_20_err, expt_20 + expt_20_err, color=OKABE_ITO["blue"], alpha=0.10)
    ax.axhspan(expt_200 - expt_200_err, expt_200 + expt_200_err, color=OKABE_ITO["vermillion"], alpha=0.10)
    ax.errorbar(df["tau"], df["eps20_isothermal_released_total_mean"], yerr=df["eps20_isothermal_released_total_ci95"], fmt="o-", color=OKABE_ITO["blue"], ms=4.8, lw=1.5, capsize=2, label="20 K")
    ax.errorbar(df["tau"], df["eps200_isothermal_mean"], yerr=df["eps200_isothermal_ci95"], fmt="s-", color=OKABE_ITO["vermillion"], ms=4.6, lw=1.5, capsize=2, label="200 K")
    style_axes(ax)
    ax.set_xlim(df["tau"].min() - 0.05, df["tau"].max() + 0.05)
    ax.set_ylim(0.15, 0.34)
    ax.set_xlabel(r"Dissociation fraction $\tau$")
    ax.set_ylabel(r"Formation efficiency $\epsilon$")
    style_legend(ax, loc="upper right", title="Temperature")
    save_figure(fig, str(out_dir / "fig06_tau_sensitivity"))
    plt.close(fig)


def make_fig14_mrn(size_summary: str, weights_path: str, comparison_path: str, out_dir: Path) -> None:
    size_df = pd.read_csv(size_summary)
    w_df = pd.read_csv(weights_path)
    comp = pd.read_csv(comparison_path)
    for c in ["surface_temperature_k", "grain_radius_um_mrn", "epsilon_mean", "h2_release_rate_cm2_s_mean"]:
        size_df[c] = pd.to_numeric(size_df[c], errors="coerce")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(SINGLE_COLUMN_WIDTH, 4.45), gridspec_kw={"hspace": 0.18})
    temps = [(20.0, OKABE_ITO["blue"], "o"), (100.0, OKABE_ITO["purple"], "^"), (200.0, OKABE_ITO["orange"], "D")]
    for temp, color, marker in temps:
        sub = size_df[np.isclose(size_df["surface_temperature_k"], temp)].sort_values("grain_radius_um_mrn")
        ax1.plot(sub["grain_radius_um_mrn"], sub["epsilon_mean"], color=color, marker=marker, ms=4.2, lw=1.5, label=f"{int(temp)} K")
    style_axes(ax1)
    ax1.set_xscale("log")
    ax1.set_xlim(0.005, 0.25)
    ax1.set_ylim(0.0, 0.4)
    ax1.set_ylabel(r"Efficiency $\epsilon$")
    style_legend(ax1, loc="upper right", title="Surface T")

    warm = size_df[np.isclose(size_df["surface_temperature_k"], 100.0)].sort_values("grain_radius_um_mrn")
    merged = warm.merge(w_df, on="grain_radius_um_mrn", how="left").sort_values("grain_radius_um_mrn")
    merged["cum_weight"] = merged["mrn_weight_area"].cumsum()
    merged["cum_eps"] = (merged["epsilon_mean"] * merged["mrn_weight_area"]).cumsum() / merged["cum_weight"]
    ax2.plot(merged["grain_radius_um_mrn"], merged["mrn_weight_area"], color="black", lw=1.2, ls="--", label="MRN weight")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim(0.005, 0.25)
    ax2.set_ylabel("MRN weight")
    style_axes(ax2)
    ax2b = ax2.twinx()
    ax2b.plot(merged["grain_radius_um_mrn"], merged["cum_eps"], color=OKABE_ITO["blue"], lw=1.8, label="Cumulative $\epsilon$")
    single_100 = float(comp.loc[np.isclose(comp["surface_temperature_k"], 100.0), "epsilon_single"].iloc[0])
    ax2b.axhline(single_100, color=OKABE_ITO["vermillion"], lw=0.9, ls=":", label="Single-size $\epsilon$")
    ax2b.set_ylim(0.18, 0.30)
    ax2b.set_ylabel(r"Cumulative integrated $\epsilon$")
    ax2.set_xlabel("Grain radius (μm)")
    lines = ax2.get_lines() + ax2b.get_lines()
    labels = [line.get_label() for line in lines]
    ax2.legend(lines, labels, loc="lower right", frameon=False, fontsize=6.9, borderpad=0.2, handletextpad=0.5)
    add_panel_labels([ax1, ax2], labels="AB")
    save_figure(fig, str(out_dir / "fig14_mrn_integration"))
    plt.close(fig)


def make_fig16_convergence(raw_path: str, out_dir: Path) -> None:
    df = pd.read_csv(raw_path)
    sub = df[np.isclose(df["surface_temperature_k"], 100.0) & np.isclose(df["h_gas_density_cm3"], 10000.0)].copy()
    eps = sub["epsilon"].to_numpy(dtype=float)
    n = np.arange(1, len(eps) + 1)
    running_mean = np.cumsum(eps) / n
    running_std = np.array([np.std(eps[:i], ddof=1) if i > 1 else 0.0 for i in n])
    sem = running_std / np.sqrt(n)

    fig, ax = plt.subplots(figsize=(SINGLE_COLUMN_WIDTH, 2.45))
    ax.plot(n, running_mean, color=OKABE_ITO["blue"], lw=1.7)
    ax.fill_between(n, running_mean - sem, running_mean + sem, color=OKABE_ITO["blue"], alpha=0.18, lw=0)
    ax.axvline(20, color="#666666", lw=0.8, ls="--")
    ax.axhline(running_mean[-1], color="#666666", lw=0.8, ls=":")
    ax.text(22, running_mean[30] if len(running_mean) > 30 else running_mean[-1], "N = 20", fontsize=6.8, color="#555555", va="bottom")
    style_axes(ax)
    ax.set_xscale("log")
    ax.set_xlim(1, len(n))
    ax.set_xlabel("Number of realizations")
    ax.set_ylabel(r"Running mean of $\epsilon$")
    inset = inset_axes(ax, width="36%", height="42%", loc="lower left", borderpad=1.2)
    inset.plot(n[1:], sem[1:], color="black", lw=1.1)
    inset.set_xscale("log")
    inset.set_yscale("log")
    style_axes(inset)
    inset.tick_params(labelsize=6)
    inset.set_title("SEM", fontsize=6.4, pad=1.5)
    save_figure(fig, str(out_dir / "fig16_ensemble_convergence"))
    plt.close(fig)


def make_fig17_sensitivity(knobs_path: str, out_dir: Path) -> None:
    df = pd.read_csv(knobs_path)
    rate_col = _pick_col(df, "h2_release_rate_cm2_s")
    df["k_eff"] = 4.0e-21 * pd.to_numeric(df[rate_col], errors="coerce") / pd.to_numeric(df["h_gas_density_cm3"], errors="coerce")
    fig, ax = plt.subplots(figsize=(SINGLE_COLUMN_WIDTH, 2.8))
    for (fchem, per), sub in df.groupby(["chemisorption_fraction", "er_reaction_probability"], dropna=False):
        sub = sub.sort_values("surface_temperature_k")
        is_base = np.isclose(float(fchem), 0.4) and np.isclose(float(per), 0.9)
        ax.plot(
            sub["surface_temperature_k"],
            sub["k_eff"],
            color=OKABE_ITO["blue"] if is_base else "#b6b6b6",
            lw=1.9 if is_base else 0.9,
            alpha=1.0 if is_base else 0.9,
            zorder=3 if is_base else 1,
        )
    env = df.groupby("surface_temperature_k")["k_eff"].agg(["min", "max"]).reset_index()
    ax.fill_between(env["surface_temperature_k"], env["min"], env["max"], color=OKABE_ITO["sky"], alpha=0.18, lw=0, zorder=0)
    style_axes(ax)
    ax.set_yscale("log")
    ax.set_xlabel("Surface temperature (K)")
    ax.set_ylabel(r"Effective $k_{\rm eff}$ (cm$^3$ s$^{-1}$)")
    save_figure(fig, str(out_dir / "fig17_sensitivity_envelope"))
    plt.close(fig)


def make_fig18_boxplots(raw_path: str, out_dir: Path) -> None:
    df = pd.read_csv(raw_path)
    temps = sorted(df["surface_temperature_k"].dropna().unique())
    fig, ax = plt.subplots(figsize=(SINGLE_COLUMN_WIDTH, 2.85))
    pos_low = np.arange(len(temps)) * 1.25 - 0.18
    pos_high = np.arange(len(temps)) * 1.25 + 0.18
    low = [df[np.isclose(df["surface_temperature_k"], t) & np.isclose(df["h_gas_density_cm3"], 10.0)]["epsilon"].to_numpy(dtype=float) for t in temps]
    high = [df[np.isclose(df["surface_temperature_k"], t) & np.isclose(df["h_gas_density_cm3"], 10000.0)]["epsilon"].to_numpy(dtype=float) for t in temps]
    common = dict(widths=0.28, patch_artist=True, showfliers=False, notch=True, medianprops={"color": "black", "lw": 1.0}, whiskerprops={"lw": 0.8}, capprops={"lw": 0.8})
    b1 = ax.boxplot(low, positions=pos_low, **common)
    b2 = ax.boxplot(high, positions=pos_high, **common)
    for patch in b1["boxes"]:
        patch.set_facecolor(OKABE_ITO["sky"])
        patch.set_edgecolor(OKABE_ITO["sky"])
        patch.set_alpha(0.9)
    for patch in b2["boxes"]:
        patch.set_facecolor(OKABE_ITO["vermillion"])
        patch.set_edgecolor(OKABE_ITO["vermillion"])
        patch.set_alpha(0.9)
    style_axes(ax)
    ax.set_xticks(np.arange(len(temps)) * 1.25)
    ax.set_xticklabels([f"{int(t)}" for t in temps])
    ax.set_xlabel("Surface temperature (K)")
    ax.set_ylabel(r"Per-run efficiency $\epsilon$")
    ax.legend([b1["boxes"][0], b2["boxes"][0]], [r"$n_{\rm H}=10$", r"$n_{\rm H}=10^4$"], frameon=False, loc="upper right", fontsize=7.0)
    save_figure(fig, str(out_dir / "fig18_transition_boxplots"))
    plt.close(fig)


def make_fig19_lh_mode(path: str, out_dir: Path) -> None:
    df = pd.read_csv(path)
    eps_col = _pick_col(df, "epsilon")
    ci_col = _maybe_ci(df, "epsilon")
    pivot = df.pivot_table(index="surface_temperature_k", columns="lh_formation_mode", values=eps_col, aggfunc="first").sort_index()
    pivot_ci = df.pivot_table(index="surface_temperature_k", columns="lh_formation_mode", values=ci_col, aggfunc="first").sort_index() if ci_col else None
    temps = pivot.index.to_numpy(dtype=float)
    x = np.arange(len(temps))
    w = 0.32
    fig, ax = plt.subplots(figsize=(SINGLE_COLUMN_WIDTH, 2.65))
    ax.bar(x - w / 2, pivot["pairs"], width=w, color=OKABE_ITO["blue"], label="pairs", yerr=pivot_ci["pairs"] if pivot_ci is not None else None, capsize=2)
    ax.bar(x + w / 2, pivot["diffusion_limited"], width=w, color=OKABE_ITO["vermillion"], label="diffusion-limited", yerr=pivot_ci["diffusion_limited"] if pivot_ci is not None else None, capsize=2)
    for i, temp in enumerate(temps):
        pct = 100.0 * (pivot.loc[temp, "diffusion_limited"] / pivot.loc[temp, "pairs"] - 1.0)
        ax.text(x[i], max(pivot.loc[temp, "pairs"], pivot.loc[temp, "diffusion_limited"]) + 0.01, f"{pct:.0f}%", ha="center", va="bottom", fontsize=6.7, color="#555555")
    style_axes(ax)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{int(t)}" for t in temps])
    ax.set_xlabel("Surface temperature (K)")
    ax.set_ylabel(r"Formation efficiency $\epsilon$")
    style_legend(ax, loc="upper right")
    save_figure(fig, str(out_dir / "fig19_lh_mode_consistency"))
    plt.close(fig)


def make_fig20_grain_size(path: str, out_dir: Path) -> None:
    df = pd.read_csv(path)
    eps_col = _pick_col(df, "epsilon")
    ci_col = _maybe_ci(df, "epsilon")
    fig, ax = plt.subplots(figsize=(SINGLE_COLUMN_WIDTH, 2.45))
    styles = {20.0: (OKABE_ITO["blue"], "o"), 150.0: (OKABE_ITO["vermillion"], "D")}
    for temp, sub in df.groupby("surface_temperature_k", dropna=False):
        sub = sub.sort_values("grain_radius_um")
        color, marker = styles.get(float(temp), (OKABE_ITO["black"], "o"))
        ax.errorbar(sub["grain_radius_um"], sub[eps_col], yerr=sub[ci_col] if ci_col else None, fmt=f"{marker}-", color=color, lw=1.5, ms=4.6, capsize=2, label=f"{int(temp)} K")
    style_axes(ax)
    ax.set_xscale("log")
    ax.set_xlabel("Grain radius (μm)")
    ax.set_ylabel(r"Formation efficiency $\epsilon$")
    style_legend(ax, loc="best")
    save_figure(fig, str(out_dir / "fig20_grain_size_dependence"))
    plt.close(fig)


def make_fig21_diagnostics(porosity_path: str, sticking_path: str, out_dir: Path) -> None:
    por = pd.read_csv(porosity_path)
    sti = pd.read_csv(sticking_path)
    eps_col_p = _pick_col(por, "epsilon")
    ci_col_p = _maybe_ci(por, "epsilon")
    eps_col_s = _pick_col(sti, "epsilon")
    ci_col_s = _maybe_ci(sti, "epsilon")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(SINGLE_COLUMN_WIDTH, 4.45), gridspec_kw={"hspace": 0.20})

    por = por.sort_values("porosity_fraction")
    ax1.errorbar(
        por["porosity_fraction"],
        por[eps_col_p],
        yerr=por[ci_col_p] if ci_col_p else None,
        fmt="o",
        color=OKABE_ITO["blue"],
        ms=5.2,
        lw=1.3,
        capsize=2,
    )
    ax1.plot(por["porosity_fraction"], por[eps_col_p], color=OKABE_ITO["blue"], lw=1.2)
    style_axes(ax1)
    ax1.set_xlim(-0.02, 0.22)
    ax1.set_xlabel("Porosity fraction")
    ax1.set_ylabel(r"Efficiency $\epsilon$")

    for model, sub in sti.groupby("sticking_temp_model", dropna=False):
        sub = sub.sort_values("surface_temperature_k")
        if str(model) == "constant":
            color, marker, label = OKABE_ITO["blue"], "o", "Constant $S=0.5$"
        else:
            color, marker, label = OKABE_ITO["vermillion"], "s", "Temperature-dependent $S(T)$"
        ax2.errorbar(sub["surface_temperature_k"], sub[eps_col_s], yerr=sub[ci_col_s] if ci_col_s else None, fmt=f"{marker}-", color=color, lw=1.5, ms=4.4, capsize=2, label=label)
    style_axes(ax2)
    ax2.set_xlabel("Surface temperature (K)")
    ax2.set_ylabel(r"Efficiency $\epsilon$")
    style_legend(ax2, loc="best")
    add_panel_labels([ax1, ax2], labels="AB")
    save_figure(fig, str(out_dir / "fig21_porosity_sticking"))
    plt.close(fig)


def make_fig22_timescales(baseline_path: str, out_dir: Path) -> None:
    df = pd.read_csv(baseline_path)
    rate_col = _pick_col(df, "h2_release_rate_cm2_s")
    for col in ["surface_temperature_k", "h_gas_density_cm3", "uv_flux_factor", rate_col]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["surface_temperature_k", "h_gas_density_cm3", "uv_flux_factor", rate_col])
    df = df[np.isclose(df["uv_flux_factor"], 0.0)].copy()

    sigma_H = 1.0e-21
    mu = 1.4
    G_CGS = 6.67430e-8
    M_H_G = 1.6735575e-24
    t_h2_s = 1.0 / (8.0 * sigma_H * df[rate_col].astype(float))
    rho = mu * M_H_G * df["h_gas_density_cm3"].astype(float)
    t_ff_s = np.sqrt((3.0 * np.pi) / (32.0 * G_CGS * rho))
    df["t_H2_over_t_ff"] = t_h2_s / t_ff_s
    df["log10_t_H2_over_t_ff"] = np.log10(df["t_H2_over_t_ff"])

    pivot = df.pivot_table(index="h_gas_density_cm3", columns="surface_temperature_k", values="log10_t_H2_over_t_ff", aggfunc="first").sort_index().sort_index(axis=1)
    dens = pivot.index.to_numpy(dtype=float)
    temps = pivot.columns.to_numpy(dtype=float)
    z = pivot.to_numpy(dtype=float)

    def _edges(vals: np.ndarray, log: bool = False) -> np.ndarray:
        vals = np.asarray(vals, dtype=float)
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

    x_edges = _edges(temps, log=False)
    y_edges = _edges(dens, log=True)
    temp_dense, dens_dense, z_dense = interpolate_phase_map(temps, dens, z, nx=260, ny=220, log_y=True)
    fig, ax = plt.subplots(figsize=(SINGLE_COLUMN_WIDTH, 3.05))
    norm = TwoSlopeNorm(vmin=-0.7, vcenter=0.0, vmax=0.7)
    mesh = ax.pcolormesh(temp_dense, dens_dense, z_dense, cmap="coolwarm", norm=norm, shading="nearest")
    tt, nn = np.meshgrid(temp_dense, dens_dense)
    ax.contour(tt, nn, z_dense, levels=[0.0], colors="white", linewidths=1.8)
    style_axes(ax)
    ax.set_yscale("log")
    ax.set_xlim(10.0, 250.0)
    ax.set_xlabel("Grain surface temperature (K)")
    ax.set_ylabel(r"$n_{\rm H}$ (cm$^{-3}$)")

    markers = [
        ("Diffuse cloud", 100.0, 10.0, "o"),
        ("Molecular cloud", 20.0, 1000.0, "s"),
        ("Dense clump", 30.0, 10000.0, "^"),
        ("z~6 clump", 60.0, 1000.0, "*"),
    ]
    for label, temp, nh, marker in markers:
        ax.plot(temp, nh, marker=marker, ms=6.2 if marker != "*" else 8.0, mec="black", mew=0.7, mfc="white", color="black")
    ax.text(108.0, 10.8, "Diffuse", fontsize=6.5, color="#333333", va="center")
    ax.text(31.0, 1220.0, "Molecular", fontsize=6.5, color="#333333", va="center")
    ax.text(35.0, 11500.0, "Dense", fontsize=6.5, color="#333333", va="bottom")
    ax.text(83.0, 1080.0, r"$z\!\sim\!6$", fontsize=6.5, color="#333333", va="center")
    ax.text(149.0, 1800.0, r"$t_{\rm H_2}=t_{\rm ff}$", fontsize=6.7, color="white", ha="left", va="center")
    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label(r"$\log_{10}(t_{\rm H_2}/t_{\rm ff})$")
    cbar.ax.tick_params(direction="in", length=3)
    save_figure(fig, str(out_dir / "fig22_timescale_phase_map"))
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser(description="Generate additional publication-style manuscript and appendix figures.")
    p.add_argument("--full-iso", default="results/grieco_validation_paper_iso_paperfit.csv")
    p.add_argument("--full-ded", default="results/grieco_ded_paper_ded_paperfit.csv")
    p.add_argument("--phys-only", default="results/grieco_physisorption_only_iso.csv")
    p.add_argument("--digitized", default="grieco_fig2_digitized.csv")
    p.add_argument("--tau-table", default="results/tables/grieco_tau_sensitivity_table.csv")
    p.add_argument("--mrn-size-summary", default="results/tables/astro_mrn_size_summary.csv")
    p.add_argument("--mrn-weights", default="results/tables/astro_mrn_weights.csv")
    p.add_argument("--mrn-comparison", default="results/tables/astro_mrn_comparison.csv")
    p.add_argument("--transition-raw", default="results/astro_transition_deep_raw.csv")
    p.add_argument("--knobs", default="results/astro_sensitivity_knobs.csv")
    p.add_argument("--lh-mode", default="results/astro_lh_mode_consistency.csv")
    p.add_argument("--grain-size", default="results/astro_grain_size_check.csv")
    p.add_argument("--porosity", default="results/astro_porosity_check.csv")
    p.add_argument("--sticking", default="results/astro_sticking_model_check.csv")
    p.add_argument("--baseline", default="results/jhub_full_merged.csv")
    p.add_argument("--out-dir", default="results/plots/manuscript")
    args = p.parse_args()

    apply_publication_style()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    make_fig03_schematic(out_dir)
    make_fig05_ablation(args.full_iso, args.full_ded, args.phys_only, args.digitized, out_dir)
    make_fig06_tau(args.tau_table, args.digitized, out_dir)
    make_fig14_mrn(args.mrn_size_summary, args.mrn_weights, args.mrn_comparison, out_dir)
    make_fig16_convergence(args.transition_raw, out_dir)
    make_fig17_sensitivity(args.knobs, out_dir)
    make_fig18_boxplots(args.transition_raw, out_dir)
    make_fig19_lh_mode(args.lh_mode, out_dir)
    make_fig20_grain_size(args.grain_size, out_dir)
    make_fig21_diagnostics(args.porosity, args.sticking, out_dir)
    make_fig22_timescales(args.baseline, out_dir)


if __name__ == "__main__":
    main()
