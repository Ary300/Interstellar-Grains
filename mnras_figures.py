# mnras_figures.py
# Shared style, palette, and helpers for all paper figures.
# All figure scripts MUST import from this module.

from __future__ import annotations
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


# -----------------------------------------------------------------------------
# COLOR PALETTE — Okabe–Ito, colorblind-safe. Use by name, never by hex.
# -----------------------------------------------------------------------------
COLORS = {
    "sky":        "#56B4E9",  # n_H = 10       (coldest density, coolest color)
    "blue":       "#0072B2",  # n_H = 100      (baseline)
    "purple":     "#CC79A7",  # n_H = 1000
    "vermillion": "#D55E00",  # n_H = 10000    (densest, warmest color)
    "green":      "#009E73",  # comparison / secondary
    "orange":     "#E69F00",  # tertiary
    "yellow":     "#F0E442",  # only with dark backgrounds
    "black":      "#000000",  # data points, reference lines
    "grey":       "#888888",  # annotations, de-emphasized
    "light_grey": "#CCCCCC",  # shaded regions background
}

# Canonical density → color mapping. Use this EVERYWHERE so figures are consistent.
DENSITY_COLOR = {
    10:    COLORS["sky"],
    100:   COLORS["blue"],
    1000:  COLORS["purple"],
    10000: COLORS["vermillion"],
}
DENSITY_MARKER = {10: "o", 100: "s", 1000: "^", 10000: "D"}

# Mechanism colors (fixed across figures)
MECH_COLOR = {"LH": COLORS["blue"], "ER": COLORS["vermillion"]}


# -----------------------------------------------------------------------------
# STYLE SETUP
# -----------------------------------------------------------------------------
STYLE_FILE = Path(__file__).parent / "mnras_style.mplstyle"


def setup_style() -> None:
    """Apply the MNRAS style. Call once at the top of every script."""
    if STYLE_FILE.exists():
        plt.style.use(str(STYLE_FILE))
    else:
        # Fallback: apply inline if style file missing
        plt.rcParams.update(_INLINE_STYLE)


_INLINE_STYLE = {
    # Fonts
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif"],
    "font.size": 8.0,
    "axes.labelsize": 9.0,
    "axes.titlesize": 9.0,
    "xtick.labelsize": 8.0,
    "ytick.labelsize": 8.0,
    "legend.fontsize": 7.5,
    "legend.frameon": False,
    "legend.handlelength": 1.8,
    "mathtext.fontset": "cm",
    "text.usetex": False,

    # Lines & markers
    "lines.linewidth": 1.3,
    "lines.markersize": 4.0,
    "lines.markeredgewidth": 0.8,

    # Axes
    "axes.linewidth": 0.7,
    "axes.edgecolor": "black",
    "axes.labelpad": 3.0,
    "axes.spines.top": True,
    "axes.spines.right": True,
    "axes.grid": False,

    # Ticks — inward, on all four sides (astronomy convention)
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
    "xtick.minor.size": 2.0,
    "ytick.minor.size": 2.0,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.minor.width": 0.5,
    "ytick.minor.width": 0.5,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,

    # Savefig
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.03,
    "pdf.fonttype": 42,  # editable text in PDFs
    "ps.fonttype": 42,
}


# -----------------------------------------------------------------------------
# FIGURE SIZE PRESETS (inches — MNRAS column widths)
# -----------------------------------------------------------------------------
SINGLE_COL = 3.33   # 240 pt = 10/3 inches
DOUBLE_COL = 7.00


def fig_single(height: float = 2.6) -> tuple:
    """Single-column figure size."""
    return (SINGLE_COL, height)


def fig_double(height: float = 3.2) -> tuple:
    """Double-column figure size."""
    return (DOUBLE_COL, height)


# -----------------------------------------------------------------------------
# HELPERS
# -----------------------------------------------------------------------------
def plot_with_ci(
    ax,
    x,
    y,
    y_err=None,
    y_low=None,
    y_high=None,
    color=None,
    label=None,
    linestyle="-",
    linewidth=1.5,
    marker=None,
    markersize=4.0,
    ci_alpha=0.12,
    ci_label=None,
    zorder=2,
    markevery=None,
):
    """
    Plot a line with an alpha-blended confidence band.

    Either pass y_err (symmetric ± band → ±1.96·y_err for 95% CI)
    OR pass y_low, y_high (explicit band edges).
    """
    x = np.asarray(x)
    y = np.asarray(y)

    if y_low is None or y_high is None:
        if y_err is not None:
            y_err = np.asarray(y_err)
            y_low = y - 1.96 * y_err
            y_high = y + 1.96 * y_err

    markerfacecolor = "white" if marker is not None else color
    markeredgewidth = 0.9 if marker is not None else 0.0

    if marker is not None and markevery is None and len(x) > 10:
        markevery = 3

    line, = ax.plot(
        x, y,
        color=color,
        linestyle=linestyle,
        linewidth=linewidth,
        marker=marker,
        markersize=markersize,
        markeredgecolor=color,
        markerfacecolor=markerfacecolor,
        markeredgewidth=markeredgewidth,
        markevery=markevery,
        label=label,
        zorder=zorder,
    )

    if y_low is not None and y_high is not None:
        ax.fill_between(
            x, y_low, y_high,
            color=color,
            alpha=ci_alpha,
            linewidth=0,
            zorder=zorder - 1,
            label=ci_label,
        )
    return line


def plot_data_points(ax, x, y, y_err, color="black", label=None, marker="o", zorder=5):
    """Plot experimental/reference data points with error bars — open markers."""
    return ax.errorbar(
        x, y,
        yerr=y_err,
        fmt=marker,
        markerfacecolor="white",
        markeredgecolor=color,
        markeredgewidth=1.0,
        markersize=4.5,
        ecolor=color,
        elinewidth=0.8,
        capsize=2.0,
        capthick=0.8,
        linestyle="none",
        label=label,
        zorder=zorder,
    )


def add_regime_bands(ax, regimes=None, text_y=0.95, fontsize=7):
    """
    Add shaded vertical regime bands with top-label text.

    regimes: list of (xmin, xmax, label, color, alpha)
    """
    if regimes is None:
        regimes = [
            (5, 20,    "Diffusion-\nlimited", COLORS["light_grey"], 0.25),
            (20, 100,  "Peak LH window",      "#FFF9C4",            0.35),
            (120, 300, "ER plateau",          "#FFEBEE",            0.45),
        ]
    for (xmin, xmax, label, color, alpha) in regimes:
        ax.axvspan(xmin, xmax, color=color, alpha=alpha, zorder=0, linewidth=0)
        xmid = 0.5 * (xmin + xmax)
        ax.text(
            xmid, text_y, label,
            transform=ax.get_xaxis_transform(),
            ha="center", va="top",
            fontsize=fontsize, color="#444444",
            style="italic",
        )


def panel_label(ax, letter, x=0.03, y=0.94, fontsize=10):
    """Add bold (a), (b), ... panel labels to multi-panel figures."""
    ax.text(
        x, y, f"({letter})",
        transform=ax.transAxes,
        fontsize=fontsize, fontweight="bold",
        ha="left", va="top",
    )


def annotate_value(ax, x, y, text, xytext_offset=(10, 10), fontsize=7):
    """Small italic grey annotation with a thin arrow."""
    ax.annotate(
        text, xy=(x, y),
        xytext=(x + xytext_offset[0], y + xytext_offset[1]),
        fontsize=fontsize, color="#444444", style="italic",
        arrowprops=dict(
            arrowstyle="-", color="#444444", lw=0.6, shrinkA=0, shrinkB=2,
        ),
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5),
    )


def textbox(ax, x, y, text, ha="right", va="bottom", fontsize=8):
    """Text box with semi-transparent white background (for chi2, etc.)."""
    ax.text(
        x, y, text,
        transform=ax.transAxes,
        ha=ha, va=va, fontsize=fontsize,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=2),
    )


def finalize(fig, name: str, outdir: str = "results/plots"):
    """
    Save figure as PDF + PNG at 300 dpi. Also save a greyscale version
    for colorblind/greyscale-readability verification.
    """
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    fig.savefig(out / f"{name}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(out / f"{name}.png", bbox_inches="tight", dpi=300)

    # Greyscale check
    try:
        from PIL import Image, ImageOps
        img = ImageOps.grayscale(Image.open(out / f"{name}.png"))
        img.save(out / f"{name}_grey.png")
    except ImportError:
        pass

    print(f"Saved: {out / name}.pdf")
    plt.close(fig)


# -----------------------------------------------------------------------------
# EXAMPLE: FIGURE 4 — Grieco validation overlay (full reference implementation)
# -----------------------------------------------------------------------------
def figure_grieco_validation(
    df_kmc,         # DataFrame with columns T_K, eps_mean, eps_sem
    df_expt,        # DataFrame with columns T_K, eps, eps_err (Grieco data)
    chi2_red: float = 0.025,
    outname: str = "fig04_grieco_validation",
):
    """
    Reference implementation. This is what 'good' looks like —
    use the same pattern for every other figure.
    """
    setup_style()
    fig, (ax1, ax2) = plt.subplots(
        1, 2,
        figsize=fig_double(height=3.2),
        gridspec_kw={"width_ratios": [2.2, 1.0], "wspace": 0.25},
    )

    # ---- Panel (a): overlay ----
    plot_with_ci(
        ax1,
        df_kmc["T_K"].values, df_kmc["eps_mean"].values,
        y_err=df_kmc["eps_sem"].values,
        color=COLORS["blue"],
        linewidth=1.8,
        label="This work (KMC)",
    )
    plot_data_points(
        ax1,
        df_expt["T_K"].values, df_expt["eps"].values,
        df_expt["eps_err"].values,
        color=COLORS["black"],
        label="Grieco et al. (2023)",
    )

    ax1.set_xlim(0, 270)
    ax1.set_ylim(0, 0.35)
    ax1.set_xlabel("Grain temperature (K)")
    ax1.set_ylabel(r"H$_2$ formation efficiency $\epsilon$")
    ax1.legend(loc="upper right")
    panel_label(ax1, "a")
    textbox(
        ax1, 0.96, 0.13,
        rf"$\chi^2_{{\mathrm{{red}}}} = {chi2_red:.3f}$" + "\n(100–250 K)",
        ha="right", va="bottom",
        fontsize=7,
    )

    # ---- Panel (b): residuals ----
    import pandas as pd
    merged = pd.merge_asof(
        df_expt.sort_values("T_K"),
        df_kmc.sort_values("T_K")[["T_K", "eps_mean"]].rename(columns={"eps_mean": "eps_model"}),
        on="T_K", direction="nearest",
    )
    residuals = merged["eps_model"] - merged["eps"]

    ax2.axhline(0, color="black", lw=0.8, zorder=0)
    for _, row in merged.iterrows():
        r = row["eps_model"] - row["eps"]
        c = COLORS["blue"] if r > 0 else COLORS["vermillion"]
        ax2.vlines(row["T_K"], 0, r, color=c, lw=1.3, alpha=0.7, zorder=1)

    plot_data_points(
        ax2,
        merged["T_K"].values, residuals.values,
        merged["eps_err"].values,
        color=COLORS["black"],
    )

    ax2.set_xlim(0, 270)
    ax2.set_xlabel("Grain temperature (K)")
    ax2.set_ylabel(r"$\epsilon_{\mathrm{model}} - \epsilon_{\mathrm{expt}}$")
    panel_label(ax2, "b")

    finalize(fig, outname)


# -----------------------------------------------------------------------------
# EXAMPLE: FIGURE 7 — ε(T) at all four densities (the money plot)
# -----------------------------------------------------------------------------
def figure_epsilon_all_densities(
    df,             # DataFrame with T_K, nH, eps_mean, eps_sem (G0=0)
    analytic_limit: float = 0.18,
    outname: str = "fig07_epsilon_all_densities",
):
    setup_style()
    fig, ax = plt.subplots(figsize=fig_double(height=3.5))

    label_positions = {}
    # Four density curves in canonical order
    for n_H in [10, 100, 1000, 10000]:
        sub = df[df["nH"] == n_H].sort_values("T_K")
        if len(sub) == 0:
            continue
        plot_with_ci(
            ax,
            sub["T_K"].values, sub["eps_mean"].values,
            y_err=sub["eps_sem"].values,
            color=DENSITY_COLOR[n_H],
            marker=DENSITY_MARKER[n_H],
            linewidth=1.5,
            markersize=4.5,
        )
        label_positions[n_H] = (float(sub["T_K"].values[-1]), float(sub["eps_mean"].values[-1]))

    # Analytic high-T limit
    ax.axhline(analytic_limit, xmin=0.48, xmax=1.0,
               color="black", linestyle="--", linewidth=0.8, zorder=1)

    ax.set_xlim(5, 270)
    ax.set_ylim(0, 0.35)
    ax.set_xlabel("Grain surface temperature (K)")
    ax.set_ylabel(r"H$_2$ formation efficiency $\epsilon$")
    offsets = {10: -0.004, 100: -0.0015, 1000: 0.0015, 10000: 0.004}
    for n_H, (_, y_end) in label_positions.items():
        ax.text(
            252,
            y_end + offsets[n_H],
            rf"$10^{{{int(np.log10(n_H))}}}$",
            color=DENSITY_COLOR[n_H],
            fontsize=7,
            ha="right",
            va="center",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=0.15),
        )
    ax.text(
        245,
        analytic_limit + 0.004,
        "analytic limit",
        color=COLORS["black"],
        fontsize=7,
        ha="right",
        va="bottom",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=0.15),
    )

    finalize(fig, outname)


# -----------------------------------------------------------------------------
# EXAMPLE: FIGURE 11 — 2D phase map ε(T, n_H)
# -----------------------------------------------------------------------------
def figure_phase_map(
    df,                     # T_K, nH, eps_mean
    outname: str = "fig11_phase_map",
):
    setup_style()
    from scipy.interpolate import griddata

    T_grid = np.linspace(10, 250, 200)
    n_grid = np.logspace(1, 4, 150)
    TT, NN = np.meshgrid(T_grid, n_grid)

    eps_grid = griddata(
        (df["T_K"].values, np.log10(df["nH"].values)),
        df["eps_mean"].values,
        (TT, np.log10(NN)),
        method="linear",
    )

    fig, ax = plt.subplots(figsize=fig_single(height=3.0))
    im = ax.pcolormesh(
        T_grid, n_grid, eps_grid,
        cmap="viridis", vmin=0, vmax=0.32,
        shading="gouraud",
    )
    ax.set_yscale("log")

    cbar = fig.colorbar(im, ax=ax, pad=0.02, aspect=25)
    cbar.set_label(r"$\epsilon$", rotation=0, labelpad=10)
    cbar.ax.tick_params(direction="in", width=0.6)

    ax.set_xlabel("Grain temperature (K)")
    ax.set_ylabel(r"Gas density $n_{\rm H}$ (cm$^{-3}$)")
    ax.set_xlim(10, 250)
    ax.set_ylim(10, 1e4)

    finalize(fig, outname)


# -----------------------------------------------------------------------------
# EXAMPLE: FIGURE 16 — Ensemble convergence (appendix)
# -----------------------------------------------------------------------------
def figure_ensemble_convergence(
    eps_raw,                # 1D array of eps from 1000 independent realisations
    outname: str = "fig16_ensemble_convergence",
    production_N: int = 20,
):
    setup_style()
    N = np.arange(1, len(eps_raw) + 1)
    running_mean = np.cumsum(eps_raw) / N
    running_sem = np.array([
        eps_raw[:n].std(ddof=1) / np.sqrt(n) if n > 1 else np.nan
        for n in N
    ])

    fig, ax = plt.subplots(figsize=fig_single(height=2.6))

    plot_with_ci(
        ax, N, running_mean,
        y_err=running_sem,
        color=COLORS["blue"],
        linewidth=1.5,
        ci_alpha=0.25,
    )

    asymptotic = running_mean[-1]
    ax.axhline(asymptotic, color="black", linestyle=":", linewidth=0.7, zorder=1)

    ax.axvline(production_N, color=COLORS["vermillion"],
               linestyle="--", linewidth=0.8, zorder=1)
    ax.text(
        production_N * 1.15, ax.get_ylim()[0] + 0.02,
        f"Production ensemble\n(N = {production_N})",
        fontsize=7, color=COLORS["vermillion"], style="italic",
    )

    ax.set_xscale("log")
    ax.set_xlabel("Number of ensemble realisations")
    ax.set_ylabel(r"Running mean $\langle\epsilon\rangle$")

    finalize(fig, outname)
