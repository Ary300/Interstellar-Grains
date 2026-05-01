from __future__ import annotations

import argparse
import csv
import os
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from paper_plot_style import DOUBLE_COLUMN_WIDTH, add_panel_labels, apply_publication_style, save_figure, style_axes, style_legend


def _read_xy_ci(
    path: str,
    x_col: str,
    y_col: str,
    ci_col: str | None,
) -> Tuple[List[float], List[float], List[float]]:
    xs: List[float] = []
    ys: List[float] = []
    cis: List[float] = []
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            try:
                xs.append(float(row[x_col]))
                ys.append(float(row[y_col]))
                if ci_col and ci_col in row and str(row[ci_col]).strip() != "":
                    cis.append(abs(float(row[ci_col])))
                else:
                    cis.append(0.0)
            except Exception:
                continue
    return xs, ys, cis


def plot(
    iso_csv: str,
    ded_csv: str,
    out_png: str,
    y_lines: List[float],
    digitized_csv: str | None = None,
) -> None:
    apply_publication_style()
    fig, (ax, ax_res) = plt.subplots(
        1,
        2,
        figsize=(DOUBLE_COLUMN_WIDTH, 2.95),
        gridspec_kw={"width_ratios": [2.2, 1.0], "wspace": 0.27},
    )

    residual_sets: list[tuple[pd.Series, np.ndarray, np.ndarray, str]] = []
    if digitized_csv and os.path.exists(digitized_csv):
        digitized = pd.read_csv(digitized_csv)
        digitized["T_K"] = pd.to_numeric(digitized["T_K"], errors="coerce")
        digitized["epsilon"] = pd.to_numeric(digitized["epsilon"], errors="coerce")
        digitized["epsilon_err_low"] = pd.to_numeric(digitized["epsilon_err_low"], errors="coerce")
        digitized["epsilon_err_high"] = pd.to_numeric(digitized["epsilon_err_high"], errors="coerce")
        for method, marker, color, label in [
            ("TPDED", "o", "#c05640", "Grieco TPDED"),
            ("isothermal", "s", "#1d3557", "Grieco isothermal"),
        ]:
            sub = digitized[digitized["method"] == method]
            if sub.empty:
                continue
            ax.errorbar(
                sub["T_K"],
                sub["epsilon"],
                yerr=[sub["epsilon_err_low"], sub["epsilon_err_high"]],
                fmt=marker,
                ms=5.2,
                color=color,
                mfc="white",
                mec=color,
                mew=1.0,
                elinewidth=0.85,
                capsize=2,
                linestyle="none",
                label=label,
                zorder=4,
            )

    if os.path.exists(ded_csv):
        # Prefer released_total ε if available (TPDED ramps can include delayed desorption of formed H2).
        # Fall back to the legacy `epsilon_mean` column for older outputs.
        try:
            with open(ded_csv, newline="") as f:
                cols = (csv.DictReader(f).fieldnames or [])
        except Exception:
            cols = []
        if "epsilon_released_total_mean" in cols:
            x, y, ci = _read_xy_ci(
                ded_csv,
                x_col="temperature_k",
                y_col="epsilon_released_total_mean",
                ci_col="epsilon_released_total_ci95",
            )
            label = "KMC TPDED"
        else:
            x, y, ci = _read_xy_ci(ded_csv, x_col="temperature_k", y_col="epsilon_mean", ci_col="epsilon_ci95")
            label = "KMC TPDED"
        if x:
            ax.plot(x, y, color="#c05640", lw=1.9, alpha=0.95, zorder=3)
            ax.fill_between(x, [yy - cc for yy, cc in zip(y, ci)], [yy + cc for yy, cc in zip(y, ci)], color="#c05640", alpha=0.16, zorder=2)
            ax.scatter(x, y, s=18, color="#c05640", edgecolor="white", linewidth=0.45, zorder=4, label=label)
            if digitized_csv and os.path.exists(digitized_csv):
                sub = digitized[digitized["method"] == "TPDED"].dropna(subset=["T_K", "epsilon"])
                if not sub.empty:
                    model = np.interp(sub["T_K"], np.asarray(x, dtype=float), np.asarray(y, dtype=float))
                    residual_sets.append((sub["T_K"], model - sub["epsilon"].to_numpy(dtype=float), sub["epsilon_err_high"].to_numpy(dtype=float), "#c05640"))

    if os.path.exists(iso_csv):
        x, y, ci = _read_xy_ci(iso_csv, x_col="surface_temperature_k", y_col="epsilon_mean", ci_col="epsilon_ci95")
        if x:
            ax.plot(x, y, color="#1d3557", lw=2.0, alpha=0.95, zorder=3)
            ax.fill_between(x, [yy - cc for yy, cc in zip(y, ci)], [yy + cc for yy, cc in zip(y, ci)], color="#1d3557", alpha=0.16, zorder=2)
            ax.scatter(x, y, s=24, marker="s", color="#1d3557", edgecolor="white", linewidth=0.55, zorder=4, label="KMC isothermal")
            if digitized_csv and os.path.exists(digitized_csv):
                sub = digitized[digitized["method"] == "isothermal"].dropna(subset=["T_K", "epsilon"])
                if not sub.empty:
                    model = np.interp(sub["T_K"], np.asarray(x, dtype=float), np.asarray(y, dtype=float))
                    residual_sets.append((sub["T_K"], model - sub["epsilon"].to_numpy(dtype=float), sub["epsilon_err_high"].to_numpy(dtype=float), "#1d3557"))

    for yl in y_lines:
        ax.axhline(float(yl), color="#84a98c", lw=1.1, alpha=0.45, linestyle="--", zorder=1)

    style_axes(ax)
    ax.grid(False)
    ax.set_xlabel("Surface temperature (K)")
    ax.set_ylabel(r"Formation efficiency $\epsilon$")
    ax.set_xlim(0, 260)
    ax.set_ylim(0, 0.55)
    style_legend(ax, loc="upper center", bbox_to_anchor=(0.56, 0.99), ncol=2)

    ax_res.axhline(0.0, color="black", lw=0.8, zorder=0)
    for temps, residuals, errs, color in residual_sets:
        for temp, resid, err in zip(temps, residuals, errs):
            ax_res.vlines(float(temp), 0.0, float(resid), color=color, lw=1.0)
            ax_res.errorbar(
                [float(temp)],
                [float(resid)],
                yerr=[[float(err)], [float(err)]],
                fmt="o",
                mfc="white",
                mec="black",
                mew=0.8,
                ms=4.0,
                capsize=2,
                elinewidth=0.8,
                color="black",
            )
    style_axes(ax_res)
    ax_res.grid(False)
    ax_res.set_xlabel("Surface temperature (K)")
    ax_res.set_ylabel(r"$\epsilon_{\rm model}-\epsilon_{\rm expt}$")
    ax_res.set_xlim(0, 260)
    ax_res.set_ylim(-0.08, 0.18)

    add_panel_labels([ax, ax_res], labels="AB")

    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.16, top=0.98, wspace=0.28)
    save_figure(fig, out_png)
    plt.close(fig)
    print(f"Wrote {os.path.splitext(out_png)[0]}.png/.pdf")


def main() -> None:
    p = argparse.ArgumentParser(description="Create a Figure-2-like plot from Grieco validation CSVs.")
    p.add_argument("--iso-csv", default="results/grieco_validation.csv", help="CSV from grieco_validation.py")
    p.add_argument("--ded-csv", default="results/grieco_ded_validation.csv", help="CSV from grieco_ded_validation.py")
    p.add_argument("--out", default="results/plots/grieco_figure2_like.png", help="Output PNG path")
    p.add_argument("--digitized-csv", default="grieco_fig2_digitized.csv", help="Digitized Grieco Figure 2 CSV")
    p.add_argument(
        "--y-lines",
        nargs="*",
        type=float,
        default=[],
        help="Horizontal reference lines (fractions).",
    )
    args = p.parse_args()
    plot(
        iso_csv=str(args.iso_csv),
        ded_csv=str(args.ded_csv),
        out_png=str(args.out),
        y_lines=list(args.y_lines),
        digitized_csv=str(args.digitized_csv),
    )


if __name__ == "__main__":
    main()
