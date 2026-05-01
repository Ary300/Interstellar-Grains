#!/usr/bin/env python3
"""
Summarize sensitivity to lh_diffusion_factor and diffusion_rate_cap_s.
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from paper_plot_style import DOUBLE_COLUMN_WIDTH, SINGLE_COLUMN_WIDTH, apply_publication_style, save_figure, style_axes, style_legend


def _pick_col(df: pd.DataFrame, base: str) -> str:
    if base in df.columns:
        return base
    mean_col = f"{base}_mean"
    if mean_col in df.columns:
        return mean_col
    raise KeyError(f"Missing column: {base} or {mean_col}")


def main() -> None:
    p = argparse.ArgumentParser(description="Plot sensitivity to lh_diffusion_factor and diffusion_rate_cap_s.")
    p.add_argument("--input", default="results/astro_transport_sensitivity.csv")
    p.add_argument("--out-dir", default="results/plots/astro_transport_sensitivity")
    p.add_argument("--out-table", default="results/tables/astro_transport_sensitivity_summary.csv")
    p.add_argument("--baseline-lh", type=float, default=0.5)
    p.add_argument("--baseline-cap", type=float, default=200.0)
    args = p.parse_args()

    df = pd.read_csv(args.input)
    eps_col = _pick_col(df, "epsilon")
    rate_col = _pick_col(df, "h2_release_rate_cm2_s")
    need = {"surface_temperature_k", "lh_diffusion_factor", "diffusion_rate_cap_s", eps_col, rate_col}
    missing = sorted(need - set(df.columns))
    if missing:
        raise SystemExit(f"Missing columns in {args.input}: {missing}")

    for c in ["surface_temperature_k", "lh_diffusion_factor", "diffusion_rate_cap_s", eps_col, rate_col]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["surface_temperature_k", "lh_diffusion_factor", "diffusion_rate_cap_s", eps_col, rate_col])

    base = df[
        np.isclose(df["lh_diffusion_factor"], float(args.baseline_lh))
        & np.isclose(df["diffusion_rate_cap_s"], float(args.baseline_cap))
    ][["surface_temperature_k", eps_col, rate_col]].rename(columns={eps_col: "epsilon_baseline", rate_col: "rate_baseline"})
    if base.empty:
        raise SystemExit("Baseline transport parameters not found in the input CSV.")

    merged = df.merge(base, on="surface_temperature_k", how="left")
    merged["epsilon_ratio_to_baseline"] = merged[eps_col] / merged["epsilon_baseline"]
    merged["rate_ratio_to_baseline"] = merged[rate_col] / merged["rate_baseline"]

    os.makedirs(os.path.dirname(args.out_table) or ".", exist_ok=True)
    merged.sort_values(["surface_temperature_k", "lh_diffusion_factor", "diffusion_rate_cap_s"]).to_csv(args.out_table, index=False)
    print(f"Wrote {args.out_table}")

    os.makedirs(args.out_dir, exist_ok=True)

    apply_publication_style()

    temps = sorted(merged["surface_temperature_k"].unique())
    for ratio_col, title, out_name in [
        ("epsilon_ratio_to_baseline", "Transport sensitivity: ε / baseline", "transport_sensitivity_epsilon_ratio.png"),
        ("rate_ratio_to_baseline", "Transport sensitivity: rate / baseline", "transport_sensitivity_rate_ratio.png"),
    ]:
        ratio_vals = pd.to_numeric(merged[ratio_col], errors="coerce").to_numpy(dtype=float)
        if not np.any(np.isfinite(ratio_vals)):
            fig, ax = plt.subplots(figsize=(SINGLE_COLUMN_WIDTH, 2.55))
            style_axes(ax)
            ax.text(
                0.5,
                0.58,
                "No finite ratio values were produced\nfor this transport-sensitivity sweep.",
                ha="center",
                va="center",
                fontsize=8.4,
                transform=ax.transAxes,
            )
            ax.text(
                0.5,
                0.38,
                "In the current sweep, the selected conditions yielded zero formation,\nso no transport-response curve can be measured here.",
                ha="center",
                va="center",
                fontsize=7.0,
                color="#444444",
                transform=ax.transAxes,
            )
            ax.set_xticks([])
            ax.set_yticks([])
            out_png = os.path.join(args.out_dir, os.path.splitext(out_name)[0])
            fig.subplots_adjust(left=0.1, right=0.98, bottom=0.12, top=0.96)
            save_figure(fig, out_png)
            plt.close(fig)
            print(f"Wrote {out_png}.png/.pdf")
            continue

        fig, axes = plt.subplots(2, 2, figsize=(DOUBLE_COLUMN_WIDTH, 4.8), sharex=True, sharey=True)
        axes = axes.flatten()
        for ax, temp in zip(axes, temps):
            sub = merged[np.isclose(merged["surface_temperature_k"], float(temp))].sort_values(["lh_diffusion_factor", "diffusion_rate_cap_s"])
            for lh_fac in sorted(sub["lh_diffusion_factor"].unique()):
                grp = sub[np.isclose(sub["lh_diffusion_factor"], float(lh_fac))].sort_values("diffusion_rate_cap_s")
                ax.plot(grp["diffusion_rate_cap_s"], grp[ratio_col], marker="o", ms=4.5, lw=2.0, label=fr"$f_{{\rm LH}}={lh_fac:g}$")
            ax.axhline(1.0, color="k", linewidth=1.0, alpha=0.5)
            ax.text(0.05, 0.92, fr"$T={temp:g}$ K", transform=ax.transAxes, fontsize=7.2, color="#444444")
            style_axes(ax)
        style_legend(axes[0], loc="best")
        fig.supxlabel("diffusion_rate_cap_s")
        fig.supylabel("ratio to baseline")
        fig.subplots_adjust(left=0.09, right=0.985, bottom=0.12, top=0.98, wspace=0.14, hspace=0.18)
        out_png = os.path.join(args.out_dir, os.path.splitext(out_name)[0])
        save_figure(fig, out_png)
        plt.close(fig)
        print(f"Wrote {out_png}.png/.pdf")


if __name__ == "__main__":
    main()
