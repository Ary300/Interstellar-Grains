#!/usr/bin/env python3
"""
Summarize the consistency between explicit-pairs LH and diffusion-limited LH.
"""

from __future__ import annotations

import argparse
import os

import pandas as pd
import matplotlib.pyplot as plt

from paper_plot_style import SINGLE_COLUMN_WIDTH, apply_publication_style, save_figure, style_axes, style_legend


def _pick_metric(df: pd.DataFrame, base: str) -> str:
    if base in df.columns:
        return base
    mean_col = f"{base}_mean"
    if mean_col in df.columns:
        return mean_col
    raise KeyError(f"Missing column: {base} or {base}_mean")


def main() -> None:
    p = argparse.ArgumentParser(description="Plot LH-mode consistency check.")
    p.add_argument("--input", default="results/astro_lh_mode_consistency.csv")
    p.add_argument("--out-dir", default="results/plots/astro_lh_mode_consistency")
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    df = pd.read_csv(args.input)

    eps_col = _pick_metric(df, "epsilon")
    rate_col = _pick_metric(df, "h2_release_rate_cm2_s")
    need = {"surface_temperature_k", "diffusion_mode", "lh_formation_mode", eps_col, rate_col}
    missing = sorted(need - set(df.columns))
    if missing:
        raise SystemExit(f"Missing columns in {args.input}: {missing}")

    df = df.copy()
    df["mode_label"] = df["diffusion_mode"].astype(str) + " / " + df["lh_formation_mode"].astype(str)

    pivot_eps = df.pivot_table(index="surface_temperature_k", columns="mode_label", values=eps_col, aggfunc="first")
    pivot_rate = df.pivot_table(index="surface_temperature_k", columns="mode_label", values=rate_col, aggfunc="first")

    pairs_col = "explicit / pairs"
    dl_col = "rate_only / diffusion_limited"
    for col in [pairs_col, dl_col]:
        if col not in pivot_eps.columns:
            raise SystemExit(f"Expected LH mode column '{col}' in {args.input}")

    summary = pivot_eps[[pairs_col, dl_col]].rename(columns={pairs_col: "epsilon_pairs", dl_col: "epsilon_diffusion_limited"})
    summary["epsilon_frac_diff"] = summary["epsilon_diffusion_limited"] / summary["epsilon_pairs"] - 1.0
    summary = summary.join(
        pivot_rate[[pairs_col, dl_col]].rename(columns={pairs_col: "rate_pairs", dl_col: "rate_diffusion_limited"})
    )
    summary["rate_frac_diff"] = summary["rate_diffusion_limited"] / summary["rate_pairs"] - 1.0
    out_csv = os.path.join(args.out_dir, "lh_mode_consistency_summary.csv")
    summary.reset_index().to_csv(out_csv, index=False)
    print(f"Wrote {out_csv}")

    apply_publication_style()
    fig, ax = plt.subplots(figsize=(SINGLE_COLUMN_WIDTH, 2.55))
    ax.plot(summary.index, 100.0 * summary["epsilon_frac_diff"], marker="o", ms=5, lw=2.2, color="#1d3557", label=r"$\epsilon$")
    ax.plot(summary.index, 100.0 * summary["rate_frac_diff"], marker="s", ms=5, lw=2.2, color="#c05640", label="Rate")
    ax.axhline(0.0, color="k", linewidth=1.0)
    ax.axhline(10.0, color="gray", linestyle="--", linewidth=1.0)
    ax.axhline(-10.0, color="gray", linestyle="--", linewidth=1.0)
    style_axes(ax)
    ax.set_xlabel("Surface temperature (K)")
    ax.set_ylabel("Difference vs pairs (%)")
    style_legend(ax, loc="best")
    fig.subplots_adjust(left=0.2, right=0.98, bottom=0.18, top=0.98)
    out_png = os.path.join(args.out_dir, "lh_mode_consistency")
    save_figure(fig, out_png)
    plt.close(fig)
    print(f"Wrote {out_png}.png/.pdf")


if __name__ == "__main__":
    main()
