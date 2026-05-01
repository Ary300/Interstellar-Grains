#!/usr/bin/env python3
"""
Plot UV response diagnostics from a merged ISM sweep CSV.

Outputs:
  - heatmaps of uv/high-to-uv0 rate ratio
  - a CSV summary table with min/median/max ratios by UV level
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd


def _ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def main() -> None:
    p = argparse.ArgumentParser(description="Plot UV response diagnostics from a merged sweep CSV.")
    p.add_argument("--input", default="jhub_full_merged.csv", help="Merged CSV")
    p.add_argument("--out-dir", default="results/plots/uv_response", help="Output directory")
    p.add_argument("--out-table", default="results/tables/uv_response_summary.csv", help="Output CSV")
    p.add_argument("--metric", default="h2_release_rate_cm2_s_mean", help="Metric column to compare across UV")
    args = p.parse_args()

    df = pd.read_csv(args.input)
    metric = args.metric
    if metric not in df.columns:
        raise SystemExit(f"Missing metric column: {metric}")
    for c in ["surface_temperature_k", "h_gas_density_cm3", "uv_flux_factor"]:
        if c not in df.columns:
            raise SystemExit(f"Missing required column: {c}")

    base = df[df["uv_flux_factor"] == 0.0][["surface_temperature_k", "h_gas_density_cm3", metric]].rename(
        columns={metric: "metric_uv0"}
    )

    rows = []
    merged_by_uv = {}
    for uv in sorted(float(x) for x in df["uv_flux_factor"].unique() if float(x) != 0.0):
        sub = df[df["uv_flux_factor"] == uv][["surface_temperature_k", "h_gas_density_cm3", metric]]
        m = sub.merge(base, on=["surface_temperature_k", "h_gas_density_cm3"], how="inner")
        m["ratio_to_uv0"] = m[metric] / m["metric_uv0"]
        merged_by_uv[uv] = m
        rows.append(
            {
                "uv_flux_factor": uv,
                "ratio_min": float(m["ratio_to_uv0"].min()),
                "ratio_median": float(m["ratio_to_uv0"].median()),
                "ratio_max": float(m["ratio_to_uv0"].max()),
            }
        )

    summary = pd.DataFrame(rows).sort_values("uv_flux_factor")
    _ensure_dir(os.path.dirname(args.out_table) or ".")
    summary.to_csv(args.out_table, index=False)
    print(f"Wrote {args.out_table}")

    _ensure_dir(args.out_dir)
    import matplotlib.pyplot as plt
    import seaborn as sns

    for uv, m in merged_by_uv.items():
        pivot = m.pivot_table(index="h_gas_density_cm3", columns="surface_temperature_k", values="ratio_to_uv0", aggfunc="mean")
        pivot = pivot.sort_index().sort_index(axis=1)
        plt.figure(figsize=(10, 4.5))
        sns.heatmap(pivot, cmap="coolwarm", center=1.0, cbar_kws={"label": f"{metric}(UV={uv}) / {metric}(UV=0)"})
        plt.xlabel("Surface temperature (K)")
        plt.ylabel(r"n(H) [cm$^{-3}$]")
        plt.title(f"UV response heatmap (UV={uv})")
        plt.tight_layout()
        out_png = os.path.join(args.out_dir, f"uv_ratio_heatmap_{int(uv)}.png")
        plt.savefig(out_png, dpi=200)
        plt.close()
        print(f"Wrote {out_png}")


if __name__ == "__main__":
    main()

