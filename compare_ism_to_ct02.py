#!/usr/bin/env python3
"""
Compare an ISM sweep (KMC) to the Cazaux & Tielens (2002) analytic prescription (CT02).

This script is designed for merged aggregated outputs (one row per condition), e.g.:
  - results/jhub_full_merged.csv
  - results/anvil_near_full_merged.csv

Outputs:
  - results/tables/ct02_comparison_<label>.csv
  - results/plots/ct02_<label>/ratio_heatmap_uv0.png
  - results/plots/ct02_<label>/ratio_heatmap_uv<uv_high>.png
  - results/plots/ct02_<label>/rate_overlay_uv0.png
  - results/plots/ct02_<label>/rate_overlay_uv<uv_high>.png
"""

from __future__ import annotations

import argparse
import os
from typing import List

import numpy as np
import pandas as pd

from cazaux_tielens_2002 import CT02Params, epsilon_ct02, rate_ct02_per_area_cm2_s


def _ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def _pick_col(df: pd.DataFrame, base: str) -> str:
    if base in df.columns:
        return base
    m = f"{base}_mean"
    if m in df.columns:
        return m
    raise KeyError(f"Missing column: {base} or {base}_mean")


def _finite(df: pd.DataFrame, cols: List[str]) -> pd.DataFrame:
    out = df.copy()
    for c in cols:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    return out.dropna(subset=cols)


def _heatmap(
    *,
    df: pd.DataFrame,
    value_col: str,
    title: str,
    out_png: str,
) -> None:
    import matplotlib.pyplot as plt
    import seaborn as sns

    pivot = df.pivot_table(index="h_gas_density_cm3", columns="surface_temperature_k", values=value_col, aggfunc="mean")
    pivot = pivot.sort_index().sort_index(axis=1)

    plt.figure(figsize=(12, 4.5))
    sns.heatmap(pivot, cmap="coolwarm", center=0.0, cbar_kws={"label": value_col})
    plt.xlabel("Surface temperature (K)")
    plt.ylabel(r"n(H) [cm$^{-3}$]")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def _rate_overlay(
    *,
    df: pd.DataFrame,
    uv: float,
    out_png: str,
) -> None:
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 6))
    variant_label = str(df["ct_model_variant"].iloc[0]) if "ct_model_variant" in df.columns and not df.empty else "ct02"
    for n in sorted(df["h_gas_density_cm3"].unique()):
        sub = df[(df["uv_flux_factor"] == uv) & (df["h_gas_density_cm3"] == n)].sort_values("surface_temperature_k")
        if sub.empty:
            continue
        plt.plot(sub["surface_temperature_k"], sub["h2_release_rate_cm2_s_mean"], marker="o", label=f"KMC n={int(n)}")
        plt.plot(sub["surface_temperature_k"], sub["ct_h2_release_rate_cm2_s"], linestyle="--", label=f"{variant_label} n={int(n)}")

    plt.yscale("log")
    plt.xlabel("Surface temperature (K)")
    plt.ylabel(r"H$_2$ release rate (cm$^{-2}$ s$^{-1}$)")
    plt.title(f"ISM H2 formation/release rate: KMC vs {variant_label} (uv_flux_factor={uv})")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend(ncol=2, fontsize=8)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def main() -> None:
    p = argparse.ArgumentParser(description="Compare ISM sweep outputs to CT02 analytic rates.")
    p.add_argument("--input", default="results/jhub_full_merged.csv", help="Merged aggregated CSV from run_sweep.py/anvil merge")
    p.add_argument("--label", default="jhub_full", help="Label for output files")
    p.add_argument(
        "--ct-variant",
        default="erratum2010_approx",
        choices=["original", "erratum2010_approx"],
        help="Which analytic Cazaux-Tielens variant to compare against",
    )
    p.add_argument("--uv-high", type=float, default=100.0, help="Second UV level to plot alongside UV=0")
    p.add_argument("--out-tables", default="results/tables", help="Output directory for CSV tables")
    p.add_argument("--out-plots", default="", help="Output directory for plots (default: results/plots/ct02_<label>)")
    args = p.parse_args()

    df = pd.read_csv(args.input)

    # Required descriptor columns
    for c in ["surface_temperature_k", "h_gas_density_cm3", "uv_flux_factor"]:
        if c not in df.columns:
            raise SystemExit(f"Missing required column: {c}")

    h2_rate_col = _pick_col(df, "h2_release_rate_cm2_s")
    stick_col = "sticking_probability" if "sticking_probability" in df.columns else None
    if stick_col is None:
        # Fall back to config defaults used in the ISM sweeps if missing.
        df["sticking_probability"] = 0.5
        stick_col = "sticking_probability"

    # Prefer gas_flux_h_cm2_s / arrival_rate_per_site_s if present; otherwise reconstruct.
    if "gas_flux_h_cm2_s" not in df.columns or "arrival_rate_per_site_s" not in df.columns:
        # Reconstruct using the same formula as run_sweep.py.
        from scientific_data import K_B_ERG, M_H

        df["gas_temperature_k"] = pd.to_numeric(df.get("gas_temperature_k", 100.0), errors="coerce").fillna(100.0)
        df["h_gas_density_cm3"] = pd.to_numeric(df["h_gas_density_cm3"], errors="coerce")
        v_th = np.sqrt(8.0 * float(K_B_ERG) * df["gas_temperature_k"] / (np.pi * float(M_H)))
        df["gas_flux_h_cm2_s"] = 0.25 * df["h_gas_density_cm3"] * v_th
        site_area_cm2 = pd.to_numeric(df.get("site_area_angstroms_sq", 25.0), errors="coerce").fillna(25.0) * 1e-16
        df["arrival_rate_per_site_s"] = df["gas_flux_h_cm2_s"] * site_area_cm2

    df = _finite(
        df,
        [
            "surface_temperature_k",
            "h_gas_density_cm3",
            "uv_flux_factor",
            h2_rate_col,
            "gas_flux_h_cm2_s",
            "arrival_rate_per_site_s",
            stick_col,
        ],
    )

    params = CT02Params()

    # Compute CT02 epsilon and rate per area for each condition.
    df["ct_model_variant"] = str(args.ct_variant)
    df["ct_epsilon"] = df.apply(
        lambda r: epsilon_ct02(
            float(r["surface_temperature_k"]),
            float(r["arrival_rate_per_site_s"]),
            params,
            variant=args.ct_variant,
        ),
        axis=1,
    )
    df["ct_h2_release_rate_cm2_s"] = df.apply(
        lambda r: rate_ct02_per_area_cm2_s(
            gas_flux_cm2_s=float(r["gas_flux_h_cm2_s"]),
            sticking=float(r[stick_col]),
            T_surf_K=float(r["surface_temperature_k"]),
            site_area_angstrom2=float(r.get("site_area_angstroms_sq", 25.0)),
            params=params,
            variant=args.ct_variant,
        ),
        axis=1,
    )
    df["ct02_epsilon"] = df["ct_epsilon"]
    df["ct02_h2_release_rate_cm2_s"] = df["ct_h2_release_rate_cm2_s"]

    df["h2_release_rate_cm2_s_mean"] = pd.to_numeric(df[h2_rate_col], errors="coerce")
    df["ratio_kmc_over_ct02"] = df["h2_release_rate_cm2_s_mean"] / df["ct_h2_release_rate_cm2_s"].replace(0.0, np.nan)
    df["log10_ratio_kmc_over_ct02"] = np.log10(df["ratio_kmc_over_ct02"])

    out_tables = args.out_tables
    _ensure_dir(out_tables)
    out_csv = os.path.join(out_tables, f"ct02_comparison_{args.label}.csv")
    keep = [
        "surface_temperature_k",
        "h_gas_density_cm3",
        "uv_flux_factor",
        "gas_temperature_k" if "gas_temperature_k" in df.columns else None,
        "ct_model_variant",
        "gas_flux_h_cm2_s",
        "arrival_rate_per_site_s",
        stick_col,
        "ct_epsilon",
        "ct_h2_release_rate_cm2_s",
        "ct02_epsilon",
        "ct02_h2_release_rate_cm2_s",
        "h2_release_rate_cm2_s_mean",
        "ratio_kmc_over_ct02",
        "log10_ratio_kmc_over_ct02",
    ]
    keep = [c for c in keep if c is not None and c in df.columns]
    df[keep].to_csv(out_csv, index=False)
    print(f"Wrote {out_csv}")

    out_plots = args.out_plots or os.path.join("results/plots", f"ct02_{args.label}")
    _ensure_dir(out_plots)

    # Heatmaps for UV=0 and UV=uv_high (log10 ratio).
    for uv in [0.0, float(args.uv_high)]:
        sub = df[df["uv_flux_factor"] == uv]
        if sub.empty:
            continue
        hm_png = os.path.join(out_plots, f"ratio_heatmap_uv{int(uv)}.png")
        _heatmap(
            df=sub,
            value_col="log10_ratio_kmc_over_ct02",
            title=f"log10(KMC/CT02) enhancement (uv_flux_factor={uv})",
            out_png=hm_png,
        )
        print(f"Wrote {hm_png}")

        overlay_png = os.path.join(out_plots, f"rate_overlay_uv{int(uv)}.png")
        _rate_overlay(df=df, uv=uv, out_png=overlay_png)
        print(f"Wrote {overlay_png}")


if __name__ == "__main__":
    main()
