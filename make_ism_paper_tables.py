#!/usr/bin/env python3
"""
Generate paper-ready tables for an ISM (gas_kinetic arrival-mode) sweep.

This mirrors the "Anvil campaign" tables but targets the ISM-focused metrics:
  - h2_release_rate_cm2_s_mean (preferred for ISM)
  - epsilon_mean (dimensionless, still useful for microphysics comparisons)

Example:
  python make_ism_paper_tables.py \
    --input results/jhub_full_merged.csv \
    --label jhub_full \
    --config config_astro_full_paperfit.yaml

Outputs (under results/tables/ by default):
  - table_campaign_overview_<label>.csv
  - table_uv_summary_<label>.csv
  - table_uv_suppression_<label>.csv
  - table_top25_<label>.csv
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

import numpy as np
import pandas as pd
import yaml


def _read_yaml(path: str) -> Dict[str, Any]:
    if not path:
        return {}
    with open(path, "r") as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise TypeError(f"Expected mapping YAML in {path}, got {type(data)}")
    return data


def _maybe_col(df: pd.DataFrame, name: str) -> Optional[str]:
    if name in df.columns:
        return name
    if f"{name}_mean" in df.columns:
        return f"{name}_mean"
    return None


def _required(df: pd.DataFrame, cols: Iterable[str]) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise SystemExit(f"Missing required columns: {missing}")


def _safe_div(a: float, b: float) -> float:
    return float(a) / float(b) if float(b) != 0 else float("nan")


def _ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def main() -> None:
    p = argparse.ArgumentParser(description="Make ISM sweep summary tables from a merged CSV.")
    p.add_argument("--input", default="results/jhub_full_merged.csv", help="Merged aggregated CSV (one row per condition)")
    p.add_argument("--label", default="jhub_full", help="Label to suffix output table filenames")
    p.add_argument("--config", default="", help="Optional YAML config used to generate the sweep (for metadata in overview)")
    p.add_argument("--out-dir", default="results/tables", help="Output directory for tables")
    p.add_argument("--uv-high", type=float, default=100.0, help="UV level used for 'suppression' comparisons (vs UV=0)")
    args = p.parse_args()

    df = pd.read_csv(args.input)
    _required(df, ["surface_temperature_k", "h_gas_density_cm3", "uv_flux_factor"])

    metric_release = _maybe_col(df, "h2_release_rate_cm2_s")
    metric_total_h2 = _maybe_col(df, "total_h2_formed")
    metric_eps = _maybe_col(df, "epsilon")

    metric_er = _maybe_col(df, "h2_formed_ER")
    metric_lh = _maybe_col(df, "h2_formed_LH")
    metric_uv = _maybe_col(df, "h2_formed_UV")

    if metric_release is None and metric_total_h2 is None:
        raise SystemExit("Expected at least one of: h2_release_rate_cm2_s(_mean) or total_h2_formed(_mean)")

    cfg = _read_yaml(args.config) if args.config else {}

    _ensure_dir(args.out_dir)

    # --- Campaign overview ---
    rows = int(len(df))
    nT = int(df["surface_temperature_k"].nunique(dropna=True))
    nN = int(df["h_gas_density_cm3"].nunique(dropna=True))
    nU = int(df["uv_flux_factor"].nunique(dropna=True))

    def _frac(numer_col: Optional[str], denom_col: Optional[str]) -> float:
        if numer_col is None or denom_col is None:
            return float("nan")
        denom = pd.to_numeric(df[denom_col], errors="coerce")
        numer = pd.to_numeric(df[numer_col], errors="coerce")
        with np.errstate(divide="ignore", invalid="ignore"):
            fr = (numer / denom).replace([np.inf, -np.inf], np.nan)
        return float(np.nanmean(fr))

    er_frac_mean = _frac(metric_er, metric_total_h2)
    lh_frac_mean = _frac(metric_lh, metric_total_h2)

    def _stat(col: Optional[str], fn: str) -> float:
        if col is None:
            return float("nan")
        s = pd.to_numeric(df[col], errors="coerce")
        if fn == "min":
            return float(np.nanmin(s))
        if fn == "max":
            return float(np.nanmax(s))
        if fn == "median":
            return float(np.nanmedian(s))
        if fn == "mean":
            return float(np.nanmean(s))
        raise ValueError(fn)

    overview = {
        "dataset": args.label,
        "rows": rows,
        "n_unique_T": nT,
        "n_unique_nH": nN,
        "n_unique_uv": nU,
        "ensemble_runs": int(cfg.get("ensemble_runs", np.nan)) if cfg else float("nan"),
        "burnin_arrivals": int(cfg.get("burnin_arrivals", np.nan)) if cfg else float("nan"),
        "measure_arrivals": int(cfg.get("measure_arrivals", np.nan)) if cfg else float("nan"),
        "er_frac_mean": er_frac_mean,
        "lh_frac_mean": lh_frac_mean,
        "h2_release_rate_cm2_s_min": _stat(metric_release, "min"),
        "h2_release_rate_cm2_s_median": _stat(metric_release, "median"),
        "h2_release_rate_cm2_s_max": _stat(metric_release, "max"),
        "total_h2_formed_min": _stat(metric_total_h2, "min"),
        "total_h2_formed_median": _stat(metric_total_h2, "median"),
        "total_h2_formed_max": _stat(metric_total_h2, "max"),
        "epsilon_min": _stat(metric_eps, "min"),
        "epsilon_median": _stat(metric_eps, "median"),
        "epsilon_max": _stat(metric_eps, "max"),
    }

    out_overview = os.path.join(args.out_dir, f"table_campaign_overview_{args.label}.csv")
    pd.DataFrame([overview]).to_csv(out_overview, index=False)
    print(f"Wrote {out_overview}")

    # --- UV summary (grouped by UV level) ---
    group = df.groupby("uv_flux_factor", dropna=False)
    uv_summary = pd.DataFrame(
        {
            "uv_flux_factor": group.size().index,
            "n_conditions": group.size().values,
            "h2_release_rate_cm2_s_mean": group[metric_release].mean().values if metric_release else np.nan,
            "h2_release_rate_cm2_s_median": group[metric_release].median().values if metric_release else np.nan,
            "total_h2_formed_mean": group[metric_total_h2].mean().values if metric_total_h2 else np.nan,
            "total_h2_formed_median": group[metric_total_h2].median().values if metric_total_h2 else np.nan,
            "epsilon_mean": group[metric_eps].mean().values if metric_eps else np.nan,
            "epsilon_median": group[metric_eps].median().values if metric_eps else np.nan,
            "er_frac_mean": (group[metric_er].mean() / group[metric_total_h2].mean()).values
            if (metric_er and metric_total_h2)
            else np.nan,
            "lh_frac_mean": (group[metric_lh].mean() / group[metric_total_h2].mean()).values
            if (metric_lh and metric_total_h2)
            else np.nan,
        }
    ).sort_values("uv_flux_factor")

    out_uv_summary = os.path.join(args.out_dir, f"table_uv_summary_{args.label}.csv")
    uv_summary.to_csv(out_uv_summary, index=False)
    print(f"Wrote {out_uv_summary}")

    # --- UV suppression table (UV=0 vs UV=uv_high) per (T, nH) ---
    key_cols = ["surface_temperature_k", "h_gas_density_cm3"]
    metric_for_supp = metric_release or metric_total_h2
    if metric_for_supp is not None:
        low = df[df["uv_flux_factor"] == 0.0][key_cols + [metric_for_supp]].rename(columns={metric_for_supp: "metric_uv0"})
        high = df[df["uv_flux_factor"] == float(args.uv_high)][key_cols + [metric_for_supp]].rename(
            columns={metric_for_supp: f"metric_uv{int(args.uv_high)}"}
        )
        merged = pd.merge(low, high, on=key_cols, how="inner")
        hi_col = f"metric_uv{int(args.uv_high)}"
        merged["ratio_uvhigh_over_uv0"] = merged.apply(lambda r: _safe_div(r[hi_col], r["metric_uv0"]), axis=1)
        merged["percent_drop"] = (1.0 - merged["ratio_uvhigh_over_uv0"]) * 100.0
        merged = merged.sort_values(key_cols)

        out_uv_supp = os.path.join(args.out_dir, f"table_uv_suppression_{args.label}.csv")
        merged.to_csv(out_uv_supp, index=False)
        print(f"Wrote {out_uv_supp}")

    # --- Top 25 conditions (by ISM-relevant metric) ---
    sort_col = metric_release or metric_total_h2
    top = df.sort_values(sort_col, ascending=False).head(25).copy()
    keep = [
        "surface_temperature_k",
        "h_gas_density_cm3",
        "uv_flux_factor",
        metric_release,
        metric_total_h2,
        metric_eps,
        metric_lh,
        metric_er,
        metric_uv,
    ]
    keep = [c for c in keep if c and c in top.columns]
    top = top[keep]
    out_top = os.path.join(args.out_dir, f"table_top25_{args.label}.csv")
    top.to_csv(out_top, index=False)
    print(f"Wrote {out_top}")


if __name__ == "__main__":
    main()

