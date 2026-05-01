#!/usr/bin/env python3
"""
Compare two Grieco-validation model variants side-by-side.

Primary use:
- baseline paperfit validation vs "prediction-model" validation
- quantify how much the fit moves when porosity / LH transport assumptions are
  changed to match the ISM production model more closely
"""

from __future__ import annotations

import argparse
import os
from typing import Dict, Iterable, List

import numpy as np
import pandas as pd

from compare_grieco_fig2_digitized import compute as compute_digitized_fit


def _read_csv(path: str) -> pd.DataFrame:
    return pd.read_csv(path)


def _pick_ded_cols(df: pd.DataFrame) -> tuple[str, str]:
    if "epsilon_released_total_mean" in df.columns:
        ci = "epsilon_released_total_ci95" if "epsilon_released_total_ci95" in df.columns else "epsilon_prompt_ci95"
        return "epsilon_released_total_mean", ci
    if "epsilon_mean" in df.columns:
        return "epsilon_mean", "epsilon_ci95"
    raise KeyError("Could not find a DED epsilon column")


def _mean_between(df: pd.DataFrame, xcol: str, ycol: str, lo: float, hi: float) -> float:
    sub = df[(pd.to_numeric(df[xcol], errors="coerce") >= float(lo)) & (pd.to_numeric(df[xcol], errors="coerce") <= float(hi))]
    if sub.empty:
        return float("nan")
    return float(pd.to_numeric(sub[ycol], errors="coerce").mean())


def _stats_rows_to_dict(rows: Iterable[Dict[str, object]]) -> Dict[str, Dict[str, float]]:
    out: Dict[str, Dict[str, float]] = {}
    for row in rows:
        out[str(row["subset"])] = {k: float(v) if isinstance(v, (int, float)) else v for k, v in row.items()}
    return out


def _summarize_variant(label: str, iso_csv: str, ded_csv: str, digitized_csv: str, n_params: int) -> Dict[str, object]:
    iso = _read_csv(iso_csv)
    ded = _read_csv(ded_csv)
    ded_y, ded_ci = _pick_ded_cols(ded)

    residual_rows, stats_rows = compute_digitized_fit(
        digitized_csv=digitized_csv,
        iso_sim_csv=iso_csv,
        ded_sim_csv=ded_csv,
        ded_y=ded_y,
        ded_ci=ded_ci,
        include_model_sigma=False,
        n_params=int(n_params),
    )
    stats = _stats_rows_to_dict(stats_rows)

    row: Dict[str, object] = {
        "variant": label,
        "iso_csv": iso_csv,
        "ded_csv": ded_csv,
        "highT_plateau_mean_100_250K": _mean_between(iso, "surface_temperature_k", "epsilon_mean", 100.0, 250.0),
        "ded_eps10_mean": _mean_between(ded, "temperature_k", ded_y, 9.5, 11.5),
        "ded_eps20_mean": _mean_between(ded, "temperature_k", ded_y, 19.5, 21.5),
        "ded_eps30_80_mean": _mean_between(ded, "temperature_k", ded_y, 30.0, 80.0),
    }
    row["ded_ratio10_over20"] = float(row["ded_eps10_mean"]) / float(row["ded_eps20_mean"]) if float(row["ded_eps20_mean"]) > 0 else float("nan")
    for subset in ["all", "tpded", "isothermal"]:
        row[f"chi2_reduced_{subset}"] = stats.get(subset, {}).get("chi2_reduced", float("nan"))
        row[f"rmse_{subset}"] = stats.get(subset, {}).get("rmse", float("nan"))
    return row


def main() -> None:
    p = argparse.ArgumentParser(description="Compare baseline vs prediction-model Grieco validation outputs.")
    p.add_argument("--digitized-csv", default="grieco_fig2_digitized.csv", help="Digitized Grieco Fig. 2 CSV")
    p.add_argument("--baseline-label", default="paperfit_baseline")
    p.add_argument("--baseline-iso", default="results/grieco_validation_paper_iso_paperfit.csv")
    p.add_argument("--baseline-ded", default="results/grieco_ded_paper_ded_paperfit.csv")
    p.add_argument("--candidate-label", default="prediction_model")
    p.add_argument("--candidate-iso", default="results/grieco_validation_paper_iso_predmodel.csv")
    p.add_argument("--candidate-ded", default="results/grieco_ded_paper_ded_predmodel.csv")
    p.add_argument("--n-params", type=int, default=0, help="Number of fitted parameters to subtract in reduced chi^2")
    p.add_argument("--out-csv", default="results/tables/grieco_model_variant_comparison.csv")
    args = p.parse_args()

    rows = [
        _summarize_variant(args.baseline_label, args.baseline_iso, args.baseline_ded, args.digitized_csv, args.n_params),
        _summarize_variant(args.candidate_label, args.candidate_iso, args.candidate_ded, args.digitized_csv, args.n_params),
    ]

    out_dir = os.path.dirname(args.out_csv) or "."
    os.makedirs(out_dir, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.out_csv, index=False)
    print(f"Wrote {args.out_csv}")


if __name__ == "__main__":
    main()

