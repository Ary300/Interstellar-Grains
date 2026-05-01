#!/usr/bin/env python3
import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


def _weighted_combine(means: np.ndarray, stds: np.ndarray, weights: np.ndarray) -> tuple[float, float, float]:
    weights = np.asarray(weights, dtype=float)
    weights = weights / np.sum(weights)
    means = np.asarray(means, dtype=float)
    stds = np.asarray(stds, dtype=float)
    wmean = float(np.sum(weights * means))
    between_var = float(np.sum(weights * (means - wmean) ** 2))
    within_var = float(np.sum((weights ** 2) * (stds ** 2)))
    wstd = float(np.sqrt(max(0.0, between_var + within_var)))
    ci95 = float(1.96 * wstd)
    return wmean, wstd, ci95


def main() -> None:
    p = argparse.ArgumentParser(description="Aggregate chunked MRN raw outputs into a paper-ready MRN CSV.")
    p.add_argument("--input-glob", required=True, help="Example: 'results/anvil_mrn_r5/raw_chunk_*.csv'")
    p.add_argument("--output", required=True, help="Merged MRN-integrated CSV output")
    p.add_argument("--size-summary", default=None, help="Optional per-size summary CSV")
    args = p.parse_args()

    paths = sorted(Path(".").glob(args.input_glob))
    if not paths:
        raise SystemExit(f"No files matched: {args.input_glob}")

    frames = [pd.read_csv(path) for path in paths]
    raw = pd.concat(frames, ignore_index=True)

    required = {
        "surface_temperature_k",
        "h_gas_density_cm3",
        "uv_flux_factor",
        "grain_radius_um_mrn",
        "_mrn_weight",
        "epsilon",
        "h2_release_rate_cm2_s",
    }
    missing = sorted(required - set(raw.columns))
    if missing:
        raise SystemExit(f"Missing required columns: {missing}")

    for col in required:
        raw[col] = pd.to_numeric(raw[col], errors="coerce")
    raw = raw.dropna(subset=list(required))

    by_size = (
        raw.groupby(
            [
                "surface_temperature_k",
                "h_gas_density_cm3",
                "uv_flux_factor",
                "grain_radius_um_mrn",
                "_mrn_weight",
            ],
            dropna=False,
        )
        .agg(
            n_runs=("epsilon", "count"),
            epsilon_mean=("epsilon", "mean"),
            epsilon_std=("epsilon", "std"),
            h2_release_rate_cm2_s_mean=("h2_release_rate_cm2_s", "mean"),
            h2_release_rate_cm2_s_std=("h2_release_rate_cm2_s", "std"),
        )
        .reset_index()
        .sort_values(["surface_temperature_k", "grain_radius_um_mrn"])
    )
    by_size["epsilon_std"] = by_size["epsilon_std"].fillna(0.0)
    by_size["h2_release_rate_cm2_s_std"] = by_size["h2_release_rate_cm2_s_std"].fillna(0.0)

    merged_rows = []
    for (temp, nh, uv), grp in by_size.groupby(
        ["surface_temperature_k", "h_gas_density_cm3", "uv_flux_factor"], dropna=False
    ):
        eps_mean, eps_std, eps_ci95 = _weighted_combine(
            grp["epsilon_mean"].to_numpy(),
            grp["epsilon_std"].to_numpy(),
            grp["_mrn_weight"].to_numpy(),
        )
        rate_mean, rate_std, rate_ci95 = _weighted_combine(
            grp["h2_release_rate_cm2_s_mean"].to_numpy(),
            grp["h2_release_rate_cm2_s_std"].to_numpy(),
            grp["_mrn_weight"].to_numpy(),
        )
        merged_rows.append(
            {
                "surface_temperature_k": float(temp),
                "h_gas_density_cm3": float(nh),
                "uv_flux_factor": float(uv),
                "mrn_aggregated": True,
                "mrn_bins": int(len(grp)),
                "rows_used": int(grp["n_runs"].sum()),
                "epsilon_mean": eps_mean,
                "epsilon_std": eps_std,
                "epsilon_ci95": eps_ci95,
                "h2_release_rate_cm2_s_mean": rate_mean,
                "h2_release_rate_cm2_s_std": rate_std,
                "h2_release_rate_cm2_s_ci95": rate_ci95,
            }
        )

    out = pd.DataFrame(merged_rows).sort_values(["surface_temperature_k", "h_gas_density_cm3", "uv_flux_factor"])
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    print(f"Wrote {out_path} ({len(out)} rows)")

    if args.size_summary:
        size_path = Path(args.size_summary)
        size_path.parent.mkdir(parents=True, exist_ok=True)
        by_size.to_csv(size_path, index=False)
        print(f"Wrote {size_path} ({len(by_size)} rows)")


if __name__ == "__main__":
    main()
