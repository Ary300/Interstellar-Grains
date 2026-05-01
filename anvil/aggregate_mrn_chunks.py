#!/usr/bin/env python3
import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd


def _combine_chunks(means: np.ndarray, stds: np.ndarray, counts: np.ndarray) -> tuple[float, float]:
    means = np.asarray(means, dtype=float)
    stds = np.asarray(stds, dtype=float)
    counts = np.asarray(counts, dtype=float)
    total = int(np.sum(counts))
    if total <= 0:
        return 0.0, 0.0
    mean = float(np.sum(counts * means) / total)
    if total <= 1:
        return mean, 0.0
    pooled_ss = float(np.sum((counts - 1.0) * (stds ** 2) + counts * ((means - mean) ** 2)))
    std = float(np.sqrt(max(0.0, pooled_ss / (total - 1.0))))
    return mean, std


def _weighted_combine(means: np.ndarray, stds: np.ndarray, weights: np.ndarray) -> tuple[float, float, float]:
    weights = np.asarray(weights, dtype=float)
    weights = weights / np.sum(weights)
    means = np.asarray(means, dtype=float)
    stds = np.asarray(stds, dtype=float)
    wmean = float(np.sum(weights * means))
    between_var = float(np.sum(weights * (means - wmean) ** 2))
    within_var = float(np.sum((weights ** 2) * (stds ** 2)))
    wstd = float(np.sqrt(max(0.0, between_var + within_var)))
    return wmean, wstd, float(1.96 * wstd)


def main() -> None:
    p = argparse.ArgumentParser(description="Aggregate per-size MRN chunk outputs into a final MRN-integrated CSV.")
    p.add_argument("--input-glob", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--size-summary", default=None)
    args = p.parse_args()

    paths = [Path(p) for p in sorted(glob.glob(args.input_glob))]
    if not paths:
        raise SystemExit(f"No files matched: {args.input_glob}")

    frames = [pd.read_csv(path) for path in paths]
    df = pd.concat(frames, ignore_index=True)

    required = {
        "surface_temperature_k",
        "h_gas_density_cm3",
        "uv_flux_factor",
        "grain_radius_um_mrn",
        "_mrn_weight",
        "mrn_chunk_n_runs",
        "epsilon_mean",
        "epsilon_std",
        "h2_release_rate_cm2_s_mean",
        "h2_release_rate_cm2_s_std",
    }
    missing = sorted(required - set(df.columns))
    if missing:
        raise SystemExit(f"Missing required columns: {missing}")

    for col in required:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=list(required))

    size_rows = []
    size_group_cols = [
        "surface_temperature_k",
        "h_gas_density_cm3",
        "uv_flux_factor",
        "grain_radius_um_mrn",
        "_mrn_weight",
    ]
    for keys, grp in df.groupby(size_group_cols, dropna=False):
        eps_mean, eps_std = _combine_chunks(
            grp["epsilon_mean"].to_numpy(),
            grp["epsilon_std"].to_numpy(),
            grp["mrn_chunk_n_runs"].to_numpy(),
        )
        rate_mean, rate_std = _combine_chunks(
            grp["h2_release_rate_cm2_s_mean"].to_numpy(),
            grp["h2_release_rate_cm2_s_std"].to_numpy(),
            grp["mrn_chunk_n_runs"].to_numpy(),
        )
        temp, nh, uv, size_um, weight = keys
        size_rows.append(
            {
                "surface_temperature_k": float(temp),
                "h_gas_density_cm3": float(nh),
                "uv_flux_factor": float(uv),
                "grain_radius_um_mrn": float(size_um),
                "_mrn_weight": float(weight),
                "n_runs": int(np.sum(grp["mrn_chunk_n_runs"].to_numpy())),
                "epsilon_mean": eps_mean,
                "epsilon_std": eps_std,
                "epsilon_ci95": float(1.96 * eps_std / np.sqrt(max(1, int(np.sum(grp["mrn_chunk_n_runs"].to_numpy()))))),
                "h2_release_rate_cm2_s_mean": rate_mean,
                "h2_release_rate_cm2_s_std": rate_std,
                "h2_release_rate_cm2_s_ci95": float(1.96 * rate_std / np.sqrt(max(1, int(np.sum(grp["mrn_chunk_n_runs"].to_numpy()))))),
            }
        )

    by_size = pd.DataFrame(size_rows).sort_values(["surface_temperature_k", "grain_radius_um_mrn"])

    final_rows = []
    for (temp, nh, uv), grp in by_size.groupby(["surface_temperature_k", "h_gas_density_cm3", "uv_flux_factor"], dropna=False):
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
        final_rows.append(
            {
                "surface_temperature_k": float(temp),
                "h_gas_density_cm3": float(nh),
                "uv_flux_factor": float(uv),
                "mrn_aggregated": True,
                "mrn_bins": int(len(grp)),
                "rows_used": int(np.sum(grp["n_runs"].to_numpy())),
                "epsilon_mean": eps_mean,
                "epsilon_std": eps_std,
                "epsilon_ci95": eps_ci95,
                "h2_release_rate_cm2_s_mean": rate_mean,
                "h2_release_rate_cm2_s_std": rate_std,
                "h2_release_rate_cm2_s_ci95": rate_ci95,
            }
        )

    out = pd.DataFrame(final_rows).sort_values(["surface_temperature_k", "h_gas_density_cm3", "uv_flux_factor"])
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
