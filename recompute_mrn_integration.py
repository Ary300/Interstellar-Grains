#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


def _logspace_edges(sizes: np.ndarray) -> np.ndarray:
    sizes = np.asarray(sizes, dtype=float)
    if sizes.ndim != 1 or sizes.size < 2:
        raise ValueError("Need at least two size bins to compute edges.")
    ratios = sizes[1:] / sizes[:-1]
    if not np.allclose(ratios, ratios[0], rtol=1e-4, atol=0):
        raise ValueError("Size bins are not log-uniform; cannot infer edges reliably.")
    r = float(ratios[0])
    edges = np.sqrt(sizes[:-1] * sizes[1:])
    first = sizes[0] / math.sqrt(r)
    last = sizes[-1] * math.sqrt(r)
    return np.concatenate([[first], edges, [last]])


def _area_weighted_mrn_weights(sizes: np.ndarray) -> np.ndarray:
    edges = _logspace_edges(sizes)
    # MRN number density ~ a^-3.5, surface area ~ a^2 => weight ~ a^-1.5
    # Integrate a^-1.5 over each bin [a_lo, a_hi].
    weights = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        w = 2.0 * (lo ** -0.5 - hi ** -0.5)
        weights.append(w)
    weights = np.asarray(weights, dtype=float)
    weights = weights / np.sum(weights)
    return weights


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
    p = argparse.ArgumentParser(
        description="Recompute MRN-integrated results using cross-section weighted MRN bins."
    )
    p.add_argument("--size-summary", default="results/tables/astro_mrn_size_summary.csv")
    p.add_argument("--output", default="results/astro_mrn_integration.csv")
    p.add_argument("--out-weights", default="results/tables/astro_mrn_weights.csv")
    args = p.parse_args()

    df = pd.read_csv(args.size_summary)
    required = {
        "surface_temperature_k",
        "h_gas_density_cm3",
        "uv_flux_factor",
        "grain_radius_um_mrn",
        "epsilon_mean",
        "epsilon_std",
        "h2_release_rate_cm2_s_mean",
        "h2_release_rate_cm2_s_std",
        "n_runs",
    }
    missing = sorted(required - set(df.columns))
    if missing:
        raise SystemExit(f"Missing required columns: {missing}")

    for col in required:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=list(required))

    sizes = np.sort(df["grain_radius_um_mrn"].unique())
    weights = _area_weighted_mrn_weights(sizes)
    weights_df = pd.DataFrame({"grain_radius_um_mrn": sizes, "mrn_weight_area": weights})
    weights_df.to_csv(args.out_weights, index=False)

    out_rows = []
    for (temp, nh, uv), grp in df.groupby(
        ["surface_temperature_k", "h_gas_density_cm3", "uv_flux_factor"], dropna=False
    ):
        grp = grp.sort_values("grain_radius_um_mrn")
        if len(grp) != len(sizes):
            raise SystemExit(
                f"Missing size bins for T={temp}, nH={nh}, UV={uv}: expected {len(sizes)}, got {len(grp)}"
            )
        eps_mean, eps_std, eps_ci95 = _weighted_combine(
            grp["epsilon_mean"].to_numpy(),
            grp["epsilon_std"].to_numpy(),
            weights,
        )
        rate_mean, rate_std, rate_ci95 = _weighted_combine(
            grp["h2_release_rate_cm2_s_mean"].to_numpy(),
            grp["h2_release_rate_cm2_s_std"].to_numpy(),
            weights,
        )
        out_rows.append(
            {
                "surface_temperature_k": float(temp),
                "h_gas_density_cm3": float(nh),
                "uv_flux_factor": float(uv),
                "mrn_aggregated": True,
                "mrn_bins": int(len(sizes)),
                "rows_used": int(grp["n_runs"].sum()),
                "epsilon_mean": eps_mean,
                "epsilon_std": eps_std,
                "epsilon_ci95": eps_ci95,
                "h2_release_rate_cm2_s_mean": rate_mean,
                "h2_release_rate_cm2_s_std": rate_std,
                "h2_release_rate_cm2_s_ci95": rate_ci95,
            }
        )

    out = pd.DataFrame(out_rows).sort_values(
        ["surface_temperature_k", "h_gas_density_cm3", "uv_flux_factor"]
    )
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    print(f"Wrote {out_path} ({len(out)} rows)")
    print(f"Wrote {args.out_weights} ({len(weights_df)} rows)")


if __name__ == "__main__":
    main()
