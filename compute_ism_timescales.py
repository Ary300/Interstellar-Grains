#!/usr/bin/env python3
"""
Compute simple H2 formation and free-fall timescales from an ISM sweep CSV.

This is intended for the "gas_kinetic" arrival-mode sweeps where the primary
observable is the *per-grain-surface-area* formation/release rate:
  h2_release_rate_cm2_s(_mean)

We convert to a volumetric formation timescale using an assumed total grain
cross-section per H nucleus, σ_H [cm^2/H]. For spherical grains:
  surface_area_per_volume = 4 * (n_d * π r^2) = 4 * σ_H * n_H

Thus:
  R_vol(H2) = R_area * (surface_area_per_volume)
            = 4 * σ_H * n_H * R_area

and the (e-folding) atomic→molecular conversion timescale is:
  t_H2 ≈ 1 / (8 * σ_H * R_area)

Free-fall time is computed from:
  t_ff = sqrt(3π / (32 G ρ)), with ρ = μ m_H n_H.
"""

from __future__ import annotations

import argparse
import math
import os
from typing import List, Optional

import numpy as np
import pandas as pd


G_CGS = 6.67430e-8  # cm^3 g^-1 s^-2
M_H_G = 1.6735575e-24  # g
CANONICAL_H2_FORMATION_RATE_CM3_S = 3.0e-17


def _pick_col(df: pd.DataFrame, base: str) -> str:
    if base in df.columns:
        return base
    m = f"{base}_mean"
    if m in df.columns:
        return m
    raise KeyError(f"Missing column: {base} or {base}_mean")


def _float_list(xs: Optional[List[str]]) -> List[float]:
    if not xs:
        return []
    out: List[float] = []
    for x in xs:
        out.append(float(x))
    return out


def main() -> None:
    p = argparse.ArgumentParser(description="Compute simple ISM H2 formation timescales from a merged sweep CSV.")
    p.add_argument("--input", default="results/jhub_full_merged.csv", help="Merged aggregated sweep CSV")
    p.add_argument("--out", default="results/tables/ism_timescales.csv", help="Output CSV with timescales per row")
    p.add_argument("--T", type=float, default=None, help="Optional filter: surface_temperature_k")
    p.add_argument("--nH", type=float, default=None, help="Optional filter: h_gas_density_cm3")
    p.add_argument("--uv", type=float, default=None, help="Optional filter: uv_flux_factor")
    p.add_argument(
        "--sigma-H",
        dest="sigma_H",
        action="append",
        default=None,
        help="Total grain geometric cross-section per H nucleus [cm^2/H]. Repeatable. Default: 1e-21",
    )
    p.add_argument("--mu", type=float, default=1.4, help="Mean mass per H nucleus in units of m_H. Default: 1.4")
    p.add_argument(
        "--benchmark-kf",
        type=float,
        default=CANONICAL_H2_FORMATION_RATE_CM3_S,
        help="Canonical observational H2 formation rate coefficient [cm^3 s^-1]. Default: 3e-17",
    )
    args = p.parse_args()

    sigma_H_list = _float_list(args.sigma_H) if args.sigma_H else [1e-21]

    df = pd.read_csv(args.input)
    for c in ["surface_temperature_k", "h_gas_density_cm3", "uv_flux_factor"]:
        if c not in df.columns:
            raise SystemExit(f"Missing required column: {c}")

    rate_col = _pick_col(df, "h2_release_rate_cm2_s")

    # Optional filtering for a single "scenario" row.
    if args.T is not None:
        df = df[np.isclose(pd.to_numeric(df["surface_temperature_k"], errors="coerce"), float(args.T))]
    if args.nH is not None:
        df = df[np.isclose(pd.to_numeric(df["h_gas_density_cm3"], errors="coerce"), float(args.nH))]
    if args.uv is not None:
        df = df[np.isclose(pd.to_numeric(df["uv_flux_factor"], errors="coerce"), float(args.uv))]

    df = df.copy()
    df["surface_temperature_k"] = pd.to_numeric(df["surface_temperature_k"], errors="coerce")
    df["h_gas_density_cm3"] = pd.to_numeric(df["h_gas_density_cm3"], errors="coerce")
    df["uv_flux_factor"] = pd.to_numeric(df["uv_flux_factor"], errors="coerce")
    df[rate_col] = pd.to_numeric(df[rate_col], errors="coerce")
    df = df.dropna(subset=["surface_temperature_k", "h_gas_density_cm3", "uv_flux_factor", rate_col])

    if df.empty:
        raise SystemExit("No rows matched the requested filters.")

    out_rows = []
    for sigma_H in sigma_H_list:
        if sigma_H <= 0:
            continue
        # t_H2 ≈ 1 / (8 σ_H R_area)
        t_h2_s = 1.0 / (8.0 * float(sigma_H) * df[rate_col].astype(float))
        t_h2_myr = t_h2_s / (3600.0 * 24.0 * 365.25 * 1.0e6)

        rho = float(args.mu) * float(M_H_G) * df["h_gas_density_cm3"].astype(float)
        t_ff_s = np.sqrt((3.0 * math.pi) / (32.0 * float(G_CGS) * rho))
        t_ff_myr = t_ff_s / (3600.0 * 24.0 * 365.25 * 1.0e6)

        block = df[["surface_temperature_k", "h_gas_density_cm3", "uv_flux_factor"]].copy()
        block["sigma_H_cm2_per_H"] = float(sigma_H)
        block["h2_release_rate_cm2_s"] = df[rate_col].astype(float)
        # Effective dust formation-rate coefficient for mostly-atomic gas:
        #   R_vol(H2) = 4 σ_H n_H R_area
        #   k_eff ≡ R_vol / n_H^2 = 4 σ_H R_area / n_H
        block["k_eff_cm3_s"] = (4.0 * float(sigma_H) * block["h2_release_rate_cm2_s"] / block["h_gas_density_cm3"]).astype(float)
        block["benchmark_kf_cm3_s"] = float(args.benchmark_kf)
        block["k_eff_over_benchmark"] = (block["k_eff_cm3_s"] / float(args.benchmark_kf)).astype(float)
        block["t_H2_Myr"] = t_h2_myr.astype(float)
        block["t_ff_Myr"] = t_ff_myr.astype(float)
        block["t_H2_over_t_ff"] = (block["t_H2_Myr"] / block["t_ff_Myr"]).astype(float)
        out_rows.append(block)

    out = pd.concat(out_rows, ignore_index=True)
    out_dir = os.path.dirname(args.out) or "."
    os.makedirs(out_dir, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"Wrote {args.out}")

    # Print a compact scenario summary if the user filtered down to a single condition.
    if out[["surface_temperature_k", "h_gas_density_cm3", "uv_flux_factor"]].drop_duplicates().shape[0] == 1:
        key = out.iloc[0][["surface_temperature_k", "h_gas_density_cm3", "uv_flux_factor"]].to_dict()
        print(f"Scenario: T={key['surface_temperature_k']} K, nH={key['h_gas_density_cm3']} cm^-3, uv={key['uv_flux_factor']}")
        for sigma_H in sigma_H_list:
            s = out[out["sigma_H_cm2_per_H"] == float(sigma_H)].iloc[0]
            print(
                f"  sigma_H={sigma_H:.1e}: k_eff={s['k_eff_cm3_s']:.2e} cm^3 s^-1 "
                f"({s['k_eff_over_benchmark']:.2f} x benchmark), "
                f"t_H2={s['t_H2_Myr']:.2f} Myr, t_ff={s['t_ff_Myr']:.2f} Myr, "
                f"t_H2/t_ff={s['t_H2_over_t_ff']:.2f}"
            )


if __name__ == "__main__":
    main()
