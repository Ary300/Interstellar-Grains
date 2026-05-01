#!/usr/bin/env python3
"""
Create a small (paper-ready) timescale table for a few representative ISM conditions.

This is a convenience wrapper around the same assumptions used in compute_ism_timescales.py,
but it selects a few scenarios and writes a compact CSV suitable for Overleaf.
"""

from __future__ import annotations

import argparse
import math
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


G_CGS = 6.67430e-8  # cm^3 g^-1 s^-2
M_H_G = 1.6735575e-24  # g
SEC_PER_MYR = 3600.0 * 24.0 * 365.25 * 1.0e6


def _pick_col(df: pd.DataFrame, base: str) -> str:
    if base in df.columns:
        return base
    m = f"{base}_mean"
    if m in df.columns:
        return m
    raise KeyError(f"Missing column: {base} or {base}_mean")


def _select_row(df: pd.DataFrame, *, T: float, nH: float, uv: float) -> pd.Series:
    sub = df[
        np.isclose(df["surface_temperature_k"], float(T))
        & np.isclose(df["h_gas_density_cm3"], float(nH))
        & np.isclose(df["uv_flux_factor"], float(uv))
    ]
    if sub.empty:
        raise SystemExit(f"No row matched T={T}, nH={nH}, uv={uv}.")
    if len(sub) != 1:
        sub = sub.iloc[[0]]
    return sub.iloc[0]


def main() -> None:
    p = argparse.ArgumentParser(description="Make a compact timescale table for representative ISM conditions.")
    p.add_argument("--input", default="results/jhub_full_merged.csv", help="Merged aggregated sweep CSV")
    p.add_argument("--out", default="results/tables/table_timescales_representative.csv", help="Output CSV path")
    p.add_argument("--sigma-H", dest="sigma_H", type=float, default=1e-21, help="Total grain cross-section per H [cm^2/H]")
    p.add_argument("--mu", type=float, default=1.4, help="Mean mass per H nucleus in units of m_H. Default: 1.4")
    args = p.parse_args()

    df = pd.read_csv(args.input)
    for c in ["surface_temperature_k", "h_gas_density_cm3", "uv_flux_factor"]:
        if c not in df.columns:
            raise SystemExit(f"Missing required column: {c}")

    for c in ["surface_temperature_k", "h_gas_density_cm3", "uv_flux_factor"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    rate_col = _pick_col(df, "h2_release_rate_cm2_s")
    df[rate_col] = pd.to_numeric(df[rate_col], errors="coerce")
    df = df.dropna(subset=["surface_temperature_k", "h_gas_density_cm3", "uv_flux_factor", rate_col])

    scenarios: List[Tuple[str, Dict[str, float]]] = [
        ("high_z_core_like", {"T": 60.0, "nH": 1000.0, "uv": 0.0}),
        ("warm_diffuse_like", {"T": 100.0, "nH": 100.0, "uv": 0.0}),
        ("dense_hot_like", {"T": 40.0, "nH": 10000.0, "uv": 0.0}),
    ]

    out_rows = []
    for name, sc in scenarios:
        row = _select_row(df, T=sc["T"], nH=sc["nH"], uv=sc["uv"])
        R_area = float(row[rate_col])
        nH = float(row["h_gas_density_cm3"])

        # Effective k_eff for mostly-atomic gas: R_vol = k_eff nH^2 with k_eff = 4 σ_H R_area / nH
        k_eff = 4.0 * float(args.sigma_H) * float(R_area) / float(nH)

        # Atomic->molecular timescale (e-folding, atomic gas): t_H2 ≈ 1 / (2 k_eff nH) = 1/(8 σ_H R_area)
        t_h2_s = 1.0 / (8.0 * float(args.sigma_H) * float(R_area))
        t_h2_myr = t_h2_s / float(SEC_PER_MYR)

        rho = float(args.mu) * float(M_H_G) * float(nH)
        t_ff_s = math.sqrt((3.0 * math.pi) / (32.0 * float(G_CGS) * float(rho)))
        t_ff_myr = t_ff_s / float(SEC_PER_MYR)

        out_rows.append(
            {
                "scenario": name,
                "surface_temperature_k": float(row["surface_temperature_k"]),
                "h_gas_density_cm3": nH,
                "uv_flux_factor": float(row["uv_flux_factor"]),
                "sigma_H_cm2_per_H": float(args.sigma_H),
                "h2_release_rate_cm2_s": float(R_area),
                "k_eff_cm3_s": float(k_eff),
                "t_H2_Myr": float(t_h2_myr),
                "t_ff_Myr": float(t_ff_myr),
                "t_H2_over_t_ff": float(t_h2_s / t_ff_s),
            }
        )

    out = pd.DataFrame(out_rows)
    out.to_csv(args.out, index=False)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()

