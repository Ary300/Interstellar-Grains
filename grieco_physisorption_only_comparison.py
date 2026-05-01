import argparse
import csv
import os
from dataclasses import dataclass
from typing import Any, Dict, List, Tuple

import numpy as np
import yaml

import matplotlib.pyplot as plt

from kmc_simulation import KineticMonteCarlo


@dataclass(frozen=True)
class Stat:
    mean: float
    std: float
    ci95: float


def _stat(x: List[float]) -> Stat:
    arr = np.array(x, dtype=float)
    n = int(arr.size)
    if n <= 0:
        return Stat(mean=0.0, std=0.0, ci95=0.0)
    mean = float(np.mean(arr))
    std = float(np.std(arr, ddof=1)) if n > 1 else 0.0
    ci95 = float(1.96 * std / float(np.sqrt(n))) if n > 1 else 0.0
    return Stat(mean=mean, std=std, ci95=ci95)


def _load_yaml(path: str) -> Dict[str, Any]:
    with open(path, "r") as f:
        obj = yaml.safe_load(f) or {}
    if not isinstance(obj, dict):
        raise ValueError(f"Expected mapping in {path}")
    return obj


def _reset_measurement_counters(kmc: KineticMonteCarlo) -> None:
    kmc.total_impinging_h_atoms = 0
    kmc.total_impinging_h2_molecules = 0
    kmc.total_adsorbed_h_atoms = 0
    kmc.total_desorbed_h_atoms = 0
    kmc.h2_molecules_formed = 0
    kmc.h2_molecules_desorbed = 0  # prompt formed-origin
    kmc.h2_molecules_desorbed_LH = 0
    kmc.h2_molecules_desorbed_ER = 0
    kmc.h2_molecules_desorbed_UV = 0
    kmc.h2_molecules_desorbed_beam = 0
    kmc.h2_molecules_released_formed = 0
    kmc.h2_molecules_released_beam = 0
    kmc.h2_molecules_formed_LH = 0
    kmc.h2_molecules_formed_ER = 0
    kmc.h2_molecules_formed_UV = 0


def _run_isothermal(
    base_params: Dict[str, Any],
    temperature_k: float,
    burnin_arrivals: int,
    measure_arrivals: int,
    seed: int,
    max_steps: int | None,
) -> Tuple[float, float]:
    """
    Returns (epsilon_prompt, epsilon_released_total) for a single realization.
    """
    p = dict(base_params)
    p["surface_temperature_k"] = float(temperature_k)
    p["rng_seed"] = int(seed)

    # Burn-in to a quasi-steady state (by arrivals).
    if burnin_arrivals > 0:
        p["max_arrivals"] = int(burnin_arrivals)
    else:
        p.pop("max_arrivals", None)

    kmc = KineticMonteCarlo(p)
    if burnin_arrivals > 0:
        kmc.run_gillespie(max_time=1e30, max_steps=max_steps)

    # Measurement window.
    _reset_measurement_counters(kmc)
    kmc.simulation_parameters["max_arrivals"] = int(measure_arrivals)
    kmc.run_gillespie(max_time=1e30, max_steps=max_steps)

    imp = int(getattr(kmc, "total_impinging_h_atoms", 0))
    if imp <= 0:
        return 0.0, 0.0

    prompt = float(getattr(kmc, "h2_molecules_desorbed", 0))
    released_total = float(getattr(kmc, "h2_molecules_desorbed", 0)) + float(getattr(kmc, "h2_molecules_released_formed", 0))
    eps_prompt = float(2.0 * prompt / float(imp))
    eps_released = float(2.0 * released_total / float(imp))
    return eps_prompt, eps_released


def _write_csv(path: str, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        return
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def _read_csv_series(path: str, x_col: str, y_col: str, ci_col: str | None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    xs: List[float] = []
    ys: List[float] = []
    cis: List[float] = []
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            try:
                xs.append(float(row[x_col]))
                ys.append(float(row[y_col]))
                if ci_col and ci_col in row and str(row[ci_col]).strip() != "":
                    cis.append(abs(float(row[ci_col])))
                else:
                    cis.append(0.0)
            except Exception:
                continue
    order = np.argsort(np.array(xs, dtype=float))
    return (
        np.array(xs, dtype=float)[order],
        np.array(ys, dtype=float)[order],
        np.array(cis, dtype=float)[order],
    )


def main() -> None:
    p = argparse.ArgumentParser(description="Physisorption-only comparison run (chemisorption OFF) vs paperfit baseline.")
    p.add_argument("--paperfit-iso", default="results/grieco_validation_paper_iso_paperfit.csv")
    p.add_argument("--paperfit-ded", default="results/grieco_ded_paper_ded_paperfit.csv")
    p.add_argument("--base-config", default="config_grieco_paper_ded_paperfit.yaml")

    p.add_argument("--temps", nargs="+", type=float, default=[5, 10, 15, 20, 25, 30, 35])
    p.add_argument("--replicates", type=int, default=100)
    p.add_argument("--burnin-arrivals", type=int, default=2000)
    p.add_argument("--measure-arrivals", type=int, default=5000)
    p.add_argument("--max-steps", type=int, default=800000)

    p.add_argument("--arrival-rate-per-site", type=float, default=0.01)
    p.add_argument("--tau", type=float, default=1.0, help="Set beam_dissociation_fraction (use 1.0 for all-atom arrivals).")

    p.add_argument("--out-csv", default="results/grieco_physisorption_only_iso.csv")
    p.add_argument("--out-plot", default="results/plots/grieco_physisorption_only_vs_paperfit.png")
    args = p.parse_args()

    base = _load_yaml(str(args.base_config))

    # Force an isothermal, arrival-mode LH-only regime.
    base.update(
        {
            "arrival_rate_per_site_s": float(args.arrival_rate_per_site),
            "beam_dissociation_fraction": float(args.tau),
            "enable_LH": True,
            "enable_diffusion": True,
            "diffusion_mode": "rate_only",
            "lh_formation_mode": "diffusion_limited",
            # Physisorption-only: no chem reservoir.
            "chemisorption_fraction": 0.0,
            # Also remove deep "defect trap" sites; Satonkin-style baselines are physisorption-only.
            "surface_defect_fraction": 0.0,
            # Keep UV off for this comparison.
            "uv_flux_factor": 0.0,
            "uv_pulse_enabled": False,
            # Speed: reuse the same grain topology/energetics across ensemble runs.
            "enable_grain_cache": True,
            "grain_cache_dir": "grain_cache",
            "grain_cache_include_rng_seed": False,
        }
    )

    # Ensure blocking is off (Satonkin-style physisorption-only baseline).
    base.update(
        {
            "enable_h2_blocking": False,
            "sticking_blocking_strength": 0.0,
            "er_blocking_strength": 0.0,
        }
    )

    rows: List[Dict[str, Any]] = []
    for T in [float(t) for t in args.temps]:
        eps_prompt_runs: List[float] = []
        eps_rel_runs: List[float] = []
        for i in range(int(args.replicates)):
            ep, er = _run_isothermal(
                base_params=base,
                temperature_k=float(T),
                burnin_arrivals=int(args.burnin_arrivals),
                measure_arrivals=int(args.measure_arrivals),
                seed=4000 + i,
                max_steps=int(args.max_steps) if args.max_steps else None,
            )
            eps_prompt_runs.append(float(ep))
            eps_rel_runs.append(float(er))
        sp = _stat(eps_prompt_runs)
        sr = _stat(eps_rel_runs)
        rows.append(
            {
                "surface_temperature_k": float(T),
                "epsilon_prompt_mean": float(sp.mean),
                "epsilon_prompt_ci95": float(sp.ci95),
                "epsilon_released_total_mean": float(sr.mean),
                "epsilon_released_total_ci95": float(sr.ci95),
                "replicates": int(args.replicates),
                "burnin_arrivals": int(args.burnin_arrivals),
                "measure_arrivals": int(args.measure_arrivals),
                "arrival_rate_per_site_s": float(args.arrival_rate_per_site),
                "beam_dissociation_fraction": float(args.tau),
                "chemisorption_fraction": 0.0,
            }
        )

    rows = sorted(rows, key=lambda r: float(r["surface_temperature_k"]))
    _write_csv(str(args.out_csv), rows)
    print(f"Wrote {args.out_csv}")

    # Build the comparison plot: paperfit curve vs physisorption-only.
    paperfit_ded_x, paperfit_ded_y, paperfit_ded_ci = _read_csv_series(
        str(args.paperfit_ded),
        x_col="temperature_k",
        y_col="epsilon_released_total_mean",
        ci_col="epsilon_released_total_ci95",
    )
    paperfit_iso_x, paperfit_iso_y, paperfit_iso_ci = _read_csv_series(
        str(args.paperfit_iso),
        x_col="surface_temperature_k",
        y_col="epsilon_mean",
        ci_col="epsilon_ci95",
    )
    phy_x = np.array([float(r["surface_temperature_k"]) for r in rows], dtype=float)
    phy_y = np.array([float(r["epsilon_released_total_mean"]) for r in rows], dtype=float)
    phy_ci = np.array([float(r["epsilon_released_total_ci95"]) for r in rows], dtype=float)

    plt.figure(figsize=(8.5, 6.0))
    if paperfit_ded_x.size:
        plt.errorbar(paperfit_ded_x, paperfit_ded_y, yerr=paperfit_ded_ci, fmt="o", ms=5, capsize=2, label="paperfit TPDED (released_total)")
    if paperfit_iso_x.size:
        plt.errorbar(paperfit_iso_x, paperfit_iso_y, yerr=paperfit_iso_ci, fmt="^", ms=7, capsize=2, label="paperfit isothermal (prompt)")

    plt.errorbar(phy_x, phy_y, yerr=phy_ci, fmt="s--", ms=6, capsize=2, label="chemisorption OFF (isothermal, LH-only)")

    plt.xlabel("T (K)")
    plt.ylabel("Recombination efficiency ε (fraction)")
    plt.title("Physisorption-only collapse vs chemisorption-enabled paperfit curve")
    plt.grid(True, alpha=0.3)
    plt.xlim(0, 260)
    plt.ylim(0, 0.55)
    plt.legend(loc="upper right", frameon=False)
    os.makedirs(os.path.dirname(str(args.out_plot)) or ".", exist_ok=True)
    plt.tight_layout()
    plt.savefig(str(args.out_plot), dpi=200)
    plt.close()
    print(f"Wrote {args.out_plot}")


if __name__ == "__main__":
    main()
