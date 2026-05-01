import argparse
import csv
import json
import os
from dataclasses import dataclass
from typing import Any, Dict, List

import numpy as np
import yaml

import matplotlib.pyplot as plt

from grieco_validation import _run_once
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


def _write_csv(path: str, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        return
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def _iso_eps_at_T(
    base_params: Dict[str, Any],
    temperature_k: float,
    tau: float,
    replicates: int,
    burnin_arrivals: int,
    measure_arrivals: int,
    burnin_exposure_atoms_cm2: float | None,
    measure_exposure_atoms_cm2: float | None,
    max_steps: int | None,
) -> Stat:
    p = dict(base_params)
    p["beam_dissociation_fraction"] = float(tau)

    eps_runs: List[float] = []
    for i in range(int(replicates)):
        r = _run_once(
            temperature_k=float(temperature_k),
            base_params=p,
            burnin_arrivals=int(burnin_arrivals),
            measure_arrivals=int(measure_arrivals),
            burnin_exposure_atoms_cm2=burnin_exposure_atoms_cm2,
            measure_exposure_atoms_cm2=measure_exposure_atoms_cm2,
            seed=1000 + i,
            max_steps=max_steps,
        )
        eps_runs.append(float(r.epsilon))
    return _stat(eps_runs)


def _iso_eps_prompt_at_T_arrivals(
    base_params: Dict[str, Any],
    temperature_k: float,
    tau: float,
    replicates: int,
    burnin_arrivals: int,
    measure_arrivals: int,
    max_steps: int | None,
) -> Stat:
    """
    Isothermal ε (prompt, formed-origin) measured using arrival-limited windows.

    This is much faster than exposure-limited windows (few × 1e15 atoms/cm²) and
    is sufficient for a τ sensitivity check where we only need the observable to
    be stable vs τ, not to exactly mirror the lab dose.
    """
    p = dict(base_params)
    p["beam_dissociation_fraction"] = float(tau)
    p["surface_temperature_k"] = float(temperature_k)

    eps_runs: List[float] = []
    for i in range(int(replicates)):
        sim_params = dict(p)
        sim_params["rng_seed"] = 2000 + i

        kmc = KineticMonteCarlo({**sim_params, "max_arrivals": int(burnin_arrivals)})
        if burnin_arrivals > 0:
            kmc.run_gillespie(max_time=1e30, max_steps=max_steps)

        _reset_measurement_counters(kmc)
        kmc.simulation_parameters["max_arrivals"] = int(measure_arrivals)
        kmc.run_gillespie(max_time=1e30, max_steps=max_steps)

        imp = int(getattr(kmc, "total_impinging_h_atoms", 0))
        if imp <= 0:
            eps_runs.append(0.0)
            continue
        prompt = float(getattr(kmc, "h2_molecules_desorbed", 0))
        eps_runs.append(float(2.0 * prompt / float(imp)))

    return _stat(eps_runs)


def _reset_measurement_counters(kmc: KineticMonteCarlo) -> None:
    kmc.total_impinging_h_atoms = 0
    kmc.total_impinging_h2_molecules = 0
    kmc.total_adsorbed_h_atoms = 0
    kmc.total_desorbed_h_atoms = 0
    kmc.h2_molecules_formed = 0
    kmc.h2_molecules_desorbed = 0
    kmc.h2_molecules_desorbed_LH = 0
    kmc.h2_molecules_desorbed_ER = 0
    kmc.h2_molecules_desorbed_UV = 0
    kmc.h2_molecules_desorbed_beam = 0
    kmc.h2_molecules_released_formed = 0
    kmc.h2_molecules_released_beam = 0
    kmc.h2_molecules_formed_LH = 0
    kmc.h2_molecules_formed_ER = 0
    kmc.h2_molecules_formed_UV = 0


def _iso_eps_released_total_at_T(
    base_params: Dict[str, Any],
    temperature_k: float,
    tau: float,
    replicates: int,
    burnin_arrivals: int,
    measure_arrivals: int,
    max_steps: int | None,
    fast_diffusion: bool,
) -> Stat:
    p = dict(base_params)
    p["beam_dissociation_fraction"] = float(tau)
    p["surface_temperature_k"] = float(temperature_k)
    if fast_diffusion:
        # For τ sensitivity we only need the *observable* to be stable vs τ; we don't need the
        # expensive explicit-diffusion TPDED microphysics. Use a fast-mixing LH approximation.
        p["diffusion_mode"] = "rate_only"
        p["lh_formation_mode"] = "diffusion_limited"

    eps_runs: List[float] = []
    for i in range(int(replicates)):
        sim_params = dict(p)
        sim_params["rng_seed"] = 3000 + i

        # Burn-in by arrivals (isothermal).
        kmc = KineticMonteCarlo({**sim_params, "max_arrivals": int(burnin_arrivals)})
        if burnin_arrivals > 0:
            kmc.run_gillespie(max_time=1e30, max_steps=max_steps)

        _reset_measurement_counters(kmc)
        kmc.simulation_parameters["max_arrivals"] = int(measure_arrivals)
        kmc.run_gillespie(max_time=1e30, max_steps=max_steps)

        imp = int(getattr(kmc, "total_impinging_h_atoms", 0))
        if imp <= 0:
            eps_runs.append(0.0)
            continue
        released_total = float(getattr(kmc, "h2_molecules_desorbed", 0)) + float(
            getattr(kmc, "h2_molecules_released_formed", 0)
        )
        eps_runs.append(float(2.0 * released_total / float(imp)))

    return _stat(eps_runs)


def _plot(out_png: str, rows: List[Dict[str, Any]]) -> None:
    taus = [float(r["tau"]) for r in rows]

    eps20 = [float(r["eps20_isothermal_released_total_mean"]) for r in rows]
    eps20_ci = [float(r["eps20_isothermal_released_total_ci95"]) for r in rows]

    eps200 = [float(r["eps200_isothermal_mean"]) for r in rows]
    eps200_ci = [float(r["eps200_isothermal_ci95"]) for r in rows]

    plt.figure(figsize=(7.5, 4.8))
    plt.errorbar(taus, eps20, yerr=eps20_ci, fmt="o-", capsize=3, label="ε(20 K) isothermal released_total")
    plt.errorbar(taus, eps200, yerr=eps200_ci, fmt="^-", capsize=3, label="ε(200 K) isothermal")
    plt.xlabel("Dissociation fraction τ")
    plt.ylabel("Recombination efficiency ε (fraction)")
    plt.title("τ sensitivity check (paperfit baseline)")
    plt.grid(True, alpha=0.3)
    plt.ylim(0.0, 0.6)
    plt.legend(frameon=False)
    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def main() -> None:
    p = argparse.ArgumentParser(description="τ sensitivity check for the Grieco paperfit baseline.")
    p.add_argument("--taus", nargs="+", type=float, default=[0.5, 0.9])

    p.add_argument("--iso-config", default="config_grieco_paper_iso_paperfit.yaml")
    p.add_argument("--ded-config", default="config_grieco_paper_ded_paperfit.yaml")

    p.add_argument("--iso-temp", type=float, default=200.0)
    p.add_argument("--iso-replicates", type=int, default=200)
    p.add_argument("--iso-burnin-arrivals", type=int, default=2000)
    p.add_argument("--iso-measure-arrivals", type=int, default=5000)
    p.add_argument("--iso-max-steps", type=int, default=500000)
    p.add_argument(
        "--iso-use-exposure-stops",
        action="store_true",
        help="Use burnin_exposure_atoms_cm2 / measure_exposure_atoms_cm2 from the isothermal config (slow).",
    )

    p.add_argument("--lowT-temp", type=float, default=20.0)
    p.add_argument("--lowT-replicates", type=int, default=200)
    p.add_argument("--lowT-burnin-arrivals", type=int, default=2000)
    p.add_argument("--lowT-measure-arrivals", type=int, default=5000)
    p.add_argument("--lowT-max-steps", type=int, default=800000)
    p.add_argument(
        "--lowT-explicit-diffusion",
        action="store_true",
        help="Use explicit diffusion + adjacent-pairs LH for the low-T 20 K τ check (slow).",
    )

    p.add_argument("--out-table", default="results/tables/grieco_tau_sensitivity_table.csv")
    p.add_argument("--out-plot", default="results/plots/grieco_tau_sensitivity.png")
    p.add_argument("--out-prefix", default="results/grieco_tau_sensitivity")
    args = p.parse_args()

    iso_cfg = _load_yaml(str(args.iso_config))
    ded_cfg = _load_yaml(str(args.ded_config))

    # Speed: reuse the same grain topology/energetics across ensemble runs.
    iso_cfg.setdefault("enable_grain_cache", True)
    iso_cfg.setdefault("grain_cache_dir", "grain_cache")
    iso_cfg.setdefault("grain_cache_include_rng_seed", False)
    ded_cfg.setdefault("enable_grain_cache", True)
    ded_cfg.setdefault("grain_cache_dir", "grain_cache")
    ded_cfg.setdefault("grain_cache_include_rng_seed", False)

    burnin_exposure = iso_cfg.get("burnin_exposure_atoms_cm2", None)
    measure_exposure = iso_cfg.get("measure_exposure_atoms_cm2", None)

    rows: List[Dict[str, Any]] = []
    for tau in [float(x) for x in args.taus]:
        if bool(args.iso_use_exposure_stops) and (burnin_exposure is not None or measure_exposure is not None):
            iso_stat = _iso_eps_at_T(
                base_params=iso_cfg,
                temperature_k=float(args.iso_temp),
                tau=float(tau),
                replicates=int(args.iso_replicates),
                burnin_arrivals=0,
                measure_arrivals=0,
                burnin_exposure_atoms_cm2=float(burnin_exposure) if burnin_exposure is not None else None,
                measure_exposure_atoms_cm2=float(measure_exposure) if measure_exposure is not None else None,
                max_steps=int(args.iso_max_steps) if args.iso_max_steps else None,
            )
        else:
            iso_stat = _iso_eps_prompt_at_T_arrivals(
                base_params=iso_cfg,
                temperature_k=float(args.iso_temp),
                tau=float(tau),
                replicates=int(args.iso_replicates),
                burnin_arrivals=int(args.iso_burnin_arrivals),
                measure_arrivals=int(args.iso_measure_arrivals),
                max_steps=int(args.iso_max_steps) if args.iso_max_steps else None,
            )

        lowT_stat = _iso_eps_released_total_at_T(
            base_params=ded_cfg,
            temperature_k=float(args.lowT_temp),
            tau=float(tau),
            replicates=int(args.lowT_replicates),
            burnin_arrivals=int(args.lowT_burnin_arrivals),
            measure_arrivals=int(args.lowT_measure_arrivals),
            max_steps=int(args.lowT_max_steps) if args.lowT_max_steps else None,
            fast_diffusion=not bool(args.lowT_explicit_diffusion),
        )

        rows.append(
            {
                "tau": float(tau),
                "eps20_isothermal_released_total_mean": float(lowT_stat.mean),
                "eps20_isothermal_released_total_ci95": float(lowT_stat.ci95),
                "eps200_isothermal_mean": float(iso_stat.mean),
                "eps200_isothermal_ci95": float(iso_stat.ci95),
                "iso_temp_k": float(args.iso_temp),
                "iso_replicates": int(args.iso_replicates),
                "iso_burnin_arrivals": int(args.iso_burnin_arrivals),
                "iso_measure_arrivals": int(args.iso_measure_arrivals),
                "iso_use_exposure_stops": bool(args.iso_use_exposure_stops),
                "lowT_temp_k": float(args.lowT_temp),
                "lowT_replicates": int(args.lowT_replicates),
                "lowT_burnin_arrivals": int(args.lowT_burnin_arrivals),
                "lowT_measure_arrivals": int(args.lowT_measure_arrivals),
                "lowT_fast_diffusion": not bool(args.lowT_explicit_diffusion),
            }
        )

    rows = sorted(rows, key=lambda r: float(r["tau"]))
    _write_csv(str(args.out_table), rows)
    _plot(str(args.out_plot), rows)
    print(f"Wrote {args.out_table}")
    print(f"Wrote {args.out_plot}")


if __name__ == "__main__":
    main()
