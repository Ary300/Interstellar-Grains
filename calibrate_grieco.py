import argparse
import csv
import json
import os
from dataclasses import dataclass
from typing import Any, Dict, List, Tuple

import numpy as np

from kmc_simulation import KineticMonteCarlo


@dataclass(frozen=True)
class IsoSummary:
    temperatures_k: List[float]
    epsilon_mean_by_T: Dict[float, float]
    epsilon_mean_plateau: float


@dataclass(frozen=True)
class DedSummary:
    bin_centers_k: List[float]
    epsilon_by_bin: List[float]
    epsilon_10k: float
    epsilon_20k: float
    epsilon_30_80k: float


def _reset_measurement_counters(kmc: KineticMonteCarlo) -> None:
    kmc.total_impinging_h_atoms = 0
    kmc.total_impinging_h2_molecules = 0
    kmc.total_adsorbed_h_atoms = 0
    kmc.total_desorbed_h_atoms = 0
    kmc.h2_molecules_formed = 0
    kmc.h2_molecules_desorbed = 0  # prompt recombination signal
    kmc.h2_molecules_desorbed_LH = 0
    kmc.h2_molecules_desorbed_ER = 0
    kmc.h2_molecules_desorbed_UV = 0
    kmc.h2_molecules_desorbed_beam = 0
    kmc.h2_molecules_released_formed = 0
    kmc.h2_molecules_released_beam = 0
    kmc.h2_molecules_formed_LH = 0
    kmc.h2_molecules_formed_ER = 0
    kmc.h2_molecules_formed_UV = 0


def _isothermal_epsilon(
    params: Dict[str, Any],
    temperature_k: float,
    burnin_arrivals: int,
    measure_arrivals: int,
    seed: int,
    max_steps: int | None,
) -> float:
    sim_params = dict(params)
    sim_params["surface_temperature_k"] = float(temperature_k)
    sim_params["rng_seed"] = int(seed)

    sim_params["max_arrivals"] = int(burnin_arrivals)
    kmc = KineticMonteCarlo(sim_params)
    if burnin_arrivals > 0:
        kmc.run_gillespie(max_time=1e30, max_steps=max_steps)

    _reset_measurement_counters(kmc)
    kmc.simulation_parameters["max_arrivals"] = int(measure_arrivals)
    kmc.run_gillespie(max_time=1e30, max_steps=max_steps)

    if kmc.total_impinging_h_atoms <= 0:
        return 0.0
    return float(2.0 * kmc.h2_molecules_desorbed / float(kmc.total_impinging_h_atoms))


def evaluate_highT_plateau(
    base_params: Dict[str, Any],
    temperatures_k: List[float],
    replicates: int,
    burnin_arrivals: int,
    measure_arrivals: int,
    max_steps: int | None,
) -> IsoSummary:
    eps_by_T: Dict[float, float] = {}
    for T in temperatures_k:
        eps_runs = [
            _isothermal_epsilon(
                params=base_params,
                temperature_k=float(T),
                burnin_arrivals=burnin_arrivals,
                measure_arrivals=measure_arrivals,
                seed=1000 + i,
                max_steps=max_steps,
            )
            for i in range(replicates)
        ]
        eps_by_T[float(T)] = float(np.mean(np.array(eps_runs, dtype=float)))

    plateau = float(np.mean(np.array(list(eps_by_T.values()), dtype=float))) if eps_by_T else 0.0
    return IsoSummary(
        temperatures_k=[float(t) for t in temperatures_k],
        epsilon_mean_by_T=eps_by_T,
        epsilon_mean_plateau=plateau,
    )


def _ded_binned_epsilon(
    params: Dict[str, Any],
    t_start_k: float,
    t_end_k: float,
    rate_k_per_min: float,
    bin_width_k: float,
    burnin_arrivals: int,
    seed: int,
    max_steps: int | None,
) -> Tuple[List[float], List[float]]:
    sim_params = dict(params)
    sim_params["rng_seed"] = int(seed)
    sim_params["surface_temperature_k"] = float(t_start_k)

    # Burn-in is isothermal at T_start_k.
    if burnin_arrivals > 0:
        sim_params["max_arrivals"] = int(burnin_arrivals)
        kmc = KineticMonteCarlo(sim_params)
        kmc.run_gillespie(max_time=1e30, max_steps=max_steps)
    else:
        kmc = KineticMonteCarlo(sim_params)

    # Measurement window: ramp.
    kmc.time = 0.0
    _reset_measurement_counters(kmc)
    kmc.simulation_parameters.pop("max_arrivals", None)
    kmc.simulation_parameters["temp_ramp"] = {
        "enabled": True,
        "T_start_K": float(t_start_k),
        "T_end_K": float(t_end_k),
        "rate_K_per_min": float(rate_k_per_min),
        "t0_s": 0.0,
    }

    rate_k_per_s = float(rate_k_per_min) / 60.0
    duration_s = float(t_end_k - t_start_k) / max(rate_k_per_s, 1e-30)

    n_bins = int(np.ceil((t_end_k - t_start_k) / float(bin_width_k)))
    atoms = np.zeros(n_bins, dtype=int)
    h2_prompt = np.zeros(n_bins, dtype=int)

    prev_atoms = 0
    prev_h2_prompt = 0

    def cb(sim: KineticMonteCarlo, _event: str) -> None:
        nonlocal prev_atoms, prev_h2_prompt
        t_now = float(sim.time)
        T_now = float(t_start_k) + rate_k_per_s * t_now
        T_now = min(float(t_end_k), max(float(t_start_k), T_now))
        idx = int((T_now - float(t_start_k)) / float(bin_width_k))
        idx = max(0, min(n_bins - 1, idx))

        cur_atoms = int(sim.total_impinging_h_atoms)
        cur_h2_prompt = int(sim.h2_molecules_desorbed)
        d_atoms = cur_atoms - prev_atoms
        d_h2 = cur_h2_prompt - prev_h2_prompt
        if d_atoms:
            atoms[idx] += int(d_atoms)
        if d_h2:
            h2_prompt[idx] += int(d_h2)
        prev_atoms = cur_atoms
        prev_h2_prompt = cur_h2_prompt

    kmc.run_gillespie(max_time=duration_s, max_steps=max_steps, callback=cb)

    bin_centers = [float(t_start_k) + (i + 0.5) * float(bin_width_k) for i in range(n_bins)]
    eps_bins: List[float] = []
    for i in range(n_bins):
        denom = float(atoms[i]) if atoms[i] > 0 else 1.0
        eps_bins.append(float(2.0 * float(h2_prompt[i]) / denom))
    return bin_centers, eps_bins


def evaluate_ded(
    base_params: Dict[str, Any],
    replicates: int,
    burnin_arrivals: int,
    t_start_k: float,
    t_end_k: float,
    rate_k_per_min: float,
    bin_width_k: float,
    max_steps: int | None,
) -> DedSummary:
    all_eps: List[List[float]] = []
    bin_centers: List[float] = []
    for i in range(replicates):
        centers, eps = _ded_binned_epsilon(
            params=base_params,
            t_start_k=float(t_start_k),
            t_end_k=float(t_end_k),
            rate_k_per_min=float(rate_k_per_min),
            bin_width_k=float(bin_width_k),
            burnin_arrivals=int(burnin_arrivals),
            seed=2000 + i,
            max_steps=max_steps,
        )
        bin_centers = centers
        all_eps.append(eps)

    eps_mean = np.mean(np.array(all_eps, dtype=float), axis=0)

    def _mean_in_range(lo: float, hi: float) -> float:
        vals = [float(eps_mean[i]) for i, T in enumerate(bin_centers) if (T >= lo and T <= hi)]
        return float(np.mean(np.array(vals, dtype=float))) if vals else 0.0

    return DedSummary(
        bin_centers_k=bin_centers,
        epsilon_by_bin=[float(x) for x in eps_mean.tolist()],
        epsilon_10k=_mean_in_range(9.5, 11.5),
        epsilon_20k=_mean_in_range(19.5, 21.5),
        epsilon_30_80k=_mean_in_range(30.0, 80.0),
    )


def _default_highT_params() -> Dict[str, Any]:
    return {
        "grain_radius_um": 0.005,
        "site_area_angstroms_sq": 25,
        "use_3d_lattice": True,
        "porosity_fraction": 0.0,
        "chemisorption_fraction": 0.5,
        "surface_defect_fraction": 0.15,
        "E_phys_mean_meV": 45.0,
        "heterogeneity_E_bind_sigma_meV": 5.0,
        "E_chem_mean_eV": 1.75,
        "heterogeneity_E_chem_sigma_eV": 0.25,
        # High-T plateau: focus on arrival-driven ER/abstraction against a chemisorbed reservoir.
        "enable_LH": False,
        "enable_diffusion": False,
        "uv_flux_factor": 0.0,
        "uv_pulse_enabled": False,
        "arrival_rate_per_site_s": 0.01,
        "sticking_probability": 0.5,
        "sticking_temp_model": "constant",
        "er_cross_section_cm2": 1e-15,
        "er_reaction_probability": 0.5,
        "beam_dissociation_fraction": 1.0,
        "gas_temperature_k": 300.0,
        "h_gas_density_cm3": 0.0,
        # Disable blocking at high T.
        "enable_h2_blocking": False,
    }


def _default_ded_params() -> Dict[str, Any]:
    return {
        "grain_radius_um": 0.005,
        "site_area_angstroms_sq": 25,
        "use_3d_lattice": True,
        "porosity_fraction": 0.0,
        "chemisorption_fraction": 0.5,
        "surface_defect_fraction": 0.15,
        "E_phys_mean_meV": 45.0,
        "heterogeneity_E_bind_sigma_meV": 5.0,
        "E_chem_mean_eV": 1.75,
        "heterogeneity_E_chem_sigma_eV": 0.25,
        "enable_LH": True,
        "enable_diffusion": True,
        "uv_flux_factor": 0.0,
        "uv_pulse_enabled": False,
        "arrival_rate_per_site_s": 0.01,
        "sticking_probability": 0.5,
        "sticking_temp_model": "constant",
        "er_cross_section_cm2": 1e-15,
        "er_reaction_probability": 0.5,
        "beam_dissociation_fraction": 1.0,
        "gas_temperature_k": 300.0,
        "h_gas_density_cm3": 0.0,
        # Blocking knobs
        "enable_h2_blocking": True,
        "E_h2_bind_eV": 0.06,
        "h2_desorption_prefactor_s": 1e12,
        "h2_stick_transition_K": 20.0,
        "h2_stick_prob_lowT": 0.9,
        "sticking_blocking_strength": 1.0,
        "er_blocking_strength": 1.0,
    }


def _objective_plateau(summary: IsoSummary, target: float) -> float:
    return float((summary.epsilon_mean_plateau - float(target)) ** 2)


def _objective_ded(summary: DedSummary, mid_target: float, drop_target_ratio: float) -> float:
    # Encourage 30–80 K band near target and a sharp drop at 10 K relative to ~20 K.
    mid_err = float((summary.epsilon_30_80k - float(mid_target)) ** 2)
    # Penalize if epsilon_10k is not sufficiently below epsilon_20k.
    desired_max_10 = float(summary.epsilon_20k) * float(drop_target_ratio)
    drop_pen = 0.0
    if summary.epsilon_10k > desired_max_10:
        drop_pen = float((summary.epsilon_10k - desired_max_10) ** 2)
    return mid_err + 5.0 * drop_pen


def _write_csv(path: str, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        return
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fieldnames = sorted({k for r in rows for k in r.keys()})
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    p = argparse.ArgumentParser(description="Two-stage calibration harness against Grieco-style ε(T) regimes.")
    p.add_argument("--out", default="results/grieco_calibration.json", help="Output JSON summary")
    p.add_argument("--log-csv", default="results/grieco_calibration_log.csv", help="CSV log of trial objectives")
    p.add_argument("--max-steps", type=int, default=500000, help="Max KMC steps per phase")

    p.add_argument("--plateau-target", type=float, default=0.20)
    p.add_argument("--plateau-temps", nargs="+", type=float, default=[100, 150, 200, 250])
    p.add_argument("--plateau-replicates", type=int, default=2)
    p.add_argument("--plateau-burnin", type=int, default=2000)
    p.add_argument("--plateau-measure", type=int, default=5000)
    p.add_argument("--sticking-grid", nargs="+", type=float, default=[0.5])
    p.add_argument("--er-cross-grid", nargs="+", type=float, default=[1e-15])
    p.add_argument("--er-prob-grid", nargs="+", type=float, default=[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])

    p.add_argument("--mid-target", type=float, default=0.30)
    p.add_argument("--drop-ratio", type=float, default=0.3, help="Require ε(10K) ≤ drop_ratio * ε(20K)")
    p.add_argument("--ded-replicates", type=int, default=1)
    p.add_argument("--ded-burnin", type=int, default=2000)
    p.add_argument("--ded-start", type=float, default=10.0)
    p.add_argument("--ded-end", type=float, default=100.0)
    p.add_argument("--ded-rate-k-per-min", type=float, default=1.0)
    p.add_argument("--ded-bin-width", type=float, default=1.0)

    p.add_argument("--E-h2-grid", nargs="+", type=float, default=[0.045, 0.055, 0.065, 0.075])
    p.add_argument("--stick-block-grid", nargs="+", type=float, default=[0.5, 1.0, 1.5, 2.0])
    p.add_argument("--er-block-grid", nargs="+", type=float, default=[0.5, 1.0, 1.5])

    args = p.parse_args()

    max_steps = int(args.max_steps) if args.max_steps else None

    # Stage 1: fit high-T plateau by scanning ER reaction probability.
    base_highT = _default_highT_params()
    best_plateau: IsoSummary | None = None
    best_er_prob = None
    best_er_cross = None
    best_sticking = None
    best_obj = None

    log_rows: List[Dict[str, Any]] = []
    for sticking in args.sticking_grid:
        for er_cross in args.er_cross_grid:
            for er_prob in args.er_prob_grid:
                base_highT["sticking_probability"] = float(sticking)
                base_highT["er_cross_section_cm2"] = float(er_cross)
                base_highT["er_reaction_probability"] = float(er_prob)
                summ = evaluate_highT_plateau(
                    base_params=base_highT,
                    temperatures_k=list(args.plateau_temps),
                    replicates=int(args.plateau_replicates),
                    burnin_arrivals=int(args.plateau_burnin),
                    measure_arrivals=int(args.plateau_measure),
                    max_steps=max_steps,
                )
                obj = _objective_plateau(summ, float(args.plateau_target))
                log_rows.append(
                    {
                        "stage": "highT",
                        "sticking_probability": float(sticking),
                        "er_cross_section_cm2": float(er_cross),
                        "er_reaction_probability": float(er_prob),
                        "plateau_eps": summ.epsilon_mean_plateau,
                        "objective": obj,
                    }
                )
                if best_obj is None or obj < best_obj:
                    best_obj = obj
                    best_plateau = summ
                    best_er_prob = float(er_prob)
                    best_er_cross = float(er_cross)
                    best_sticking = float(sticking)

    assert best_plateau is not None and best_er_prob is not None and best_er_cross is not None and best_sticking is not None

    # Stage 2: freeze high-T knobs; search blocking knobs for low-T ramp structure.
    base_ded = _default_ded_params()
    base_ded["er_reaction_probability"] = float(best_er_prob)
    base_ded["er_cross_section_cm2"] = float(best_er_cross)
    base_ded["sticking_probability"] = float(best_sticking)

    best_ded: DedSummary | None = None
    best_ded_params: Dict[str, Any] | None = None
    best_ded_obj = None

    for E_h2 in args.E_h2_grid:
        for sb in args.stick_block_grid:
            for eb in args.er_block_grid:
                base_ded["E_h2_bind_eV"] = float(E_h2)
                base_ded["sticking_blocking_strength"] = float(sb)
                base_ded["er_blocking_strength"] = float(eb)
                ded = evaluate_ded(
                    base_params=base_ded,
                    replicates=int(args.ded_replicates),
                    burnin_arrivals=int(args.ded_burnin),
                    t_start_k=float(args.ded_start),
                    t_end_k=float(args.ded_end),
                    rate_k_per_min=float(args.ded_rate_k_per_min),
                    bin_width_k=float(args.ded_bin_width),
                    max_steps=max_steps,
                )
                obj = _objective_ded(ded, float(args.mid_target), float(args.drop_ratio))
                log_rows.append(
                    {
                        "stage": "ded",
                        "E_h2_bind_eV": float(E_h2),
                        "sticking_blocking_strength": float(sb),
                        "er_blocking_strength": float(eb),
                        "er_reaction_probability": float(best_er_prob),
                        "er_cross_section_cm2": float(best_er_cross),
                        "sticking_probability": float(best_sticking),
                        "eps_10k": float(ded.epsilon_10k),
                        "eps_20k": float(ded.epsilon_20k),
                        "eps_30_80k": float(ded.epsilon_30_80k),
                        "objective": float(obj),
                    }
                )
                if best_ded_obj is None or obj < best_ded_obj:
                    best_ded_obj = obj
                    best_ded = ded
                    best_ded_params = {
                        "E_h2_bind_eV": float(E_h2),
                        "sticking_blocking_strength": float(sb),
                        "er_blocking_strength": float(eb),
                    }

    assert best_ded is not None and best_ded_params is not None

    out = {
        "stage1": {
            "target_plateau": float(args.plateau_target),
            "best_er_reaction_probability": float(best_er_prob),
            "best_er_cross_section_cm2": float(best_er_cross),
            "best_sticking_probability": float(best_sticking),
            "epsilon_by_T": {str(k): float(v) for k, v in best_plateau.epsilon_mean_by_T.items()},
            "plateau_epsilon_mean": float(best_plateau.epsilon_mean_plateau),
        },
        "stage2": {
            "targets": {"mid_target": float(args.mid_target), "drop_ratio": float(args.drop_ratio)},
            "best_params": best_ded_params,
            "epsilon_10k": float(best_ded.epsilon_10k),
            "epsilon_20k": float(best_ded.epsilon_20k),
            "epsilon_30_80k": float(best_ded.epsilon_30_80k),
        },
    }

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2, sort_keys=True)

    _write_csv(str(args.log_csv), log_rows)
    print(f"Wrote summary to {args.out}")
    print(f"Wrote trial log to {args.log_csv}")
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
