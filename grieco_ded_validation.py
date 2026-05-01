from __future__ import annotations

import argparse
import csv
import json
import os
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import yaml

from kmc_simulation import KineticMonteCarlo


@dataclass(frozen=True)
class DedBin:
    temperature_k: float
    impinging_atoms: int
    h2_desorbed: int
    h2_released_total: int
    h2_desorbed_LH: int
    h2_desorbed_ER: int
    h2_desorbed_UV: int
    epsilon_prompt: float
    epsilon_released_total: float
    theta_h2_mean: float


def _default_ded_params() -> Dict[str, float]:
    """
    A pragmatic parameter set for *protocol-shape* validation:
    - <100 K: DED-style ramp + allow LH/diffusion + include H2 blocking at low T
    - >100 K: use the existing isothermal validation harness (grieco_validation.py)
    """
    return {
        "grain_radius_um": 0.005,
        "site_area_angstroms_sq": 25,
        "use_3d_lattice": True,
        "porosity_fraction": 0.0,
        "chemisorption_fraction": 0.5,
        "surface_defect_fraction": 0.15,
        # Energetics
        "E_phys_mean_meV": 45.0,
        "heterogeneity_E_bind_sigma_meV": 5.0,
        "E_chem_mean_eV": 1.75,
        "heterogeneity_E_chem_sigma_eV": 0.25,
        # Enable surface mechanisms for the low-T part of the curve
        "enable_LH": True,
        "enable_diffusion": True,
        "uv_flux_factor": 0.0,
        "uv_pulse_enabled": False,
        # Beam arrivals (~1 landing per site per 100 s)
        "arrival_rate_per_site_s": 0.01,
        # Sticking: calibrated effective accretion probability for this model/protocol harness.
        "sticking_probability": 0.3,
        "sticking_temp_model": "constant",
        # Arrival-driven ER/abstraction at chemisorbed reservoir (tunable)
        "er_cross_section_cm2": 1.5e-15,
        "er_reaction_probability": 0.9,
        # Low-T blocking by adsorbed molecules
        "enable_h2_blocking": True,
        # Use a stronger effective binding at 10–20 K so coverage can build and block accretion,
        # reproducing the sharp low-T efficiency reduction reported by Grieco et al.
        "E_h2_bind_eV": 0.06,
        "h2_desorption_prefactor_s": 1e12,
        "h2_stick_transition_K": 20.0,
        "h2_stick_prob_lowT": 0.9,
        "er_blocking_strength": 1.0,
        "sticking_blocking_strength": 1.0,
        # Imperfect dissociation: each arrival is an atom with prob tau, else a molecule that may stick/block
        "beam_dissociation_fraction": 1.0,
        "h2_beam_stick_probability": 1.0,
        # Gas values are unused in arrival mode, but keep defined.
        "gas_temperature_k": 300.0,
        "h_gas_density_cm3": 0.0,
    }


def _read_yaml_params(path: str) -> Dict[str, float]:
    with open(path, "r") as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise TypeError(f"Expected mapping YAML in {path}, got {type(data)}")
    return data


def _run_ded_once(
    base_params: Dict[str, float],
    t_start_k: float,
    t_end_k: float,
    rate_k_per_min: float,
    bin_width_k: float,
    burnin_arrivals: int,
    seed: int,
    max_steps: int | None,
) -> List[DedBin]:
    params = dict(base_params)
    params["rng_seed"] = int(seed)
    params["surface_temperature_k"] = float(t_start_k)

    # Optional burn-in at fixed low temperature to reach a quasi steady state before ramping.
    if burnin_arrivals > 0:
        params["max_arrivals"] = int(burnin_arrivals)
        kmc = KineticMonteCarlo(params)
        kmc.run_gillespie(max_time=1e30, max_steps=max_steps)
    else:
        kmc = KineticMonteCarlo(params)

    # Reset counters for the ramp measurement window (keep surface state).
    kmc.time = 0.0
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

    # Configure temperature ramp.
    rate_k_per_s = float(rate_k_per_min) / 60.0
    duration_s = float(t_end_k - t_start_k) / max(rate_k_per_s, 1e-30)
    kmc.simulation_parameters.pop("max_arrivals", None)
    kmc.simulation_parameters["temp_ramp"] = {
        "enabled": True,
        "T_start_K": float(t_start_k),
        "T_end_K": float(t_end_k),
        "rate_K_per_min": float(rate_k_per_min),
        "t0_s": 0.0,
    }

    n_bins = int(np.ceil((t_end_k - t_start_k) / bin_width_k))
    atoms = np.zeros(n_bins, dtype=int)
    h2_des = np.zeros(n_bins, dtype=int)
    h2_released_total = np.zeros(n_bins, dtype=int)
    h2_des_lh = np.zeros(n_bins, dtype=int)
    h2_des_er = np.zeros(n_bins, dtype=int)
    h2_des_uv = np.zeros(n_bins, dtype=int)

    # Time-weighted coverage tracking: accumulate ∫ theta_h2 dt per temperature bin.
    theta_dt = np.zeros(n_bins, dtype=float)
    time_dt = np.zeros(n_bins, dtype=float)

    prev_atoms = 0
    prev_h2_des = 0
    prev_h2_released_total = 0
    prev_h2_lh = 0
    prev_h2_er = 0
    prev_h2_uv = 0
    prev_time = 0.0
    denom_sites = float(getattr(kmc, "total_accessible_surface_sites", 0) or 0)
    prev_theta = 0.0
    if denom_sites > 0:
        prev_theta = float(getattr(kmc, "h2_molecules_on_surface", 0) or 0) / denom_sites

    def callback(sim: KineticMonteCarlo, _event: str) -> None:
        nonlocal prev_atoms, prev_h2_des, prev_h2_released_total, prev_h2_lh, prev_h2_er, prev_h2_uv, prev_time, prev_theta
        # Temperature at the event time.
        t_now = float(sim.time)
        T_now = float(t_start_k) + rate_k_per_s * t_now
        T_now = min(float(t_end_k), max(float(t_start_k), T_now))

        idx = int((T_now - float(t_start_k)) / float(bin_width_k))
        if idx < 0:
            idx = 0
        if idx >= n_bins:
            idx = n_bins - 1

        # Attribute the *interval since the previous event* to the previous coverage.
        dt = float(getattr(sim, "last_delta_t", 0.0) or 0.0)
        if dt > 0:
            t_mid = float(prev_time) + 0.5 * dt
            T_mid = float(t_start_k) + rate_k_per_s * t_mid
            T_mid = min(float(t_end_k), max(float(t_start_k), T_mid))
            idx_mid = int((T_mid - float(t_start_k)) / float(bin_width_k))
            idx_mid = max(0, min(n_bins - 1, idx_mid))
            theta_dt[idx_mid] += dt * float(prev_theta)
            time_dt[idx_mid] += dt

        cur_atoms = int(sim.total_impinging_h_atoms)
        cur_h2_des = int(sim.h2_molecules_desorbed)
        cur_h2_released_total = int(sim.h2_molecules_desorbed + sim.h2_molecules_released_formed)
        cur_h2_lh = int(getattr(sim, "h2_molecules_desorbed_LH", 0))
        cur_h2_er = int(getattr(sim, "h2_molecules_desorbed_ER", 0))
        cur_h2_uv = int(getattr(sim, "h2_molecules_desorbed_UV", 0))
        d_atoms = cur_atoms - prev_atoms
        d_h2_des = cur_h2_des - prev_h2_des
        d_h2_released_total = cur_h2_released_total - prev_h2_released_total
        if d_atoms:
            atoms[idx] += int(d_atoms)
        if d_h2_des:
            h2_des[idx] += int(d_h2_des)
        if d_h2_released_total:
            h2_released_total[idx] += int(d_h2_released_total)
        d_lh = cur_h2_lh - prev_h2_lh
        d_er = cur_h2_er - prev_h2_er
        d_uv = cur_h2_uv - prev_h2_uv
        if d_lh:
            h2_des_lh[idx] += int(d_lh)
        if d_er:
            h2_des_er[idx] += int(d_er)
        if d_uv:
            h2_des_uv[idx] += int(d_uv)
        prev_atoms = cur_atoms
        prev_h2_des = cur_h2_des
        prev_h2_released_total = cur_h2_released_total
        prev_h2_lh = cur_h2_lh
        prev_h2_er = cur_h2_er
        prev_h2_uv = cur_h2_uv
        prev_time = t_now

        denom = float(getattr(sim, "total_accessible_surface_sites", 0) or 0)
        if denom > 0:
            prev_theta = float(getattr(sim, "h2_molecules_on_surface", 0) or 0) / denom
        else:
            prev_theta = 0.0

    kmc.run_gillespie(max_time=duration_s, max_steps=max_steps, callback=callback)

    # Account for the final interval from the last event to the end of the ramp.
    tail_dt = float(duration_s) - float(prev_time)
    if tail_dt > 0:
        t_mid = float(prev_time) + 0.5 * tail_dt
        T_mid = float(t_start_k) + rate_k_per_s * t_mid
        T_mid = min(float(t_end_k), max(float(t_start_k), T_mid))
        idx_mid = int((T_mid - float(t_start_k)) / float(bin_width_k))
        idx_mid = max(0, min(n_bins - 1, idx_mid))
        theta_dt[idx_mid] += tail_dt * float(prev_theta)
        time_dt[idx_mid] += tail_dt

    out: List[DedBin] = []
    for i in range(n_bins):
        t_center = float(t_start_k) + (i + 0.5) * float(bin_width_k)
        denom = float(max(atoms[i], 1))
        eps_prompt = float(2.0 * float(h2_des[i]) / denom)
        eps_released_total = float(2.0 * float(h2_released_total[i]) / denom)
        theta_mean = float(theta_dt[i] / time_dt[i]) if float(time_dt[i]) > 0 else 0.0
        out.append(
            DedBin(
                temperature_k=t_center,
                impinging_atoms=int(atoms[i]),
                h2_desorbed=int(h2_des[i]),
                h2_released_total=int(h2_released_total[i]),
                h2_desorbed_LH=int(h2_des_lh[i]),
                h2_desorbed_ER=int(h2_des_er[i]),
                h2_desorbed_UV=int(h2_des_uv[i]),
                epsilon_prompt=eps_prompt,
                epsilon_released_total=eps_released_total,
                theta_h2_mean=theta_mean,
            )
        )
    return out


def run_ded_validation(
    output_csv: str,
    summary_json: str | None,
    replicates: int,
    t_start_k: float,
    t_end_k: float,
    rate_k_per_min: float,
    bin_width_k: float,
    burnin_arrivals: int,
    max_steps: int | None,
    base_params: Dict[str, float] | None = None,
) -> None:
    if base_params is None:
        base_params = _default_ded_params()

    all_bins: List[List[DedBin]] = [
        _run_ded_once(
            base_params=base_params,
            t_start_k=float(t_start_k),
            t_end_k=float(t_end_k),
            rate_k_per_min=float(rate_k_per_min),
            bin_width_k=float(bin_width_k),
            burnin_arrivals=int(burnin_arrivals),
            seed=2000 + i,
            max_steps=max_steps,
        )
        for i in range(replicates)
    ]

    def _summary_for_run(run_bins: List[DedBin]) -> Dict[str, float]:
        def _mean_eps(lo: float, hi: float) -> float:
            vals = [
                float(b.epsilon_released_total)
                for b in run_bins
                if float(b.temperature_k) >= lo and float(b.temperature_k) <= hi
            ]
            return float(np.mean(np.array(vals, dtype=float))) if vals else 0.0

        eps10 = float(_mean_eps(9.5, 11.5))
        eps20 = float(_mean_eps(19.5, 21.5))
        eps30_80 = float(_mean_eps(30.0, 80.0))
        ratio = float(eps10 / eps20) if eps20 > 0 else 0.0
        return {
            "eps10": eps10,
            "eps20": eps20,
            "eps30_80": eps30_80,
            "ratio10_over_20": ratio,
        }

    per_run_summaries = [_summary_for_run(bins) for bins in all_bins]
    eps10_runs = np.array([s["eps10"] for s in per_run_summaries], dtype=float)
    eps20_runs = np.array([s["eps20"] for s in per_run_summaries], dtype=float)
    eps30_runs = np.array([s["eps30_80"] for s in per_run_summaries], dtype=float)
    ratio_runs = np.array([s["ratio10_over_20"] for s in per_run_summaries], dtype=float)

    # Aggregate by bin index.
    rows: List[dict] = []
    for b in range(len(all_bins[0])):
        temps = np.array([run[b].temperature_k for run in all_bins], dtype=float)
        eps_prompt = np.array([run[b].epsilon_prompt for run in all_bins], dtype=float)
        eps_released_total = np.array([run[b].epsilon_released_total for run in all_bins], dtype=float)
        atoms = np.array([run[b].impinging_atoms for run in all_bins], dtype=float)
        h2d = np.array([run[b].h2_desorbed for run in all_bins], dtype=float)
        h2rt = np.array([run[b].h2_released_total for run in all_bins], dtype=float)
        theta = np.array([run[b].theta_h2_mean for run in all_bins], dtype=float)
        h2_lh = np.array([run[b].h2_desorbed_LH for run in all_bins], dtype=float)
        h2_er = np.array([run[b].h2_desorbed_ER for run in all_bins], dtype=float)
        h2_uv = np.array([run[b].h2_desorbed_UV for run in all_bins], dtype=float)

        # Mechanism fractions (computed from aggregated counts for stability).
        tot_h2 = float(np.sum(h2d))
        frac_lh = float(np.sum(h2_lh) / tot_h2) if tot_h2 > 0 else 0.0
        frac_er = float(np.sum(h2_er) / tot_h2) if tot_h2 > 0 else 0.0
        frac_uv = float(np.sum(h2_uv) / tot_h2) if tot_h2 > 0 else 0.0
        rows.append(
            {
                "temperature_k": float(np.mean(temps)),
                "epsilon_prompt_mean": float(np.mean(eps_prompt)),
                "epsilon_prompt_std": float(np.std(eps_prompt, ddof=1)) if replicates > 1 else 0.0,
                "epsilon_prompt_ci95": float(1.96 * float(np.std(eps_prompt, ddof=1)) / float(np.sqrt(replicates))) if replicates > 1 else 0.0,
                "epsilon_released_total_mean": float(np.mean(eps_released_total)),
                "epsilon_released_total_std": float(np.std(eps_released_total, ddof=1)) if replicates > 1 else 0.0,
                "epsilon_released_total_ci95": float(1.96 * float(np.std(eps_released_total, ddof=1)) / float(np.sqrt(replicates))) if replicates > 1 else 0.0,
                "impinging_atoms_mean": float(np.mean(atoms)),
                "h2_desorbed_mean": float(np.mean(h2d)),
                "h2_released_total_mean": float(np.mean(h2rt)),
                "theta_h2_mean": float(np.mean(theta)),
                "theta_h2_std": float(np.std(theta, ddof=1)) if replicates > 1 else 0.0,
                "theta_h2_ci95": float(1.96 * float(np.std(theta, ddof=1)) / float(np.sqrt(replicates))) if replicates > 1 else 0.0,
                "h2_desorbed_LH_mean": float(np.mean(h2_lh)),
                "h2_desorbed_ER_mean": float(np.mean(h2_er)),
                "h2_desorbed_UV_mean": float(np.mean(h2_uv)),
                "frac_desorbed_LH": frac_lh,
                "frac_desorbed_ER": frac_er,
                "frac_desorbed_UV": frac_uv,
            }
        )

    os.makedirs(os.path.dirname(output_csv) or ".", exist_ok=True)
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    def _mean_eps(lo: float, hi: float) -> float:
        vals = [
            float(r["epsilon_released_total_mean"])
            for r in rows
            if float(r["temperature_k"]) >= lo and float(r["temperature_k"]) <= hi
        ]
        return float(np.mean(np.array(vals, dtype=float))) if vals else 0.0

    eps10 = _mean_eps(9.5, 11.5)
    eps20 = _mean_eps(19.5, 21.5)
    eps30_80 = _mean_eps(30.0, 80.0)
    ratio = float(eps10 / eps20) if eps20 > 0 else 0.0
    print(f"Wrote {len(rows)} DED bins to {output_csv}")
    print(f"DED summary: eps10={eps10:.4g} eps20={eps20:.4g} eps30_80={eps30_80:.4g} ratio10/20={ratio:.4g}")

    if replicates > 1:
        eps10_ci = float(1.96 * float(np.std(eps10_runs, ddof=1)) / float(np.sqrt(replicates)))
        eps20_ci = float(1.96 * float(np.std(eps20_runs, ddof=1)) / float(np.sqrt(replicates)))
        eps30_ci = float(1.96 * float(np.std(eps30_runs, ddof=1)) / float(np.sqrt(replicates)))
        ratio_ci = float(1.96 * float(np.std(ratio_runs, ddof=1)) / float(np.sqrt(replicates)))
        print(
            "DED summary (mean ± 95% CI over replicates): "
            f"eps10={float(np.mean(eps10_runs)):.4g}±{eps10_ci:.3g} "
            f"eps20={float(np.mean(eps20_runs)):.4g}±{eps20_ci:.3g} "
            f"eps30_80={float(np.mean(eps30_runs)):.4g}±{eps30_ci:.3g} "
            f"ratio10/20={float(np.mean(ratio_runs)):.4g}±{ratio_ci:.3g}"
        )

    if summary_json:
        summary = {
            "replicates": int(replicates),
            "t_start_k": float(t_start_k),
            "t_end_k": float(t_end_k),
            "rate_k_per_min": float(rate_k_per_min),
            "bin_width_k": float(bin_width_k),
            "eps10_mean": float(np.mean(eps10_runs)),
            "eps20_mean": float(np.mean(eps20_runs)),
            "eps30_80_mean": float(np.mean(eps30_runs)),
            "ratio10_over_20_mean": float(np.mean(ratio_runs)),
        }
        if replicates > 1:
            summary.update(
                {
                    "eps10_ci95": float(1.96 * float(np.std(eps10_runs, ddof=1)) / float(np.sqrt(replicates))),
                    "eps20_ci95": float(1.96 * float(np.std(eps20_runs, ddof=1)) / float(np.sqrt(replicates))),
                    "eps30_80_ci95": float(1.96 * float(np.std(eps30_runs, ddof=1)) / float(np.sqrt(replicates))),
                    "ratio10_over_20_ci95": float(1.96 * float(np.std(ratio_runs, ddof=1)) / float(np.sqrt(replicates))),
                }
            )
        os.makedirs(os.path.dirname(str(summary_json)) or ".", exist_ok=True)
        with open(str(summary_json), "w") as f:
            json.dump(summary, f, indent=2, sort_keys=True)
        print(f"Wrote DED summary JSON to {summary_json}")


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="DED-style (<100 K) validation with temperature ramp and H2 blocking.")
    p.add_argument("--base-config", default=None, help="Optional YAML file providing a base parameter set")
    p.add_argument("--output", default="results/grieco_ded_validation.csv", help="Output CSV path")
    p.add_argument("--summary-json", default=None, help="Optional output JSON with key DED metrics and CI")
    p.add_argument("--replicates", type=int, default=3, help="Independent ramps to average")
    p.add_argument("--t-start", type=float, default=10.0, help="Ramp start temperature (K)")
    p.add_argument("--t-end", type=float, default=100.0, help="Ramp end temperature (K)")
    p.add_argument("--rate-k-per-min", type=float, default=1.0, help="Ramp rate (K/min)")
    p.add_argument("--bin-width", type=float, default=1.0, help="Bin width (K)")
    p.add_argument("--burnin-arrivals", type=int, default=2000, help="Arrivals discarded for burn-in at T_start")
    p.add_argument("--max-steps", type=int, default=800000, help="Hard cap on KMC steps per replicate")
    p.add_argument("--arrival-rate-per-site", type=float, default=None, help="Override arrival_rate_per_site_s")
    p.add_argument("--uv-flux-factor", type=float, default=None, help="Override uv_flux_factor")
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    base = _default_ded_params()
    if args.base_config:
        base.update(_read_yaml_params(str(args.base_config)))
    if args.arrival_rate_per_site is not None:
        base["arrival_rate_per_site_s"] = float(args.arrival_rate_per_site)
    if args.uv_flux_factor is not None:
        base["uv_flux_factor"] = float(args.uv_flux_factor)
    run_ded_validation(
        output_csv=str(args.output),
        summary_json=str(args.summary_json) if args.summary_json else None,
        replicates=int(args.replicates),
        t_start_k=float(args.t_start),
        t_end_k=float(args.t_end),
        rate_k_per_min=float(args.rate_k_per_min),
        bin_width_k=float(args.bin_width),
        burnin_arrivals=int(args.burnin_arrivals),
        max_steps=int(args.max_steps) if args.max_steps else None,
        base_params=base,
    )
