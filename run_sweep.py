import os
import json
import yaml
import math
from itertools import product
from typing import Dict, Any, List, Tuple
import numpy as np
import pandas as pd
from kmc_simulation import KineticMonteCarlo

def _ensure_numeric(d: Dict[str, Any]) -> Dict[str, Any]:
    out = {}
    for k, v in d.items():
        if isinstance(v, str):
            try:
                if 'e' in v or '.' in v:
                    out[k] = float(v)
                elif v.isdigit():
                    out[k] = int(v)
                else:
                    out[k] = float(v)
            except ValueError:
                out[k] = v
        else:
            out[k] = v
    return out

def _aggregate_runs(rows: List[Dict[str, Any]], keys: List[str]) -> Dict[str, Any]:
    n = len(rows)
    agg = {}
    for k in keys:
        arr = np.array([float(r.get(k, 0.0)) for r in rows], dtype=float)
        mean = float(np.mean(arr))
        std = float(np.std(arr, ddof=1)) if n > 1 else 0.0
        ci95 = float(1.96 * std / math.sqrt(n)) if n > 1 else 0.0
        agg[f"{k}_mean"] = mean
        agg[f"{k}_std"] = std
        agg[f"{k}_ci95"] = ci95
    return agg

def _mrn_weights(min_um: float, max_um: float, bins: int) -> Tuple[np.ndarray, np.ndarray]:
    a = np.logspace(np.log10(min_um), np.log10(max_um), bins)
    w = a**(-3.5)
    w = w / np.sum(w)
    return a, w

def _weighted_aggregate_across_sizes(size_blocks: List[Dict[str, Any]], metrics: List[str]) -> Dict[str, Any]:
    weights = np.array([blk["_mrn_weight"] for blk in size_blocks], dtype=float)
    weights = weights / np.sum(weights)
    out = {}
    for m in metrics:
        means = np.array([blk[f"{m}_mean"] for blk in size_blocks], dtype=float)
        stds = np.array([blk[f"{m}_std"] for blk in size_blocks], dtype=float)
        wmean = float(np.sum(weights * means))
        between_var = float(np.sum(weights * (means - wmean)**2))
        within_var = float(np.sum((weights**2) * (stds**2)))
        wstd = float(np.sqrt(max(0.0, between_var + within_var)))
        ci95 = float(1.96 * wstd)
        out[f"{m}_mean"] = wmean
        out[f"{m}_std"] = wstd
        out[f"{m}_ci95"] = ci95
    return out

def run_sweep(config_file="config.yaml"):
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    output_filename = config.get("output_filename", "results/parameter_sweep_results.csv")
    raw_runs_enabled = bool(config.get("save_raw_runs", True))
    raw_runs_filename = config.get("raw_runs_output", "results/raw_runs.csv")

    mrn_output_filename = config.get("mrn_output_filename", "results/science_campaign_v1_mrn.csv")

    os.makedirs(os.path.dirname(output_filename), exist_ok=True)
    if raw_runs_enabled:
        os.makedirs(os.path.dirname(raw_runs_filename), exist_ok=True)
    os.makedirs(os.path.dirname(mrn_output_filename), exist_ok=True)

    base_params = {
        "grain_radius_um": config.get("grain_radius_um", 0.1),
        "site_area_angstroms_sq": config.get("site_area_angstroms_sq", 9),
        "surface_temperature_k": config.get("surface_temperature_k", 10.0),
        "gas_temperature_k": config.get("gas_temperature_k", 100.0),
        "h_gas_density_cm3": config.get("h_gas_density_cm3", 1e2),
        "sticking_probability": config.get("sticking_probability", 0.3),
        "initial_h_coverage": config.get("initial_h_coverage", 0.0),
        "uv_flux_factor": config.get("uv_flux_factor", 1.0),
        "uv_stimulated_diffusion_factor": config.get("uv_stimulated_diffusion_factor", 2.0),
        "use_3d_lattice": config.get("use_3d_lattice", True),
        "porosity_fraction": config.get("porosity_fraction", 0.2),
        "E_phys_mean_meV": config.get("E_phys_mean_meV", 45.0),
        "heterogeneity_E_bind_sigma_meV": config.get("heterogeneity_E_bind_sigma_meV", 5.0),
        "uv_pulse_enabled": config.get("uv_pulse_enabled", True),
        "uv_defect_creation_rate": config.get("uv_defect_creation_rate", 0.5),
        "uv_pulse_duration": config.get("uv_pulse_duration", 1e-6),
    }
    base_params = {k: v for k, v in base_params.items() if v is not None}

    max_time_s = float(config.get("max_time_s", 3.154e7))
    ensemble_runs = int(config.get("ensemble_runs", 1))
    if ensemble_runs < 20:
        ensemble_runs = 20

    parameter_sweeps = config.get("parameter_sweeps", {})

    use_mrn = bool(config.get("use_mrn", False))
    aggregate_across_sizes = bool(config.get("aggregate_across_sizes", True))
    mrn_min_um = float(config.get("mrn_min_um", 0.005))
    mrn_max_um = float(config.get("mrn_max_um", 0.25))
    mrn_bins = int(config.get("mrn_bins", 20))

    aggregated_rows = []
    raw_rows = []

    def run_one_condition(sim_params: Dict[str, Any], run_id: int) -> Dict[str, Any]:
        kmc_sim = KineticMonteCarlo(sim_params)
        kmc_sim.run_gillespie(max_time=max_time_s)
        return {
            **sim_params,
            "run_id": run_id,
            "final_time": kmc_sim.time,
            "final_h_atoms_on_surface": kmc_sim.h_atoms_on_surface,
            "h2_formed_LH": kmc_sim.h2_molecules_formed_LH,
            "h2_formed_ER": kmc_sim.h2_molecules_formed_ER,
            "h2_formed_UV": kmc_sim.h2_molecules_formed_UV,
            "total_h2_formed": kmc_sim.h2_molecules_formed
        }

    metrics = ["h2_formed_LH", "h2_formed_ER", "h2_formed_UV", "total_h2_formed"]

    if not parameter_sweeps:
        if not use_mrn:
            rows = []
            for run_id in range(ensemble_runs):
                sim_params = _ensure_numeric(base_params.copy())
                result = run_one_condition(sim_params, run_id)
                rows.append(result)
                if raw_runs_enabled:
                    raw_rows.append(result)
            agg = _aggregate_runs(rows, metrics)
            descriptor = {k: v for k, v in rows[0].items() if k not in ["run_id", "final_time",
                                                                        "final_h_atoms_on_surface",
                                                                        "h2_formed_LH", "h2_formed_ER",
                                                                        "h2_formed_UV", "total_h2_formed"]}
            aggregated_rows.append({**descriptor, **agg})
        else:
            sizes_um, weights = _mrn_weights(mrn_min_um, mrn_max_um, mrn_bins)
            size_blocks = []
            for a_um, w in zip(sizes_um, weights):
                rows = []
                for run_id in range(ensemble_runs):
                    sim_params = _ensure_numeric({**base_params, "grain_radius_um": float(a_um)})
                    result = run_one_condition(sim_params, run_id)
                    rows.append(result)
                    if raw_runs_enabled:
                        raw_rows.append({**result, "grain_radius_um_mrn": float(a_um), "_mrn_weight": float(w)})
                agg = _aggregate_runs(rows, metrics)
                block_descriptor = {k: v for k, v in rows[0].items() if k not in ["run_id", "final_time",
                                                                                  "final_h_atoms_on_surface",
                                                                                  "h2_formed_LH", "h2_formed_ER",
                                                                                  "h2_formed_UV", "total_h2_formed"]}
                block = {**block_descriptor, **agg, "grain_radius_um_mrn": float(a_um), "_mrn_weight": float(w)}
                size_blocks.append(block)
            if aggregate_across_sizes:
                merged = _weighted_aggregate_across_sizes(size_blocks, metrics)
                descriptor = {k: v for k, v in base_params.items() if k not in ["grain_radius_um"]}
                aggregated_rows.append({**descriptor, "mrn_min_um": mrn_min_um, "mrn_max_um": mrn_max_um,
                                        "mrn_bins": mrn_bins, "mrn_aggregated": True, **merged})
            else:
                aggregated_rows.extend(size_blocks)
    else:
        sweep_keys = list(parameter_sweeps.keys())
        sweep_values = [parameter_sweeps[key] for key in sweep_keys]
        for i, combo in enumerate(product(*sweep_values)):
            current_sweep_params = dict(zip(sweep_keys, combo))
            if not use_mrn:
                rows = []
                for run_id in range(ensemble_runs):
                    sim_params = base_params.copy()
                    sim_params.update(current_sweep_params)
                    sim_params = _ensure_numeric(sim_params)
                    result = run_one_condition(sim_params, run_id)
                    rows.append(result)
                    if raw_runs_enabled:
                        raw_rows.append(result)
                agg = _aggregate_runs(rows, metrics)
                descriptor = {k: v for k, v in rows[0].items() if k not in ["run_id", "final_time",
                                                                            "final_h_atoms_on_surface",
                                                                            "h2_formed_LH", "h2_formed_ER",
                                                                            "h2_formed_UV", "total_h2_formed"]}
                aggregated_rows.append({**descriptor, **agg})
            else:
                sizes_um, weights = _mrn_weights(mrn_min_um, mrn_max_um, mrn_bins)
                size_blocks = []
                for a_um, w in zip(sizes_um, weights):
                    rows = []
                    for run_id in range(ensemble_runs):
                        sim_params = base_params.copy()
                        sim_params.update(current_sweep_params)
                        sim_params["grain_radius_um"] = float(a_um)
                        sim_params = _ensure_numeric(sim_params)
                        result = run_one_condition(sim_params, run_id)
                        rows.append(result)
                        if raw_runs_enabled:
                            raw_rows.append({**result, "grain_radius_um_mrn": float(a_um), "_mrn_weight": float(w)})
                    agg = _aggregate_runs(rows, metrics)
                    block_descriptor = {k: v for k, v in rows[0].items() if k not in ["run_id", "final_time",
                                                                                      "final_h_atoms_on_surface",
                                                                                      "h2_formed_LH", "h2_formed_ER",
                                                                                      "h2_formed_UV", "total_h2_formed"]}
                    block = {**block_descriptor, **agg, "grain_radius_um_mrn": float(a_um), "_mrn_weight": float(w)}
                    size_blocks.append(block)
                if aggregate_across_sizes:
                    merged = _weighted_aggregate_across_sizes(size_blocks, metrics)
                    descriptor = {k: v for k, v in size_blocks[0].items() if k not in ["run_id", "final_time",
                                                                                       "final_h_atoms_on_surface",
                                                                                       "h2_formed_LH", "h2_formed_ER",
                                                                                       "h2_formed_UV", "total_h2_formed",
                                                                                       "grain_radius_um_mrn", "_mrn_weight"]}
                    aggregated_rows.append({**descriptor, "mrn_min_um": mrn_min_um, "mrn_max_um": mrn_max_um,
                                            "mrn_bins": mrn_bins, "mrn_aggregated": True, **merged})
                else:
                    aggregated_rows.extend(size_blocks)

    df_agg = pd.DataFrame(aggregated_rows)
    target_output = mrn_output_filename if use_mrn else output_filename
    df_agg.to_csv(target_output, index=False)
    print(f"Aggregated simulation results saved to {target_output}")

    if raw_runs_enabled and raw_rows:
        df_raw = pd.DataFrame(raw_rows)
        df_raw.to_csv(raw_runs_filename, index=False)
        print(f"Raw per-run results saved to {raw_runs_filename}")

if __name__ == "__main__":
    run_sweep()
