import os
import sys
import json
import yaml
import math
import csv
from itertools import product
from typing import Dict, Any, List, Tuple
import numpy as np
try:
    import pandas as pd
except ImportError:  # pragma: no cover
    pd = None
from kmc_simulation import KineticMonteCarlo
from scientific_data import K_B_ERG, M_H

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

def _write_csv(path: str, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        return
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    # Use a stable union of keys so missing keys become blank cells.
    fieldnames: List[str] = sorted({k for row in rows for k in row.keys()})
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

def _apply_arrival_rate_mode(sim_params: Dict[str, Any], config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Populate explicit arrival knobs from higher-level sweep intent.

    In particular, `arrival_rate_mode: gas_kinetic` should compute an
    `arrival_rate_per_site_s` from gas density, thermal velocity, and site area.
    This keeps the KMC on the arrival/impingement path used by the ISM campaigns,
    rather than silently falling back to the adsorption-only path.
    """
    arm = str(
        sim_params.get("arrival_rate_mode",
                       config.get("arrival_rate_mode", ""))
    ).strip().lower()

    if arm != "gas_kinetic":
        return sim_params

    existing = sim_params.get("arrival_rate_per_site_s", None)
    try:
        if existing is not None and float(existing) > 0.0:
            return sim_params
    except (TypeError, ValueError):
        pass

    n_h = float(sim_params.get("h_gas_density_cm3", 0.0) or 0.0)
    t_gas = float(sim_params.get("gas_temperature_k", 100.0) or 100.0)
    site_area_ang2 = float(sim_params.get("site_area_angstroms_sq", 25.0) or 25.0)
    site_area_cm2 = site_area_ang2 * 1e-16

    if n_h <= 0.0 or t_gas <= 0.0 or site_area_cm2 <= 0.0:
        return sim_params

    v_th = math.sqrt(8.0 * float(K_B_ERG) * float(t_gas) / (math.pi * float(M_H)))
    flux_cm2_s = 0.25 * float(n_h) * float(v_th)
    sim_params["arrival_rate_per_site_s"] = float(flux_cm2_s) * float(site_area_cm2)
    sim_params["arrival_rate_mode"] = "gas_kinetic"
    return sim_params

def run_sweep(config_file="config.yaml"):
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    mode = str(config.get("mode", "default")).strip().lower()
    grieco_mode = mode in {"grieco", "grieco_validation", "experiment", "verification"}

    output_filename = config.get("output_filename", "results/parameter_sweep_results.csv")
    raw_runs_enabled = bool(config.get("save_raw_runs", True))
    raw_runs_filename = config.get("raw_runs_output", "results/raw_runs.csv")

    mrn_output_filename = config.get("mrn_output_filename", "results/science_campaign_v1_mrn.csv")

    os.makedirs(os.path.dirname(output_filename), exist_ok=True)
    if raw_runs_enabled:
        os.makedirs(os.path.dirname(raw_runs_filename), exist_ok=True)
    os.makedirs(os.path.dirname(mrn_output_filename), exist_ok=True)

    base_params = {
        "rng_seed": config.get("rng_seed", None),
        "grain_radius_um": config.get("grain_radius_um", 0.005 if grieco_mode else 0.1),
        "site_area_angstroms_sq": config.get("site_area_angstroms_sq", 25 if grieco_mode else 9),
        "surface_temperature_k": config.get("surface_temperature_k", 250.0 if grieco_mode else 10.0),
        "gas_temperature_k": config.get("gas_temperature_k", 300.0 if grieco_mode else 100.0),
        "h_gas_density_cm3": config.get("h_gas_density_cm3", 0.0 if grieco_mode else 1e2),
        "sticking_probability": config.get("sticking_probability", 0.5 if grieco_mode else 0.3),
        "sticking_temp_model": config.get("sticking_temp_model", "constant" if grieco_mode else None),
        "initial_h_count": config.get("initial_h_count", None),
        "initial_h_coverage": config.get("initial_h_coverage", 0.0),
        "initial_h_chemisorption_only": config.get("initial_h_chemisorption_only", None),
        "uv_flux_factor": config.get("uv_flux_factor", 0.0 if grieco_mode else 1.0),
        "uv_mode": config.get("uv_mode", "pulse"),
        "uv_h2_mode": config.get("uv_h2_mode", None),
        "uv_photofrag_min_chemisorbed_h": config.get("uv_photofrag_min_chemisorbed_h", None),
        "uv_photofrag_h_per_event": config.get("uv_photofrag_h_per_event", None),
        "uv_photofrag_cross_section_cm2": config.get("uv_photofrag_cross_section_cm2", None),
        "uv_photofrag_branching_ratio": config.get("uv_photofrag_branching_ratio", None),
        "uv_stimulated_diffusion_factor": config.get("uv_stimulated_diffusion_factor", 2.0),
        "enable_LH": config.get("enable_LH", False if grieco_mode else None),
        "enable_diffusion": config.get("enable_diffusion", False if grieco_mode else None),
        "use_3d_lattice": config.get("use_3d_lattice", True),
        "porosity_fraction": config.get("porosity_fraction", 0.0 if grieco_mode else 0.2),
        "chemisorption_fraction": config.get("chemisorption_fraction", 0.5 if grieco_mode else None),
        "surface_defect_fraction": config.get("surface_defect_fraction", 0.15 if grieco_mode else None),
        "E_phys_mean_meV": config.get("E_phys_mean_meV", 45.0),
        "heterogeneity_E_bind_sigma_meV": config.get("heterogeneity_E_bind_sigma_meV", 5.0),
        "E_chem_mean_eV": config.get("E_chem_mean_eV", 1.75 if grieco_mode else None),
        "heterogeneity_E_chem_sigma_eV": config.get("heterogeneity_E_chem_sigma_eV", 0.25 if grieco_mode else None),
        "E_defect_mean_eV": config.get("E_defect_mean_eV", None),
        "heterogeneity_E_defect_sigma_eV": config.get("heterogeneity_E_defect_sigma_eV", None),
        "diffusion_to_binding_ratio_physisorption": config.get("diffusion_to_binding_ratio_physisorption", None),
        "diffusion_to_binding_ratio_chemisorption": config.get("diffusion_to_binding_ratio_chemisorption", None),
        "diffusion_to_binding_ratio_defect": config.get("diffusion_to_binding_ratio_defect", None),
        "uv_pulse_enabled": config.get("uv_pulse_enabled", False if grieco_mode else True),
        "uv_defect_creation_rate": config.get("uv_defect_creation_rate", 0.5),
        "uv_pulse_duration": config.get("uv_pulse_duration", 1e-6),
        "uv_pulse_start_rate_s": config.get("uv_pulse_start_rate_s", None),
        # Experiment/beam-style arrivals
        "arrival_rate_s": config.get("arrival_rate_s", None),
        "arrival_rate_per_site_s": config.get("arrival_rate_per_site_s", 0.01 if grieco_mode else None),
        "max_arrivals": config.get("max_arrivals", None),
        "er_cross_section_cm2": config.get("er_cross_section_cm2", 1e-15 if grieco_mode else None),
        "er_reaction_probability": config.get("er_reaction_probability", 0.5 if grieco_mode else None),
        "diffusion_mode": config.get("diffusion_mode", None),
        "diffusion_rate_cap_s": config.get("diffusion_rate_cap_s", None),
        "lh_formation_mode": config.get("lh_formation_mode", None),
        "lh_exclude_chemisorption": config.get("lh_exclude_chemisorption", None),
        "enable_chemisorption_diffusion": config.get("enable_chemisorption_diffusion", None),
        "lh_diffusion_factor": config.get("lh_diffusion_factor", None),
        # Low-T protocol pieces (optional)
        "enable_h2_blocking": config.get("enable_h2_blocking", None),
        "E_h2_bind_eV": config.get("E_h2_bind_eV", None),
        "h2_desorption_prefactor_s": config.get("h2_desorption_prefactor_s", None),
        "h2_stick_transition_K": config.get("h2_stick_transition_K", None),
        "h2_stick_prob_lowT": config.get("h2_stick_prob_lowT", None),
        "sticking_blocking_strength": config.get("sticking_blocking_strength", None),
        "er_blocking_strength": config.get("er_blocking_strength", None),
        "beam_dissociation_fraction": config.get("beam_dissociation_fraction", None),
        "h2_beam_stick_probability": config.get("h2_beam_stick_probability", None),
        "temp_ramp": config.get("temp_ramp", None),
        # Optional grain structure caching
        "enable_grain_cache": config.get("enable_grain_cache", None),
        "grain_cache_dir": config.get("grain_cache_dir", None),
        "grain_cache_include_rng_seed": config.get("grain_cache_include_rng_seed", None),
    }
    base_params = {k: v for k, v in base_params.items() if v is not None}

    max_time_s = float(config.get("max_time_s", 3.154e7))
    max_steps = config.get("max_steps", None)
    max_steps = int(max_steps) if isinstance(max_steps, (int, float, str)) and str(max_steps).strip() != "" else None

    burnin_arrivals = int(config.get("burnin_arrivals", 2000 if grieco_mode else 0))
    measure_arrivals = config.get("measure_arrivals", 5000 if grieco_mode else None)
    measure_arrivals = int(measure_arrivals) if isinstance(measure_arrivals, (int, float, str)) and str(measure_arrivals).strip() != "" else None

    ensemble_runs = int(config.get("ensemble_runs", 1))
    min_ensemble_default = 1 if grieco_mode else 20
    min_ensemble_runs = int(config.get("min_ensemble_runs", min_ensemble_default))
    if ensemble_runs < min_ensemble_runs:
        ensemble_runs = min_ensemble_runs

    parameter_sweeps = config.get("parameter_sweeps", {})
    explicit_conditions = config.get("explicit_conditions", None)
    if isinstance(explicit_conditions, list) and len(explicit_conditions) > 0:
        # explicit_conditions takes precedence over cartesian parameter_sweeps.
        parameter_sweeps = {}
    else:
        explicit_conditions = None

    use_mrn = bool(config.get("use_mrn", False))
    aggregate_across_sizes = bool(config.get("aggregate_across_sizes", True))
    mrn_min_um = float(config.get("mrn_min_um", 0.005))
    mrn_max_um = float(config.get("mrn_max_um", 0.25))
    mrn_bins = int(config.get("mrn_bins", 20))

    aggregated_rows = []
    raw_rows = []

    def run_one_condition(sim_params: Dict[str, Any], run_id: int) -> Dict[str, Any]:
        if isinstance(sim_params.get("rng_seed", None), (int, float, str)) and str(sim_params.get("rng_seed")).strip() != "":
            try:
                sim_params["rng_seed"] = int(float(sim_params["rng_seed"])) + int(run_id)
            except (TypeError, ValueError):
                pass
        kmc_sim = KineticMonteCarlo(sim_params)
        arrival_mode = (
            (isinstance(sim_params.get("arrival_rate_s", None), (int, float)) and float(sim_params["arrival_rate_s"]) > 0)
            or (isinstance(sim_params.get("arrival_rate_per_site_s", None), (int, float)) and float(sim_params["arrival_rate_per_site_s"]) > 0)
        )

        if arrival_mode and (burnin_arrivals > 0 or measure_arrivals is not None):
            temp_ramp_cfg = kmc_sim.simulation_parameters.get("temp_ramp", None)
            has_enabled_ramp = isinstance(temp_ramp_cfg, dict) and bool(temp_ramp_cfg.get("enabled", False))
            if burnin_arrivals > 0:
                # Burn-in should be isothermal; do not advance a temperature ramp during burn-in.
                if has_enabled_ramp:
                    kmc_sim.simulation_parameters["temp_ramp"] = {**temp_ramp_cfg, "enabled": False}
                kmc_sim.simulation_parameters["max_arrivals"] = int(burnin_arrivals)
                kmc_sim.run_gillespie(max_time=1e30, max_steps=max_steps)
                if has_enabled_ramp:
                    kmc_sim.simulation_parameters["temp_ramp"] = temp_ramp_cfg

            # Reset counters for the measurement window while keeping surface state.
            kmc_sim.total_impinging_h_atoms = 0
            kmc_sim.total_impinging_h2_molecules = 0
            kmc_sim.total_adsorbed_h_atoms = 0
            kmc_sim.total_desorbed_h_atoms = 0
            kmc_sim.h2_molecules_formed = 0
            kmc_sim.h2_molecules_desorbed = 0
            kmc_sim.h2_molecules_desorbed_LH = 0
            kmc_sim.h2_molecules_desorbed_ER = 0
            kmc_sim.h2_molecules_desorbed_UV = 0
            kmc_sim.h2_molecules_desorbed_beam = 0
            kmc_sim.h2_molecules_released_formed = 0
            kmc_sim.h2_molecules_released_beam = 0
            kmc_sim.h2_molecules_formed_LH = 0
            kmc_sim.h2_molecules_formed_ER = 0
            kmc_sim.h2_molecules_formed_UV = 0
            if isinstance(sim_params.get("temp_ramp", None), dict) and bool(sim_params["temp_ramp"].get("enabled", False)):
                kmc_sim.time = 0.0

            if measure_arrivals is not None:
                kmc_sim.simulation_parameters["max_arrivals"] = int(measure_arrivals)
            kmc_sim.run_gillespie(max_time=1e30, max_steps=max_steps)
        else:
            kmc_sim.run_gillespie(max_time=max_time_s, max_steps=max_steps)

        epsilon = 0.0
        if getattr(kmc_sim, "total_impinging_h_atoms", 0) > 0:
            epsilon = float(2.0 * kmc_sim.h2_molecules_desorbed / kmc_sim.total_impinging_h_atoms)

        surface_area_cm2 = (
            float(sim_params.get("site_area_angstroms_sq", 25.0) or 25.0)
            * 1e-16
            * float(getattr(kmc_sim, "total_accessible_surface_sites", 0) or 0)
        )
        h2_release_rate_cm2_s = 0.0
        if float(kmc_sim.time) > 0.0 and surface_area_cm2 > 0.0:
            h2_release_rate_cm2_s = float(kmc_sim.h2_molecules_desorbed) / (float(kmc_sim.time) * float(surface_area_cm2))
        final_h_atoms_on_surface = int(kmc_sim.h_atoms_on_surface)
        final_h_surface_coverage = 0.0
        if float(getattr(kmc_sim, "total_accessible_surface_sites", 0) or 0) > 0.0:
            final_h_surface_coverage = float(final_h_atoms_on_surface) / float(kmc_sim.total_accessible_surface_sites)
        return {
            **sim_params,
            "run_id": run_id,
            "final_time": kmc_sim.time,
            "final_h_atoms_on_surface": final_h_atoms_on_surface,
            "final_h_surface_coverage": final_h_surface_coverage,
            "h2_formed_LH": kmc_sim.h2_molecules_formed_LH,
            "h2_formed_ER": kmc_sim.h2_molecules_formed_ER,
            "h2_formed_UV": kmc_sim.h2_molecules_formed_UV,
            "total_h2_formed": kmc_sim.h2_molecules_formed,
            "total_h2_desorbed": kmc_sim.h2_molecules_desorbed,
            "total_h2_desorbed_LH": getattr(kmc_sim, "h2_molecules_desorbed_LH", 0),
            "total_h2_desorbed_ER": getattr(kmc_sim, "h2_molecules_desorbed_ER", 0),
            "total_h2_desorbed_UV": getattr(kmc_sim, "h2_molecules_desorbed_UV", 0),
            "total_h2_desorbed_beam": getattr(kmc_sim, "h2_molecules_desorbed_beam", 0),
            "total_h2_released_formed": getattr(kmc_sim, "h2_molecules_released_formed", 0),
            "total_h2_released_beam": getattr(kmc_sim, "h2_molecules_released_beam", 0),
            "h2_on_surface": getattr(kmc_sim, "h2_molecules_on_surface", 0),
            "impinging_h_atoms": getattr(kmc_sim, "total_impinging_h_atoms", 0),
            "impinging_h2_molecules": getattr(kmc_sim, "total_impinging_h2_molecules", 0),
            "epsilon": epsilon,
            "h2_release_rate_cm2_s": h2_release_rate_cm2_s,
        }

    metrics = [
        "final_time",
        "final_h_atoms_on_surface",
        "final_h_surface_coverage",
        "h2_formed_LH",
        "h2_formed_ER",
        "h2_formed_UV",
        "total_h2_formed",
        "total_h2_desorbed",
        "total_h2_desorbed_LH",
        "total_h2_desorbed_ER",
        "total_h2_desorbed_UV",
        "total_h2_desorbed_beam",
        "total_h2_released_formed",
        "total_h2_released_beam",
        "h2_on_surface",
        "epsilon",
        "impinging_h_atoms",
        "impinging_h2_molecules",
        "h2_release_rate_cm2_s",
    ]
    exclude_descriptor_keys = ["run_id", *metrics]

    if explicit_conditions:
        for i, condition in enumerate(explicit_conditions):
            current_sweep_params = _ensure_numeric(dict(condition))
            if not use_mrn:
                rows = []
                for run_id in range(ensemble_runs):
                    sim_params = base_params.copy()
                    sim_params.update(current_sweep_params)
                    sim_params = _ensure_numeric(sim_params)
                    sim_params = _apply_arrival_rate_mode(sim_params, config)
                    result = run_one_condition(sim_params, run_id)
                    rows.append(result)
                    if raw_runs_enabled:
                        raw_rows.append(result)
                agg = _aggregate_runs(rows, metrics)
                descriptor = {k: v for k, v in rows[0].items() if k not in exclude_descriptor_keys}
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
                        sim_params = _apply_arrival_rate_mode(sim_params, config)
                        result = run_one_condition(sim_params, run_id)
                        rows.append(result)
                        if raw_runs_enabled:
                            raw_rows.append({**result, "grain_radius_um_mrn": float(a_um), "_mrn_weight": float(w)})
                    agg = _aggregate_runs(rows, metrics)
                    block_descriptor = {k: v for k, v in rows[0].items() if k not in exclude_descriptor_keys}
                    block = {**block_descriptor, **agg, "grain_radius_um_mrn": float(a_um), "_mrn_weight": float(w)}
                    size_blocks.append(block)
                if aggregate_across_sizes:
                    merged = _weighted_aggregate_across_sizes(size_blocks, metrics)
                    mrn_exclude = [*exclude_descriptor_keys, "grain_radius_um_mrn", "_mrn_weight"]
                    descriptor = {k: v for k, v in size_blocks[0].items() if k not in mrn_exclude}
                    aggregated_rows.append({**descriptor, "mrn_min_um": mrn_min_um, "mrn_max_um": mrn_max_um,
                                            "mrn_bins": mrn_bins, "mrn_aggregated": True, **merged})
                else:
                    aggregated_rows.extend(size_blocks)
    elif not parameter_sweeps:
        if not use_mrn:
            rows = []
            for run_id in range(ensemble_runs):
                sim_params = _ensure_numeric(base_params.copy())
                sim_params = _apply_arrival_rate_mode(sim_params, config)
                result = run_one_condition(sim_params, run_id)
                rows.append(result)
                if raw_runs_enabled:
                    raw_rows.append(result)
            agg = _aggregate_runs(rows, metrics)
            descriptor = {k: v for k, v in rows[0].items() if k not in exclude_descriptor_keys}
            aggregated_rows.append({**descriptor, **agg})
        else:
            sizes_um, weights = _mrn_weights(mrn_min_um, mrn_max_um, mrn_bins)
            size_blocks = []
            for a_um, w in zip(sizes_um, weights):
                rows = []
                for run_id in range(ensemble_runs):
                    sim_params = _ensure_numeric({**base_params, "grain_radius_um": float(a_um)})
                    sim_params = _apply_arrival_rate_mode(sim_params, config)
                    result = run_one_condition(sim_params, run_id)
                    rows.append(result)
                    if raw_runs_enabled:
                        raw_rows.append({**result, "grain_radius_um_mrn": float(a_um), "_mrn_weight": float(w)})
                agg = _aggregate_runs(rows, metrics)
                block_descriptor = {k: v for k, v in rows[0].items() if k not in exclude_descriptor_keys}
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
                    sim_params = _apply_arrival_rate_mode(sim_params, config)
                    result = run_one_condition(sim_params, run_id)
                    rows.append(result)
                    if raw_runs_enabled:
                        raw_rows.append(result)
                agg = _aggregate_runs(rows, metrics)
                descriptor = {k: v for k, v in rows[0].items() if k not in exclude_descriptor_keys}
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
                        sim_params = _apply_arrival_rate_mode(sim_params, config)
                        result = run_one_condition(sim_params, run_id)
                        rows.append(result)
                        if raw_runs_enabled:
                            raw_rows.append({**result, "grain_radius_um_mrn": float(a_um), "_mrn_weight": float(w)})
                    agg = _aggregate_runs(rows, metrics)
                    block_descriptor = {k: v for k, v in rows[0].items() if k not in exclude_descriptor_keys}
                    block = {**block_descriptor, **agg, "grain_radius_um_mrn": float(a_um), "_mrn_weight": float(w)}
                    size_blocks.append(block)
                if aggregate_across_sizes:
                    merged = _weighted_aggregate_across_sizes(size_blocks, metrics)
                    mrn_exclude = [*exclude_descriptor_keys, "grain_radius_um_mrn", "_mrn_weight"]
                    descriptor = {k: v for k, v in size_blocks[0].items() if k not in mrn_exclude}
                    aggregated_rows.append({**descriptor, "mrn_min_um": mrn_min_um, "mrn_max_um": mrn_max_um,
                                            "mrn_bins": mrn_bins, "mrn_aggregated": True, **merged})
                else:
                    aggregated_rows.extend(size_blocks)

    target_output = mrn_output_filename if use_mrn else output_filename
    if pd is not None:
        pd.DataFrame(aggregated_rows).to_csv(target_output, index=False)
    else:
        _write_csv(target_output, aggregated_rows)
    print(f"Aggregated simulation results saved to {target_output}")

    if raw_runs_enabled and raw_rows:
        if pd is not None:
            pd.DataFrame(raw_rows).to_csv(raw_runs_filename, index=False)
        else:
            _write_csv(raw_runs_filename, raw_rows)
        print(f"Raw per-run results saved to {raw_runs_filename}")

if __name__ == "__main__":
    config_file = sys.argv[1] if len(sys.argv) > 1 else "config.yaml"
    run_sweep(config_file=config_file)
