import argparse
import csv
import os
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np

from kmc_simulation import KineticMonteCarlo


@dataclass(frozen=True)
class GriecoRunResult:
    temperature_k: float
    epsilon: float
    h2_total: int
    h2_er: int
    chemisorbed_surface_h: int


def _default_coronene_like_params() -> Dict[str, float]:
    """
    A pragmatic, coronene-like parameter set intended to reproduce the *qualitative*
    high-temperature behavior reported by Grieco et al. (Nature Astronomy 2023):
    efficient H2 formation up to ~250 K enabled by a chemisorbed H reservoir.

    This is not a faithful reproduction of the FORMOLISM setup or coronene film
    microphysics; it is a lightweight validation harness for this codebase.
    """
    return {
        # Geometry: 4–200 Å grains ~ 0.0004–0.02 µm; use upper end to get enough sites.
        "grain_radius_um": 0.005,
        "site_area_angstroms_sq": 25,
        "use_3d_lattice": True,
        "porosity_fraction": 0.0,
        # Site population: coronene films are rich in C–H termination / reactive sites in steady state.
        "chemisorption_fraction": 0.5,
        "surface_defect_fraction": 0.15,
        # Energetics
        "E_phys_mean_meV": 45.0,
        "heterogeneity_E_bind_sigma_meV": 5.0,
        "E_chem_mean_eV": 1.75,
        "heterogeneity_E_chem_sigma_eV": 0.25,
        # Dynamics / mechanisms
        "enable_LH": False,  # focus on chemisorption-driven ER/abstraction at high T
        "enable_diffusion": False,  # speeds up arrival-mode runs; not needed for ER-only plateau
        "uv_flux_factor": 0.0,
        "uv_pulse_enabled": False,
        # Beam-like arrivals: ~1 landing per adsorption site per 100 s (main-text statement).
        "arrival_rate_per_site_s": 0.01,
        # Sticking: calibrated effective accretion probability for this model/protocol harness.
        "sticking_probability": 0.3,
        "sticking_temp_model": "constant",
        # ER/abstraction: effective cross section & reaction probability (tunable)
        "er_cross_section_cm2": 1.5e-15,
        "er_reaction_probability": 0.9,
        # Gas values are unused in arrival mode, but keep defined.
        "gas_temperature_k": 300.0,
        "h_gas_density_cm3": 0.0,
    }


def _run_once(
    temperature_k: float,
    base_params: Dict[str, float],
    burnin_arrivals: int,
    measure_arrivals: int,
    seed: int,
    max_steps: int | None,
) -> GriecoRunResult:
    sim_params = dict(base_params)
    sim_params["surface_temperature_k"] = float(temperature_k)
    sim_params["rng_seed"] = int(seed)

    # Burn-in to reach a quasi steady-state chemisorbed reservoir.
    if burnin_arrivals > 0:
        sim_params["max_arrivals"] = int(burnin_arrivals)
        kmc = KineticMonteCarlo(sim_params)
        kmc.run_gillespie(max_time=1e30, max_steps=max_steps)
    else:
        kmc = KineticMonteCarlo(sim_params)

    # Reset counters for measurement window (keep surface state).
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

    kmc.simulation_parameters["max_arrivals"] = int(measure_arrivals)
    kmc.run_gillespie(max_time=1e30, max_steps=max_steps)

    eps = 0.0
    if kmc.total_impinging_h_atoms > 0:
        # ε ≡ 2·N(H2 prompt-desorbed from formation)/N(atoms impinging)
        eps = float(2.0 * kmc.h2_molecules_desorbed / kmc.total_impinging_h_atoms)

    return GriecoRunResult(
        temperature_k=float(temperature_k),
        epsilon=eps,
        h2_total=int(kmc.h2_molecules_formed),
        h2_er=int(kmc.h2_molecules_formed_ER),
        chemisorbed_surface_h=int(len(kmc.occupied_chemisorption_surface)),
    )


def run_grieco_validation(
    temperatures_k: List[float],
    output_csv: str,
    replicates: int,
    burnin_arrivals: int,
    measure_arrivals: int,
    max_steps: int | None,
    base_params: Dict[str, float] | None = None,
) -> None:
    if base_params is None:
        base_params = _default_coronene_like_params()

    rows: List[Dict[str, float]] = []
    for t in temperatures_k:
        results = [
            _run_once(
                temperature_k=float(t),
                base_params=base_params,
                burnin_arrivals=burnin_arrivals,
                measure_arrivals=measure_arrivals,
                seed=1000 + i,
                max_steps=max_steps,
            )
            for i in range(replicates)
        ]

        eps_arr = np.array([r.epsilon for r in results], dtype=float)
        h2_arr = np.array([r.h2_total for r in results], dtype=float)
        chem_arr = np.array([r.chemisorbed_surface_h for r in results], dtype=float)

        rows.append(
            {
                "surface_temperature_k": float(t),
                "epsilon_mean": float(np.mean(eps_arr)),
                "epsilon_std": float(np.std(eps_arr, ddof=1)) if replicates > 1 else 0.0,
                "epsilon_ci95": float(1.96 * float(np.std(eps_arr, ddof=1)) / float(np.sqrt(replicates))) if replicates > 1 else 0.0,
                "h2_total_mean": float(np.mean(h2_arr)),
                "chemisorbed_surface_h_mean": float(np.mean(chem_arr)),
            }
        )

    os.makedirs(os.path.dirname(output_csv) or ".", exist_ok=True)
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} temperature points to {output_csv}")


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Validate high-T H2 formation efficiency (Grieco et al. 2023) in this KMC model.")
    p.add_argument("--output", default="results/grieco_validation.csv", help="Output CSV path")
    p.add_argument("--temps", nargs="+", type=float, default=[100, 150, 200, 250], help="Surface temperatures (K)")
    p.add_argument("--replicates", type=int, default=5, help="Runs per temperature")
    p.add_argument("--burnin-arrivals", type=int, default=2000, help="Arrivals discarded for burn-in")
    p.add_argument("--measure-arrivals", type=int, default=5000, help="Arrivals used for epsilon measurement")
    p.add_argument("--max-steps", type=int, default=500000, help="Hard cap on KMC steps per phase")
    p.add_argument("--arrival-rate-per-site", type=float, default=None, help="Override arrival_rate_per_site_s")
    p.add_argument("--sticking-probability", type=float, default=None, help="Override sticking_probability")
    p.add_argument("--er-cross-section-cm2", type=float, default=None, help="Override er_cross_section_cm2")
    p.add_argument("--er-reaction-probability", type=float, default=None, help="Override er_reaction_probability")
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    base = _default_coronene_like_params()
    if args.arrival_rate_per_site is not None:
        base["arrival_rate_per_site_s"] = float(args.arrival_rate_per_site)
    if args.sticking_probability is not None:
        base["sticking_probability"] = float(args.sticking_probability)
    if args.er_cross_section_cm2 is not None:
        base["er_cross_section_cm2"] = float(args.er_cross_section_cm2)
    if args.er_reaction_probability is not None:
        base["er_reaction_probability"] = float(args.er_reaction_probability)
    run_grieco_validation(
        temperatures_k=list(args.temps),
        output_csv=str(args.output),
        replicates=int(args.replicates),
        burnin_arrivals=int(args.burnin_arrivals),
        measure_arrivals=int(args.measure_arrivals),
        max_steps=int(args.max_steps) if args.max_steps else None,
        base_params=base,
    )
