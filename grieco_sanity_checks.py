import argparse
from dataclasses import dataclass
from typing import Dict, Any, List

from kmc_simulation import KineticMonteCarlo


@dataclass(frozen=True)
class SanityResult:
    name: str
    impinging_atoms: int
    impinging_h2: int
    h2_formed: int
    h2_desorbed_formed: int
    h2_desorbed_beam: int
    h2_released_formed: int
    h2_released_beam: int
    h2_on_surface: int
    epsilon: float


def _run_case(name: str, params: Dict[str, Any], max_steps: int | None) -> SanityResult:
    kmc = KineticMonteCarlo(params)
    kmc.run_gillespie(max_time=1e30, max_steps=max_steps)
    denom = float(max(getattr(kmc, "total_impinging_h_atoms", 0), 0))
    eps = 0.0
    if denom > 0:
        eps = float(2.0 * getattr(kmc, "h2_molecules_desorbed", 0) / denom)
    return SanityResult(
        name=name,
        impinging_atoms=int(getattr(kmc, "total_impinging_h_atoms", 0)),
        impinging_h2=int(getattr(kmc, "total_impinging_h2_molecules", 0)),
        h2_formed=int(getattr(kmc, "h2_molecules_formed", 0)),
        h2_desorbed_formed=int(getattr(kmc, "h2_molecules_desorbed", 0)),
        h2_desorbed_beam=int(getattr(kmc, "h2_molecules_desorbed_beam", 0)),
        h2_released_formed=int(getattr(kmc, "h2_molecules_released_formed", 0)),
        h2_released_beam=int(getattr(kmc, "h2_molecules_released_beam", 0)),
        h2_on_surface=int(getattr(kmc, "h2_molecules_on_surface", 0)),
        epsilon=float(eps),
    )


def run_sanity(max_steps: int | None = 600000) -> List[SanityResult]:
    # Small grain for fast checks.
    base: Dict[str, Any] = {
        "rng_seed": 1234,
        "grain_radius_um": 0.005,
        "site_area_angstroms_sq": 25,
        "use_3d_lattice": True,
        "porosity_fraction": 0.0,
        "surface_defect_fraction": 0.15,
        "chemisorption_fraction": 0.5,
        "E_phys_mean_meV": 45.0,
        "heterogeneity_E_bind_sigma_meV": 5.0,
        "E_chem_mean_eV": 1.75,
        "heterogeneity_E_chem_sigma_eV": 0.25,
        "uv_flux_factor": 0.0,
        "uv_pulse_enabled": False,
        # Arrival-mode (per-site arrival rate: ~1 per 100 s)
        "arrival_rate_per_site_s": 0.01,
        "sticking_probability": 0.5,
        "sticking_temp_model": "constant",
        "max_arrivals": 2000,
        # ER/abstraction params (tunable)
        "er_cross_section_cm2": 1e-15,
        "er_reaction_probability": 0.5,
        # Blocking params
        "enable_h2_blocking": True,
        "E_h2_bind_eV": 0.03,
        "h2_stick_transition_K": 20.0,
        "h2_stick_prob_lowT": 0.9,
        "h2_beam_stick_probability": 1.0,
    }

    results: List[SanityResult] = []

    # 1) Dissociation off: no impinging atoms => epsilon must be ~0 even if lots of beam H2 cycles.
    p1 = dict(base)
    p1.update(
        {
            "surface_temperature_k": 30.0,
            "beam_dissociation_fraction": 0.0,
            "E_h2_bind_eV": 0.01,  # make beam H2 desorb readily
        }
    )
    results.append(_run_case("dissociation_off", p1, max_steps))

    # 2) Strong H2 binding: lots of formation can happen, but detection should drop (H2 remains stuck).
    p2 = dict(base)
    p2.update(
        {
            "surface_temperature_k": 10.0,
            "beam_dissociation_fraction": 1.0,
            "enable_LH": True,
            "enable_diffusion": False,
            "h2_stick_transition_K": 200.0,
            "h2_stick_prob_lowT": 1.0,
            "E_h2_bind_eV": 1.0,  # very strong binding => negligible desorption
            "max_arrivals": 500,
        }
    )
    results.append(_run_case("strong_h2_binding", p2, max_steps))

    # 3) No surface chemistry: atoms impinge but never stick/recombine; beam H2 can still desorb as baseline.
    p3 = dict(base)
    p3.update(
        {
            "surface_temperature_k": 30.0,
            "beam_dissociation_fraction": 0.5,
            "sticking_probability": 0.0,  # no adsorption => no LH reservoir, no ER via chem reservoir
            "E_h2_bind_eV": 0.01,
            "max_arrivals": 3000,
        }
    )
    results.append(_run_case("no_surface_chemistry", p3, max_steps))

    return results


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Quick sanity checks for Grieco-style epsilon accounting.")
    p.add_argument("--max-steps", type=int, default=600000, help="Hard cap on KMC steps per case")
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    rows = run_sanity(max_steps=int(args.max_steps) if args.max_steps else None)
    for r in rows:
        print(
            f"{r.name}: eps={r.epsilon:.6g} "
            f"imp_atoms={r.impinging_atoms} imp_h2={r.impinging_h2} "
            f"h2_formed={r.h2_formed} h2_des_formed={r.h2_desorbed_formed} "
            f"h2_des_beam={r.h2_desorbed_beam} "
            f"h2_rel_formed={r.h2_released_formed} h2_rel_beam={r.h2_released_beam} "
            f"h2_surf={r.h2_on_surface}"
        )
