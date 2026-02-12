import argparse
import csv
import os
from typing import Any, Dict, List

import numpy as np

from calibrate_grieco import _default_ded_params, _default_highT_params, _isothermal_epsilon, evaluate_ded


def _estimate_surface_sites(grain_radius_um: float, site_area_angstroms_sq: float) -> float:
    # Surface area (cm^2) ~ 4πr^2. Convert µm -> cm.
    r_cm = float(grain_radius_um) * 1e-4
    area_cm2 = 4.0 * np.pi * (r_cm**2)
    site_area_cm2 = float(site_area_angstroms_sq) * 1e-16
    if site_area_cm2 <= 0:
        return 0.0
    return float(area_cm2 / site_area_cm2)


def _calibrated_params() -> tuple[Dict[str, Any], Dict[str, Any]]:
    highT = _default_highT_params()
    ded = _default_ded_params()
    calibrated = {
        "sticking_probability": 0.3,
        "sticking_temp_model": "constant",
        "er_cross_section_cm2": 1.5e-15,
        "er_reaction_probability": 0.9,
        "beam_dissociation_fraction": 1.0,
        "enable_h2_blocking": True,
        "E_h2_bind_eV": 0.06,
        "h2_desorption_prefactor_s": 1e12,
        "h2_stick_transition_K": 20.0,
        "h2_stick_prob_lowT": 0.9,
        "sticking_blocking_strength": 1.0,
        "er_blocking_strength": 1.0,
    }
    # Keep high-T harness in its intended regime.
    highT.update({k: v for k, v in calibrated.items() if k in highT})
    highT["enable_h2_blocking"] = False
    ded.update(calibrated)
    return highT, ded


def main() -> None:
    p = argparse.ArgumentParser(description="Lightweight convergence/regression report for the calibrated Grieco harness.")
    p.add_argument("--out", default="results/grieco_convergence.csv", help="Output CSV path")
    p.add_argument("--radii-um", nargs="+", type=float, default=[0.003, 0.005, 0.01])
    p.add_argument("--site-area-a2", type=float, default=25.0)
    p.add_argument(
        "--reference-radius-um",
        type=float,
        default=0.005,
        help="Reference radius used to interpret burn-in/measure arrivals (scaled by estimated surface sites).",
    )

    p.add_argument("--plateau-T", type=float, default=150.0)
    p.add_argument("--plateau-burnin", type=int, default=2000)
    p.add_argument("--plateau-measure", type=int, default=5000)
    p.add_argument("--plateau-replicates", type=int, default=2)

    p.add_argument("--ded-start", type=float, default=10.0)
    p.add_argument("--ded-end", type=float, default=80.0)
    p.add_argument("--ded-rate-k-per-min", type=float, default=1.0)
    p.add_argument("--ded-bin-width", type=float, default=1.0)
    p.add_argument("--ded-burnin", type=int, default=2000)
    p.add_argument("--ded-replicates", type=int, default=1)

    p.add_argument("--max-steps", type=int, default=800000)
    args = p.parse_args()

    max_steps = int(args.max_steps) if args.max_steps else None
    radii = [float(x) for x in args.radii_um]
    site_area = float(args.site_area_a2)

    base_highT, base_ded = _calibrated_params()

    # Scale burn-in/measurement arrivals by estimated surface sites so that runs are comparable across radii.
    # The provided plateau_burnin/measure values are interpreted for the reference radius, and then scaled.
    ref_sites = _estimate_surface_sites(float(args.reference_radius_um), float(site_area))
    burnin_per_site = float(args.plateau_burnin) / ref_sites if ref_sites > 0 else float(args.plateau_burnin)
    measure_per_site = float(args.plateau_measure) / ref_sites if ref_sites > 0 else float(args.plateau_measure)
    ded_burnin_per_site = float(args.ded_burnin) / ref_sites if ref_sites > 0 else float(args.ded_burnin)

    rows: List[Dict[str, Any]] = []
    for r_um in radii:
        highT = dict(base_highT)
        ded = dict(base_ded)
        highT["grain_radius_um"] = float(r_um)
        ded["grain_radius_um"] = float(r_um)
        highT["site_area_angstroms_sq"] = float(site_area)
        ded["site_area_angstroms_sq"] = float(site_area)

        sites_est = _estimate_surface_sites(float(r_um), float(site_area))
        plateau_burnin = max(1, int(round(burnin_per_site * sites_est)))
        plateau_measure = max(1, int(round(measure_per_site * sites_est)))
        ded_burnin = max(1, int(round(ded_burnin_per_site * sites_est)))

        eps_runs = [
            _isothermal_epsilon(
                params=highT,
                temperature_k=float(args.plateau_T),
                burnin_arrivals=int(plateau_burnin),
                measure_arrivals=int(plateau_measure),
                seed=7000 + i,
                max_steps=max_steps,
            )
            for i in range(int(args.plateau_replicates))
        ]
        eps_mean = float(np.mean(np.array(eps_runs, dtype=float)))
        eps_ci95 = (
            float(1.96 * float(np.std(np.array(eps_runs, dtype=float), ddof=1)) / float(np.sqrt(int(args.plateau_replicates))))
            if int(args.plateau_replicates) > 1
            else 0.0
        )

        ded_summary = evaluate_ded(
            base_params=ded,
            replicates=int(args.ded_replicates),
            burnin_arrivals=int(ded_burnin),
            t_start_k=float(args.ded_start),
            t_end_k=float(args.ded_end),
            rate_k_per_min=float(args.ded_rate_k_per_min),
            bin_width_k=float(args.ded_bin_width),
            max_steps=max_steps,
        )
        eps10 = float(ded_summary.epsilon_10k)
        eps20 = float(ded_summary.epsilon_20k)
        eps30_80 = float(ded_summary.epsilon_30_80k)
        ratio = float(eps10 / eps20) if eps20 > 0 else 0.0

        rows.append(
            {
                "grain_radius_um": float(r_um),
                "site_area_angstroms_sq": float(site_area),
                "estimated_surface_sites": sites_est,
                "reference_radius_um": float(args.reference_radius_um),
                "plateau_burnin_arrivals": int(plateau_burnin),
                "plateau_measure_arrivals": int(plateau_measure),
                "ded_burnin_arrivals": int(ded_burnin),
                "plateau_T_k": float(args.plateau_T),
                "plateau_epsilon_mean": eps_mean,
                "plateau_epsilon_ci95": eps_ci95,
                "ded_start_k": float(args.ded_start),
                "ded_end_k": float(args.ded_end),
                "ded_eps10": eps10,
                "ded_eps20": eps20,
                "ded_eps30_80": eps30_80,
                "ded_ratio10_over_20": ratio,
            }
        )

    os.makedirs(os.path.dirname(str(args.out)) or ".", exist_ok=True)
    with open(str(args.out), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows)} rows to {args.out}")


if __name__ == "__main__":
    main()
