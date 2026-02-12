import argparse
import json
import os
from typing import Any, Dict

from calibrate_grieco import _default_ded_params, _default_highT_params, _isothermal_epsilon, evaluate_ded
from kmc_simulation import KineticMonteCarlo


def _calibrated_highT_params() -> Dict[str, Any]:
    p = _default_highT_params()
    p.update(
        {
            "sticking_probability": 0.3,
            "sticking_temp_model": "constant",
            "er_cross_section_cm2": 1.5e-15,
            "er_reaction_probability": 0.9,
            "beam_dissociation_fraction": 1.0,
            "enable_h2_blocking": False,
        }
    )
    return p


def _calibrated_ded_params() -> Dict[str, Any]:
    p = _default_ded_params()
    p.update(
        {
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
    )
    return p


def _beam_baseline_separation_check(max_steps: int | None) -> Dict[str, Any]:
    """
    Demonstrate that undissociated beam H2 can adsorb/desorb without contaminating
    the ε numerator (which counts only formed, prompt-desorbed H2).
    """
    p: Dict[str, Any] = _calibrated_ded_params()
    p.update(
        {
            "surface_temperature_k": 100.0,
            "beam_dissociation_fraction": 0.0,  # no atoms impinging => ε should be 0
            "E_h2_bind_eV": 0.01,  # ensure beam H2 desorbs readily on this timescale
            "enable_LH": False,
            "enable_diffusion": False,
        }
    )
    p.pop("max_arrivals", None)
    kmc = KineticMonteCarlo(p)
    kmc.run_gillespie(max_time=1.0, max_steps=max_steps)
    return {
        "impinging_atoms": int(getattr(kmc, "total_impinging_h_atoms", 0)),
        "impinging_h2": int(getattr(kmc, "total_impinging_h2_molecules", 0)),
        "h2_prompt_desorbed_formed": int(getattr(kmc, "h2_molecules_desorbed", 0)),
        "h2_released_beam": int(getattr(kmc, "h2_molecules_released_beam", 0)),
        "h2_released_formed": int(getattr(kmc, "h2_molecules_released_formed", 0)),
    }


def run_checks(burnin: int, measure: int, max_steps: int | None) -> Dict[str, Any]:
    out: Dict[str, Any] = {}

    # 1) High-T plateau depends on chemisorption reservoir.
    base_highT = _calibrated_highT_params()
    eps_baseline_150 = _isothermal_epsilon(base_highT, 150.0, burnin, measure, seed=1234, max_steps=max_steps)
    out["highT_baseline_150K"] = {"epsilon": eps_baseline_150}

    no_chem = dict(base_highT)
    no_chem["chemisorption_fraction"] = 0.0
    eps_no_chem = _isothermal_epsilon(no_chem, 150.0, burnin, measure, seed=1234, max_steps=max_steps)
    out["highT_no_chemisorption_150K"] = {"epsilon": eps_no_chem}

    # 2) Low-T collapse depends on blocking (use the DED-style ramp, which is the experimental protocol).
    base_ded = _calibrated_ded_params()
    ded_block = evaluate_ded(
        base_params=base_ded,
        replicates=1,
        burnin_arrivals=int(burnin),
        t_start_k=10.0,
        t_end_k=80.0,
        rate_k_per_min=1.0,
        bin_width_k=1.0,
        max_steps=max_steps,
    )
    ratio_block = float(ded_block.epsilon_10k / ded_block.epsilon_20k) if float(ded_block.epsilon_20k) > 0 else 0.0
    out["ded_blocking_on"] = {
        "epsilon_10K": float(ded_block.epsilon_10k),
        "epsilon_20K": float(ded_block.epsilon_20k),
        "epsilon_30_80K": float(ded_block.epsilon_30_80k),
        "ratio_10_over_20": ratio_block,
    }

    no_block = dict(base_ded)
    no_block["enable_h2_blocking"] = False
    no_block["sticking_blocking_strength"] = 0.0
    no_block["er_blocking_strength"] = 0.0
    ded_noblock = evaluate_ded(
        base_params=no_block,
        replicates=1,
        burnin_arrivals=int(burnin),
        t_start_k=10.0,
        t_end_k=80.0,
        rate_k_per_min=1.0,
        bin_width_k=1.0,
        max_steps=max_steps,
    )
    ratio_noblock = float(ded_noblock.epsilon_10k / ded_noblock.epsilon_20k) if float(ded_noblock.epsilon_20k) > 0 else 0.0
    out["ded_blocking_off"] = {
        "epsilon_10K": float(ded_noblock.epsilon_10k),
        "epsilon_20K": float(ded_noblock.epsilon_20k),
        "epsilon_30_80K": float(ded_noblock.epsilon_30_80k),
        "ratio_10_over_20": ratio_noblock,
    }

    # 3) Beam baseline separation: beam-origin H2 can desorb without inflating ε.
    out["beam_baseline_separation"] = _beam_baseline_separation_check(max_steps=max_steps)

    return out


def main() -> None:
    p = argparse.ArgumentParser(description="Mechanism turn-off checks for Grieco-style ε(T) validation.")
    p.add_argument("--out", default="results/grieco_mechanism_checks.json")
    p.add_argument("--burnin", type=int, default=2000)
    p.add_argument("--measure", type=int, default=5000)
    p.add_argument("--max-steps", type=int, default=500000)
    p.add_argument("--assert", dest="do_assert", action="store_true", help="Exit nonzero if checks fail thresholds")
    p.add_argument("--min-highT-eps", type=float, default=0.05)
    p.add_argument("--max-no-chem-eps", type=float, default=0.01)
    p.add_argument("--max-blocked-ratio10-20", type=float, default=0.30)
    p.add_argument("--min-unblocked-ratio10-20", type=float, default=0.80)
    p.add_argument("--min-beam-h2-released", type=int, default=1)
    args = p.parse_args()

    max_steps = int(args.max_steps) if args.max_steps else None
    out = run_checks(burnin=int(args.burnin), measure=int(args.measure), max_steps=max_steps)
    os.makedirs(os.path.dirname(str(args.out)) or ".", exist_ok=True)
    with open(str(args.out), "w") as f:
        json.dump(out, f, indent=2, sort_keys=True)
    print(json.dumps(out, indent=2, sort_keys=True))
    print(f"Wrote {args.out}")

    if args.do_assert:
        errors: list[str] = []

        eps_highT = float(out["highT_baseline_150K"]["epsilon"])
        eps_no_chem = float(out["highT_no_chemisorption_150K"]["epsilon"])
        if eps_highT < float(args.min_highT_eps):
            errors.append(f"highT_baseline_150K epsilon {eps_highT:.3g} < {float(args.min_highT_eps):.3g}")
        if eps_no_chem > float(args.max_no_chem_eps):
            errors.append(f"highT_no_chemisorption_150K epsilon {eps_no_chem:.3g} > {float(args.max_no_chem_eps):.3g}")

        ratio_block = float(out["ded_blocking_on"]["ratio_10_over_20"])
        ratio_noblock = float(out["ded_blocking_off"]["ratio_10_over_20"])
        if ratio_block > float(args.max_blocked_ratio10_20):
            errors.append(f"ded_blocking_on ratio10/20 {ratio_block:.3g} > {float(args.max_blocked_ratio10_20):.3g}")
        if ratio_noblock < float(args.min_unblocked_ratio10_20):
            errors.append(f"ded_blocking_off ratio10/20 {ratio_noblock:.3g} < {float(args.min_unblocked_ratio10_20):.3g}")

        beam = out["beam_baseline_separation"]
        imp_atoms = int(beam["impinging_atoms"])
        imp_h2 = int(beam["impinging_h2"])
        h2_prompt = int(beam["h2_prompt_desorbed_formed"])
        h2_rel_beam = int(beam["h2_released_beam"])
        if imp_atoms != 0:
            errors.append(f"beam baseline: expected 0 impinging atoms, got {imp_atoms}")
        if imp_h2 <= 0:
            errors.append("beam baseline: expected >0 impinging H2 molecules")
        if h2_prompt != 0:
            errors.append(f"beam baseline: expected 0 prompt-desorbed formed H2, got {h2_prompt}")
        if h2_rel_beam < int(args.min_beam_h2_released):
            errors.append(
                f"beam baseline: expected >= {int(args.min_beam_h2_released)} beam H2 releases, got {h2_rel_beam}"
            )

        if errors:
            print("ASSERTION FAILURES:")
            for e in errors:
                print(f"- {e}")
            raise SystemExit(1)


if __name__ == "__main__":
    main()
