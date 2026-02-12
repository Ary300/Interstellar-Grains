import argparse
import json
import os
from dataclasses import dataclass
from typing import Any, Dict, List, Tuple

import numpy as np

from calibrate_grieco import _default_ded_params, _default_highT_params, _isothermal_epsilon, evaluate_ded, evaluate_highT_plateau


@dataclass(frozen=True)
class HoldoutPlateauResult:
    best_params: Dict[str, float]
    train_temps_k: List[float]
    test_temps_k: List[float]
    train_plateau_eps: float
    test_eps_by_T: Dict[float, float]
    test_plateau_eps: float


@dataclass(frozen=True)
class HoldoutDedResult:
    best_params: Dict[str, float]
    objective: str
    eps10: float
    eps20: float
    eps30_80: float
    ratio10_over_20: float


def _fit_plateau(
    *,
    train_temps_k: List[float],
    target_plateau: float,
    sticking_grid: List[float],
    er_cross_grid: List[float],
    er_prob_grid: List[float],
    replicates: int,
    burnin_arrivals: int,
    measure_arrivals: int,
    max_steps: int | None,
) -> HoldoutPlateauResult:
    base = _default_highT_params()
    best_obj: float | None = None
    best_params: Dict[str, float] | None = None
    best_train: float | None = None

    for sticking in sticking_grid:
        for er_cross in er_cross_grid:
            for er_prob in er_prob_grid:
                base["sticking_probability"] = float(sticking)
                base["er_cross_section_cm2"] = float(er_cross)
                base["er_reaction_probability"] = float(er_prob)
                summ = evaluate_highT_plateau(
                    base_params=base,
                    temperatures_k=train_temps_k,
                    replicates=replicates,
                    burnin_arrivals=burnin_arrivals,
                    measure_arrivals=measure_arrivals,
                    max_steps=max_steps,
                )
                obj = float((summ.epsilon_mean_plateau - float(target_plateau)) ** 2)
                if best_obj is None or obj < best_obj:
                    best_obj = obj
                    best_params = {
                        "sticking_probability": float(sticking),
                        "er_cross_section_cm2": float(er_cross),
                        "er_reaction_probability": float(er_prob),
                    }
                    best_train = float(summ.epsilon_mean_plateau)

    assert best_params is not None and best_train is not None
    return HoldoutPlateauResult(
        best_params=best_params,
        train_temps_k=[float(x) for x in train_temps_k],
        test_temps_k=[],
        train_plateau_eps=float(best_train),
        test_eps_by_T={},
        test_plateau_eps=0.0,
    )


def _eval_plateau_holdout(
    *,
    plateau_fit: HoldoutPlateauResult,
    test_temps_k: List[float],
    burnin_arrivals: int,
    measure_arrivals: int,
    replicates: int,
    max_steps: int | None,
) -> HoldoutPlateauResult:
    base = _default_highT_params()
    base.update(plateau_fit.best_params)

    eps_by_T: Dict[float, float] = {}
    for t in test_temps_k:
        eps_runs = [
            _isothermal_epsilon(
                params=base,
                temperature_k=float(t),
                burnin_arrivals=burnin_arrivals,
                measure_arrivals=measure_arrivals,
                seed=9000 + i,
                max_steps=max_steps,
            )
            for i in range(replicates)
        ]
        eps_by_T[float(t)] = float(np.mean(np.array(eps_runs, dtype=float)))

    test_plateau = float(np.mean(np.array(list(eps_by_T.values()), dtype=float))) if eps_by_T else 0.0
    return HoldoutPlateauResult(
        best_params=dict(plateau_fit.best_params),
        train_temps_k=list(plateau_fit.train_temps_k),
        test_temps_k=[float(x) for x in test_temps_k],
        train_plateau_eps=float(plateau_fit.train_plateau_eps),
        test_eps_by_T=eps_by_T,
        test_plateau_eps=float(test_plateau),
    )


def _fit_ded_blocking(
    *,
    objective: str,
    target_mid: float,
    target_drop_ratio: float,
    E_h2_grid: List[float],
    sticking_block_grid: List[float],
    er_block_grid: List[float],
    base_params: Dict[str, Any],
    replicates: int,
    burnin_arrivals: int,
    t_start_k: float,
    t_end_k: float,
    rate_k_per_min: float,
    bin_width_k: float,
    max_steps: int | None,
) -> HoldoutDedResult:
    best_obj: float | None = None
    best_params: Dict[str, float] | None = None
    best_summary = None

    for E_h2 in E_h2_grid:
        for sb in sticking_block_grid:
            for eb in er_block_grid:
                params = dict(base_params)
                params["E_h2_bind_eV"] = float(E_h2)
                params["sticking_blocking_strength"] = float(sb)
                params["er_blocking_strength"] = float(eb)
                summ = evaluate_ded(
                    base_params=params,
                    replicates=replicates,
                    burnin_arrivals=burnin_arrivals,
                    t_start_k=t_start_k,
                    t_end_k=t_end_k,
                    rate_k_per_min=rate_k_per_min,
                    bin_width_k=bin_width_k,
                    max_steps=max_steps,
                )

                eps10 = float(summ.epsilon_10k)
                eps20 = float(summ.epsilon_20k)
                eps30_80 = float(summ.epsilon_30_80k)
                ratio = float(eps10 / eps20) if eps20 > 0 else 0.0

                if objective == "drop_only":
                    desired_max_10 = float(eps20) * float(target_drop_ratio)
                    obj = float((eps10 - desired_max_10) ** 2) if eps10 > desired_max_10 else 0.0
                elif objective == "mid_only":
                    obj = float((eps30_80 - float(target_mid)) ** 2)
                else:
                    # both: mid + drop penalty (same spirit as calibrate_grieco.py)
                    mid_err = float((eps30_80 - float(target_mid)) ** 2)
                    desired_max_10 = float(eps20) * float(target_drop_ratio)
                    drop_pen = float((eps10 - desired_max_10) ** 2) if eps10 > desired_max_10 else 0.0
                    obj = float(mid_err + 5.0 * drop_pen)

                if best_obj is None or obj < best_obj:
                    best_obj = obj
                    best_summary = summ
                    best_params = {
                        "E_h2_bind_eV": float(E_h2),
                        "sticking_blocking_strength": float(sb),
                        "er_blocking_strength": float(eb),
                    }

    assert best_params is not None and best_summary is not None
    eps10 = float(best_summary.epsilon_10k)
    eps20 = float(best_summary.epsilon_20k)
    eps30_80 = float(best_summary.epsilon_30_80k)
    ratio = float(eps10 / eps20) if eps20 > 0 else 0.0
    return HoldoutDedResult(
        best_params=best_params,
        objective=str(objective),
        eps10=eps10,
        eps20=eps20,
        eps30_80=eps30_80,
        ratio10_over_20=ratio,
    )


def main() -> None:
    p = argparse.ArgumentParser(description="Hold-out style checks to guard against one-off curve fitting.")
    p.add_argument("--out", default="results/grieco_holdout.json", help="Output JSON path")
    p.add_argument("--max-steps", type=int, default=500000)

    p.add_argument("--plateau-target", type=float, default=0.20)
    p.add_argument("--plateau-train-temps", nargs="+", type=float, default=[100, 200])
    p.add_argument("--plateau-test-temps", nargs="+", type=float, default=[150, 250])
    p.add_argument("--plateau-replicates", type=int, default=2)
    p.add_argument("--plateau-burnin", type=int, default=2000)
    p.add_argument("--plateau-measure", type=int, default=5000)
    p.add_argument("--sticking-grid", nargs="+", type=float, default=[0.2, 0.3, 0.4])
    p.add_argument("--er-cross-grid", nargs="+", type=float, default=[1e-15, 1.5e-15, 2e-15])
    p.add_argument("--er-prob-grid", nargs="+", type=float, default=[0.6, 0.75, 0.9])

    p.add_argument("--ded-objective", choices=["drop_only", "mid_only", "both"], default="drop_only")
    p.add_argument("--ded-mid-target", type=float, default=0.30)
    p.add_argument("--ded-drop-ratio", type=float, default=0.30, help="Require eps10 <= drop_ratio * eps20")
    p.add_argument("--ded-replicates", type=int, default=1)
    p.add_argument("--ded-burnin", type=int, default=2000)
    p.add_argument("--ded-start", type=float, default=10.0)
    p.add_argument("--ded-end", type=float, default=80.0)
    p.add_argument("--ded-rate-k-per-min", type=float, default=1.0)
    p.add_argument("--ded-bin-width", type=float, default=1.0)
    p.add_argument("--E-h2-grid", nargs="+", type=float, default=[0.045, 0.055, 0.065, 0.075])
    p.add_argument("--stick-block-grid", nargs="+", type=float, default=[0.5, 1.0, 1.5])
    p.add_argument("--er-block-grid", nargs="+", type=float, default=[0.5, 1.0, 1.5])

    args = p.parse_args()

    max_steps = int(args.max_steps) if args.max_steps else None

    plateau_fit = _fit_plateau(
        train_temps_k=[float(x) for x in args.plateau_train_temps],
        target_plateau=float(args.plateau_target),
        sticking_grid=[float(x) for x in args.sticking_grid],
        er_cross_grid=[float(x) for x in args.er_cross_grid],
        er_prob_grid=[float(x) for x in args.er_prob_grid],
        replicates=int(args.plateau_replicates),
        burnin_arrivals=int(args.plateau_burnin),
        measure_arrivals=int(args.plateau_measure),
        max_steps=max_steps,
    )
    plateau_holdout = _eval_plateau_holdout(
        plateau_fit=plateau_fit,
        test_temps_k=[float(x) for x in args.plateau_test_temps],
        burnin_arrivals=int(args.plateau_burnin),
        measure_arrivals=int(args.plateau_measure),
        replicates=int(args.plateau_replicates),
        max_steps=max_steps,
    )

    # DED hold-out (by default, train on drop-only and report the mid-band as "free").
    base_ded = _default_ded_params()
    base_ded.update(
        {
            "sticking_probability": float(plateau_holdout.best_params["sticking_probability"]),
            "er_cross_section_cm2": float(plateau_holdout.best_params["er_cross_section_cm2"]),
            "er_reaction_probability": float(plateau_holdout.best_params["er_reaction_probability"]),
        }
    )
    ded_fit = _fit_ded_blocking(
        objective=str(args.ded_objective),
        target_mid=float(args.ded_mid_target),
        target_drop_ratio=float(args.ded_drop_ratio),
        E_h2_grid=[float(x) for x in args.E_h2_grid],
        sticking_block_grid=[float(x) for x in args.stick_block_grid],
        er_block_grid=[float(x) for x in args.er_block_grid],
        base_params=base_ded,
        replicates=int(args.ded_replicates),
        burnin_arrivals=int(args.ded_burnin),
        t_start_k=float(args.ded_start),
        t_end_k=float(args.ded_end),
        rate_k_per_min=float(args.ded_rate_k_per_min),
        bin_width_k=float(args.ded_bin_width),
        max_steps=max_steps,
    )

    out = {
        "plateau_holdout": {
            "target_plateau": float(args.plateau_target),
            "train_temps_k": plateau_holdout.train_temps_k,
            "test_temps_k": plateau_holdout.test_temps_k,
            "best_params": plateau_holdout.best_params,
            "train_plateau_eps": plateau_holdout.train_plateau_eps,
            "test_eps_by_T": {str(k): float(v) for k, v in plateau_holdout.test_eps_by_T.items()},
            "test_plateau_eps": plateau_holdout.test_plateau_eps,
        },
        "ded_holdout": {
            "objective": ded_fit.objective,
            "targets": {"mid_target": float(args.ded_mid_target), "drop_ratio": float(args.ded_drop_ratio)},
            "best_params": ded_fit.best_params,
            "eps10": ded_fit.eps10,
            "eps20": ded_fit.eps20,
            "eps30_80": ded_fit.eps30_80,
            "ratio10_over_20": ded_fit.ratio10_over_20,
        },
    }

    os.makedirs(os.path.dirname(str(args.out)) or ".", exist_ok=True)
    with open(str(args.out), "w") as f:
        json.dump(out, f, indent=2, sort_keys=True)
    print(json.dumps(out, indent=2, sort_keys=True))
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()

