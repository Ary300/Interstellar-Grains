import argparse
import csv
import os
from dataclasses import dataclass
from typing import Any, Dict, List, Tuple

import numpy as np

from calibrate_grieco import _default_ded_params, _default_highT_params, evaluate_ded, evaluate_highT_plateau


@dataclass(frozen=True)
class Metrics:
    plateau_eps: float
    ded_eps10: float
    ded_eps20: float
    ded_eps30_80: float

    @property
    def drop_ratio(self) -> float:
        if self.ded_eps20 <= 0:
            return 0.0
        return float(self.ded_eps10) / float(self.ded_eps20)


def _calibrated_base_params() -> Tuple[Dict[str, Any], Dict[str, Any]]:
    # Defaults match the calibration harness geometry/energetics.
    highT = _default_highT_params()
    ded = _default_ded_params()

    # Calibrated knobs from results/grieco_calibration_local_v3.json.
    calibrated = {
        "sticking_probability": 0.3,
        "er_cross_section_cm2": 1.5e-15,
        "er_reaction_probability": 0.9,
        "enable_h2_blocking": True,
        "E_h2_bind_eV": 0.06,
        "h2_desorption_prefactor_s": 1e12,
        "h2_stick_transition_K": 20.0,
        "h2_stick_prob_lowT": 0.9,
        "sticking_blocking_strength": 1.0,
        "er_blocking_strength": 1.0,
        "beam_dissociation_fraction": 1.0,
    }

    highT.update({k: v for k, v in calibrated.items() if k in highT})
    # Ensure high-T harness stays in its intended regime.
    highT.update(
        {
            "sticking_probability": calibrated["sticking_probability"],
            "er_cross_section_cm2": calibrated["er_cross_section_cm2"],
            "er_reaction_probability": calibrated["er_reaction_probability"],
            "enable_h2_blocking": False,
        }
    )

    ded.update(calibrated)
    return highT, ded


def _evaluate(
    highT_params: Dict[str, Any],
    ded_params: Dict[str, Any],
    plateau_temps: List[float],
    plateau_replicates: int,
    plateau_burnin: int,
    plateau_measure: int,
    ded_replicates: int,
    ded_burnin: int,
    ded_start: float,
    ded_end: float,
    ded_rate_k_per_min: float,
    ded_bin_width: float,
    max_steps: int | None,
) -> Metrics:
    iso = evaluate_highT_plateau(
        base_params=highT_params,
        temperatures_k=plateau_temps,
        replicates=plateau_replicates,
        burnin_arrivals=plateau_burnin,
        measure_arrivals=plateau_measure,
        max_steps=max_steps,
    )
    ded = evaluate_ded(
        base_params=ded_params,
        replicates=ded_replicates,
        burnin_arrivals=ded_burnin,
        t_start_k=ded_start,
        t_end_k=ded_end,
        rate_k_per_min=ded_rate_k_per_min,
        bin_width_k=ded_bin_width,
        max_steps=max_steps,
    )
    return Metrics(
        plateau_eps=float(iso.epsilon_mean_plateau),
        ded_eps10=float(ded.epsilon_10k),
        ded_eps20=float(ded.epsilon_20k),
        ded_eps30_80=float(ded.epsilon_30_80k),
    )


def _rank_columns(x: np.ndarray) -> np.ndarray:
    # Simple rank transform (ties are unlikely for continuous sampling).
    out = np.zeros_like(x, dtype=float)
    for j in range(x.shape[1]):
        order = np.argsort(x[:, j], kind="mergesort")
        ranks = np.empty(x.shape[0], dtype=float)
        ranks[order] = np.arange(1, x.shape[0] + 1, dtype=float)
        out[:, j] = ranks
    return out


def _prcc(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """
    Partial rank correlation coefficients between columns of x and scalar y,
    controlling for other x columns.
    """
    x_r = _rank_columns(x)
    y_r = _rank_columns(y.reshape(-1, 1)).reshape(-1)
    z = np.column_stack([x_r, y_r])
    c = np.corrcoef(z, rowvar=False)
    p = np.linalg.inv(c)
    k = x.shape[1]
    out = np.zeros(k, dtype=float)
    for i in range(k):
        out[i] = float(-p[i, k] / np.sqrt(p[i, i] * p[k, k]))
    return out


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
    p = argparse.ArgumentParser(description="Sensitivity/identifiability helpers for calibrated Grieco-style ε(T).")
    p.add_argument("--mode", choices=["oat", "prcc", "both"], default="oat")
    p.add_argument("--variation-frac", type=float, default=0.2, help="±fraction for one-at-a-time sensitivity")

    p.add_argument("--oat-out", default="results/grieco_sensitivity_oat.csv")
    p.add_argument("--prcc-out", default="results/grieco_sensitivity_prcc.csv")

    p.add_argument("--plateau-temps", nargs="+", type=float, default=[100, 150, 200, 250])
    p.add_argument("--plateau-replicates", type=int, default=2)
    p.add_argument("--plateau-burnin", type=int, default=2000)
    p.add_argument("--plateau-measure", type=int, default=5000)

    p.add_argument("--ded-replicates", type=int, default=1)
    p.add_argument("--ded-burnin", type=int, default=2000)
    p.add_argument("--ded-start", type=float, default=10.0)
    p.add_argument("--ded-end", type=float, default=80.0)
    p.add_argument("--ded-rate-k-per-min", type=float, default=1.0)
    p.add_argument("--ded-bin-width", type=float, default=1.0)

    p.add_argument("--max-steps", type=int, default=500000)

    p.add_argument("--prcc-samples", type=int, default=30, help="Random samples for PRCC (can be slow)")
    p.add_argument("--seed", type=int, default=1234)
    args = p.parse_args()

    max_steps = int(args.max_steps) if args.max_steps else None
    plateau_temps = [float(t) for t in args.plateau_temps]

    base_highT, base_ded = _calibrated_base_params()

    knobs = [
        ("sticking_probability", "prob"),
        ("er_cross_section_cm2", "pos"),
        ("er_reaction_probability", "prob"),
        ("E_h2_bind_eV", "pos"),
        ("sticking_blocking_strength", "pos"),
        ("er_blocking_strength", "pos"),
        ("diffusion_to_binding_ratio_physisorption", "ratio"),
    ]

    baseline = _evaluate(
        highT_params=base_highT,
        ded_params=base_ded,
        plateau_temps=plateau_temps,
        plateau_replicates=int(args.plateau_replicates),
        plateau_burnin=int(args.plateau_burnin),
        plateau_measure=int(args.plateau_measure),
        ded_replicates=int(args.ded_replicates),
        ded_burnin=int(args.ded_burnin),
        ded_start=float(args.ded_start),
        ded_end=float(args.ded_end),
        ded_rate_k_per_min=float(args.ded_rate_k_per_min),
        ded_bin_width=float(args.ded_bin_width),
        max_steps=max_steps,
    )

    if args.mode in {"oat", "both"}:
        rows: List[Dict[str, Any]] = []
        rows.append(
            {
                "knob": "baseline",
                "variant": "baseline",
                "value": "",
                "plateau_eps": baseline.plateau_eps,
                "ded_eps10": baseline.ded_eps10,
                "ded_eps20": baseline.ded_eps20,
                "ded_eps30_80": baseline.ded_eps30_80,
                "drop_ratio": baseline.drop_ratio,
            }
        )

        for name, kind in knobs:
            base_val = float(base_ded.get(name, base_highT.get(name, np.nan)))
            if not np.isfinite(base_val):
                continue
            for tag, factor in [("low", 1.0 - float(args.variation_frac)), ("high", 1.0 + float(args.variation_frac))]:
                val = float(base_val) * float(factor)
                if kind == "prob":
                    val = min(1.0, max(0.0, val))
                elif kind == "ratio":
                    val = min(1.0, max(0.0, val))
                else:
                    val = max(0.0, val)

                hT = dict(base_highT)
                dd = dict(base_ded)
                if name in hT:
                    hT[name] = val
                if name in dd:
                    dd[name] = val

                m = _evaluate(
                    highT_params=hT,
                    ded_params=dd,
                    plateau_temps=plateau_temps,
                    plateau_replicates=int(args.plateau_replicates),
                    plateau_burnin=int(args.plateau_burnin),
                    plateau_measure=int(args.plateau_measure),
                    ded_replicates=int(args.ded_replicates),
                    ded_burnin=int(args.ded_burnin),
                    ded_start=float(args.ded_start),
                    ded_end=float(args.ded_end),
                    ded_rate_k_per_min=float(args.ded_rate_k_per_min),
                    ded_bin_width=float(args.ded_bin_width),
                    max_steps=max_steps,
                )
                rows.append(
                    {
                        "knob": name,
                        "variant": tag,
                        "value": val,
                        "plateau_eps": m.plateau_eps,
                        "ded_eps10": m.ded_eps10,
                        "ded_eps20": m.ded_eps20,
                        "ded_eps30_80": m.ded_eps30_80,
                        "drop_ratio": m.drop_ratio,
                        "delta_plateau": m.plateau_eps - baseline.plateau_eps,
                        "delta_mid": m.ded_eps30_80 - baseline.ded_eps30_80,
                        "delta_drop_ratio": m.drop_ratio - baseline.drop_ratio,
                    }
                )

        _write_csv(str(args.oat_out), rows)
        print(f"Wrote OAT sensitivity to {args.oat_out}")

    if args.mode in {"prcc", "both"}:
        rng = np.random.default_rng(int(args.seed))
        k = len(knobs)
        n = int(args.prcc_samples)
        x = np.zeros((n, k), dtype=float)
        y_plateau = np.zeros(n, dtype=float)
        y_mid = np.zeros(n, dtype=float)
        y_drop = np.zeros(n, dtype=float)

        base_vals = []
        for name, _kind in knobs:
            base_vals.append(float(base_ded.get(name, base_highT.get(name, 0.0))))
        base_vals = np.array(base_vals, dtype=float)

        for i in range(n):
            factors = rng.uniform(1.0 - float(args.variation_frac), 1.0 + float(args.variation_frac), size=k)
            vals = base_vals * factors
            for j, (_name, kind) in enumerate(knobs):
                if kind == "prob":
                    vals[j] = min(1.0, max(0.0, float(vals[j])))
                elif kind == "ratio":
                    vals[j] = min(1.0, max(0.0, float(vals[j])))
                else:
                    vals[j] = max(0.0, float(vals[j]))

            hT = dict(base_highT)
            dd = dict(base_ded)
            for (name, _kind), v in zip(knobs, vals.tolist()):
                if name in hT:
                    hT[name] = float(v)
                if name in dd:
                    dd[name] = float(v)

            m = _evaluate(
                highT_params=hT,
                ded_params=dd,
                plateau_temps=plateau_temps,
                plateau_replicates=1,
                plateau_burnin=int(args.plateau_burnin),
                plateau_measure=int(args.plateau_measure),
                ded_replicates=1,
                ded_burnin=int(args.ded_burnin),
                ded_start=float(args.ded_start),
                ded_end=float(args.ded_end),
                ded_rate_k_per_min=float(args.ded_rate_k_per_min),
                ded_bin_width=float(args.ded_bin_width),
                max_steps=max_steps,
            )

            x[i, :] = vals
            y_plateau[i] = float(m.plateau_eps)
            y_mid[i] = float(m.ded_eps30_80)
            y_drop[i] = float(m.drop_ratio)

        prcc_plateau = _prcc(x, y_plateau)
        prcc_mid = _prcc(x, y_mid)
        prcc_drop = _prcc(x, y_drop)

        rows = []
        for j, (name, _kind) in enumerate(knobs):
            rows.append(
                {
                    "knob": name,
                    "prcc_plateau": float(prcc_plateau[j]),
                    "prcc_mid_30_80": float(prcc_mid[j]),
                    "prcc_drop_ratio": float(prcc_drop[j]),
                }
            )
        _write_csv(str(args.prcc_out), rows)
        print(f"Wrote PRCC table to {args.prcc_out}")


if __name__ == "__main__":
    main()

