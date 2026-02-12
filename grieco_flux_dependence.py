import argparse
import csv
import json
import os
from typing import Any, Dict, List

import numpy as np

from calibrate_grieco import _isothermal_epsilon
from grieco_ded_validation import _default_ded_params, run_ded_validation
from grieco_validation import _default_coronene_like_params


def _write_csv(path: str, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        return
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def _mean_field_in_range(csv_path: str, lo: float, hi: float, field: str, temp_field: str = "temperature_k") -> float:
    rows = list(csv.DictReader(open(csv_path)))
    vals = [float(r[field]) for r in rows if float(lo) <= float(r[temp_field]) <= float(hi)]
    if not vals:
        return 0.0
    return float(np.mean(np.array(vals, dtype=float)))


def main() -> None:
    p = argparse.ArgumentParser(description="Flux dependence check for calibrated Grieco-style ε(T) harness.")
    p.add_argument("--rates", nargs="+", type=float, default=[0.005, 0.01, 0.02], help="arrival_rate_per_site_s values")
    p.add_argument("--out-csv", default="results/grieco_flux_dependence.csv", help="Aggregated CSV output")
    p.add_argument("--out-json", default="results/grieco_flux_dependence.json", help="Aggregated JSON output")

    # High-T (isothermal) settings
    p.add_argument("--highT-temp", type=float, default=150.0)
    p.add_argument("--highT-replicates", type=int, default=20)
    p.add_argument("--highT-burnin-arrivals", type=int, default=2000)
    p.add_argument("--highT-measure-arrivals", type=int, default=5000)
    p.add_argument("--highT-max-steps", type=int, default=500000)

    # DED settings
    p.add_argument("--ded-replicates", type=int, default=10)
    p.add_argument("--ded-t-start", type=float, default=10.0)
    p.add_argument("--ded-t-end", type=float, default=80.0)
    p.add_argument("--ded-rate-k-per-min", type=float, default=1.0)
    p.add_argument("--ded-bin-width", type=float, default=1.0)
    p.add_argument("--ded-burnin-arrivals", type=int, default=2000)
    p.add_argument("--ded-max-steps", type=int, default=800000)
    args = p.parse_args()

    rates = [float(x) for x in args.rates]

    out_rows: List[Dict[str, Any]] = []
    json_out: Dict[str, Any] = {"rates": rates, "runs": []}

    for rate in rates:
        # High-T: same calibrated harness, vary flux only.
        highT_params = _default_coronene_like_params()
        highT_params["arrival_rate_per_site_s"] = float(rate)
        eps_runs = [
            _isothermal_epsilon(
                params=highT_params,
                temperature_k=float(args.highT_temp),
                burnin_arrivals=int(args.highT_burnin_arrivals),
                measure_arrivals=int(args.highT_measure_arrivals),
                seed=1000 + i,
                max_steps=int(args.highT_max_steps) if args.highT_max_steps else None,
            )
            for i in range(int(args.highT_replicates))
        ]
        eps_arr = np.array(eps_runs, dtype=float)
        highT_mean = float(np.mean(eps_arr))
        highT_ci95 = (
            float(1.96 * float(np.std(eps_arr, ddof=1)) / float(np.sqrt(int(args.highT_replicates))))
            if int(args.highT_replicates) > 1
            else 0.0
        )

        # DED: run the calibrated DED harness and use its summary JSON for ε metrics/CI.
        ded_params = _default_ded_params()
        ded_params["arrival_rate_per_site_s"] = float(rate)
        tag = str(rate).replace(".", "p")
        ded_csv = f"results/grieco_ded_validation_flux_{tag}.csv"
        ded_summary = f"results/grieco_ded_validation_flux_{tag}_summary.json"
        run_ded_validation(
            output_csv=ded_csv,
            summary_json=ded_summary,
            replicates=int(args.ded_replicates),
            t_start_k=float(args.ded_t_start),
            t_end_k=float(args.ded_t_end),
            rate_k_per_min=float(args.ded_rate_k_per_min),
            bin_width_k=float(args.ded_bin_width),
            burnin_arrivals=int(args.ded_burnin_arrivals),
            max_steps=int(args.ded_max_steps) if args.ded_max_steps else None,
            base_params=ded_params,
        )
        ded_s = json.load(open(ded_summary))

        # Coverage diagnostics (time-weighted theta per temperature bin; averaged over bins near 10K/20K).
        theta10 = _mean_field_in_range(ded_csv, 9.5, 11.5, "theta_h2_mean")
        theta20 = _mean_field_in_range(ded_csv, 19.5, 21.5, "theta_h2_mean")
        theta30_80 = _mean_field_in_range(ded_csv, 30.0, 80.0, "theta_h2_mean")

        row = {
            "arrival_rate_per_site_s": float(rate),
            "highT_T_k": float(args.highT_temp),
            "highT_replicates": int(args.highT_replicates),
            "highT_eps_mean": highT_mean,
            "highT_eps_ci95": highT_ci95,
            "ded_t_start_k": float(args.ded_t_start),
            "ded_t_end_k": float(args.ded_t_end),
            "ded_replicates": int(args.ded_replicates),
            "ded_eps10_mean": float(ded_s["eps10_mean"]),
            "ded_eps10_ci95": float(ded_s.get("eps10_ci95", 0.0)),
            "ded_eps20_mean": float(ded_s["eps20_mean"]),
            "ded_eps20_ci95": float(ded_s.get("eps20_ci95", 0.0)),
            "ded_eps30_80_mean": float(ded_s["eps30_80_mean"]),
            "ded_eps30_80_ci95": float(ded_s.get("eps30_80_ci95", 0.0)),
            "ded_ratio10_over_20_mean": float(ded_s["ratio10_over_20_mean"]),
            "ded_ratio10_over_20_ci95": float(ded_s.get("ratio10_over_20_ci95", 0.0)),
            "ded_theta10_mean": theta10,
            "ded_theta20_mean": theta20,
            "ded_theta30_80_mean": theta30_80,
            "ded_bins_csv": ded_csv,
            "ded_summary_json": ded_summary,
        }
        out_rows.append(row)
        json_out["runs"].append(dict(row))

    _write_csv(str(args.out_csv), out_rows)
    os.makedirs(os.path.dirname(str(args.out_json)) or ".", exist_ok=True)
    with open(str(args.out_json), "w") as f:
        json.dump(json_out, f, indent=2, sort_keys=True)

    print(f"Wrote {args.out_csv}")
    print(f"Wrote {args.out_json}")


if __name__ == "__main__":
    main()

