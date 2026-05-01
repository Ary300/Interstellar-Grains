from __future__ import annotations

import argparse
import csv
import math
import os
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np

import matplotlib.pyplot as plt


@dataclass(frozen=True)
class Series:
    x: np.ndarray
    y: np.ndarray
    sigma: np.ndarray  # 1-sigma uncertainty (same units as y)


def _read_csv_rows(path: str) -> List[Dict[str, str]]:
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        return list(r)


def _float(v: object, default: float | None = None) -> float | None:
    if v is None:
        return default
    s = str(v).strip()
    if s == "":
        return default
    try:
        return float(s)
    except ValueError:
        return default


def _norm_method(v: str) -> str:
    s = str(v or "").strip().lower()
    if s in {"tpded", "t-pded", "temp-programmed", "temperature-programmed", "ramp", "circle", "circles"}:
        return "tpded"
    if s in {"isothermal", "iso", "triangle", "triangles", "ded"}:
        return "isothermal"
    return s


def _read_digitized_points(path: str) -> List[Dict[str, object]]:
    rows = _read_csv_rows(path)
    out: List[Dict[str, object]] = []
    for row in rows:
        # Support either T_K or temperature_k
        t = _float(row.get("T_K", None), None)
        if t is None:
            t = _float(row.get("temperature_k", None), None)
        if t is None:
            continue

        eps = _float(row.get("epsilon", None), None)
        if eps is None:
            # Support epsilon_percent -> epsilon fraction
            eps_pct = _float(row.get("epsilon_percent", None), None)
            if eps_pct is not None:
                eps = float(eps_pct) / 100.0
        if eps is None:
            continue

        err = _float(row.get("epsilon_err", None), None)
        err_lo = _float(row.get("epsilon_err_low", None), None)
        err_hi = _float(row.get("epsilon_err_high", None), None)
        if err_lo is None and err_hi is None and err is not None:
            err_lo = float(err)
            err_hi = float(err)
        if err_lo is None:
            err_lo = 0.0
        if err_hi is None:
            err_hi = 0.0

        method = _norm_method(str(row.get("method", "")))
        out.append(
            {
                "T_K": float(t),
                "epsilon": float(eps),
                "epsilon_err_low": float(abs(err_lo)),
                "epsilon_err_high": float(abs(err_hi)),
                "method": method,
                "notes": str(row.get("notes", "")),
            }
        )
    return out


def _read_iso_sim(path: str) -> Series:
    rows = _read_csv_rows(path)
    xs: List[float] = []
    ys: List[float] = []
    sig: List[float] = []
    for row in rows:
        t = _float(row.get("surface_temperature_k", None), None)
        eps = _float(row.get("epsilon_mean", None), None)
        if t is None or eps is None:
            continue
        # Simulation uncertainty is stored as CI95 in these CSVs; convert to sigma.
        ci95 = _float(row.get("epsilon_ci95", None), 0.0) or 0.0
        sigma = float(ci95) / 1.96 if ci95 > 0 else 0.0
        xs.append(float(t))
        ys.append(float(eps))
        sig.append(float(sigma))
    order = np.argsort(np.array(xs, dtype=float))
    return Series(
        x=np.array(xs, dtype=float)[order],
        y=np.array(ys, dtype=float)[order],
        sigma=np.array(sig, dtype=float)[order],
    )


def _read_ded_sim(path: str, y_col: str, ci_col: str) -> Series:
    rows = _read_csv_rows(path)
    xs: List[float] = []
    ys: List[float] = []
    sig: List[float] = []
    for row in rows:
        t = _float(row.get("temperature_k", None), None)
        eps = _float(row.get(y_col, None), None)
        if t is None or eps is None:
            continue
        ci95 = _float(row.get(ci_col, None), 0.0) or 0.0
        sigma = float(ci95) / 1.96 if ci95 > 0 else 0.0
        xs.append(float(t))
        ys.append(float(eps))
        sig.append(float(sigma))
    order = np.argsort(np.array(xs, dtype=float))
    return Series(
        x=np.array(xs, dtype=float)[order],
        y=np.array(ys, dtype=float)[order],
        sigma=np.array(sig, dtype=float)[order],
    )


def _interp(series: Series, x: float) -> Tuple[float, float]:
    y = float(np.interp(float(x), series.x, series.y))
    s = float(np.interp(float(x), series.x, series.sigma)) if series.sigma.size else 0.0
    return y, s


def _write_csv(path: str, rows: List[Dict[str, object]]) -> None:
    if not rows:
        return
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fieldnames = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _rmse(arr: np.ndarray) -> float:
    if arr.size == 0:
        return 0.0
    return float(np.sqrt(float(np.mean(arr**2))))


def _mae(arr: np.ndarray) -> float:
    if arr.size == 0:
        return 0.0
    return float(np.mean(np.abs(arr)))


def compute(
    digitized_csv: str,
    iso_sim_csv: str,
    ded_sim_csv: str,
    ded_y: str,
    ded_ci: str,
    include_model_sigma: bool,
    n_params: int,
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    pts = _read_digitized_points(digitized_csv)
    iso = _read_iso_sim(iso_sim_csv)
    ded = _read_ded_sim(ded_sim_csv, y_col=ded_y, ci_col=ded_ci)

    residual_rows: List[Dict[str, object]] = []

    for p in pts:
        t = float(p["T_K"])
        method = str(p["method"])
        exp_eps = float(p["epsilon"])
        exp_sig = 0.5 * (float(p["epsilon_err_low"]) + float(p["epsilon_err_high"]))

        if method == "tpded":
            model_eps, model_sig = _interp(ded, t)
            model_source = "TPDED sim"
        else:
            model_eps, model_sig = _interp(iso, t)
            model_source = "isothermal sim"

        sigma = float(exp_sig)
        if include_model_sigma:
            sigma = float(np.sqrt(float(exp_sig) ** 2 + float(model_sig) ** 2))

        used = bool(sigma > 0)
        resid = float(model_eps - exp_eps)
        z = float(resid / sigma) if used else 0.0
        residual_rows.append(
            {
                "T_K": t,
                "method": method,
                "epsilon_exp": exp_eps,
                "sigma_exp": exp_sig,
                "epsilon_model": model_eps,
                "sigma_model": model_sig,
                "sigma_used": sigma,
                "residual": resid,
                "z": z,
                "used_in_chi2": used,
                "model_source": model_source,
                "notes": str(p.get("notes", "")),
            }
        )

    stats_rows: List[Dict[str, object]] = []

    def add_stats(label: str, rows: List[Dict[str, object]]) -> None:
        used = [r for r in rows if bool(r["used_in_chi2"])]
        zs = np.array([float(r["z"]) for r in used], dtype=float)
        res = np.array([float(r["residual"]) for r in used], dtype=float)
        chi2 = float(np.sum(zs**2)) if zs.size else 0.0
        n = int(zs.size)
        dof = int(n - int(n_params))
        chi2_red = float(chi2 / dof) if dof > 0 else float("nan")
        stats_rows.append(
            {
                "subset": label,
                "n_points_used": n,
                "n_params_assumed": int(n_params),
                "dof": dof,
                "chi2": chi2,
                "chi2_reduced": chi2_red,
                "rmse": _rmse(res),
                "mae": _mae(res),
                "bias_mean_residual": float(np.mean(res)) if res.size else 0.0,
            }
        )

    add_stats("all", residual_rows)
    add_stats("tpded", [r for r in residual_rows if str(r["method"]) == "tpded"])
    add_stats("isothermal", [r for r in residual_rows if str(r["method"]) != "tpded"])

    return residual_rows, stats_rows


def plot_overlay(
    residual_rows: List[Dict[str, object]],
    iso_sim_csv: str,
    ded_sim_csv: str,
    ded_y: str,
    ded_ci: str,
    out_png: str,
) -> None:
    iso = _read_iso_sim(iso_sim_csv)
    ded = _read_ded_sim(ded_sim_csv, y_col=ded_y, ci_col=ded_ci)

    plt.figure(figsize=(8.5, 6.0))

    # Experimental points (digitized)
    tpded = [r for r in residual_rows if str(r["method"]) == "tpded"]
    iso_pts = [r for r in residual_rows if str(r["method"]) != "tpded"]

    def scatter(rows: List[Dict[str, object]], marker: str, label: str) -> None:
        if not rows:
            return
        xs = np.array([float(r["T_K"]) for r in rows], dtype=float)
        ys = np.array([float(r["epsilon_exp"]) for r in rows], dtype=float)
        yerr_lo = np.array([float(r["sigma_exp"]) for r in rows], dtype=float)
        yerr = np.vstack([yerr_lo, yerr_lo])
        plt.errorbar(xs, ys, yerr=yerr, fmt=marker, ms=7, capsize=3, label=label)

    scatter(tpded, "o", "Grieco Fig. 2 (digitized TPDED)")
    scatter(iso_pts, "^", "Grieco Fig. 2 (digitized isothermal)")

    # Simulation series
    plt.plot(ded.x, ded.y, "-", lw=2, label=f"TPDED sim ({ded_y})")
    plt.plot(iso.x, iso.y, "-", lw=2, label="Isothermal sim")

    plt.xlabel("T (K)")
    plt.ylabel("Recombination efficiency ε (fraction)")
    plt.title("Grieco Fig. 2: digitized experiment vs paperfit simulation")
    plt.grid(True, alpha=0.3)
    plt.xlim(0, 260)
    plt.ylim(0, 0.55)
    plt.legend(loc="upper right", frameon=False)

    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def plot_residuals(residual_rows: List[Dict[str, object]], out_png: str) -> None:
    plt.figure(figsize=(8.5, 4.5))
    tpded = [r for r in residual_rows if bool(r["used_in_chi2"]) and str(r["method"]) == "tpded"]
    iso = [r for r in residual_rows if bool(r["used_in_chi2"]) and str(r["method"]) != "tpded"]

    def pts(rows: List[Dict[str, object]], marker: str, label: str) -> None:
        if not rows:
            return
        xs = np.array([float(r["T_K"]) for r in rows], dtype=float)
        ys = np.array([float(r["residual"]) for r in rows], dtype=float)
        sig = np.array([float(r["sigma_used"]) for r in rows], dtype=float)
        yerr = np.vstack([sig, sig])
        plt.errorbar(xs, ys, yerr=yerr, fmt=marker, ms=7, capsize=3, label=label)

    pts(tpded, "o", "TPDED residuals")
    pts(iso, "^", "Isothermal residuals")
    plt.axhline(0.0, color="black", lw=1, alpha=0.7)
    plt.xlabel("T (K)")
    plt.ylabel("ε_model − ε_exp")
    plt.title("Residuals vs temperature (digitized Grieco Fig. 2)")
    plt.grid(True, alpha=0.3)
    plt.xlim(0, 260)
    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def _default_digitized_path() -> str:
    # Prefer repo-root path, but also support results/ path.
    if os.path.exists("grieco_fig2_digitized.csv"):
        return "grieco_fig2_digitized.csv"
    return "results/grieco_fig2_digitized.csv"


def main() -> None:
    p = argparse.ArgumentParser(description="Compare digitized Grieco Fig. 2 points to paperfit simulation outputs.")
    p.add_argument("--digitized-csv", default=_default_digitized_path())
    p.add_argument("--iso-sim-csv", default="results/grieco_validation_paper_iso_paperfit.csv")
    p.add_argument("--ded-sim-csv", default="results/grieco_ded_paper_ded_paperfit.csv")
    p.add_argument(
        "--ded-y",
        default="epsilon_released_total_mean",
        help="DED simulation y-column to compare against TPDED points (recommended: epsilon_released_total_mean).",
    )
    p.add_argument(
        "--ded-ci",
        default="epsilon_released_total_ci95",
        help="DED simulation CI95 column matching --ded-y (converted to sigma if --include-model-sigma is set).",
    )
    p.add_argument("--include-model-sigma", action="store_true", help="Combine experimental σ with simulation σ in χ².")
    p.add_argument(
        "--n-params",
        type=int,
        default=0,
        help="Number of calibrated parameters to use for reduced χ² (dof = N - n_params).",
    )
    p.add_argument("--out-residuals-csv", default="results/tables/grieco_fig2_digitized_residuals.csv")
    p.add_argument("--out-stats-csv", default="results/tables/grieco_fig2_digitized_fit_stats.csv")
    p.add_argument("--out-overlay-png", default="results/plots/grieco_fig2_overlay_digitized.png")
    p.add_argument("--out-residuals-png", default="results/plots/grieco_fig2_residuals_digitized.png")
    args = p.parse_args()

    if not os.path.exists(str(args.digitized_csv)):
        raise SystemExit(
            f"Digitized CSV not found: {args.digitized_csv}. Create it from digitization/grieco_fig2_digitized_template.csv."
        )

    residual_rows, stats_rows = compute(
        digitized_csv=str(args.digitized_csv),
        iso_sim_csv=str(args.iso_sim_csv),
        ded_sim_csv=str(args.ded_sim_csv),
        ded_y=str(args.ded_y),
        ded_ci=str(args.ded_ci),
        include_model_sigma=bool(args.include_model_sigma),
        n_params=int(args.n_params),
    )

    _write_csv(str(args.out_residuals_csv), residual_rows)
    _write_csv(str(args.out_stats_csv), stats_rows)
    plot_overlay(
        residual_rows=residual_rows,
        iso_sim_csv=str(args.iso_sim_csv),
        ded_sim_csv=str(args.ded_sim_csv),
        ded_y=str(args.ded_y),
        ded_ci=str(args.ded_ci),
        out_png=str(args.out_overlay_png),
    )
    plot_residuals(residual_rows=residual_rows, out_png=str(args.out_residuals_png))

    # Print a small console summary (no pandas required).
    for r in stats_rows:
        if r["subset"] == "all":
            print(
                f"χ²={r['chi2']:.3g}  dof={r['dof']}  χ²_red={r['chi2_reduced']:.3g}  "
                f"RMSE={r['rmse']:.3g}  N={r['n_points_used']}"
            )


if __name__ == "__main__":
    main()
