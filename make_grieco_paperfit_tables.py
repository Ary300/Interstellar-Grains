#!/usr/bin/env python3
"""
Generate "paper-ready" CSV tables for the Grieco paperfit validation state.

Outputs:
- results/tables/grieco_paperfit_parameter_table.csv
- results/tables/grieco_paperfit_validation_summary_table.csv
"""

from __future__ import annotations

import csv
import os
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import yaml


def _read_yaml(path: str) -> Dict[str, Any]:
    with open(path, "r") as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise TypeError(f"Expected dict YAML in {path}, got {type(data)}")
    return data


def _read_csv_rows(path: str) -> List[Dict[str, Any]]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def _write_csv(path: str, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("No rows to write")
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fieldnames = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def _float_or_none(v: Any) -> float | None:
    if v is None:
        return None
    try:
        return float(v)
    except Exception:
        return None


def _describe_param(
    *,
    name: str,
    value: Any,
    units: str,
    role: str,
    status: str,
    constraint_or_range: str,
    reference: str,
    notes: str = "",
) -> Dict[str, Any]:
    return {
        "parameter": name,
        "value": value,
        "units": units,
        "role": role,
        "status": status,
        "constraint_or_range": constraint_or_range,
        "reference": reference,
        "notes": notes,
    }


def make_parameter_table(iso_cfg: Dict[str, Any], ded_cfg: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Table is intentionally explicit about what is:
    - fixed from the Grieco protocol/paper,
    - literature-constrained,
    - assumed,
    - calibrated as an effective parameter,
    - numerical (performance) choice.
    """
    rows: List[Dict[str, Any]] = []

    # Protocol / beam geometry
    rows.append(
        _describe_param(
            name="beam_incidence_angle_deg",
            value=_float_or_none(iso_cfg.get("beam_incidence_angle_deg")),
            units="deg",
            role="FORMOLISM beam incidence angle (projection factor cosθ for flux→arrival rate).",
            status="fixed (paper protocol)",
            constraint_or_range="from Grieco et al. (2023): 40°",
            reference="Grieco et al. 2023, Nat Astron, doi:10.1038/s41550-023-01902-4",
        )
    )
    rows.append(
        _describe_param(
            name="beam_flux_total_cm2_s",
            value=_float_or_none(iso_cfg.get("beam_flux_total_cm2_s")),
            units="cm^-2 s^-1",
            role="Directed beam flux at the surface (converted to arrivals via surface area and cosθ).",
            status="derived (paper conditions)",
            constraint_or_range="order 1e12–1e13 cm^-2 s^-1 (from partial pressure + kinetic theory)",
            reference="Grieco et al. 2023 (Methods; partial pressure estimate); kinetic theory",
            notes="Used only for lab/protocol-matching configs, not for ISM runs.",
        )
    )
    rows.append(
        _describe_param(
            name="beam_dissociation_fraction",
            value=_float_or_none(iso_cfg.get("beam_dissociation_fraction")),
            units="fraction",
            role="Fraction of arrivals that are atoms (τ in the paper’s ε definition).",
            status="assumed (with sensitivity recommended)",
            constraint_or_range="0–1 (paper measures τ before/after each experiment)",
            reference="Grieco et al. 2023, Nat Astron, doi:10.1038/s41550-023-01902-4",
        )
    )
    rows.append(
        _describe_param(
            name="burnin_exposure_atoms_cm2",
            value=_float_or_none(iso_cfg.get("burnin_exposure_atoms_cm2")),
            units="atoms cm^-2",
            role="Dose to reach steady state ('slightly superhydrogenated').",
            status="fixed (paper protocol)",
            constraint_or_range="few × 1e15 atoms/cm^2",
            reference="Grieco et al. 2023, Nat Astron, doi:10.1038/s41550-023-01902-4",
        )
    )
    rows.append(
        _describe_param(
            name="measure_exposure_atoms_cm2",
            value=_float_or_none(iso_cfg.get("measure_exposure_atoms_cm2")),
            units="atoms cm^-2",
            role="Measurement window dose for ε at fixed T.",
            status="fixed (paper protocol-inspired)",
            constraint_or_range="order 1e15 atoms/cm^2",
            reference="Grieco et al. 2023, Nat Astron, doi:10.1038/s41550-023-01902-4",
        )
    )

    # Microphysics (energetics + reaction knobs)
    rows.append(
        _describe_param(
            name="E_phys_mean_meV",
            value=_float_or_none(iso_cfg.get("E_phys_mean_meV")),
            units="meV",
            role="Mean binding energy for physisorbed H (sets residence time at moderate T).",
            status="fixed (paper-aligned)",
            constraint_or_range="typical ~40–45 meV used in literature",
            reference="Grieco et al. 2023 (typical 45 meV); Cuppen & Hornekær 2008 (40 meV), arXiv:0807.0108",
        )
    )
    rows.append(
        _describe_param(
            name="sticking_probability",
            value=_float_or_none(iso_cfg.get("sticking_probability")),
            units="probability",
            role="Accretion/sticking probability (enters arrivals→adsorption).",
            status="fixed (paper statement)",
            constraint_or_range="~0.5 (weak T dependence over this range)",
            reference="Grieco et al. 2023, Nat Astron, doi:10.1038/s41550-023-01902-4",
        )
    )
    rows.append(
        _describe_param(
            name="chemisorption_fraction",
            value=_float_or_none(iso_cfg.get("chemisorption_fraction")),
            units="fraction",
            role="Effective fraction of surface sites that act as chemisorption reservoirs (enables high-T plateau).",
            status="calibrated (effective parameter)",
            constraint_or_range="0–1 (bounded; represents reactive-site density)",
            reference="Model parameter (effective); calibrated to reproduce Grieco ε(T) regimes",
        )
    )
    rows.append(
        _describe_param(
            name="E_chem_mean_eV",
            value=_float_or_none(iso_cfg.get("E_chem_mean_eV")),
            units="eV",
            role="Mean binding energy for chemisorbed H reservoir sites.",
            status="literature-constrained (chosen)",
            constraint_or_range="~0.8–1.9 eV (DFT range on graphene/graphite configurations)",
            reference="Casolo et al. 2008, arXiv:0808.1312",
        )
    )
    rows.append(
        _describe_param(
            name="er_cross_section_cm2",
            value=_float_or_none(iso_cfg.get("er_cross_section_cm2")),
            units="cm^2",
            role="Effective ER/abstraction cross section per chemisorbed H target (sets plateau amplitude).",
            status="literature-constrained (chosen)",
            constraint_or_range="4–17 Å^2 (≈4e-16–1.7e-15 cm^2) reported for abstraction on graphite",
            reference="Cuppen & Hornekær 2008, arXiv:0807.0108 (cites Zecho et al. 2002)",
        )
    )
    rows.append(
        _describe_param(
            name="er_reaction_probability",
            value=_float_or_none(iso_cfg.get("er_reaction_probability")),
            units="probability",
            role="Reaction probability given an ER encounter (absorbs unresolved microphysics).",
            status="calibrated (effective parameter)",
            constraint_or_range="0–1",
            reference="Model parameter (effective); calibrated to reproduce high-T plateau",
        )
    )

    # Low-T blocking / TPDED mapping (DED config)
    rows.append(
        _describe_param(
            name="E_h2_bind_eV",
            value=_float_or_none(ded_cfg.get("E_h2_bind_eV")),
            units="eV",
            role="Effective H2 binding energy for low-T blocking (controls D2 coverage build-up).",
            status="calibrated within literature range",
            constraint_or_range="H2 physisorption binding energies on graphene/graphite reported ~0.01–0.06 eV",
            reference="Tozzini & Pellegrini 2013, PCCP, doi:10.1039/C2CP42538F",
        )
    )
    rows.append(
        _describe_param(
            name="h2_stick_prob_lowT",
            value=_float_or_none(ded_cfg.get("h2_stick_prob_lowT")),
            units="probability",
            role="Probability that newly formed (or beam-origin) H2 sticks at low T (blocking channel).",
            status="calibrated (effective parameter)",
            constraint_or_range="0–1",
            reference="Model parameter (effective); calibrated to reproduce 10–20 K suppression/peak structure",
        )
    )
    rows.append(
        _describe_param(
            name="h2_stick_transition_K",
            value=_float_or_none(ded_cfg.get("h2_stick_transition_K")),
            units="K",
            role="Temperature threshold separating low-T sticking/blocking vs negligible sticking.",
            status="fixed (regime boundary)",
            constraint_or_range="~20 K (matches the onset of the low-T regime in Grieco TPDED curve)",
            reference="Grieco et al. 2023, Fig. 2 (regime change around 20 K)",
        )
    )

    # Numerical/performance knobs (documented to avoid confusion with physics parameters)
    rows.append(
        _describe_param(
            name="diffusion_mode",
            value=str(ded_cfg.get("diffusion_mode")),
            units="enum",
            role="How diffusion is represented (explicit micro-events vs rate-only approximation).",
            status="numerical/performance choice",
            constraint_or_range="explicit | rate_only",
            reference="Code design choice (see kmc_simulation.py)",
            notes="Paperfit TPDED uses explicit diffusion; ISM configs use rate_only to avoid timestep collapse.",
        )
    )
    rows.append(
        _describe_param(
            name="diffusion_rate_cap_s",
            value=_float_or_none(ded_cfg.get("diffusion_rate_cap_s")),
            units="s^-1",
            role="Caps very fast diffusion rates to keep Gillespie step counts tractable.",
            status="numerical/performance choice",
            constraint_or_range=">=0",
            reference="Code design choice (see kmc_simulation.py)",
        )
    )

    return rows


def make_validation_summary_table(
    iso_csv: str,
    ded_summary_json: str,
) -> List[Dict[str, Any]]:
    iso_rows = _read_csv_rows(iso_csv)
    eps = np.array([float(r["epsilon_mean"]) for r in iso_rows], dtype=float)
    summary: List[Dict[str, Any]] = [
        {
            "metric": "highT_plateau_mean_100_250K",
            "value": float(np.mean(eps)) if len(eps) else 0.0,
            "note": "mean of isothermal ε points (100–250 K)",
        },
        {"metric": "highT_min", "value": float(np.min(eps)) if len(eps) else 0.0, "note": "min across isothermal points"},
        {"metric": "highT_max", "value": float(np.max(eps)) if len(eps) else 0.0, "note": "max across isothermal points"},
    ]

    import json

    with open(ded_summary_json, "r") as f:
        js = json.load(f)

    rel = js.get("released_total", js)
    summary.extend(
        [
            {
                "metric": "ded_eps10_released_total_mean",
                "value": float(rel.get("eps10_mean", 0.0)),
                "note": "mean over 9.5–11.5 K bins (released_total observable)",
            },
            {
                "metric": "ded_eps10_released_total_ci95",
                "value": float(rel.get("eps10_ci95", 0.0)),
                "note": "95% CI across TPDED replicates",
            },
            {
                "metric": "ded_eps20_released_total_mean",
                "value": float(rel.get("eps20_mean", 0.0)),
                "note": "mean over 19.5–21.5 K bins (released_total observable)",
            },
            {
                "metric": "ded_eps20_released_total_ci95",
                "value": float(rel.get("eps20_ci95", 0.0)),
                "note": "95% CI across TPDED replicates",
            },
            {
                "metric": "ded_eps30_80_released_total_mean",
                "value": float(rel.get("eps30_80_mean", 0.0)),
                "note": "mean over 30–80 K bins (released_total observable)",
            },
            {
                "metric": "ded_eps30_80_released_total_ci95",
                "value": float(rel.get("eps30_80_ci95", 0.0)),
                "note": "95% CI across TPDED replicates",
            },
            {
                "metric": "ded_ratio10_over_20_released_total_mean",
                "value": float(rel.get("ratio10_over_20_mean", 0.0)),
                "note": "eps10/eps20 (released_total observable)",
            },
            {
                "metric": "ded_ratio10_over_20_released_total_ci95",
                "value": float(rel.get("ratio10_over_20_ci95", 0.0)),
                "note": "95% CI across TPDED replicates",
            },
        ]
    )
    return summary


def main() -> None:
    iso_cfg_path = "config_grieco_paper_iso_paperfit.yaml"
    ded_cfg_path = "config_grieco_paper_ded_paperfit.yaml"

    iso_csv = "results/grieco_validation_paper_iso_paperfit.csv"
    ded_summary_json = "results/grieco_ded_paper_ded_paperfit_summary.json"

    iso_cfg = _read_yaml(iso_cfg_path)
    ded_cfg = _read_yaml(ded_cfg_path)

    param_rows = make_parameter_table(iso_cfg=iso_cfg, ded_cfg=ded_cfg)
    _write_csv("results/tables/grieco_paperfit_parameter_table.csv", param_rows)
    print("Wrote results/tables/grieco_paperfit_parameter_table.csv")

    summary_rows = make_validation_summary_table(iso_csv=iso_csv, ded_summary_json=ded_summary_json)
    _write_csv("results/tables/grieco_paperfit_validation_summary_table.csv", summary_rows)
    print("Wrote results/tables/grieco_paperfit_validation_summary_table.csv")


if __name__ == "__main__":
    main()

