#!/usr/bin/env python3
from __future__ import annotations

import argparse
import pickle
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
plt.style.use("./mnras_style.mplstyle")

import numpy as np
import pandas as pd
import yaml
from matplotlib.colors import ListedColormap, TwoSlopeNorm
from matplotlib.patches import Patch
from PIL import Image
from scipy.interpolate import griddata
from scipy.ndimage import binary_erosion, gaussian_filter, gaussian_filter1d

from kmc_simulation import KineticMonteCarlo
from mnras_figures import (
    COLORS,
    DENSITY_COLOR,
    DENSITY_MARKER,
    MECH_COLOR,
    add_regime_bands,
    annotate_value,
    fig_double,
    fig_single,
    figure_ensemble_convergence,
    figure_epsilon_all_densities,
    figure_grieco_validation,
    figure_phase_map,
    finalize,
    panel_label,
    plot_data_points,
    plot_with_ci,
    textbox,
)


ROOT = Path(__file__).resolve().parent
OUTDIR = ROOT / "results" / "plots" / "manuscript"
REPORT_PATH = OUTDIR / "BUILD_REPORT.md"
FIG03_TEX = OUTDIR / "fig03_mechanism_schematic.tex"
POVRAY_BIN = Path("/opt/homebrew/bin/povray") if Path("/opt/homebrew/bin/povray").exists() else Path("povray")

NON_SIM_KEYS = {
    "mode",
    "output_filename",
    "raw_runs_output",
    "mrn_output_filename",
    "parameter_sweeps",
    "explicit_conditions",
    "save_raw_runs",
    "ensemble_runs",
    "min_ensemble_runs",
    "max_time_s",
    "max_steps",
    "use_mrn",
    "aggregate_across_sizes",
    "mrn_min_um",
    "mrn_max_um",
    "mrn_bins",
    "burnin_arrivals",
    "measure_arrivals",
}


def _safe_sem(ci95: pd.Series | np.ndarray) -> np.ndarray:
    values = np.asarray(ci95, dtype=float)
    return values / 1.96


def _load_yaml_params(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle) or {}
    params = {}
    for key, value in raw.items():
        if key in NON_SIM_KEYS:
            continue
        if isinstance(value, list):
            continue
        if isinstance(value, dict) and key != "temp_ramp":
            continue
        params[key] = value
    return params


def _grain_pickle_path() -> Path:
    candidates = sorted(ROOT.glob("results/grain_cache*/*.pkl"))
    if not candidates:
        raise FileNotFoundError("No grain cache pickle found under results/grain_cache*")
    return candidates[0]


def _surface_shell_mask(lattice: np.ndarray) -> np.ndarray:
    occupied = lattice != None
    mask = np.zeros_like(occupied, dtype=bool)
    neighbor_offsets = (
        (-1, 0, 0),
        (1, 0, 0),
        (0, -1, 0),
        (0, 1, 0),
        (0, 0, -1),
        (0, 0, 1),
    )
    depth, rows, cols = occupied.shape
    for d, r, c in np.argwhere(occupied):
        for dd, dr, dc in neighbor_offsets:
            nd, nr, nc = d + dd, r + dr, c + dc
            if nd < 0 or nr < 0 or nc < 0 or nd >= depth or nr >= rows or nc >= cols:
                mask[d, r, c] = True
                break
            if not occupied[nd, nr, nc]:
                mask[d, r, c] = True
                break
    return mask


def _illustrative_grain_surface(
    grain_radius_um: float,
    site_area_angstroms_sq: float,
    porosity_fraction: float,
    chemisorption_fraction: float,
    surface_defect_fraction: float,
    rng_seed: int = 1000,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    spacing = float(np.sqrt(site_area_angstroms_sq))
    radius_a = float(grain_radius_um) * 1.0e4
    radius_cells = max(6, int(round(radius_a / spacing)))
    pad = 4
    n = 2 * (radius_cells + pad) + 1

    axis = np.arange(n, dtype=float) - (n - 1) / 2.0
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    rr = np.sqrt(xx**2 + yy**2 + zz**2)
    sphere = rr <= radius_cells

    rng = np.random.default_rng(rng_seed)

    # Low-frequency roughening gives the exterior a more amorphous silhouette.
    rough = gaussian_filter(rng.normal(size=sphere.shape), sigma=1.1)
    boundary_zone = np.abs(rr - radius_cells) <= 1.8
    sphere = sphere & ~(boundary_zone & (rough < np.quantile(rough[boundary_zone], 0.15)))

    # Carve a small number of spherical voids so the rendered shell shows pores/holes.
    n_voids = max(3, int(round(5 + 10 * porosity_fraction)))
    for _ in range(n_voids):
        direction = rng.normal(size=3)
        direction /= np.linalg.norm(direction)
        center_radius = rng.uniform(0.15, 0.72) * radius_cells
        cx, cy, cz = direction * center_radius
        void_radius = rng.uniform(0.14, 0.30) * radius_cells
        void = (xx - cx) ** 2 + (yy - cy) ** 2 + (zz - cz) ** 2 <= void_radius**2
        # Only keep cavities that intersect the outer shell so they are visually legible.
        if np.any(void & boundary_zone):
            sphere &= ~void

    # Keep a single coherent grain body.
    sphere &= rr <= radius_cells
    shell = sphere & ~binary_erosion(sphere, structure=np.ones((3, 3, 3), dtype=bool), border_value=0)
    coords = np.argwhere(shell)
    positions = coords[:, [2, 1, 0]].astype(float) * spacing
    positions -= positions.mean(axis=0, keepdims=True)

    # Cut away a front-side wedge so interior porosity is visible in projection.
    cut_plane = positions[:, 0] + 0.65 * positions[:, 1] - 0.15 * positions[:, 2]
    keep = cut_plane < 0.52 * radius_a
    positions = positions[keep]
    positions_clean = positions.copy()

    # Break the obvious cubic-grid rows so the render reads as an amorphous grain,
    # not a synthetic checkerboard of spheres.
    jitter = rng.normal(scale=0.18 * spacing, size=positions.shape)
    positions = np.asarray(positions + jitter, dtype=float)
    positions -= positions.mean(axis=0, keepdims=True)

    # Create spatially coherent site patches instead of salt-and-pepper random colors.
    norms = np.linalg.norm(positions, axis=1, keepdims=True)
    norms = np.where(norms == 0.0, 1.0, norms)
    unit = np.nan_to_num(np.asarray(positions / norms, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)

    def _patch_score(n_centers: int, width: float) -> np.ndarray:
        centers = np.asarray(rng.normal(size=(n_centers, 3)), dtype=float)
        center_norms = np.linalg.norm(centers, axis=1, keepdims=True)
        center_norms = np.where(center_norms == 0.0, 1.0, center_norms)
        centers /= center_norms
        dots = np.clip(np.asarray(np.einsum("ij,kj->ik", unit, centers), dtype=float), -1.0, 1.0)
        ang = np.arccos(dots)
        score = np.exp(-(ang**2) / (2.0 * width**2)).sum(axis=1)
        score += 0.12 * rng.normal(size=score.shape[0])
        return score

    n_surface = len(positions)
    n_chem = int(round(n_surface * chemisorption_fraction))
    n_def = int(round(n_surface * surface_defect_fraction))
    chem_score = _patch_score(6, 0.34)
    def_score = _patch_score(5, 0.26)

    chem_idx = np.argsort(chem_score)[-n_chem:]
    remaining = np.setdiff1d(np.arange(n_surface), chem_idx, assume_unique=False)
    def_rank = remaining[np.argsort(def_score[remaining])[-n_def:]]

    site_kind = np.full(n_surface, "regular", dtype=object)
    site_kind[chem_idx] = "chem"
    site_kind[def_rank] = "defect"
    positions = np.nan_to_num(np.asarray(positions, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    return positions, site_kind, positions_clean


def _write_grain_pov(pov_path: Path, positions: np.ndarray, site_kind: np.ndarray) -> None:
    body_positions = np.vstack([positions, 0.82 * positions, 0.62 * positions])
    visible = positions[:, 2] < np.quantile(positions[:, 2], 0.72)
    chem = positions[(site_kind == "chem") & visible]
    defect = positions[(site_kind == "defect") & visible]
    norms = np.linalg.norm(positions, axis=1, keepdims=True)
    norms = np.where(norms == 0.0, 1.0, norms)
    outward = positions / norms
    chem = chem + 2.1 * (chem / np.where(np.linalg.norm(chem, axis=1, keepdims=True) == 0.0, 1.0, np.linalg.norm(chem, axis=1, keepdims=True))) if len(chem) else chem
    defect = defect + 2.0 * (defect / np.where(np.linalg.norm(defect, axis=1, keepdims=True) == 0.0, 1.0, np.linalg.norm(defect, axis=1, keepdims=True))) if len(defect) else defect

    lines = [
        "#version 3.7;",
        "global_settings { assumed_gamma 1.0 }",
        "background { color rgb <1, 1, 1> }",
        "camera { perspective location <0, 0, -255> look_at <0, 0, 0> angle 20 }",
        "light_source { <140, 210, -260> color rgb <1.0, 1.0, 1.0> }",
        "light_source { <-180, 120, -120> color rgb <0.75, 0.75, 0.75> }",
        "light_source { <40, -180, -90> color rgb <0.35, 0.35, 0.35> }",
        "#declare GrainBody = blob {",
        "  threshold 0.58",
    ]
    for x, y, z in body_positions:
        lines.append(f"  sphere {{ <{x:.3f}, {y:.3f}, {z:.3f}>, 2.55, 1.0 }}")
    lines.extend(
        [
            "  texture {",
            "    pigment { color rgbf <0.86, 0.86, 0.86, 0.18> }",
            "    finish { ambient 0.18 diffuse 0.80 brilliance 1.10 phong 0.10 phong_size 16 specular 0.04 roughness 0.06 }",
            "  }",
            "}",
            "#declare ChemSites = union {",
        ]
    )
    for x, y, z in chem:
        lines.append(
            f"  sphere {{ <{x:.3f}, {y:.3f}, {z:.3f}>, 1.18 texture {{ pigment {{ color rgb <0.10, 0.47, 0.73> }} finish {{ ambient 0.14 diffuse 0.86 phong 0.20 phong_size 20 specular 0.08 roughness 0.03 }} }} no_shadow }}"
        )
    lines.extend(
        [
            "}",
            "#declare DefectSites = union {",
        ]
    )
    for x, y, z in defect:
        lines.append(
            f"  sphere {{ <{x:.3f}, {y:.3f}, {z:.3f}>, 1.13 texture {{ pigment {{ color rgb <0.79, 0.43, 0.10> }} finish {{ ambient 0.14 diffuse 0.84 phong 0.18 phong_size 18 specular 0.08 roughness 0.03 }} }} no_shadow }}"
        )
    lines.extend(
        [
            "}",
            "union {",
            "  object { GrainBody }",
            "  object { ChemSites }",
            "  object { DefectSites }",
            "  rotate <20, -34, 7>",
            "}",
            "",
        ]
    )
    pov_path.write_text("\n".join(lines), encoding="utf-8")


def _load_baseline() -> pd.DataFrame:
    df = pd.read_csv(ROOT / "results/jhub_full_merged.csv").dropna(
        subset=["surface_temperature_k", "h_gas_density_cm3", "uv_flux_factor", "epsilon_mean"]
    )
    df = df[np.isclose(df["uv_flux_factor"], 0.0)].copy()
    df["T_K"] = pd.to_numeric(df["surface_temperature_k"], errors="coerce")
    df["nH"] = pd.to_numeric(df["h_gas_density_cm3"], errors="coerce")
    df["eps_mean"] = pd.to_numeric(df["epsilon_mean"], errors="coerce")
    df["eps_sem"] = _safe_sem(pd.to_numeric(df["epsilon_ci95"], errors="coerce"))
    df["rate_mean"] = pd.to_numeric(df["h2_release_rate_cm2_s_mean"], errors="coerce")
    df["rate_sem"] = _safe_sem(pd.to_numeric(df["h2_release_rate_cm2_s_ci95"], errors="coerce"))
    total = pd.to_numeric(df["h2_formed_LH_mean"], errors="coerce") + pd.to_numeric(df["h2_formed_ER_mean"], errors="coerce")
    total = total.replace(0.0, np.nan)
    df["lh_fraction"] = pd.to_numeric(df["h2_formed_LH_mean"], errors="coerce") / total
    df["er_fraction"] = pd.to_numeric(df["h2_formed_ER_mean"], errors="coerce") / total
    df["surface_h_mean"] = pd.to_numeric(df["final_h_atoms_on_surface_mean"], errors="coerce")
    df["surface_h_sem"] = _safe_sem(pd.to_numeric(df["final_h_atoms_on_surface_ci95"], errors="coerce"))
    return df.sort_values(["nH", "T_K"]).dropna(subset=["T_K", "nH", "eps_mean"])


def _load_validation_curve() -> pd.DataFrame:
    df = pd.read_csv(ROOT / "results/grieco_validation_paper_iso_paperfit.csv").dropna(
        subset=["surface_temperature_k", "epsilon_mean"]
    )
    out = pd.DataFrame(
        {
            "T_K": pd.to_numeric(df["surface_temperature_k"], errors="coerce"),
            "eps_mean": pd.to_numeric(df["epsilon_mean"], errors="coerce"),
            "eps_sem": _safe_sem(pd.to_numeric(df["epsilon_ci95"], errors="coerce")),
        }
    )
    out = out.dropna().sort_values("T_K")
    out["T_K"] = out["T_K"].astype(float)
    out["eps_mean"] = out["eps_mean"].astype(float)
    out["eps_sem"] = out["eps_sem"].astype(float)
    return out


def _load_digitized(method: str | None = None, combine_errors: bool = True) -> pd.DataFrame:
    df = pd.read_csv(ROOT / "grieco_fig2_digitized.csv").dropna(subset=["T_K", "epsilon"])
    if method is not None and "method" in df.columns:
        df = df[df["method"] == method].copy()
    df["T_K"] = pd.to_numeric(df["T_K"], errors="coerce")
    df["eps"] = pd.to_numeric(df["epsilon"], errors="coerce")
    if combine_errors:
        df["eps_err"] = 0.5 * (
            pd.to_numeric(df["epsilon_err_low"], errors="coerce")
            + pd.to_numeric(df["epsilon_err_high"], errors="coerce")
        )
    df = df.dropna(subset=["T_K", "eps"]).sort_values("T_K")
    df["T_K"] = df["T_K"].astype(float)
    df["eps"] = df["eps"].astype(float)
    if "eps_err" in df.columns:
        df["eps_err"] = pd.to_numeric(df["eps_err"], errors="coerce").astype(float)
    return df


def _load_full_model_curve() -> pd.DataFrame:
    ded = pd.read_csv(ROOT / "results/grieco_ded_paper_ded_paperfit.csv")
    iso = pd.read_csv(ROOT / "results/grieco_validation_paper_iso_paperfit.csv")
    ded_curve = pd.DataFrame(
        {
            "T_K": pd.to_numeric(ded["temperature_k"], errors="coerce"),
            "eps_mean": pd.to_numeric(ded["epsilon_released_total_mean"], errors="coerce"),
            "eps_sem": _safe_sem(pd.to_numeric(ded["epsilon_released_total_ci95"], errors="coerce")),
        }
    )
    iso_curve = pd.DataFrame(
        {
            "T_K": pd.to_numeric(iso["surface_temperature_k"], errors="coerce"),
            "eps_mean": pd.to_numeric(iso["epsilon_mean"], errors="coerce"),
            "eps_sem": _safe_sem(pd.to_numeric(iso["epsilon_ci95"], errors="coerce")),
        }
    )
    combined = pd.concat([ded_curve[ded_curve["T_K"] < 100.0], iso_curve], ignore_index=True)
    return combined.dropna().sort_values("T_K")


def _write_fig03_tikz() -> None:
    tikz = r"""\begin{tikzpicture}[x=1cm,y=1cm,>=stealth]
\tikzstyle{surface}=[line width=0.9pt, draw=gray!70]
\tikzstyle{atom}=[circle, draw=white, line width=0.5pt, minimum size=0.34cm, inner sep=0pt]
\tikzstyle{labeltext}=[font=\footnotesize]
\tikzstyle{paneltitle}=[font=\small]

% Panel A
\begin{scope}[shift={(0,0)}]
  \node[anchor=west,font=\bfseries] at (0.0,2.25) {(a)};
  \node[paneltitle] at (2.1,2.22) {Langmuir-Hinshelwood};
  \draw[surface] (0.3,0.7)--(1.0,0.7)--(1.2,0.45)--(1.4,0.7)--(2.1,0.7)--(2.3,0.45)--(2.5,0.7)--(3.2,0.7)--(3.4,0.45)--(3.6,0.7)--(4.3,0.7);
  \node[atom, fill=blue!70] (h1) at (1.2,0.48) {\tiny H};
  \node[atom, fill=blue!70] (h2) at (3.4,0.48) {\tiny H};
  \draw[->, line width=0.8pt] (1.45,1.35) .. controls (1.8,1.0) and (2.0,0.9) .. (2.05,0.73);
  \draw[->, line width=0.8pt] (3.15,1.35) .. controls (2.8,1.0) and (2.6,0.9) .. (2.55,0.73);
  \draw[->, line width=0.9pt] (2.25,0.95) -- (2.25,1.55);
  \node[atom, fill=blue!70, minimum size=0.26cm] at (2.1,1.72) {};
  \node[atom, fill=blue!70, minimum size=0.26cm] at (2.4,1.72) {};
  \draw[line width=0.7pt] (2.2,1.72)--(2.3,1.72);
  \node[labeltext] at (2.25,1.96) {H$_2$};
\end{scope}

% Panel B
\begin{scope}[shift={(5.0,0)}]
  \node[anchor=west,font=\bfseries] at (0.0,2.25) {(b)};
  \node[paneltitle] at (2.1,2.22) {Eley-Rideal};
  \draw[surface] (0.3,0.7)--(1.8,0.7)--(2.1,0.32)--(2.4,0.7)--(4.3,0.7);
  \node[atom, fill=orange!80!black] (hs) at (2.1,0.42) {\tiny H};
  \node[atom, fill=blue!70] (hg) at (2.1,1.7) {\tiny H};
  \draw[->, line width=0.85pt] (2.1,1.5)--(2.1,0.86);
  \draw[->, line width=0.9pt] (2.35,0.95) -- (3.3,1.55);
  \node[atom, fill=orange!80!black, minimum size=0.26cm] at (3.55,1.7) {};
  \node[atom, fill=blue!70, minimum size=0.26cm] at (3.85,1.7) {};
  \draw[line width=0.7pt] (3.65,1.7)--(3.75,1.7);
  \node[labeltext] at (3.7,1.96) {H$_2$};
\end{scope}

% Panel C
\begin{scope}[shift={(10.0,0)}, opacity=0.62]
  \node[anchor=west,font=\bfseries] at (0.0,2.25) {(c)};
  \node[paneltitle] at (2.1,2.22) {UV (future work)};
  \draw[rounded corners=0.08cm, dashed, draw=gray!70] (0.15,0.15) rectangle (4.35,2.05);
  \foreach \x/\y in {2.1/1.1,1.7/1.45,1.25/1.3,1.25/0.9,1.7/0.75,2.5/0.75,2.95/0.9,2.95/1.3,2.5/1.45}
    \node[atom, fill=gray!70, minimum size=0.24cm] at (\x,\y) {};
  \node[atom, fill=blue!70, minimum size=0.24cm] at (1.0,1.55) {};
  \node[atom, fill=blue!70, minimum size=0.24cm] at (3.2,1.55) {};
  \draw[->, dashed, line width=0.8pt] (0.45,1.95) -- (1.15,1.45);
  \node[labeltext] at (0.55,2.1) {$h\nu$};
  \draw[->, line width=0.85pt] (3.0,1.45) -- (4.0,1.8);
  \node[labeltext] at (4.0,2.0) {H$_2$};
\end{scope}
\end{tikzpicture}
"""
    FIG03_TEX.write_text(tikz + "\n", encoding="utf-8")


def _crop_white_border(image_path: Path, padding: int = 12) -> Image.Image:
    image = Image.open(image_path).convert("RGB")
    arr = np.asarray(image)
    mask = np.any(arr < 248, axis=2)
    if not mask.any():
        return image
    ys, xs = np.where(mask)
    left = max(int(xs.min()) - padding, 0)
    right = min(int(xs.max()) + padding + 1, image.size[0])
    top = max(int(ys.min()) - padding, 0)
    bottom = min(int(ys.max()) + padding + 1, image.size[1])
    return image.crop((left, top, right, bottom))


def _cleanup_legacy_outputs(outdir: Path) -> None:
    legacy_stems = [
        "fig01_grain_structure",
        "fig01_grieco_validation",
        "fig02_binding_energy_distribution",
        "fig02_ism_overview",
        "fig03_mechanism_schematic",
        "fig05_physisorption_ablation",
        "fig07_efficiency_density",
        "fig08_mechanism_decomposition",
        "fig09_surface_inventory",
        "fig18_transition_boxplots",
        "fig19_lh_mode_consistency",
        "fig20_grain_size_dependence",
        "fig22_timescale_phase_map",
    ]
    suffixes = [".pdf", ".png", "_grey.png"]
    for stem in legacy_stems:
        for suffix in suffixes:
            path = outdir / f"{stem}{suffix}"
            if path.exists():
                path.unlink()


def make_fig01_render(outdir: str) -> None:
    params = _load_yaml_params(ROOT / "config_astro_full_paperfit.yaml")
    spacing = float(np.sqrt(params.get("site_area_angstroms_sq", 25.0)))
    positions, site_kind, _ = _illustrative_grain_surface(
        grain_radius_um=float(params.get("grain_radius_um", 0.005)),
        site_area_angstroms_sq=float(params.get("site_area_angstroms_sq", 25.0)),
        porosity_fraction=float(params.get("porosity_fraction", 0.2)),
        chemisorption_fraction=float(params.get("chemisorption_fraction", 0.4)),
        surface_defect_fraction=float(params.get("surface_defect_fraction", 0.15)),
        rng_seed=1000,
    )
    with _grain_pickle_path().open("rb") as handle:
        grain = pickle.load(handle)

    lattice = np.asarray(grain["lattice"], dtype=object)
    site_types = np.asarray(grain["site_types"], dtype=int)
    depth_layers, rows, cols = lattice.shape

    x_coords = (np.arange(cols, dtype=float) - 0.5 * (cols - 1)) * spacing
    y_coords = (np.arange(rows, dtype=float) - 0.5 * (rows - 1)) * spacing
    depth_coords = np.arange(depth_layers, dtype=float) * spacing

    top_occ = lattice[0] != None
    row_center = int(np.argmax(top_occ.sum(axis=1)))
    row_lo = max(0, row_center - 1)
    row_hi = min(rows, row_center + 2)
    section_map = np.full((depth_layers, cols), np.nan)
    for d in range(depth_layers):
        occ_strip = lattice[d, row_lo:row_hi, :] != None
        type_strip = site_types[d, row_lo:row_hi, :]
        for c in range(cols):
            vals = type_strip[:, c][occ_strip[:, c]]
            if vals.size == 0:
                continue
            counts = np.bincount(vals, minlength=4)
            if counts[2] >= max(counts[1], counts[3]):
                section_map[d, c] = 2.0
            elif counts[3] >= counts[1]:
                section_map[d, c] = 1.0
            else:
                section_map[d, c] = 0.0

    cmap = ListedColormap([COLORS["light_grey"], COLORS["vermillion"], COLORS["blue"]])
    cmap.set_bad(color="white", alpha=0.0)
    extent_section = [
        x_coords.min() - 0.5 * spacing,
        x_coords.max() + 0.5 * spacing,
        depth_coords.max() + 0.5 * spacing,
        depth_coords.min() - 0.5 * spacing,
    ]

    fig, (ax1, ax2) = plt.subplots(
        2,
        1,
        figsize=fig_single(3.5),
        gridspec_kw={"hspace": 0.16, "height_ratios": [1.0, 0.72]},
    )

    rx, ry, rz = np.deg2rad([25.0, -38.0, 8.0])
    rot_x = np.array([[1, 0, 0], [0, np.cos(rx), -np.sin(rx)], [0, np.sin(rx), np.cos(rx)]])
    rot_y = np.array([[np.cos(ry), 0, np.sin(ry)], [0, 1, 0], [-np.sin(ry), 0, np.cos(ry)]])
    rot_z = np.array([[np.cos(rz), -np.sin(rz), 0], [np.sin(rz), np.cos(rz), 0], [0, 0, 1]])
    proj = np.dot(np.nan_to_num(np.asarray(positions, dtype=float), nan=0.0, posinf=0.0, neginf=0.0), (rot_z @ rot_y @ rot_x).T)
    px, py, pz = proj[:, 0], proj[:, 1], proj[:, 2]
    bins = 180
    pad = 5.0
    x_edges = np.linspace(px.min() - pad, px.max() + pad, bins)
    y_edges = np.linspace(py.min() - pad, py.max() + pad, bins)
    density, _, _ = np.histogram2d(px, py, bins=[x_edges, y_edges])
    density = gaussian_filter(density.T, sigma=3.6)
    vmax = float(np.nanmax(density)) if np.isfinite(density).any() else 1.0
    levels = [0.12 * vmax, 0.25 * vmax, 0.45 * vmax, 0.72 * vmax, 1.01 * vmax]
    ax1.contourf(
        0.5 * (x_edges[:-1] + x_edges[1:]),
        0.5 * (y_edges[:-1] + y_edges[1:]),
        density,
        levels=levels,
        colors=["#f2f2f2", "#e4e4e4", "#d1d1d1", "#bdbdbd"],
        antialiased=True,
    )
    shell_r = np.sqrt(px**2 + py**2)
    rim = shell_r >= np.quantile(shell_r, 0.72)
    front = pz <= np.quantile(pz, 0.58)
    ax1.scatter(px[rim & front & (site_kind == "regular")], py[rim & front & (site_kind == "regular")], s=11, c=COLORS["light_grey"], alpha=0.9, linewidths=0, rasterized=True)
    ax1.scatter(px[rim & front & (site_kind == "defect")], py[rim & front & (site_kind == "defect")], s=17, c=COLORS["vermillion"], alpha=0.92, linewidths=0, rasterized=True)
    ax1.scatter(px[rim & front & (site_kind == "chem")], py[rim & front & (site_kind == "chem")], s=17, c=COLORS["blue"], alpha=0.92, linewidths=0, rasterized=True)
    ax1.set_aspect("equal")
    ax1.set_xticks([])
    ax1.set_yticks([])
    for spine in ax1.spines.values():
        spine.set_visible(False)
    ax1.text(
        0.84,
        0.08,
        r"radius $\approx$ 50 A",
        transform=ax1.transAxes,
        fontsize=7,
        color=COLORS["grey"],
        ha="right",
        va="bottom",
        style="italic",
        bbox=dict(fc="white", ec="none", alpha=0.82, pad=0.15),
    )
    ax1.text(0.04, 0.88, "shell cutaway", transform=ax1.transAxes, fontsize=7, color=COLORS["grey"], ha="left", va="top", style="italic")

    ax2.imshow(section_map, origin="upper", interpolation="nearest", cmap=cmap, vmin=0, vmax=2, extent=extent_section, aspect="auto")
    ax2.set_xlabel("x (A)")
    ax2.set_ylabel("Depth (A)")
    ax2.set_xlim(x_coords.min() - 2.0 * spacing, x_coords.max() + 2.0 * spacing)
    ax2.set_yticks(depth_coords)
    ax2.text(0.15, 0.92, "three-row meridional slice", transform=ax2.transAxes, fontsize=7, color=COLORS["grey"], ha="left", va="top", style="italic", bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.2))
    panel_label(ax1, "a")
    panel_label(ax2, "b")
    handles = [
        Patch(facecolor=COLORS["light_grey"], edgecolor="none", label="Regular site"),
        Patch(facecolor=COLORS["vermillion"], edgecolor="none", label="Surface defect"),
        Patch(facecolor=COLORS["blue"], edgecolor="none", label="Chemisorption site"),
    ]
    fig.legend(handles=handles, frameon=False, loc="lower center", bbox_to_anchor=(0.5, 0.01), ncol=3, fontsize=7)
    finalize(fig, "fig01_grain_lattice", outdir=outdir)


def make_fig24_surface_energy_map(outdir: str) -> None:
    spacing = 5.0
    with _grain_pickle_path().open("rb") as handle:
        grain = pickle.load(handle)
    lattice = np.asarray(grain["lattice"], dtype=object)
    site_types = np.asarray(grain["site_types"], dtype=int)
    ebind = np.asarray(grain["E_bind_eV_map"], dtype=float)
    top_occ = lattice[0] != None
    rows, cols = top_occ.shape
    x_coords = (np.arange(cols, dtype=float) - 0.5 * (cols - 1)) * spacing
    y_coords = (np.arange(rows, dtype=float) - 0.5 * (rows - 1)) * spacing
    extent = [
        x_coords.min() - 0.5 * spacing,
        x_coords.max() + 0.5 * spacing,
        y_coords.min() - 0.5 * spacing,
        y_coords.max() + 0.5 * spacing,
    ]

    site_map = np.full((rows, cols), np.nan)
    site_map[top_occ & (site_types[0] == 1)] = 0.0
    site_map[top_occ & (site_types[0] == 3)] = 1.0
    site_map[top_occ & (site_types[0] == 2)] = 2.0
    class_cmap = ListedColormap([COLORS["light_grey"], COLORS["vermillion"], COLORS["blue"]])
    class_cmap.set_bad(color="white", alpha=0.0)

    energy_map = np.full((rows, cols), np.nan)
    energy_map[top_occ] = ebind[0][top_occ]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=fig_double(2.8), gridspec_kw={"wspace": 0.15})
    ax1.imshow(site_map, origin="lower", interpolation="nearest", cmap=class_cmap, vmin=0, vmax=2, extent=extent)
    ax1.set_aspect("equal")
    ax1.set_xlabel("x (A)")
    ax1.set_ylabel("y (A)")
    panel_label(ax1, "a")
    ax1.text(0.04, 0.96, "top-layer site classes", transform=ax1.transAxes, fontsize=7, color=COLORS["grey"], ha="left", va="top", style="italic", bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.2))

    im = ax2.imshow(energy_map, origin="lower", interpolation="nearest", cmap="cividis", extent=extent)
    ax2.set_aspect("equal")
    ax2.set_xlabel("x (A)")
    ax2.set_ylabel("y (A)")
    panel_label(ax2, "b")
    ax2.text(0.04, 0.96, r"$E_{\rm bind}$ top layer", transform=ax2.transAxes, fontsize=7, color=COLORS["grey"], ha="left", va="top", style="italic", bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.2))
    cbar = fig.colorbar(im, ax=ax2, pad=0.02, fraction=0.046)
    cbar.set_label(r"$E_{\rm bind}$ (eV)")

    handles = [
        Patch(facecolor=COLORS["light_grey"], edgecolor="none", label="Regular"),
        Patch(facecolor=COLORS["vermillion"], edgecolor="none", label="Defect"),
        Patch(facecolor=COLORS["blue"], edgecolor="none", label="Chemisorption"),
    ]
    ax1.legend(handles=handles, frameon=False, loc="lower left", fontsize=6.6)
    finalize(fig, "fig24_surface_energy_map", outdir=outdir)


def make_fig25_layer_gallery(outdir: str) -> None:
    spacing = 5.0
    with _grain_pickle_path().open("rb") as handle:
        grain = pickle.load(handle)
    lattice = np.asarray(grain["lattice"], dtype=object)
    site_types = np.asarray(grain["site_types"], dtype=int)
    depth_layers, rows, cols = lattice.shape
    x_coords = (np.arange(cols, dtype=float) - 0.5 * (cols - 1)) * spacing
    y_coords = (np.arange(rows, dtype=float) - 0.5 * (rows - 1)) * spacing
    extent = [
        x_coords.min() - 0.5 * spacing,
        x_coords.max() + 0.5 * spacing,
        y_coords.min() - 0.5 * spacing,
        y_coords.max() + 0.5 * spacing,
    ]

    cmap = ListedColormap([COLORS["light_grey"], COLORS["vermillion"], COLORS["blue"]])
    cmap.set_bad(color="white", alpha=0.0)

    layer_indices = [0, depth_layers // 2, depth_layers - 1]
    layer_titles = ["Surface layer", "Mid-depth layer", "Deep layer"]
    fig, axes = plt.subplots(1, 3, figsize=fig_double(2.45), gridspec_kw={"wspace": 0.08})
    for ax, layer, title, letter in zip(axes, layer_indices, layer_titles, ["a", "b", "c"]):
        occ = lattice[layer] != None
        site_map = np.full((rows, cols), np.nan)
        site_map[occ & (site_types[layer] == 1)] = 0.0
        site_map[occ & (site_types[layer] == 3)] = 1.0
        site_map[occ & (site_types[layer] == 2)] = 2.0
        ax.imshow(site_map, origin="lower", interpolation="nearest", cmap=cmap, vmin=0, vmax=2, extent=extent)
        ax.set_aspect("equal")
        ax.set_xlabel("x (A)")
        if ax is axes[0]:
            ax.set_ylabel("y (A)")
        else:
            ax.set_yticklabels([])
        panel_label(ax, letter)
        ax.text(0.04, 0.96, title, transform=ax.transAxes, fontsize=6.9, color=COLORS["grey"], ha="left", va="top", style="italic", bbox=dict(fc="white", ec="none", alpha=0.82, pad=0.15))
    handles = [
        Patch(facecolor=COLORS["light_grey"], edgecolor="none", label="Regular"),
        Patch(facecolor=COLORS["vermillion"], edgecolor="none", label="Defect"),
        Patch(facecolor=COLORS["blue"], edgecolor="none", label="Chemisorption"),
    ]
    fig.legend(handles=handles, frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.02), ncol=3, fontsize=6.8)
    finalize(fig, "fig25_layer_gallery", outdir=outdir)


def make_fig02_binding(outdir: str) -> None:
    with _grain_pickle_path().open("rb") as handle:
        grain = pickle.load(handle)
    surface = np.asarray(grain["lattice"][0], dtype=object) != None
    site_types = np.asarray(grain["site_types"][0][surface], dtype=int)
    ebind = np.asarray(grain["E_bind_eV_map"][0][surface], dtype=float)
    phys = ebind[site_types == 1]
    chem = ebind[site_types == 2]

    fig, (ax1, ax2) = plt.subplots(
        1,
        2,
        figsize=fig_single(2.6),
        sharey=True,
        gridspec_kw={"width_ratios": [3.2, 1.4], "wspace": 0.05},
    )
    bins = np.logspace(-2, 0.5, 41)
    ax1.hist(phys, bins=bins, color=COLORS["blue"], alpha=0.78, edgecolor="none")
    ax2.hist(chem, bins=bins, color=COLORS["vermillion"], alpha=0.78, edgecolor="none")
    for ax in (ax1, ax2):
        ax.set_xscale("log")
    ax1.set_xlim(8.0e-3, 8.0e-2)
    ax2.set_xlim(0.8, 2.4)
    ax1.set_ylabel("Number of sites")
    fig.supxlabel(r"Binding energy $E_{\mathrm{bind}}$ (eV)")
    ymax = max(ax1.get_ylim()[1], ax2.get_ylim()[1])
    ax1.set_ylim(0, ymax)
    for temp in (20, 100, 250):
        energy = 8.617e-5 * temp
        ax1.axvline(energy, color=COLORS["grey"], lw=0.7, ls="--")
        x_text = max(energy * 1.08, ax1.get_xlim()[0] * 1.18)
        ax1.text(
            x_text,
            0.92 * ymax,
            f"{temp} K",
            rotation=90,
            fontsize=6.5,
            color=COLORS["grey"],
            ha="left",
            va="top",
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.25),
        )
    ax1.text(
        1.15e-2,
        0.86 * ymax,
        "Physisorption\n45 +/- 5 meV",
        color=COLORS["blue"],
        fontsize=7,
        ha="left",
        va="top",
        bbox=dict(fc="white", ec="none", alpha=0.9, pad=0.3),
    )
    ax2.text(
        0.92,
        0.12 * ymax,
        "Chemisorption\n1.75 +/- 0.25 eV",
        color=COLORS["vermillion"],
        fontsize=7,
        ha="left",
        va="top",
        bbox=dict(fc="white", ec="none", alpha=0.9, pad=0.3),
    )
    ax2.tick_params(labelleft=False, left=False)
    ax1.spines["right"].set_visible(False)
    ax2.spines["left"].set_visible(False)
    d = 0.012
    kwargs = dict(transform=ax1.transAxes, color="k", clip_on=False, lw=0.7)
    ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)
    ax1.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
    kwargs = dict(transform=ax2.transAxes, color="k", clip_on=False, lw=0.7)
    ax2.plot((-d, +d), (-d, +d), **kwargs)
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)
    finalize(fig, "fig02_binding_energies", outdir=outdir)


def make_fig05_physisorption_only(outdir: str) -> None:
    full = _load_full_model_curve()
    phys = pd.read_csv(ROOT / "results/grieco_physisorption_only_iso.csv")
    phys = pd.DataFrame(
        {
            "T_K": pd.to_numeric(phys["surface_temperature_k"], errors="coerce"),
            "eps_mean": pd.to_numeric(phys["epsilon_released_total_mean"], errors="coerce"),
            "eps_sem": _safe_sem(pd.to_numeric(phys["epsilon_released_total_ci95"], errors="coerce")),
        }
    ).dropna().sort_values("T_K")
    expt = _load_digitized()

    fig, ax = plt.subplots(figsize=fig_single(2.6))
    plot_with_ci(
        ax,
        full["T_K"],
        full["eps_mean"],
        y_err=full["eps_sem"],
        color=COLORS["blue"],
        linewidth=1.8,
        label="Full model",
    )
    plot_with_ci(
        ax,
        phys["T_K"],
        phys["eps_mean"],
        y_err=phys["eps_sem"],
        color=COLORS["vermillion"],
        linestyle="--",
        linewidth=1.8,
        label="Physisorption only",
    )
    plot_data_points(ax, expt["T_K"], expt["eps"], expt["eps_err"], color=COLORS["black"], label="Grieco et al. (2023)")
    ax.set_xlim(0, 270)
    ax.set_ylim(0, 0.35)
    ax.set_xlabel("Grain temperature (K)")
    ax.set_ylabel(r"H$_2$ formation efficiency $\epsilon$")
    ax.legend(frameon=False, loc="upper right")
    finalize(fig, "fig05_physisorption_only", outdir=outdir)


def make_fig06_tau_sensitivity(outdir: str) -> None:
    df = pd.read_csv(ROOT / "results/tables/grieco_tau_sensitivity_table.csv").dropna(subset=["tau"])
    expt = _load_digitized()
    fig, ax = plt.subplots(figsize=fig_single(2.4))

    for temp, color, marker, mean_col, ci_col in (
        (20, COLORS["blue"], "o", "eps20_isothermal_released_total_mean", "eps20_isothermal_released_total_ci95"),
        (200, COLORS["vermillion"], "s", "eps200_isothermal_mean", "eps200_isothermal_ci95"),
    ):
        exp_row = expt[np.isclose(expt["T_K"], float(temp))]
        if not exp_row.empty:
            y = float(exp_row["eps"].iloc[0])
            err = float(exp_row["eps_err"].iloc[0])
            ax.axhspan(y - err, y + err, color=color, alpha=0.07, lw=0)
        plot_with_ci(
            ax,
            df["tau"],
            df[mean_col],
            y_err=_safe_sem(df[ci_col]),
            color=color,
            marker=marker,
            linewidth=1.5,
            markevery=1,
        )

    ax.set_xlim(0.45, 0.95)
    ax.set_ylim(0.15, 0.32)
    ax.set_xticks(df["tau"].to_numpy(dtype=float))
    ax.set_xlabel(r"Dissociation fraction $\tau$")
    ax.set_ylabel(r"H$_2$ formation efficiency $\epsilon$")
    ax.text(0.905, float(df["eps20_isothermal_released_total_mean"].iloc[-1]), "20 K", color=COLORS["blue"], fontsize=7, ha="left", va="center", bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.2))
    ax.text(0.905, float(df["eps200_isothermal_mean"].iloc[-1]), "200 K", color=COLORS["vermillion"], fontsize=7, ha="left", va="center", bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.2))
    finalize(fig, "fig06_tau_sensitivity", outdir=outdir)


def make_fig08_mechanism_decomp(df: pd.DataFrame, outdir: str) -> None:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=fig_double(3.0), gridspec_kw={"wspace": 0.28})

    sub = df[np.isclose(df["nH"], 100.0)].sort_values("T_K")
    x_vals = sub["T_K"].to_numpy(dtype=float)
    lh_vals = sub["lh_fraction"].fillna(0.0).to_numpy(dtype=float)
    er_vals = sub["er_fraction"].fillna(0.0).to_numpy(dtype=float)
    ax1.stackplot(
        x_vals,
        lh_vals,
        er_vals,
        colors=[MECH_COLOR["LH"], MECH_COLOR["ER"]],
        alpha=0.85,
        edgecolor="white",
        linewidth=0.5,
    )
    ax1.axvline(110, color=COLORS["grey"], lw=0.8, ls="--")
    ax1.set_xlim(10, 270)
    ax1.set_ylim(0, 1)
    ax1.set_xlabel("Grain temperature (K)")
    ax1.set_ylabel("Fraction of formed H$_2$")
    ax1.text(55, 0.78, "LH", color="white", fontsize=9, fontweight="bold", ha="center", va="center")
    ax1.text(185, 0.80, "ER", color="white", fontsize=9, fontweight="bold", ha="center", va="center")
    panel_label(ax1, "a")

    for n_h in (10, 100, 1000, 10000):
        sub = df[np.isclose(df["nH"], float(n_h))].sort_values("T_K")
        mask = sub["T_K"].between(90, 150)
        ax2.plot(
            sub.loc[mask, "T_K"],
            100.0 * sub.loc[mask, "er_fraction"],
            color=DENSITY_COLOR[n_h],
            marker=DENSITY_MARKER[n_h],
            lw=1.5,
            ms=3.8,
            label=rf"$10^{{{int(np.log10(n_h))}}}$",
        )
    ax2.set_xlim(90, 150)
    ax2.set_ylim(0, 100)
    ax2.set_xlabel("Grain temperature (K)")
    ax2.set_ylabel("ER fraction (%)")
    ax2.legend(frameon=False, loc="lower right", fontsize=7)
    panel_label(ax2, "b")
    finalize(fig, "fig08_mechanism_decomp", outdir=outdir)


def make_fig09_surface_h(df: pd.DataFrame, outdir: str) -> None:
    sub = df[np.isclose(df["nH"], 100.0)].sort_values("T_K")
    fig, ax = plt.subplots(figsize=fig_single(2.7))
    plot_with_ci(
        ax,
        sub["T_K"],
        sub["surface_h_mean"],
        y_err=sub["surface_h_sem"],
        color=COLORS["blue"],
        marker="o",
        linewidth=1.8,
    )
    ax.axhline(212.0, color=COLORS["grey"], lw=0.8, ls=":")
    ax.set_xlim(0, 270)
    ax.set_yscale("log")
    ax.set_ylim(1, 500)
    ax.set_xlabel("Grain temperature (K)")
    ax.set_ylabel(r"Mean surface H count $\langle N_{\rm H}\rangle$")
    finalize(fig, "fig09_surface_h", outdir=outdir)


def make_fig10_transition_zoom(df: pd.DataFrame, outdir: str) -> None:
    fig, ax = plt.subplots(figsize=fig_single(2.8))
    label_positions = {}
    for n_h in (10, 100, 1000, 10000):
        sub = df[np.isclose(df["nH"], float(n_h))].sort_values("T_K")
        mask = sub["T_K"].between(80, 140)
        sub_plot = sub.loc[mask]
        plot_with_ci(
            ax,
            sub_plot["T_K"],
            sub_plot["eps_mean"],
            y_err=sub_plot["eps_sem"],
            color=DENSITY_COLOR[n_h],
            marker=DENSITY_MARKER[n_h],
            linewidth=1.5,
            markersize=3.8,
        )
        label_positions[n_h] = (float(sub_plot["T_K"].iloc[-1]), float(sub_plot["eps_mean"].iloc[-1]))

    low = float(df[np.isclose(df["nH"], 10.0) & np.isclose(df["T_K"], 100.0)]["eps_mean"].iloc[0])
    high = float(df[np.isclose(df["nH"], 10000.0) & np.isclose(df["T_K"], 100.0)]["eps_mean"].iloc[0])
    ax.axvline(100.0, color=COLORS["grey"], lw=0.8, ls="--")
    ax.annotate("", xy=(100.0, high), xytext=(100.0, low), arrowprops=dict(arrowstyle="<->", lw=0.9, color=COLORS["grey"]))
    ax.text(
        102.0,
        0.5 * (low + high),
        "16%",
        fontsize=7,
        color=COLORS["grey"],
        va="center",
        bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.2),
    )
    for n_h, (x_end, y_end) in label_positions.items():
        ax.text(
            x_end + 1.6,
            y_end,
            rf"$10^{{{int(np.log10(n_h))}}}$",
            color=DENSITY_COLOR[n_h],
            fontsize=7,
            ha="left",
            va="center",
            bbox=dict(fc="white", ec="none", alpha=0.9, pad=0.2),
        )
    ax.set_xlim(80, 140)
    ax.set_ylim(0.17, 0.30)
    ax.set_xlabel("Grain temperature (K)")
    ax.set_ylabel(r"H$_2$ formation efficiency $\epsilon$")
    finalize(fig, "fig10_transition_zoom", outdir=outdir)


def make_fig12_release_rate(df: pd.DataFrame, outdir: str) -> None:
    fig, ax = plt.subplots(figsize=fig_single(2.6))
    label_positions = {}
    for n_h in (10, 100, 1000, 10000):
        sub = df[np.isclose(df["nH"], float(n_h))].sort_values("T_K")
        plot_with_ci(
            ax,
            sub["T_K"],
            sub["rate_mean"],
            y_err=sub["rate_sem"],
            color=DENSITY_COLOR[n_h],
            marker=DENSITY_MARKER[n_h],
            linewidth=1.5,
            markersize=3.8,
        )
        label_positions[n_h] = (float(sub["T_K"].iloc[-1]), float(sub["rate_mean"].iloc[-1]))
    n_ref = 100.0
    sigma_h = 1.0e-21
    k_band = np.array([1.0e-17, 6.0e-17])
    rate_band = 0.5 * k_band * n_ref / sigma_h
    ax.axhspan(rate_band[0], rate_band[1], color=COLORS["grey"], alpha=0.11, lw=0)
    ax.set_yscale("log")
    ax.set_xlim(10, 250)
    ax.set_ylim(1e3, 1e8)
    ax.set_xlabel("Grain temperature (K)")
    ax.set_ylabel(r"$R({\rm H}_2)$ (cm$^{-2}$ s$^{-1}$)")
    for n_h, (_, y_end) in label_positions.items():
        ax.text(
            246,
            y_end,
            rf"$10^{{{int(np.log10(n_h))}}}$",
            color=DENSITY_COLOR[n_h],
            fontsize=7,
            ha="right",
            va="center",
            bbox=dict(fc="white", ec="none", alpha=0.9, pad=0.2),
        )
    ax.text(
        238,
        np.sqrt(rate_band[0] * rate_band[1]),
        "Jura (1975)\ncanonical range",
        color=COLORS["grey"],
        fontsize=6.8,
        ha="right",
        va="center",
        bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.2),
    )
    finalize(fig, "fig12_release_rate", outdir=outdir)


def make_fig13_ct10(df: pd.DataFrame, outdir: str) -> None:
    ct = pd.read_csv(ROOT / "results/tables/ct02_comparison_jhub_full.csv")
    ct = ct[np.isclose(ct["uv_flux_factor"], 0.0) & np.isclose(ct["h_gas_density_cm3"], 1000.0)].copy()
    ct["T_K"] = pd.to_numeric(ct["surface_temperature_k"], errors="coerce")
    ct["ct_rate"] = pd.to_numeric(ct["ct02_h2_release_rate_cm2_s"], errors="coerce")
    ct["ratio"] = pd.to_numeric(ct["ratio_kmc_over_ct02"], errors="coerce")
    ct = ct.dropna(subset=["T_K", "ct_rate", "ratio"]).sort_values("T_K")
    sub = df[np.isclose(df["nH"], 1000.0)].sort_values("T_K")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=fig_double(3.0), gridspec_kw={"wspace": 0.28})
    plot_with_ci(ax1, sub["T_K"], sub["rate_mean"], y_err=sub["rate_sem"], color=COLORS["blue"], marker="o", linewidth=1.8, label="KMC")
    ax1.plot(ct["T_K"], ct["ct_rate"], color=COLORS["vermillion"], ls="--", lw=1.6, label="CT10")
    ax1.set_yscale("log")
    ax1.set_xlim(10, 250)
    ax1.set_xlabel("Grain temperature (K)")
    ax1.set_ylabel(r"$R({\rm H}_2)$ (cm$^{-2}$ s$^{-1}$)")
    ax1.legend(frameon=False, loc="upper right")
    panel_label(ax1, "a")

    ax2.axhline(1.0, color=COLORS["grey"], lw=0.8, ls="--")
    ax2.fill_between(ct["T_K"], ct["ratio"], 1.0, where=ct["ratio"] < 1.0, color="#fde1de", alpha=0.45, lw=0)
    ax2.plot(ct["T_K"], ct["ratio"], color=COLORS["black"], lw=1.5)
    ax2.set_xlim(10, 250)
    ax2.set_ylim(0, 2)
    ax2.set_xlabel("Grain temperature (K)")
    ax2.set_ylabel(r"$R_{\rm KMC}/R_{\rm CT10}$")
    panel_label(ax2, "b")
    finalize(fig, "fig13_ct10_comparison", outdir=outdir)


def make_fig14_mrn(outdir: str) -> None:
    size_df = pd.read_csv(ROOT / "results/tables/astro_mrn_size_summary.csv")
    w_df = pd.read_csv(ROOT / "results/tables/astro_mrn_weights.csv")
    comp = pd.read_csv(ROOT / "results/tables/astro_mrn_comparison.csv")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=fig_single(4.5), gridspec_kw={"hspace": 0.18})
    for temp, color, marker in ((20, COLORS["blue"], "o"), (100, COLORS["purple"], "^"), (200, COLORS["orange"], "D")):
        sub = size_df[np.isclose(size_df["surface_temperature_k"], float(temp))].sort_values("grain_radius_um_mrn")
        if sub.empty:
            continue
        ax1.plot(sub["grain_radius_um_mrn"], sub["epsilon_mean"], color=color, marker=marker, ms=4.0, lw=1.5, label=f"{temp} K")
        ax1.text(
            float(sub["grain_radius_um_mrn"].iloc[-1]) * 0.90,
            float(sub["epsilon_mean"].iloc[-1]),
            f"{temp} K",
            color=color,
            fontsize=7,
            ha="right",
            va="center",
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.15),
        )
    ax1.set_xscale("log")
    ax1.set_xlim(0.005, 0.25)
    ax1.set_ylim(0, 0.4)
    ax1.set_ylabel(r"Efficiency $\epsilon$")
    panel_label(ax1, "a")

    warm = size_df[np.isclose(size_df["surface_temperature_k"], 100.0)].sort_values("grain_radius_um_mrn")
    merged = warm.merge(w_df, on="grain_radius_um_mrn", how="left").sort_values("grain_radius_um_mrn")
    merged["cum_weight"] = merged["mrn_weight_area"].cumsum()
    merged["cum_eps"] = (merged["epsilon_mean"] * merged["mrn_weight_area"]).cumsum() / merged["cum_weight"]
    ax2.plot(merged["grain_radius_um_mrn"], merged["mrn_weight_area"], color=COLORS["black"], lw=1.2, ls="--", label="MRN weight")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim(0.005, 0.25)
    ax2.set_ylabel("MRN weight")
    ax2.set_xlabel("Grain radius ($\\mu$m)")
    ax2.axvline(0.005, color=COLORS["grey"], lw=0.8, ls=":")
    ax2b = ax2.twinx()
    ax2b.plot(merged["grain_radius_um_mrn"], merged["cum_eps"], color=COLORS["blue"], lw=1.7, label="Integrated $\\epsilon$")
    single_100 = float(comp.loc[np.isclose(comp["surface_temperature_k"], 100.0), "epsilon_single"].iloc[0])
    ax2b.axhline(single_100, color=COLORS["vermillion"], lw=0.9, ls=":", label="Single-grain $\\epsilon$")
    ax2b.set_ylabel(r"Cumulative integrated $\epsilon$")
    ax2.text(0.006, 0.0285, "baseline", color=COLORS["grey"], fontsize=6.8, ha="left", va="bottom", bbox=dict(fc="white", ec="none", alpha=0.8, pad=0.15))
    ax2.text(0.19, float(merged["mrn_weight_area"].iloc[-1]) * 1.05, "MRN weight", color=COLORS["black"], fontsize=6.8, ha="right", va="bottom", bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.15))
    ax2b.text(0.19, float(merged["cum_eps"].iloc[-1]) - 0.0002, r"Integrated $\epsilon$", color=COLORS["blue"], fontsize=6.8, ha="right", va="top", bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.15))
    ax2b.text(0.19, single_100 + 0.00015, r"Single-grain $\epsilon$", color=COLORS["vermillion"], fontsize=6.8, ha="right", va="bottom", bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.15))
    panel_label(ax2, "b", y=0.88)
    finalize(fig, "fig14_mrn_integration", outdir=outdir)


def make_fig17_sensitivity_envelope(outdir: str) -> None:
    df = pd.read_csv(ROOT / "results/astro_sensitivity_knobs.csv").dropna(
        subset=["surface_temperature_k", "chemisorption_fraction", "er_reaction_probability", "h2_release_rate_cm2_s_mean", "h_gas_density_cm3"]
    )
    df["T_K"] = pd.to_numeric(df["surface_temperature_k"], errors="coerce")
    df["k_eff"] = 4.0e-21 * pd.to_numeric(df["h2_release_rate_cm2_s_mean"], errors="coerce") / pd.to_numeric(df["h_gas_density_cm3"], errors="coerce")
    env = df.groupby("T_K")["k_eff"].agg(["min", "max"]).reset_index().sort_values("T_K")

    fig, ax = plt.subplots(figsize=fig_single(2.8))
    ax.fill_between(env["T_K"], env["min"], env["max"], color=COLORS["blue"], alpha=0.08, lw=0)
    ax.plot(env["T_K"], env["min"], color=COLORS["grey"], lw=0.8, alpha=0.55)
    ax.plot(env["T_K"], env["max"], color=COLORS["grey"], lw=0.8, alpha=0.55)
    base_curve = None
    for (f_chem, p_er), sub in df.groupby(["chemisorption_fraction", "er_reaction_probability"]):
        sub = sub.sort_values("T_K")
        is_base = np.isclose(float(f_chem), 0.4) and np.isclose(float(p_er), 0.9)
        line, = ax.plot(
            sub["T_K"],
            sub["k_eff"],
            color=COLORS["blue"] if is_base else COLORS["grey"],
            lw=1.8 if is_base else 0.8,
            alpha=1.0 if is_base else 0.28,
            zorder=3 if is_base else 1,
        )
        if is_base:
            base_curve = (sub["T_K"].to_numpy(dtype=float), sub["k_eff"].to_numpy(dtype=float))
    if base_curve is not None:
        ax.text(
            base_curve[0][-1] - 2.5,
            base_curve[1][-1] * 1.03,
            "baseline",
            color=COLORS["blue"],
            fontsize=7,
            ha="right",
            va="bottom",
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.15),
        )
    ax.set_yscale("log")
    ax.set_xlabel("Grain temperature (K)")
    ax.set_ylabel(r"Effective $k_{\rm eff}$ (cm$^3$ s$^{-1}$)")
    finalize(fig, "fig17_sensitivity_envelope", outdir=outdir)


def make_fig18_transition_distributions(outdir: str) -> None:
    df = pd.read_csv(ROOT / "results/astro_transition_deep_raw.csv").dropna(
        subset=["surface_temperature_k", "h_gas_density_cm3", "epsilon"]
    )
    temps = [80, 90, 100, 110, 120]
    positions_low = [t - 0.15 for t in temps]
    positions_high = [t + 0.15 for t in temps]
    low = [df[np.isclose(df["surface_temperature_k"], t) & np.isclose(df["h_gas_density_cm3"], 10.0)]["epsilon"].to_numpy(dtype=float) for t in temps]
    high = [df[np.isclose(df["surface_temperature_k"], t) & np.isclose(df["h_gas_density_cm3"], 10000.0)]["epsilon"].to_numpy(dtype=float) for t in temps]

    fig, ax = plt.subplots(figsize=fig_single(2.8))
    common = dict(
        widths=0.25,
        patch_artist=True,
        notch=True,
        showfliers=True,
        flierprops=dict(marker=".", markersize=2, markeredgecolor=COLORS["grey"]),
        medianprops=dict(color="black", lw=1.0),
        whiskerprops=dict(lw=0.8),
        capprops=dict(lw=0.8),
    )
    b1 = ax.boxplot(low, positions=positions_low, **common)
    b2 = ax.boxplot(high, positions=positions_high, **common)
    for patch in b1["boxes"]:
        patch.set_facecolor(COLORS["sky"])
        patch.set_edgecolor(COLORS["sky"])
        patch.set_alpha(0.9)
    for patch in b2["boxes"]:
        patch.set_facecolor(COLORS["vermillion"])
        patch.set_edgecolor(COLORS["vermillion"])
        patch.set_alpha(0.9)
    ax.set_xticks(temps)
    ax.set_xlim(77, 123)
    ax.set_xlabel("Grain temperature (K)")
    ax.set_ylabel(r"Per-run efficiency $\epsilon$")
    ax.legend(
        [Patch(facecolor=COLORS["sky"]), Patch(facecolor=COLORS["vermillion"])],
        [r"$n_{\rm H}=10$", r"$n_{\rm H}=10^4$"],
        frameon=False,
        loc="upper right",
    )
    finalize(fig, "fig18_transition_distributions", outdir=outdir)


def make_fig19_lh_mode(outdir: str) -> None:
    df = pd.read_csv(ROOT / "results/astro_lh_mode_consistency.csv").dropna(
        subset=["surface_temperature_k", "lh_formation_mode", "epsilon_mean"]
    )
    pivot = df.pivot_table(index="surface_temperature_k", columns="lh_formation_mode", values="epsilon_mean", aggfunc="first").sort_index()
    pivot_ci = df.pivot_table(index="surface_temperature_k", columns="lh_formation_mode", values="epsilon_ci95", aggfunc="first").sort_index()
    temps = pivot.index.to_numpy(dtype=float)
    x = np.arange(len(temps))
    width = 0.35

    fig, ax = plt.subplots(figsize=fig_single(2.6))
    ax.bar(x - 0.2, pivot["pairs"], width=width, color=COLORS["blue"], yerr=pivot_ci["pairs"], capsize=2, label="pairs")
    ax.bar(x + 0.2, pivot["diffusion_limited"], width=width, color=COLORS["vermillion"], yerr=pivot_ci["diffusion_limited"], capsize=2, label="diffusion-limited")
    for i, temp in enumerate(temps):
        pct = 100.0 * (pivot.loc[temp, "diffusion_limited"] / pivot.loc[temp, "pairs"] - 1.0)
        ax.text(x[i], max(pivot.loc[temp, "pairs"], pivot.loc[temp, "diffusion_limited"]) + 0.008, f"{pct:.0f}%", ha="center", va="bottom", fontsize=6.8, color=COLORS["grey"])
    ax.set_xticks(x)
    ax.set_xticklabels([f"{int(t)}" for t in temps])
    ax.set_xlabel("Grain temperature (K)")
    ax.set_ylabel(r"H$_2$ formation efficiency $\epsilon$")
    ax.legend(frameon=False, loc="upper right")
    finalize(fig, "fig19_lh_mode", outdir=outdir)


def make_fig20_grain_size(outdir: str) -> None:
    df = pd.read_csv(ROOT / "results/astro_grain_size_check.csv").dropna(
        subset=["grain_radius_um", "surface_temperature_k", "epsilon_mean"]
    )
    fig, ax = plt.subplots(figsize=fig_single(2.4))
    for temp, color, marker in ((20.0, COLORS["blue"], "o"), (150.0, COLORS["vermillion"], "s")):
        sub = df[np.isclose(df["surface_temperature_k"], temp)].sort_values("grain_radius_um")
        if sub.empty:
            continue
        plot_with_ci(
            ax,
            sub["grain_radius_um"],
            sub["epsilon_mean"],
            y_err=_safe_sem(sub["epsilon_ci95"]),
            color=color,
            marker=marker,
            linewidth=1.5,
            label=f"{int(temp)} K",
        )
    ax.set_xscale("log")
    ax.set_ylim(0.0, 0.35)
    ax.set_xlabel("Grain radius ($\\mu$m)")
    ax.set_ylabel(r"H$_2$ formation efficiency $\epsilon$")
    ax.legend(frameon=False, loc="upper right")
    finalize(fig, "fig20_grain_size", outdir=outdir)


def make_fig21_porosity_sticking(outdir: str) -> None:
    por = pd.read_csv(ROOT / "results/astro_porosity_check.csv").dropna(
        subset=["porosity_fraction", "surface_temperature_k", "h_gas_density_cm3", "epsilon_mean"]
    )
    por = por[np.isclose(por["surface_temperature_k"], 100.0) & np.isclose(por["h_gas_density_cm3"], 100.0)].sort_values("porosity_fraction")
    sti = pd.read_csv(ROOT / "results/astro_sticking_model_check.csv").dropna(
        subset=["surface_temperature_k", "sticking_temp_model", "epsilon_mean"]
    )

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=fig_single(4.5), gridspec_kw={"hspace": 0.22})
    x = np.arange(len(por))
    ax1.bar(x, por["epsilon_mean"], width=0.6, color=[COLORS["blue"], COLORS["purple"]][: len(por)], yerr=por["epsilon_ci95"], capsize=2)
    ax1.set_xticks(x)
    ax1.set_xticklabels([f"{p:.1f}" for p in por["porosity_fraction"]])
    ax1.set_xlabel("Porosity fraction")
    ax1.set_ylabel(r"Efficiency $\epsilon$")
    panel_label(ax1, "a")

    for model, color, linestyle, label in (
        ("constant", COLORS["blue"], "-", "Constant $S=0.5$"),
        ("empirical_exp", COLORS["vermillion"], "--", "Temperature-dependent $S(T)$"),
    ):
        sub = sti[sti["sticking_temp_model"] == model].sort_values("surface_temperature_k")
        if sub.empty:
            continue
        plot_with_ci(
            ax2,
            sub["surface_temperature_k"],
            sub["epsilon_mean"],
            y_err=_safe_sem(sub["epsilon_ci95"]),
            color=color,
            marker="o" if model == "constant" else "s",
            linestyle=linestyle,
            linewidth=1.5,
            label=label,
        )
    ax2.set_xlabel("Grain temperature (K)")
    ax2.set_ylabel(r"Efficiency $\epsilon$")
    ax2.legend(frameon=False, loc="best")
    panel_label(ax2, "b")
    finalize(fig, "fig21_porosity_sticking", outdir=outdir)


def make_fig22_timescale_map(df: pd.DataFrame, outdir: str) -> None:
    sigma_h = 1.0e-21
    mu = 1.4
    g_cgs = 6.67430e-8
    m_h = 1.6735575e-24
    work = df.copy()
    work["t_h2_over_tff"] = (1.0 / (8.0 * sigma_h * work["rate_mean"])) / np.sqrt((3.0 * np.pi) / (32.0 * g_cgs * mu * m_h * work["nH"]))
    work["log_ratio"] = np.log10(work["t_h2_over_tff"])

    T_grid = np.linspace(10, 250, 240)
    n_grid = np.logspace(1, 4, 180)
    TT, NN = np.meshgrid(T_grid, n_grid)
    Z = griddata((work["T_K"], np.log10(work["nH"])), work["log_ratio"], (TT, np.log10(NN)), method="linear")

    fig, ax = plt.subplots(figsize=fig_single(3.0))
    mesh = ax.pcolormesh(T_grid, n_grid, Z, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-1, vcenter=0, vmax=2), shading="gouraud")
    contours = ax.contour(TT, NN, Z, levels=[0.0], colors="white", linewidths=1.5)
    ax.set_yscale("log")
    ax.set_xlim(10, 250)
    ax.set_xlabel("Grain temperature (K)")
    ax.set_ylabel(r"$n_{\rm H}$ (cm$^{-3}$)")
    marker_specs = [
        ("Diffuse cloud", 100.0, 30.0, "o"),
        ("Molecular cloud", 20.0, 1000.0, "s"),
        ("Dense clump", 30.0, 10000.0, "^"),
        ("z~6 clump", 60.0, 1000.0, "*"),
    ]
    legend_handles = []
    legend_labels = []
    for label, temp, dens, marker in marker_specs:
        handle = ax.plot(temp, dens, marker=marker, ms=6 if marker != "*" else 8, mec="black", mew=0.7, mfc="white", color="black", linestyle="none")[0]
        legend_handles.append(handle)
        legend_labels.append(label)
    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label(r"$\log_{10}(t_{\rm H_2}/t_{\rm ff})$")
    ax.legend(legend_handles, legend_labels, frameon=False, loc="lower right", fontsize=6.2, handletextpad=0.5)
    finalize(fig, "fig22_timescale_map", outdir=outdir)


def make_fig23_kinetic_trajectory(outdir: str) -> None:
    params = _load_yaml_params(ROOT / "config_astro_full_paperfit.yaml")
    params.update(
        {
            "surface_temperature_k": 100.0,
            "gas_temperature_k": 100.0,
            "h_gas_density_cm3": 1000.0,
            "uv_flux_factor": 0.0,
            "uv_mode": "continuous",
            "uv_pulse_enabled": False,
            "enable_grain_cache": True,
            "grain_cache_dir": "grain_cache_test",
            "rng_seed": 1000,
        }
    )

    kmc = KineticMonteCarlo(params)
    burnin_arrivals = 400
    measure_arrivals = 1200

    times = []
    surface_h = []
    formed_h2 = []
    adsorption_cum = []
    desorption_cum = []
    lh_cum = []
    er_cum = []
    burnin_end_time = {"value": None}

    def callback(sim: KineticMonteCarlo, event_name: str) -> None:
        times.append(float(sim.time))
        surface_h.append(float(sim.h_atoms_on_surface))
        formed_h2.append(float(sim.h2_molecules_formed))
        adsorption_cum.append(float(sim.total_adsorbed_h_atoms))
        desorption_cum.append(float(sim.total_desorbed_h_atoms))
        lh_cum.append(float(sim.h2_molecules_formed_LH))
        er_cum.append(float(sim.h2_molecules_formed_ER))
        if burnin_end_time["value"] is None and sim.total_arrivals >= burnin_arrivals:
            burnin_end_time["value"] = float(sim.time)

    kmc.simulation_parameters["max_arrivals"] = burnin_arrivals + measure_arrivals
    kmc.run_gillespie(max_time=1.0e12, max_steps=250000, callback=callback)

    if not times:
        raise RuntimeError("No trajectory recorded for Figure 23")

    times_arr = np.asarray(times, dtype=float)
    t_eval = np.logspace(np.log10(max(times_arr[times_arr > 0].min(), 1e-12)), np.log10(times_arr.max()), 180)

    def _interp_and_rate(values: list[float]) -> tuple[np.ndarray, np.ndarray]:
        arr = np.asarray(values, dtype=float)
        interp = np.interp(t_eval, times_arr, arr)
        rate = np.gradient(interp, t_eval)
        return interp, np.clip(rate, 0.0, None)

    surf_interp = np.interp(t_eval, times_arr, np.asarray(surface_h, dtype=float))
    formed_interp = np.interp(t_eval, times_arr, np.asarray(formed_h2, dtype=float))
    _, ads_rate = _interp_and_rate(adsorption_cum)
    _, des_rate = _interp_and_rate(desorption_cum)
    _, lh_rate = _interp_and_rate(lh_cum)
    _, er_rate = _interp_and_rate(er_cum)
    surf_smooth = gaussian_filter1d(surf_interp, sigma=2.2)
    formed_smooth = gaussian_filter1d(formed_interp, sigma=2.2)
    ads_smooth = gaussian_filter1d(ads_rate, sigma=2.0)
    des_smooth = gaussian_filter1d(des_rate, sigma=2.0)
    lh_smooth = gaussian_filter1d(lh_rate, sigma=2.0)
    er_smooth = gaussian_filter1d(er_rate, sigma=2.0)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=fig_single(3.2), sharex=True, gridspec_kw={"hspace": 0.12})
    burnin_end = burnin_end_time["value"] if burnin_end_time["value"] is not None else t_eval[len(t_eval) // 3]
    for axis in (ax1, ax2):
        axis.axvspan(t_eval.min(), burnin_end, color=COLORS["light_grey"], alpha=0.25, lw=0)
        axis.axvspan(burnin_end, t_eval.max(), color="#f6f6f6", alpha=0.35, lw=0)

    ax1.plot(t_eval, surf_interp, color=COLORS["blue"], lw=0.8, alpha=0.22)
    ax1.plot(t_eval, surf_smooth, color=COLORS["blue"], lw=1.7, label="Surface H")
    ax1.set_ylabel(r"$N_{\rm H}$ on surface")
    ax1_t = ax1.twinx()
    ax1_t.plot(t_eval, formed_interp, color=COLORS["vermillion"], lw=0.8, alpha=0.18)
    ax1_t.plot(t_eval, formed_smooth, color=COLORS["vermillion"], lw=1.7, label=r"Cumulative H$_2$ formed")
    ax1_t.set_ylabel(r"Cumulative H$_2$ formed")
    ax1.text(0.10, 0.88, "burn-in", transform=ax1.transAxes, fontsize=7, color=COLORS["grey"], style="italic")
    ax1.text(0.66, 0.88, "measured", transform=ax1.transAxes, fontsize=7, color=COLORS["grey"], style="italic")
    panel_label(ax1, "a")

    ax2.stackplot(
        t_eval,
        ads_smooth,
        des_smooth,
        lh_smooth,
        er_smooth,
        colors=[COLORS["blue"], COLORS["orange"], COLORS["purple"], COLORS["vermillion"]],
        alpha=0.82,
        labels=["Adsorption", "Desorption", "LH", "ER"],
    )
    ax2.set_ylabel("Event rate (s$^{-1}$)")
    ax2.set_xlabel("Simulation time (s)")
    ax2.legend(frameon=False, loc="upper left", ncol=2)
    ax2.set_xscale("log")
    panel_label(ax2, "b")
    finalize(fig, "fig23_kinetic_trajectory", outdir=outdir)


def _compute_chi2_red(df_kmc: pd.DataFrame, df_expt: pd.DataFrame) -> float:
    df_kmc = df_kmc.copy()
    df_expt = df_expt.copy()
    df_kmc["T_K"] = pd.to_numeric(df_kmc["T_K"], errors="coerce").astype(float)
    df_expt["T_K"] = pd.to_numeric(df_expt["T_K"], errors="coerce").astype(float)
    merged = pd.merge_asof(
        df_expt.sort_values("T_K"),
        df_kmc.sort_values("T_K")[["T_K", "eps_mean"]],
        on="T_K",
        direction="nearest",
    )
    warm = merged["T_K"].between(100.0, 250.0)
    dof = max(int(warm.sum()) - 1, 1)
    return float(np.sum(((merged.loc[warm, "eps_mean"] - merged.loc[warm, "eps"]) / merged.loc[warm, "eps_err"]) ** 2) / dof)


def write_report(outdir: Path, generated: list[str], notes: list[str]) -> None:
    lines = [
        "# Canonical MNRAS Figure Build",
        "",
        "Generated with `mnras_style.mplstyle` + `mnras_figures.py` helpers and saved under `results/plots/manuscript/`.",
        "",
        "## Generated",
    ]
    lines.extend([f"- `{name}`" for name in generated])
    lines.extend(["", "## Notes"])
    lines.extend([f"- {note}" for note in notes])
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build canonical publication-ready MNRAS figures.")
    parser.add_argument("--out-dir", default=str(OUTDIR))
    args = parser.parse_args()

    outdir = Path(args.out_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    _cleanup_legacy_outputs(outdir)

    generated: list[str] = []
    notes: list[str] = []

    baseline = _load_baseline()
    validation = _load_validation_curve()
    expt_iso = _load_digitized("isothermal")
    chi2_red = _compute_chi2_red(validation, expt_iso)

    make_fig01_render(str(outdir))
    generated.append("fig01_grain_lattice.pdf/.png")
    notes.append("Figure 1 is a hybrid setup figure: a schematic shell cutaway consistent with the manuscript geometry parameters plus a discrete slice from the actual cached grain lattice.")

    make_fig02_binding(str(outdir))
    generated.append("fig02_binding_energies.pdf/.png")

    _write_fig03_tikz()
    generated.append("fig03_mechanism_schematic.tex")
    notes.append("Figure 3 is emitted as TikZ source for manual Inkscape/LaTeX finishing rather than as a Matplotlib plot.")

    figure_grieco_validation(validation, expt_iso, chi2_red=chi2_red, outname="fig04_grieco_validation")
    generated.append("fig04_grieco_validation.pdf/.png")
    # Move outputs into manuscript folder if finalize used default directory.
    for ext in (".pdf", ".png", "_grey.png"):
        src = ROOT / "results" / "plots" / f"fig04_grieco_validation{ext}"
        if src.exists():
            src.replace(outdir / src.name)

    make_fig05_physisorption_only(str(outdir))
    generated.append("fig05_physisorption_only.pdf/.png")

    make_fig06_tau_sensitivity(str(outdir))
    generated.append("fig06_tau_sensitivity.pdf/.png")
    notes.append("Figure 6 uses the available 20 K and 200 K tau-sensitivity table; the local repo does not currently contain a matching 10/20/50 K set.")

    figure_epsilon_all_densities(baseline[["T_K", "nH", "eps_mean", "eps_sem"]], analytic_limit=0.18, outname="fig07_epsilon_all_densities")
    generated.append("fig07_epsilon_all_densities.pdf/.png")
    for ext in (".pdf", ".png", "_grey.png"):
        src = ROOT / "results" / "plots" / f"fig07_epsilon_all_densities{ext}"
        if src.exists():
            src.replace(outdir / src.name)

    make_fig08_mechanism_decomp(baseline, str(outdir))
    generated.append("fig08_mechanism_decomp.pdf/.png")

    make_fig09_surface_h(baseline, str(outdir))
    generated.append("fig09_surface_h.pdf/.png")

    make_fig10_transition_zoom(baseline, str(outdir))
    generated.append("fig10_transition_zoom.pdf/.png")

    figure_phase_map(baseline[["T_K", "nH", "eps_mean"]], outname="fig11_phase_map")
    generated.append("fig11_phase_map.pdf/.png")
    for ext in (".pdf", ".png", "_grey.png"):
        src = ROOT / "results" / "plots" / f"fig11_phase_map{ext}"
        if src.exists():
            src.replace(outdir / src.name)

    make_fig12_release_rate(baseline, str(outdir))
    generated.append("fig12_release_rate.pdf/.png")

    make_fig13_ct10(baseline, str(outdir))
    generated.append("fig13_ct10_comparison.pdf/.png")

    make_fig14_mrn(str(outdir))
    generated.append("fig14_mrn_integration.pdf/.png")
    notes.append("MRN figure styling is final-form, but the local warm-regime MRN provenance still deserves a science-level recheck before submission.")

    raw = pd.read_csv(ROOT / "results/astro_transition_deep_raw.csv")
    eps_raw = raw[np.isclose(raw["surface_temperature_k"], 100.0) & np.isclose(raw["h_gas_density_cm3"], 10000.0)]["epsilon"].to_numpy(dtype=float)
    figure_ensemble_convergence(eps_raw, outname="fig16_ensemble_convergence", production_N=20)
    generated.append("fig16_ensemble_convergence.pdf/.png")
    for ext in (".pdf", ".png", "_grey.png"):
        src = ROOT / "results" / "plots" / f"fig16_ensemble_convergence{ext}"
        if src.exists():
            src.replace(outdir / src.name)

    make_fig17_sensitivity_envelope(str(outdir))
    generated.append("fig17_sensitivity_envelope.pdf/.png")

    make_fig18_transition_distributions(str(outdir))
    generated.append("fig18_transition_distributions.pdf/.png")

    make_fig19_lh_mode(str(outdir))
    generated.append("fig19_lh_mode.pdf/.png")

    make_fig20_grain_size(str(outdir))
    generated.append("fig20_grain_size.pdf/.png")

    make_fig21_porosity_sticking(str(outdir))
    generated.append("fig21_porosity_sticking.pdf/.png")

    make_fig22_timescale_map(baseline, str(outdir))
    generated.append("fig22_timescale_map.pdf/.png")

    make_fig23_kinetic_trajectory(str(outdir))
    generated.append("fig23_kinetic_trajectory.pdf/.png")

    make_fig24_surface_energy_map(str(outdir))
    generated.append("fig24_surface_energy_map.pdf/.png")
    notes.append("Figure 24 is an additional support-style figure showing the top-layer site classes and binding-energy map from the actual cached grain, useful for methods or supplementary placement.")

    make_fig25_layer_gallery(str(outdir))
    generated.append("fig25_layer_gallery.pdf/.png")
    notes.append("Figure 25 shows how site classes are distributed across surface, mid-depth, and deep cached grain layers; it is a non-line support figure tied directly to the cached lattice.")

    notes.append("Figure 15 (UV suppression) was intentionally not rebuilt in this pass because the current manuscript direction is non-UV.")

    write_report(outdir, generated, notes)


if __name__ == "__main__":
    main()
