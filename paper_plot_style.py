from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


SINGLE_COLUMN_WIDTH = 10.0 / 3.0
DOUBLE_COLUMN_WIDTH = 7.0

OKABE_ITO = {
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "green": "#009E73",
    "purple": "#CC79A7",
    "yellow": "#F0E442",
    "sky": "#56B4E9",
    "orange": "#E69F00",
    "black": "#000000",
}

DENSITY_STYLES = {
    10.0: {"color": OKABE_ITO["sky"], "marker": "o", "linestyle": "-"},
    100.0: {"color": OKABE_ITO["blue"], "marker": "s", "linestyle": "-"},
    1000.0: {"color": OKABE_ITO["purple"], "marker": "^", "linestyle": "-"},
    10000.0: {"color": OKABE_ITO["vermillion"], "marker": "D", "linestyle": "-"},
}

MECHANISM_COLORS = {
    "LH": OKABE_ITO["blue"],
    "ER": OKABE_ITO["vermillion"],
    "UV": "#7f7f7f",
}


def apply_publication_style() -> None:
    style_path = Path(__file__).resolve().with_name("mnras.mplstyle")
    plt.style.use(str(style_path))
    mpl.rcParams.update(
        {
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.transparent": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def style_axes(ax: plt.Axes) -> None:
    ax.set_facecolor("white")
    for spine in ("left", "bottom", "top", "right"):
        ax.spines[spine].set_color("black")
        ax.spines[spine].set_linewidth(0.7)
    ax.minorticks_on()
    ax.tick_params(which="major", direction="in", top=True, right=True, length=3.5, width=0.7, colors="black")
    ax.tick_params(which="minor", direction="in", top=True, right=True, length=2.0, width=0.5, colors="black")


def density_color(n_h: float) -> str:
    return density_style(n_h)["color"]


def density_style(n_h: float) -> dict[str, str]:
    return DENSITY_STYLES.get(float(n_h), {"color": OKABE_ITO["black"], "marker": "o", "linestyle": "-"})


def figure_size(columns: str, height: float) -> tuple[float, float]:
    width = SINGLE_COLUMN_WIDTH if columns == "single" else DOUBLE_COLUMN_WIDTH
    return (width, height)


def add_regime_band(
    ax: plt.Axes,
    x0: float,
    x1: float,
    label: str,
    *,
    color: str = "#e9c46a",
    alpha: float = 0.08,
    y: float = 0.975,
) -> None:
    ax.axvspan(x0, x1, color=color, alpha=alpha, zorder=0)
    ax.text(
        0.5 * (x0 + x1),
        y,
        label,
        transform=ax.get_xaxis_transform(),
        ha="center",
        va="top",
        fontsize=6.8,
        color="#5b5b5b",
    )


def add_panel_labels(axes: Iterable[plt.Axes], labels: str = "ABCDEFGHIJKLMNOPQRSTUVWXYZ") -> None:
    for ax, label in zip(axes, labels):
        ax.text(
            0.03,
            0.94,
            f"({label.lower()})",
            transform=ax.transAxes,
            fontsize=10,
            fontweight="bold",
            va="top",
            ha="left",
            color="#202020",
        )


def style_legend(ax: plt.Axes, **kwargs):
    defaults = {
        "frameon": False,
        "borderpad": 0.2,
        "handletextpad": 0.5,
        "labelspacing": 0.35,
        "columnspacing": 0.9,
    }
    defaults.update(kwargs)
    return ax.legend(**defaults)


def save_figure(fig: plt.Figure, out_base: str) -> None:
    root, ext = os.path.splitext(out_base)
    base = root if ext else out_base
    os.makedirs(os.path.dirname(base) or ".", exist_ok=True)
    metadata = {"Creator": "Codex manuscript figure pipeline"}
    fig.savefig(base + ".png", bbox_inches="tight", metadata=metadata)
    fig.savefig(base + ".pdf", bbox_inches="tight", metadata=metadata)


def interpolate_phase_map(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    *,
    nx: int = 240,
    ny: int = 220,
    log_y: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return a denser grid for publication-style phase maps.

    The interpolation is only for visual smoothing of a regularly sampled grid.
    If SciPy is unavailable, fall back to the original coarse grid.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    y_work = np.log10(y) if log_y else y
    try:
        from scipy.interpolate import RegularGridInterpolator
    except Exception:
        return x, y, z

    interp = RegularGridInterpolator(
        (y_work, x),
        z,
        method="linear",
        bounds_error=False,
        fill_value=np.nan,
    )
    x_dense = np.linspace(float(x.min()), float(x.max()), int(nx))
    y_dense_work = np.linspace(float(y_work.min()), float(y_work.max()), int(ny))
    xx, yy = np.meshgrid(x_dense, y_dense_work)
    z_dense = interp(np.column_stack([yy.ravel(), xx.ravel()])).reshape(yy.shape)
    y_dense = 10 ** y_dense_work if log_y else y_dense_work
    return x_dense, y_dense, z_dense
