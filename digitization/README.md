# Digitizing Grieco et al. (2023) Fig. 2 (manual workflow)

Goal: create a simple CSV of experimental points from Grieco et al. (2023) Fig. 2 so we can compute
residuals and χ² against the simulation (“paperfit”) curve.

## 1) What to digitize (recommended minimal set)

Digitize **only the points we directly simulate in the paperfit baseline**:

- **TPDED circles**: points in the **10–80 K** range.
- **Isothermal triangles**: points in the **100–250 K** range.

This avoids ambiguity from any additional isothermal points in the low‑T regime.

## 2) Use WebPlotDigitizer (fastest)

1. Open WebPlotDigitizer (browser tool).
2. Load the Fig. 2 image (screenshot or PDF extract).
3. Calibrate axes (2D, X‑Y):
   - X axis: `T (K)` (linear)
   - Y axis: recombination efficiency `ε` shown as a **percentage** in the plot.
4. Digitize points:
   - For each point, record `T_K` and `epsilon_percent`.
   - Also digitize error bars if you can (upper/lower).

Important:
- Convert to **fraction** before saving: `epsilon = epsilon_percent / 100`.
- The paper caption indicates the error bars are **1σ SEM**. Treat digitized error bars as **σ**.

## 3) Save as CSV

Copy `digitization/grieco_fig2_digitized_template.csv` to one of:

- `grieco_fig2_digitized.csv` (repo root), or
- `results/grieco_fig2_digitized.csv`

Fill one row per point, with columns:

- `T_K` (float)
- `epsilon` (fraction, 0–1)
- `epsilon_err_low` (σ, fraction)
- `epsilon_err_high` (σ, fraction)
- `method` (`TPDED` or `isothermal`)

## 4) Compute residuals + χ² + plots

Run:

```bash
python3 compare_grieco_fig2_digitized.py
```

Outputs:

- `results/tables/grieco_fig2_digitized_residuals.csv`
- `results/tables/grieco_fig2_digitized_fit_stats.csv`
- `results/plots/grieco_fig2_overlay_digitized.png`
- `results/plots/grieco_fig2_residuals_digitized.png`

