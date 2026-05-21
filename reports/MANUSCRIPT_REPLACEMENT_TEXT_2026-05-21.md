# Manuscript Replacement Text (Repo-Verified)

This note gives ready-to-paste wording that matches the current repository outputs checked on 2026-05-21.

Verified numerical anchors used below:

- Grieco isothermal calibration-set fit statistic: `chi2_red = 0.0249` from `results/tables/grieco_fig2_digitized_fit_stats.csv`
- Warm plateau mean over 100–250 K: `epsilon = 0.1946` from `results/tables/grieco_paperfit_validation_summary_table.csv`
- Representative solar-metallicity timescale case:
  - `T = 60 K`, `n_H = 10^3 cm^-3`, `t_H2 / t_ff = 0.933`
- Same case with metallicity scaling:
  - `Z' = 0.1`: `t_H2 / t_ff = 9.33`
  - `Z' = 0.01`: `t_H2 / t_ff = 93.3`
- LH mode offset from explicit-pairs vs diffusion-limited consistency check:
  - `20 K`: `-25.0%`
  - `50 K`: `-7.9%`
  - `80 K`: `-9.3%`

## Abstract replacement

Use this version if you want the safest honest framing:

> We present a kinetic Monte Carlo model for H2 formation on porous carbonaceous grains, calibrated against the Grieco et al. (2023) laboratory benchmark for D2 formation on coronene. The model reproduces the isothermal 100–250 K calibration data with `chi2_red = 0.025`, which we interpret as a calibration-set goodness-of-fit statistic rather than an independent validation. Across the full temperature range, the model yields a low-temperature classical suppression regime near 10 K, a peak Langmuir-Hinshelwood efficiency of `epsilon approx 0.28` at 20–80 K, and a warm Eley-Rideal plateau near `epsilon approx 0.19`, consistent with the analytic limit `epsilon_high-T = f_chem P_ER S` to which the simulation converges above about 150 K. A stochastic transition zone appears around 100–120 K, where the efficiency becomes density-dependent at the ~16% level. For a representative solar-metallicity dense-cloud case (`T = 60 K`, `n_H = 10^3 cm^-3`), we obtain `t_H2 / t_ff approx 0.93`; however, applying the same linear dust-metallicity scaling gives `t_H2 / t_ff approx 9.3` at `Z' = 0.1` and `approx 93` at `Z' = 0.01`, implying that grain-surface formation alone cannot keep pace with collapse in low-metallicity high-redshift gas.

## Section 3.1 calibration wording

Replace any "validation" phrasing with:

> The model reproduces the Grieco laboratory calibration dataset with `chi2_red = 0.025` over the isothermal 100–250 K points. Because key microphysical parameters were tuned against this same benchmark, we interpret `chi2_red` as a calibration-set goodness-of-fit statistic rather than an independent validation metric.

If you want one extra sentence:

> The full digitized comparison statistics are reported in `results/tables/grieco_fig2_digitized_fit_stats.csv`, for which the isothermal subset gives `chi2_red = 0.0249`.

## Warm-plateau framing

Use this in the results or discussion where the ER plateau is introduced:

> In the warm regime, the simulation approaches an intensive analytic limit rather than uncovering a new stochastic effect. Once the chemisorption reservoir is active, the efficiency converges to `epsilon_high-T = f_chem P_ER S`, giving a plateau near `epsilon approx 0.19`; the KMC therefore serves mainly as a consistency check on the approach to this limit above about 150 K.

## LH-regime systematic caveat

Use this wherever `epsilon approx 0.285` at 20 K is quoted:

> The quoted low-temperature peak should be interpreted together with the LH-mode systematic identified in the consistency test: the diffusion-limited ISM mode underestimates the explicit-pairs result by about 25% at 20 K and by about 8–9% at 50–80 K. Statistical confidence intervals therefore capture only part of the uncertainty in the LH-dominated regime.

Shorter version for a caption or footnote:

> Additional systematic uncertainty applies below 100 K from the LH encounter approximation: about 25% at 20 K and about 8–9% at 50–80 K.

## High-redshift discussion replacement

Replace any sentence claiming that the solar-metallicity `t_H2 / t_ff approx 0.93` result is directly applicable to `Z' = 0.01–0.1` high-redshift gas with:

> At solar metallicity, the representative `T = 60 K`, `n_H = 10^3 cm^-3` case gives `t_H2 / t_ff approx 0.93`, placing grain-surface formation near the collapse timescale. However, the model's own linear dust-metallicity scaling implies `t_H2 / t_ff approx 9.3` at `Z' = 0.1` and `approx 93` at `Z' = 0.01`. Dust-catalyzed H2 formation alone is therefore not fast enough to keep pace with gravitational collapse in low-metallicity high-redshift gas, and additional channels or more enriched environments are required.

## Conclusions replacement bullets

If your conclusion section is bulletized, these are safe replacements for the most exposed claims:

- The Grieco benchmark should be described as a calibration target rather than an independent validation set; the current model reproduces the isothermal 100–250 K laboratory points with `chi2_red = 0.025`.
- The warm `epsilon approx 0.19` plateau is consistent with the analytic ER limit `epsilon_high-T = f_chem P_ER S`, to which the KMC converges above about 150 K.
- The most genuinely stochastic result is the 100–120 K transition band, where the model predicts a density-dependent efficiency enhancement of order 16%.
- Low-temperature efficiencies below 100 K carry an additional LH-mode systematic uncertainty beyond the reported statistical confidence intervals, reaching about 25% at 20 K.
- The representative solar-metallicity dense-cloud case gives `t_H2 / t_ff approx 0.93`, but this rises to `9.33` at `Z' = 0.1` and `93.3` at `Z' = 0.01`, so the near-threshold result does not transfer directly to low-metallicity high-redshift gas.

## Parameter-degeneracy sentence

Use this near the analytic plateau discussion:

> The warm plateau constrains only the product `f_chem P_ER` (for fixed `S`), not the two parameters independently; the adopted values therefore represent one physically motivated choice within that degeneracy.

## 10 K caveat sentence

Use this in the abstract, conclusions, or a table footnote:

> The 10 K efficiency should be interpreted as a classical lower bound, since tunnelling is omitted and would tend to raise the cold-regime formation efficiency.

## What still needs manuscript-side editing

This repository does not contain the main paper `.tex` source, so the following still need to be applied in the manuscript itself:

- Replace `validation` with `calibration` where the Grieco dataset is the fitting target.
- Reframe any high-redshift discussion using the metallicity-scaled `t_H2 / t_ff` values above.
- Add the LH systematic caveat anywhere `epsilon approx 0.285` or low-temperature rate coefficients are quoted as headline values.
- Add the analytic-limit caveat wherever the warm plateau is presented as a major result.
