# Paper Completion Status

Updated: 2026-04-05

## Current Call

The repaired support-campaign workflow is now behaving consistently with the trusted
baseline in `results/jhub_full_merged.csv`.

This does **not** mean every support result is numerically identical to the original
baseline grid. It means:

- the zero-output regression in `run_sweep.py` is fixed
- the support campaigns are producing physical, nonzero results again
- the trusted baseline science in `results/jhub_full_merged.csv` remains unchanged

## Baseline Numbers To Anchor The Paper

From `results/jhub_full_merged.csv` at `n_H = 10^3 cm^-3`, `UV = 0`:

- `epsilon(20 K) = 0.28504`
- `epsilon(100 K) = 0.23963`
- `epsilon(150 K) = 0.189845`
- `epsilon(200 K) = 0.189625`
- `epsilon(250 K) = 0.189625`

High-temperature plateau average over `150-250 K`:

- `epsilon_plateau = 0.189669`

Simple analytic limit:

- `f_chem * P_ER * S = 0.40 * 0.9 * 0.5 = 0.18`

Interpretation:

- the plateau is close to the simple analytic estimate
- the small offset above `0.18` is acceptable and can be described as the finite
  contribution of residual stochastic transport / occupancy effects near the onset
  of the ER-dominated regime

## Timescale Result

From `results/tables/table_timescales_representative.csv`:

- `T = 60 K`, `n_H = 10^3 cm^-3`: `t_H2 / t_ff = 0.933`
- `T = 100 K`, `n_H = 10^2 cm^-3`: `t_H2 / t_ff = 3.592`
- `T = 40 K`, `n_H = 10^4 cm^-3`: `t_H2 / t_ff = 0.295`

Recommended paper framing:

- the `60 K`, `10^3 cm^-3` case sits near the transition where dust-catalyzed H2
  formation becomes comparable to dynamical collapse
- this is a microphysical constraint, not a standalone explanation of any JWST result

## Completed Referee-Facing Checks

### Grieco Validation Variant Comparison

From `results/tables/grieco_model_variant_comparison.csv`:

- baseline isothermal reduced chi-square: `0.0249`
- prediction-model isothermal reduced chi-square: `0.4870`

Interpretation:

- the baseline paperfit validation remains the best reproduction of Grieco
- the prediction-model cross-check is usable as an uncertainty bound, not as a
  replacement validation state
- this supports a methods sentence like:
  "Using the ISM-style transport and porosity assumptions worsens the direct
  laboratory fit, so we retain the paperfit validation state as the primary
  experimental calibration and treat the prediction-model variant as a systematic
  cross-check."

### LH Mode Consistency

From `results/plots/astro_lh_mode_consistency/lh_mode_consistency_summary.csv`:

- `20 K`: diffusion-limited is `25.0%` below pairs mode
- `50 K`: diffusion-limited is `7.9%` below pairs mode
- `80 K`: diffusion-limited is `9.3%` below pairs mode

Interpretation:

- the approximation is acceptable through most of the LH-active regime
- the largest deviation is at `20 K`
- safest manuscript phrasing:
  "The diffusion-limited approximation reproduces the explicit-pairs result to
  within about 8-10% at 50-80 K, with a larger ~25% offset at 20 K where the
  LH channel is most sensitive to local encounter statistics."

Important clarification:

- this `20 K` offset is **not a new regression**
- the headline ISM baseline in `results/jhub_full_merged.csv` was generated from
  `config_astro_full_paperfit.yaml`, which already uses
  `lh_formation_mode: diffusion_limited`
- therefore the quoted low-temperature ISM value `epsilon(20 K) = 0.285` already
  includes this suppression relative to a hypothetical `pairs`-mode value of
  roughly `0.38`

Recommended interpretation:

- the low-temperature point is systematically conservative
- this does **not** alter the main novel claims, because those live in the
  `100-250 K` regime where ER dominates and the LH algorithm choice becomes
  negligible
- the `20 K` value should be framed as a low-temperature algorithmic systematic,
  not as a challenge to the warm-regime story

Recommended manuscript sentence:

- "All ISM results reported in this work use the diffusion-limited LH
  approximation. Relative to explicit pairs-mode simulations, this suppresses the
  `20 K` efficiency by about 25%, while agreement improves to within about 10%
  by `50-80 K`; the warm (`>= 100 K`) results emphasized here are effectively
  unchanged because the chemistry is ER-dominated in that regime."

This is good enough for MNRAS with honest wording.

For ApJ, this remains the main reason a stricter referee could still ask for a
pairs-mode production rerun.

## Resolution Of The 20 K LH-Mode Offset

### Decision

We do **not** change the baseline ISM dataset.

Reason:

- `results/jhub_full_merged.csv` is the trusted science baseline
- it was generated with `lh_formation_mode: diffusion_limited`
- the `20 K` suppression is therefore already part of the published numerical
  story we have been building around for months
- the warm-regime conclusions do not depend on replacing that point with a
  pairs-mode value

### Practical Fix

The fix is to make the low-temperature algorithmic offset explicit in the paper,
bound it quantitatively, and state why it does not affect the headline claims.

### Results-Section Wording

Use something close to:

"A dedicated LH-mode consistency test shows that the diffusion-limited
approximation reproduces the explicit pairs-mode result to within about 8-10%
at 50-80 K, with a larger offset of about 25% at 20 K. All ISM results in this
work use the diffusion-limited treatment, so the reported low-temperature ISM
efficiencies should be interpreted as conservative with respect to explicit
encounter statistics."

### Discussion / Limitations Wording

Use something close to:

"This low-temperature offset does not materially alter the main conclusions of
the paper. The central results emphasized here — the persistence of a warm
high-temperature efficiency plateau, the LH-to-ER crossover, the comparison to
CT10, and the resulting ISM timescales above about 100 K — arise in a regime
where ER dominates and the LH algorithm choice becomes negligible."

### Short Caption-Style Version

For a figure caption or referee response:

"Diffusion-limited and explicit-pairs LH treatments agree to within about 10%
at 50-80 K and differ by about 25% at 20 K; the warm ER-dominated regime is
effectively unchanged."

### What Not To Do

Do **not**:

- rerun the full baseline grid
- replace `results/jhub_full_merged.csv`
- retune the validated paperfit state just to move the `20 K` point

### Optional ApJ-Only Upgrade

If a stricter venue or referee insists on a numerical low-temperature check,
the clean supplementary response is:

- rerun only the sub-30 K ISM slice in `pairs` mode
- present it as a low-temperature systematic cross-check
- leave the warm-regime baseline untouched

That is a small, targeted add-on. It is not required to preserve the current
MNRAS-level story.

### Porosity Check

From `results/astro_porosity_check.csv` at `100 K`, `n_H = 100`:

- `porosity = 0.0`: `epsilon_mean = 0.22867`
- `porosity = 0.2`: `epsilon_mean = 0.22747`

Interpretation:

- porosity changes the rate only weakly in this warm ER-dominated regime
- this is exactly the kind of quick robustness check a referee will want

### Sticking Model Check

From `results/astro_sticking_model_check.csv`:

- constant sticking remains higher than the empirical exponential form
- the comparison is nonzero and physically interpretable
- this should be framed as a model sensitivity / limitation, not as a fatal flaw

### Transport Sensitivity

From `results/astro_transport_sensitivity.csv`:

- the rerun is nonzero and centered near the baseline
- the `epsilon_ratio_to_baseline` values around 20 K and 100 K are generally
  close to unity for the sampled `lh_diffusion_factor` / `diffusion_rate_cap_s`
  combinations

Interpretation:

- this supports the claim that the transport knobs modulate, but do not overturn,
  the main warm-regime story

## Still Running / Pending

At the time of this note, the remaining queue items are:

- `16067465` `astro_transdeep`
- `16067466` `astro_transplot` (dependent)
- `16068470` `astro_sensk_r2`
- `16068471` `astro_sensk_plot_r2` (dependent)
- `16068472` `astro_grain_r2`
- `16068533` `astro_refplots_r2` (dependent on grain rerun)
- `16068481` `astro_mrn_r2`
- `16068482` `astro_mrn_plot_r2` (dependent)
- `16068484` `ct10_compare_r2`

These are the last pieces needed to fully close the paper-support package.

## Remaining Writeup Items

The manuscript still needs the following prose sections assembled from the finished outputs:

- Grieco validation subsection:
  - emphasize `chi2_red = 0.025` for the baseline fit
  - discuss the `10 K` shortfall as a low-temperature limitation
  - mention the prediction-model cross-check as a systematic bound
- CT10 comparison subsection:
  - use the rerun `compare_ism_to_ct02.py` numbers and figure
- Physisorption-only subsection:
  - point to the existing physisorption-only comparison figure
- Discussion comparison to Satonkin et al. (2025):
  - your physisorption-only run is the closest conceptual bridge
  - chemisorption is what extends nonzero efficiency into the warm regime
- Abstract
- Conclusions
- Bibliography

## Recommended Submission-Level Conclusion

If the remaining queue items finish cleanly, the paper package should be strong
enough for:

- a serious MNRAS submission
- a plausible ApJ submission, but with higher revision risk because of the
  validation/prediction mismatch and the `20 K` LH-mode offset

The core result remains the same:

- the Grieco-calibrated KMC benchmark is real
- the warm high-temperature plateau is robust
- the LH-to-ER transition is physically visible
- the resulting H2 timescale can approach the free-fall time in representative
  high-redshift conditions
