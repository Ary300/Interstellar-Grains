# Workspace Audit

Date: 2026-04-16

This audit summarizes what is currently strong in the workspace, what is risky,
and what appears to be missing before more science runs should be trusted.

## Executive Summary

The workspace already contains a strong paper-grade science package in saved
artifacts, especially for the Grieco validation and the main ISM baseline. The
biggest current risk is not "the science looks weak"; it is provenance and
freeze discipline.

The saved baseline artifacts support the story you described:

- Grieco paperfit validation is strong.
- The ISM baseline in `results/jhub_full_merged.csv` is coherent and nonzero.
- The LH to ER crossover is visible in the saved baseline.
- UV is structurally suppressed in the default ISM setup and should not be a
  main validated paper result in the current form.

The main workspace risk is that the repo is not in a clean, reproducible git
state:

- `git status` shows many modified and untracked files.
- `git log` shows only the initial commit, so the actual science evolution is
  effectively unversioned.
- Core files have drifted heavily from `HEAD`.
- Some support-run artifacts in the tree are clearly bad zero-output outputs,
  while others show the repaired nonzero behavior.

## What Looks Trusted

### 1. Grieco validation package

Existing saved tables support the validation story:

- `results/tables/grieco_paperfit_validation_summary_table.csv`
  - high-T plateau mean over 100-250 K: `0.194630`
  - DED released-total epsilon at 10 K: `0.211518 +/- 0.007809`
  - DED released-total epsilon at 20 K: `0.459920 +/- 0.008432`
  - DED 30-80 K mean: `0.298280 +/- 0.001295`
- `results/tables/grieco_paperfit_key_numbers.csv`
  - confirms the same paper-ready numbers
- `reports/PAPER_COMPLETION_STATUS.md`
  - cites reduced chi-square `0.0249` for the baseline isothermal fit

This is a strong validation anchor.

### 2. Main ISM baseline

The saved baseline in `results/jhub_full_merged.csv` is coherent and matches the
paper-facing reports.

For `n_H = 10^3 cm^-3`, `uv_flux_factor = 0`:

- `epsilon(20 K) = 0.28504`
- `epsilon(100 K) = 0.23963`
- `epsilon(150 K) = 0.189845`
- `epsilon(200 K) = 0.189625`
- `epsilon(250 K) = 0.189625`
- plateau mean over `150, 180, 200, 220, 250 K` is `0.189669`

Mechanism decomposition in the same saved baseline is also clean:

- `20 K`: LH dominates
  - `h2_formed_LH_mean = 2796.55`
  - `h2_formed_ER_mean = 53.85`
  - LH fraction about `98.1%`
- `100 K`: still mostly LH in the saved baseline point
  - `h2_formed_LH_mean = 2277.6`
  - `h2_formed_ER_mean = 118.7`
- `150 K`: essentially ER dominated
  - `h2_formed_LH_mean = 0.55`
  - `h2_formed_ER_mean = 1897.9`
- `200 K` and `250 K`: fully ER dominated in practice

That supports the warm-regime story and a crossover in roughly the quoted
100-120 K neighborhood.

### 3. Timescale result

`results/tables/table_timescales_representative.csv` contains the headline
timescale result:

- `T = 60 K`, `n_H = 10^3 cm^-3`: `t_H2 / t_ff = 0.933198`

This agrees with the paper-facing summary.

### 4. Other saved support pieces that already exist

The workspace already contains these useful artifacts:

- `results/tables/ct02_comparison_jhub_full.csv`
- `results/tables/grieco_tau_sensitivity_table.csv`
- `results/grieco_physisorption_only_iso.csv`
- `results/tables/astro_mrn_comparison.csv`
- `results/plots/paper/grieco_validation_overlay.png`
- `results/plots/paper/ism_main_results.png`
- `results/plots/paper_suite/` figures for epsilon, release rate, mechanism
  crossover, CT10 ratio, and UV suppression

## What Looks Risky

### 1. Provenance is the main technical problem

The workspace is not frozen:

- many tracked files are modified
- many project files are untracked
- core simulation files differ strongly from `HEAD`
- there is no useful commit history beyond the initial commit

Direct checks:

- `kmc_simulation.py`: current 1102 lines vs `HEAD` 469 lines
- `run_sweep.py`: current 492 lines vs `HEAD` 415 lines
- `physical_rates.py`: current 308 lines vs `HEAD` 275 lines
- `scientific_data.py`: current 76 lines vs `HEAD` 71 lines

That means the "trusted baseline" is represented mainly by saved CSV/plot/table
artifacts, not by a clean git snapshot.

### 2. Mixed good and bad artifacts coexist

The workspace contains explicit evidence of earlier zero-output problems:

- `results/tmp_restore_sanity_seed42.csv` gives `epsilon_mean = 0.0`
- `astro_full_paperfit_uv_only.csv` contains many zero rows
- older lightweight files like `results/quick_test.csv` are all-zero and are not
  trustworthy science outputs

But the workspace also contains repaired nonzero outputs:

- `results/tmp_restore_sanity.csv` gives
  - `epsilon_mean = 0.235`
  - `h2_formed_LH_mean = 220`
  - `h2_formed_ER_mean = 15`
- `results/astro_pilot_paperfit.csv` is nonzero and physically coherent
- `reports/PAPER_COMPLETION_STATUS.md` explicitly says the zero-output regression
  in `run_sweep.py` was fixed

Conclusion:

- the workspace is not uniformly bad
- but the artifact layer is mixed enough that careful curation is required

### 3. Some reports reference outputs that are not present locally

`reports/PAPER_COMPLETION_STATUS.md` cites files like:

- `results/plots/astro_lh_mode_consistency/lh_mode_consistency_summary.csv`
- `results/astro_porosity_check.csv`
- `results/astro_sticking_model_check.csv`
- `results/astro_transport_sensitivity.csv`

Those files are not currently present in the workspace. So the narrative report
is ahead of the local artifact set.

This does not invalidate the main baseline science, but it does mean the local
workspace is incomplete as a referee-support package.

## UV Status

The code and configs directly confirm the structural UV issue from the prior
session.

### What the ISM configs do

The main ISM configs use:

- `diffusion_mode: rate_only`
- `lh_formation_mode: diffusion_limited`

This is true in:

- `config_astro_full_paperfit.yaml`
- `config_astro_transition_deep.yaml`
- `config_astro_sensitivity_knobs.yaml`
- `config_astro_grain_size_check.yaml`
- `config_astro_porosity_check.yaml`
- `config_astro_sticking_model_check.yaml`

### What the default UV channel does

The default UV H2 logic in `kmc_simulation.py` still depends on explicit
adjacent H pairs:

- `uv_h2_mode` default is `adjacent_pair`
- UV H2 formation uses `adjacent_h_pairs_count`

So in low-coverage ISM runs with `rate_only` diffusion, `h2_formed_UV_mean = 0`
is structurally expected rather than surprising.

### Exploratory photofragmentation branch

There is now an exploratory photofragmentation path:

- `uv_h2_mode: chemisorption_photofrag`
- related parameters exist in `scientific_data.py`
- `config_astro_uv_pilot_photofrag.yaml` is present

But no trustworthy validated photofragmentation science package is present in
the local artifacts, and this branch should still be treated as exploratory.

## Run Coverage Status

The "must-run" and "quick-check" configs mostly exist, but their expected local
outputs are not present.

Missing local outputs:

- `config_astro_transition_deep.yaml`
  - `results/astro_transition_deep.csv`
  - `results/astro_transition_deep_raw.csv`
- `config_astro_sensitivity_knobs.yaml`
  - `results/astro_sensitivity_knobs.csv`
- `config_astro_lh_mode_consistency.yaml`
  - `results/astro_lh_mode_consistency.csv`
- `config_astro_porosity_check.yaml`
  - `results/astro_porosity_check.csv`
- `config_astro_grain_size_check.yaml`
  - `results/astro_grain_size_check.csv`
- `config_astro_sticking_model_check.yaml`
  - `results/astro_sticking_model_check.csv`
- `config_astro_transport_sensitivity.yaml`
  - `results/astro_transport_sensitivity.csv`
- `config_astro_uv_pilot_baseline.yaml`
  - `results/astro_uv_pilot_baseline.csv`
- `config_astro_uv_pilot_photofrag.yaml`
  - `results/astro_uv_pilot_photofrag.csv`

Partially present:

- `config_astro_mrn_integration.yaml`
  - aggregated result exists as `results/astro_mrn_integration.csv`
  - the exact configured raw/singlegrain outputs are not present

## Script Readiness

The figure-support tooling is in decent shape.

Ready-to-use scripts include:

- `plot_ensemble_convergence.py`
- `plot_transition_distributions.py`
- `plot_sensitivity_knobs.py`
- `plot_referee_quick_checks.py`

Notably, `plot_ensemble_convergence.py` already does the right kind of analysis
for the requested convergence figure:

- running mean vs number of runs
- SEM envelope
- selectable condition from raw transition-deep output

So the convergence-figure implementation gap is basically not "missing code"
anymore. The missing piece is the raw transition-deep dataset.

## Repository Health

### Tests

The current test suite passes:

- `pytest -q -k 'not slow'`
- result: `19 passed in 25.19s`

That is encouraging, but it is not enough to establish scientific provenance.

### Snapshot clue

`kmc_simulation_current_anvil.py` exists and is byte-for-byte identical to the
current `kmc_simulation.py`.

That suggests the current working copy may reflect an Anvil-era branch snapshot,
but because it is uncommitted, it still is not a reliable frozen provenance
anchor by itself.

## Best Current Interpretation

The project is in a better position scientifically than operationally.

Scientifically:

- the baseline paper story is strong
- Grieco validation is publishable
- the warm-regime ISM plateau is supported
- the timescale result is present
- the LH to ER transition is visible
- UV is not ready as a main result

Operationally:

- the repo is not clean
- artifact provenance is mixed
- some local reports are ahead of the available files
- the referee-proofing support package is incomplete in this checkout

## Recommended Next Actions

1. Freeze the science state before launching more runs.
   - At minimum, duplicate the current core simulation files plus key configs
     into a clearly named frozen snapshot directory or a real git commit/tag.

2. Treat these as the primary trusted science anchors for writing:
   - `results/jhub_full_merged.csv`
   - `results/tables/grieco_paperfit_validation_summary_table.csv`
   - `results/tables/grieco_paperfit_key_numbers.csv`
   - `results/tables/table_timescales_representative.csv`
   - `results/tables/ct02_comparison_jhub_full.csv`
   - `results/tables/astro_mrn_comparison.csv`

3. Do not promote UV into the main paper result set.
   - Keep it as negligible or exploratory unless a real validation/calibration
     campaign is done.

4. Rebuild the referee-support layer from the frozen baseline state.
   - transition deep
   - sensitivity knobs
   - LH consistency
   - porosity, grain size, and sticking quick checks

5. Use the saved baseline artifacts, not miscellaneous older zeroed files, when
   drafting text and figures.

