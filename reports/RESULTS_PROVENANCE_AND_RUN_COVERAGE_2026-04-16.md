# Results Provenance And Run Coverage

Date: 2026-04-16

This is a compact inventory of which artifacts currently exist locally, which
appear trustworthy, and which support campaigns still need reruns or recovery.

## Trusted-First Artifacts

These are the most useful current anchors for manuscript writing.

### Core science

- `results/jhub_full_merged.csv`
- `results/tables/table_campaign_overview_jhub_full.csv`
- `results/tables/table_top25_jhub_full.csv`
- `results/tables/table_timescales_representative.csv`
- `results/tables/ism_timescales_scenario.csv`
- `results/tables/ct02_comparison_jhub_full.csv`
- `results/tables/astro_mrn_comparison.csv`
- `results/astro_mrn_integration.csv`

### Grieco validation

- `results/grieco_validation_paper_iso_paperfit.csv`
- `results/grieco_ded_paper_ded_paperfit.csv`
- `results/grieco_ded_paper_ded_paperfit_summary.json`
- `results/tables/grieco_paperfit_validation_summary_table.csv`
- `results/tables/grieco_paperfit_key_numbers.csv`
- `results/tables/grieco_paperfit_parameter_table.csv`
- `results/tables/grieco_fig2_digitized_fit_stats.csv`
- `results/tables/grieco_fig2_digitized_residuals.csv`

### Supporting science already present

- `results/grieco_physisorption_only_iso.csv`
- `results/tables/grieco_tau_sensitivity_table.csv`
- `results/tables/grieco_validation_summary_table.csv`
- `results/tables/table_uv_summary_jhub_full.csv`
- `results/tables/table_uv_suppression_jhub_full.csv`

### Existing paper plots

- `results/plots/paper/grieco_validation_overlay.png`
- `results/plots/paper/ism_main_results.png`
- `results/plots/paper_suite/epsilon_vs_temperature.png`
- `results/plots/paper_suite/release_rate_vs_temperature.png`
- `results/plots/paper_suite/mechanism_crossover.png`
- `results/plots/paper_suite/ct10_ratio_vs_temperature.png`

## Existing But Use Carefully

These files are informative but should not be treated as baseline anchors without
context.

- `results/tmp_restore_sanity.csv`
  - useful as evidence that repaired code can produce nonzero outputs
- `results/astro_pilot_paperfit.csv`
  - useful as a coherent smoke or pilot check
- `astro_full_paperfit_uv_only.csv`
  - useful mainly as evidence of UV suppression or mixed exploratory behavior

## Files That Look Like Old Or Bad Outputs

These should not be used as science anchors.

- `results/quick_test.csv`
- `results/raw_runs.csv`
- `results/tmp_restore_sanity_seed42.csv`

Reason:

- they contain all-zero or obviously non-production outputs
- they are inconsistent with the trusted baseline and repaired support-story

## Support Campaign Coverage

### Configs present, expected local outputs missing

- `config_astro_transition_deep.yaml`
  - missing `results/astro_transition_deep.csv`
  - missing `results/astro_transition_deep_raw.csv`
- `config_astro_sensitivity_knobs.yaml`
  - missing `results/astro_sensitivity_knobs.csv`
  - missing `results/astro_sensitivity_knobs_raw.csv`
- `config_astro_lh_mode_consistency.yaml`
  - missing `results/astro_lh_mode_consistency.csv`
  - missing `results/astro_lh_mode_consistency_raw.csv`
- `config_astro_porosity_check.yaml`
  - missing `results/astro_porosity_check.csv`
  - missing `results/astro_porosity_check_raw.csv`
- `config_astro_grain_size_check.yaml`
  - missing `results/astro_grain_size_check.csv`
  - missing `results/astro_grain_size_check_raw.csv`
- `config_astro_sticking_model_check.yaml`
  - missing `results/astro_sticking_model_check.csv`
  - missing `results/astro_sticking_model_check_raw.csv`
- `config_astro_transport_sensitivity.yaml`
  - missing `results/astro_transport_sensitivity.csv`
  - missing `results/astro_transport_sensitivity_raw.csv`
- `config_astro_uv_pilot_baseline.yaml`
  - missing `results/astro_uv_pilot_baseline.csv`
  - missing `results/astro_uv_pilot_baseline_raw.csv`
- `config_astro_uv_pilot_photofrag.yaml`
  - missing `results/astro_uv_pilot_photofrag.csv`
  - missing `results/astro_uv_pilot_photofrag_raw.csv`

### Partially covered

- `config_astro_mrn_integration.yaml`
  - aggregate-style result exists as `results/astro_mrn_integration.csv`
  - configured raw output is missing locally

## Report Consistency Notes

`reports/PAPER_COMPLETION_STATUS.md` appears to reflect a more complete support
package than what is available in this checkout.

It references outputs that are not currently present locally, including:

- `results/astro_porosity_check.csv`
- `results/astro_sticking_model_check.csv`
- `results/astro_transport_sensitivity.csv`
- `results/plots/astro_lh_mode_consistency/lh_mode_consistency_summary.csv`

So that report should be read as a project-status memo, not a strict file-system
inventory.

## UV Provenance Note

The local code confirms the current UV mismatch:

- ISM configs mostly use `diffusion_mode: rate_only`
- default UV H2 uses explicit adjacent pairs
- therefore `h2_formed_UV_mean = 0` is structurally expected in many ISM runs

The exploratory photofragmentation path exists in code and config, but there is
not yet a trustworthy validated artifact package for using it as a main paper
result.

## Practical Use Order

If writing or briefing from this workspace today, use artifacts in this order:

1. `results/jhub_full_merged.csv`
2. `results/tables/grieco_paperfit_*`
3. `results/tables/table_timescales_representative.csv`
4. `results/tables/ct02_comparison_jhub_full.csv`
5. `results/tables/astro_mrn_comparison.csv`
6. `results/plots/paper/` and `results/plots/paper_suite/`

Avoid using miscellaneous older smoke files unless the goal is explicitly to
document regression history.

