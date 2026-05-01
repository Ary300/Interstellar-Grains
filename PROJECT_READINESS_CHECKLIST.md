# Project Readiness Checklist

This checklist maps your research scope to code artifacts in this repository.

## 1) Interstellar Dust Grain Dynamics (KMC core)

- [x] Stochastic KMC / Gillespie engine: `kmc_simulation.py`
- [x] 3D heterogeneous grain lattice + porosity + defect + chemisorption site types: `kmc_simulation.py`
- [x] LH / ER / UV-assisted channels tracked separately: `kmc_simulation.py`
- [x] Temperature ramp (DED-style) support: `kmc_simulation.py`
- [x] Grain structure caching for repeated runs: `kmc_simulation.py`

## 2) Validation discipline (Grieco benchmark)

- [x] Isothermal high-T validation: `grieco_validation.py`
- [x] DED ramp validation with blocking + binning: `grieco_ded_validation.py`
- [x] Mechanism ablation checks with assertions: `grieco_mechanism_checks.py`
- [x] Hold-out calibration/generalization check: `grieco_holdout_validation.py`
- [x] Sensitivity tooling (OAT/PRCC): `grieco_sensitivity.py`

## 3) Astrophysical campaign setup

- [x] Pilot ISM campaign config: `config_astro_pilot.yaml`
- [x] Full ISM campaign config: `config_astro_full.yaml`
- [x] MRN support in sweep runner: `run_sweep.py`
- [x] Cache passthrough in sweep runner (`enable_grain_cache`, `grain_cache_dir`): `run_sweep.py`

## 4) HPC (Anvil) operationalization

- [x] Config sharding script: `anvil/generate_combo_configs.py`
- [x] Slurm array script: `anvil/run_array.sbatch`
- [x] Submit helper: `anvil/submit_array.sh`
- [x] Result merge utility: `anvil/merge_results.py`
- [x] End-to-end Anvil instructions: `anvil/README.md`

## 5) Remaining science work (not boilerplate)

- [ ] Run transition-region deep ensembles (`config_astro_transition_deep.yaml`) and plot the stochastic distributions.
- [ ] Run sensitivity envelope checks (`config_astro_sensitivity_knobs.yaml`, `config_astro_transport_sensitivity.yaml`).
- [ ] Run the validation-vs-prediction cross-check (`config_grieco_paper_*_predmodel.yaml`) and summarize with `compare_grieco_model_variants.py`.
- [ ] Run quick referee checks (LH mode consistency, porosity, grain size, sticking model).
- [ ] Run MRN integration slice (`config_astro_mrn_integration.yaml`) and compare to the single-grain campaign.
- [ ] Produce final paper figures/tables (validation summary, mechanism decomposition, timescale table, MRN correction).
