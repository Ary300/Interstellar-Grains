# Grieco (2023) validation + ISM sweep report (paperfit baseline)

This report documents:

1) **How the KMC model works** (Gillespie time, events, observables).
2) **How we mirror the Grieco et al. (2023) protocol** (beam angle, dissociation fraction, dose).
3) **What was calibrated vs fixed** (parameterization, not unconstrained curve-fitting).
4) **What the “paperfit” validation results are** and where the artifacts live.
5) **How this carries into ISM / UV sweeps on Anvil** and how to interpret those outputs.

Key reference: *Grieco et al.*, **Nature Astronomy** (2023), doi:10.1038/s41550-023-01902-4.

---

## 0) “Parameterization” vs “curve-fitting” (how to frame this honestly)

What we are doing is **calibration of an effective, physically motivated model**:

- The code is not a generic polynomial fit; it is a mechanistic stochastic model with explicit processes:
  adsorption (accretion), diffusion, desorption, LH formation, ER formation, chemisorption reservoir, and (optionally) UV processes + H2 blocking.
- Some parameters are **protocol constraints** (beam incidence angle, dissociation fraction, dose).
- Some are **literature-constrained physics** (binding energies, ER abstraction cross-sections within measured/DFT ranges).
- A small subset are **effective parameters** that absorb unresolved microscopic detail (e.g., an ER “reaction probability” given an encounter; an effective chemisorption-site fraction; low‑T H2 sticking/blocking knobs).

The way to keep this defensible:

1) **Count degrees of freedom**: explicitly list which knobs were tuned.
2) **Show constraints**: for each tuned knob, provide a physically plausible range and cite a source (or state it is an effective/unknown quantity).
3) **Separate calibration from prediction**: matching Grieco is **calibration**; the novel science is using the *same calibrated parameters* to predict ISM/UV dependence outside lab conditions.

We generate a paper-ready parameter table for this purpose:

- `results/tables/grieco_paperfit_parameter_table.csv`

---

## 1) What “DED / TPDED” means in the Grieco paper

**DED** = *During Exposure Desorption*: you expose the surface to the beam while monitoring desorption products with the QMS.

**TPDED** = *Temperature‑Programmed DED*: the surface temperature is ramped (paper uses ~1 K/min) while exposure continues, and the QMS monitors the evolving desorption signal during the ramp.

In Fig. 2 of Grieco et al. (2023):

- **Low‑T points (≈10–80 K)** come from **TPDED** ramps.
- **High‑T points (≥100 K)** are **isothermal** measurements.

The key experimental observable is the **recombination efficiency** ε(T): the probability that an *impinging atom* recombines (directly or indirectly) into a molecule that is detected leaving the surface.

---

## 2) Observable definition: ε in the paper vs in the code

### 2.1 Paper concept (QMS + background subtraction + incomplete dissociation)

Experimentally, the QMS sees:

- a background signal (chamber + walls),
- a contribution from the *undissociated* molecular beam fraction,
- plus additional D2 formed on the surface when the plasma is ON.

The paper’s ε definition corrects for:

- background subtraction (flag ON/OFF),
- non‑perfect dissociation (τ).

### 2.2 Simulation bookkeeping (what we count)

In **arrival/beam mode**, each “arrival” event is classified as:

- **atom** with probability `beam_dissociation_fraction = τ`, or
- **molecule** with probability `1 − τ` (represents the undissociated beam fraction).

We track “formed-origin” vs “beam-origin” H2 separately:

- `h2_molecules_desorbed` = **prompt** desorption of *formed-origin* H2 (what the QMS essentially measures promptly).
- `h2_molecules_released_formed` = **delayed** gas release of *formed-origin* H2 that temporarily stuck (important at low T).
- `h2_molecules_desorbed_beam` / `h2_molecules_released_beam` = beam-origin H2 baseline channels.
- `total_impinging_h_atoms` = count of atomic arrivals in the measurement window.

From this we define (dimensionless) ε:

- **Prompt** efficiency:
  - `epsilon_prompt = 2 * h2_molecules_desorbed / total_impinging_h_atoms`
- **Released-total** efficiency (recommended for TPDED / QMS comparisons when low‑T sticking matters):
  - `epsilon_released_total = 2 * (h2_molecules_desorbed + h2_molecules_released_formed) / total_impinging_h_atoms`

The “baseline subtraction” idea is mirrored because **beam-origin H2 is explicitly excluded** from the “formed-origin” numerator.

---

## 3) KMC method: what Gillespie “time steps” mean

The simulation uses a standard **Gillespie / stochastic simulation algorithm (SSA)**:

1) At state S, compute all event rates `{k_i(S)}`.
2) Total rate `K = Σ_i k_i`.
3) Sample physical time increment: `Δt ~ Exp(K)`, i.e. `Δt = -ln(u)/K` for uniform `u∈(0,1)`.
4) Choose which event happens with probability `k_i/K`.
5) Update the surface state and repeat.

Important implications:

- There is **no fixed “time step”**. The time increments are random and represent **physical seconds**.
- If diffusion is extremely fast, it can dominate `K`, causing very small `Δt` and huge step counts.
  - This is why ISM configs often use `diffusion_mode: rate_only` + `lh_formation_mode: diffusion_limited` (see §7).

Stopping conditions:

- **Lab / paperfit validation** typically stops by **arrivals or dose**, not by wall time:
  - `max_arrivals` or `target_exposure_atoms_cm2` (dose).
- **ISM campaigns** can stop by:
  - arrivals (burn-in + measurement arrivals), and/or
  - `max_time_s` and `max_steps` as safety caps.

---

## 4) Mapping Grieco experimental conditions → simulation parameters

### 4.1 Beam incidence angle (40°)

Paper: atoms approach at **40° from the surface normal**.

In the code, beam flux to total arrival rate uses:

- `effective_flux = beam_flux_total_cm2_s * cos(beam_incidence_angle_deg)`
- `arrival_rate_total = effective_flux * surface_area_cm2`

This is implemented in `kmc_simulation.py` (arrival-mode conversion of flux → total arrivals).

### 4.2 Dose “few × 10^15 atoms/cm²” (steady state)

Paper: results are stable after exposure of **a few × 10^15 atoms/cm²** (steady “slightly superhydrogenated” state).

Paperfit configs implement this as:

- `burnin_exposure_atoms_cm2: 3.0e15` (reach steady state)
- `measure_exposure_atoms_cm2: 1.5e15` (measurement window)

So we mirror the protocol: **burn-in to the same dose scale**, then measure ε.

### 4.3 “One atom lands per site every ~100 s”

Paper: “D atoms land on a surface adsorption site roughly every 100 s.”

With our paperfit beam flux:

- `beam_flux_total_cm2_s = 7.46e12`
- `cos(40°) ≈ 0.766`
- `site_area = 25 Å² = 2.5e-15 cm²`

Then the per-site arrival rate is approximately:

- `λ_site ≈ 7.46e12 * 0.766 * 2.5e-15 ≈ 1.4e-2 s^-1` → one arrival per ~70 s.

That is in the same ballpark as the paper’s “~100 s” statement.

---

## 5) The “paperfit” validation artifacts (what files to cite in Overleaf)

### 5.1 Configs (protocol + microphysics)

- Isothermal (high‑T plateau):
  - `config_grieco_paper_iso_paperfit.yaml`
- TPDED ramp (low‑T structure + blocking):
  - `config_grieco_paper_ded_paperfit.yaml`

### 5.2 Outputs (numbers)

- Isothermal ε(T) CSV:
  - `results/grieco_validation_paper_iso_paperfit.csv`
- TPDED binned ε(T) CSV:
  - `results/grieco_ded_paper_ded_paperfit.csv`
- TPDED summary JSON (key metrics + CI):
  - `results/grieco_ded_paper_ded_paperfit_summary.json`
- Key-number table (paper-ready):
  - `results/tables/grieco_paperfit_key_numbers.csv`
- Validation summary table (paper-ready):
  - `results/tables/grieco_paperfit_validation_summary_table.csv`
- Parameter table with constraints (paper-ready):
  - `results/tables/grieco_paperfit_parameter_table.csv`

### 5.3 Figure

- Figure‑2‑style plot (simulation-only, with reference lines):
  - `results/plots/grieco_figure2_like_paperfit.png`

Regenerating it:

```bash
python3 plot_grieco_figure2_like.py \
  --iso-csv results/grieco_validation_paper_iso_paperfit.csv \
  --ded-csv results/grieco_ded_paper_ded_paperfit.csv \
  --out results/plots/grieco_figure2_like_paperfit.png
```

---

## 6) Paperfit validation results (what they mean)

From `results/tables/grieco_paperfit_validation_summary_table.csv`:

- **High‑T plateau (100–250 K)**:
  - mean ε ≈ **0.195** (flat within small CI)
  - Interpretation: dominated by **chemisorption reservoir + ER-on-arrival** sustaining formation at high T.

- **TPDED released_total ε** (QMS-like observable):
  - ε(10 K) ≈ **0.212 ± 0.0078** (95% CI)
  - ε(20 K) ≈ **0.460 ± 0.0084**
  - ε(30–80 K) ≈ **0.298 ± 0.0013**
  - ratio ε(10)/ε(20) ≈ **0.461 ± 0.0186**

Qualitative regime match to Grieco Fig. 2:

- **Sharp drop at 10 K relative to the peak near 20 K** (captured via H2 sticking + blocking).
- **~30% band** across 30–80 K.
- **~20% plateau** above 100 K.

---

## 7) What was tuned vs fixed (degrees of freedom)

A defensible way to present this is to group parameters into:

1) **Paper protocol** (fixed): incidence angle, dose scale, ramp rate.
2) **Literature physics** (constrained ranges): binding energies, ER cross section range.
3) **Effective calibrated parameters** (tuned within bounds):
   - `chemisorption_fraction` (effective reactive-site density)
   - `er_reaction_probability` (effective microphysics; bounded [0,1])
   - low‑T blocking knobs (`E_h2_bind_eV`, `h2_stick_prob_lowT`, etc.)
4) **Numerical/performance** (not physics): diffusion mode choice, diffusion rate cap.

The canonical paper-ready table is:

- `results/tables/grieco_paperfit_parameter_table.csv`

This explicitly documents constraints and references for each key parameter.

---

## 8) Transition to ISM / UV: what changes, what stays the same

### 8.1 What stays the same

- The **microphysics parameters** you validated against Grieco (chemisorption fraction, ER cross section/probability, binding energies, blocking model) should be reused unchanged.

That’s the point of the validation: the ISM predictions are not re-tuned.

### 8.2 What changes (no directed beam)

In ISM mode, there is no directed molecular beam. Instead, H atoms arrive **isotropically** from gas kinetics:

`F = 0.25 * n_H * v_th`

and per-site arrival rate is:

`arrival_rate_per_site_s = F * site_area_cm2`

This is supported in `run_sweep.py` via:

- `arrival_rate_mode: gas_kinetic`

### 8.3 Performance note: diffusion dominates in ISM

In many ISM regimes:

- diffusion can be *much faster* than gas arrivals,
- explicitly simulating every diffusion hop can cause “timestep collapse” (huge step counts, tiny Δt).

So ISM configs typically use:

- `diffusion_mode: rate_only`
- `lh_formation_mode: diffusion_limited`

This preserves the *diffusion-limited scaling* of LH formation without simulating each hop.

See:

- `config_astro_pilot_paperfit.yaml`
- `config_astro_full_paperfit.yaml`

---

## 9) Current Anvil campaign outputs (what they are and how to interpret them)

You have three completed campaign summaries in:

- `results/tables/table_campaign_overview.csv`
- `results/tables/table_uv_summary_near_full.csv`
- `results/tables/table_uv_suppression_near_full.csv`
- `results/tables/table_nearfull_vs_present150.csv`
- `results/tables/table_top25_near_full.csv`

Important interpretation note:

- Many Anvil campaign tables show `epsilon_mean = 0.0` because those runs were done in **gas-adsorption mode**, and ε is only defined/recorded cleanly in **arrival/beam mode** (where `total_impinging_h_atoms` is tracked).
- For ISM science, you should focus on **physical formation rates**:
  - `h2_release_rate_cm2_s_*` (molecules cm^-2 s^-1)
  - and mechanism breakdowns (`h2_formed_LH_*`, `h2_formed_ER_*`, `h2_formed_UV_*`)

The Anvil plots you generated live under:

- `results/plots/anvil_1h_merged/`
- `results/plots/anvil_near_full_merged/`
- `results/plots/anvil_present_150_merged/`

(These are the folders you rsynced back.)

---

## 10) Exact “paper-matching” next steps (recommended)

To make the Grieco match *bulletproof* for a manuscript:

1) **Digitize Grieco Fig. 2 points** (with error bars) and compute residuals / χ².
2) Run a **τ sensitivity** at (e.g.) 20 K and 200 K to show ε(T) shape is robust.
3) Run a **flux dependence** check (2–3 fluxes) to show:
   - high‑T plateau is relatively flux-insensitive,
   - low‑T blocking regime is coverage/flux sensitive (physically sensible).

Then pivot to novelty:

- Run `config_astro_full_paperfit.yaml` (or a reduced grid) on Anvil
- Report `h2_release_rate_cm2_s` vs (T, n_H, UV) and mechanism fractions.

### 10.1 Fig. 2 digitization workflow + χ² script

Digitization instructions + template live in:

- `digitization/README.md`
- `digitization/grieco_fig2_digitized_template.csv`

Once you create `grieco_fig2_digitized.csv` (or `results/grieco_fig2_digitized.csv`), run:

```bash
python3 compare_grieco_fig2_digitized.py
```

Outputs:

- `results/tables/grieco_fig2_digitized_fit_stats.csv`
- `results/tables/grieco_fig2_digitized_residuals.csv`
- `results/plots/grieco_fig2_overlay_digitized.png`
- `results/plots/grieco_fig2_residuals_digitized.png`
