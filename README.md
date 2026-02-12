Kinetic Monte Carlo Simulation of H₂ Formation on Carbonaceous Interstellar Grains

## Research Goals

We present a comprehensive kinetic Monte Carlo (KMC) simulation framework for studying molecular hydrogen (H₂) formation on carbonaceous interstellar dust grains under astrophysical conditions. Our research addresses a critical gap in astrochemical modeling by focusing on carbonaceous grains, which constitute a significant fraction of interstellar dust but have received less attention than silicate grains in H₂ formation studies. The primary objectives are to: (1) quantify the efficiency of H₂ formation on carbonaceous surfaces across relevant temperature (10-50 K) and density (10²-10⁴ cm⁻³) ranges, (2) elucidate the relative contributions of Langmuir-Hinshelwood (LH), Eley-Rideal (ER), and UV-assisted formation mechanisms, and (3) investigate the role of surface heterogeneity and stochastic effects in low-density interstellar environments.

## Methodology

Our approach implements a sophisticated 3D amorphous carbon lattice model with realistic surface physics. **Key methodological innovations include:**

- Novel 3D Surface Structure: Unlike traditional 2D models, we implement a 3D amorphous carbon lattice with porosity (20%), surface defects (15%), and chemisorption sites (10%), capturing the complex morphology of carbonaceous grains.

- Multi-Mechanism Chemistry: The simulation incorporates three distinct H₂ formation pathways: (i) LH mechanism requiring adjacent adsorbed H atoms, (ii) ER mechanism involving gas-phase H + surface H, and (iii) UV-assisted formation during stochastic photon pulses.

- Stochastic UV Treatment: We introduce a novel stochastic UV pulse model (1-10 photons grain⁻¹ yr⁻¹) that creates temporary surface defects and stimulates H₂ formation, addressing the role of intermittent UV irradiation in star-forming regions.

- Surface Heterogeneity: Binding energies follow realistic distributions (physisorption: 50±5 meV, chemisorption: 1.75±0.25 eV) with site-specific diffusion barriers, capturing the energy landscape complexity of carbonaceous surfaces.

The KMC algorithm uses the Gillespie n-fold way method with 10⁻⁶ second time resolution, ensuring accurate treatment of rare events in low-density environments.

## Grieco et al. (2023) ε(T) validation harness (calibrated)

To benchmark this simulator against the “efficient H₂ formation up to 250 K on carbonaceous (coronene-like) surfaces” result (Nature Astronomy, 2023), this repo includes an optional **arrival/beam mode** plus a **DED-style temperature ramp** for the <100 K part of ε(T).

**Important:** the settings here produce a *calibrated* minimal KMC that reproduces the **reported regime behavior** under a comparable observable definition. It is **not** a parameter-free recreation of the FORMOLISM setup or coronene film microphysics.

- Use `arrival_rate_per_site_s` (recommended) or `arrival_rate_s` to model a Poisson stream of particles impinging the surface (instead of computing adsorption from `h_gas_density_cm3`).
- Handle ER/abstraction at the moment of arrival using `er_cross_section_cm2` and `er_reaction_probability`.
- Use strong chemisorption energetics (`E_chem_mean_eV`, `heterogeneity_E_chem_sigma_eV`) to sustain a chemisorbed H reservoir at high T.

Run the included sweep:

```bash
python grieco_validation.py --temps 100 150 200 250 --output results/grieco_validation.csv
```

This produces an ε(T) CSV where **ε ≡ 2·N(H₂ prompt-desorbed from formation)/N(atoms impinging)** over a measurement window (i.e., “detected/desorbed recombination products”, not total formed).

You can also run the same idea through the general sweep runner by using `mode: grieco` (see `config_grieco.yaml`):

```bash
python run_sweep.py config_grieco.yaml
```

To include the two protocol-critical low-temperature features discussed in the paper (1 K/min DED ramp below 100 K + blocking by adsorbed molecules), run:

```bash
python grieco_ded_validation.py --t-start 10 --t-end 100 --rate-k-per-min 1 --output results/grieco_ded_validation.csv --summary-json results/grieco_ded_validation_summary.json
```

For a small, two-stage “publishability discipline” calibration run (fit only the high-T plateau first, then fit low-T blocking knobs), run:

```bash
python calibrate_grieco.py --out results/grieco_calibration.json
```

The repo also includes a calibrated example config matching the three regime targets (high‑T plateau, mid‑T band, low‑T collapse) at the level of this minimal KMC model:

```bash
python run_sweep.py config_grieco_calibrated.yaml
python grieco_validation.py --output results/grieco_validation_calibrated.csv
python grieco_ded_validation.py --output results/grieco_ded_validation_calibrated.csv --t-end 80 --replicates 1
```

For “defensibility” tooling:

- `python grieco_mechanism_checks.py --assert` runs “turn off a mechanism” checks and fails nonzero if expectations are violated (no chemisorption ⇒ high‑T plateau vanishes; no blocking ⇒ low‑T collapse disappears; beam-only baseline doesn’t inflate ε).
- `python grieco_holdout_validation.py` runs simple hold-out checks (fit on a subset of conditions and report withheld conditions).
- `python grieco_convergence.py --out results/grieco_convergence.csv` generates a small convergence/regression table across grain radii.
- `python grieco_flux_dependence.py --out-csv results/grieco_flux_dependence.csv` runs a small arrival-flux dependence sweep and records ε and coverage metrics.
- `python grieco_sanity_checks.py` runs quick ε-accounting sanity cases (dissociation off, strong H2 binding, etc.).
- `python grieco_sensitivity.py --mode oat` produces a one‑at‑a‑time sensitivity CSV (and `--mode prcc` can produce a small PRCC table; it can be slow).

Runtime note: the default harness runs are designed to be local-friendly; HPC is mainly helpful for dense grids and large ensembles (e.g., ≥20–50 replicates across many temperatures/UV values).

Optional performance knob: set `enable_grain_cache: true` in configs to reuse pre-generated grain topology/energy maps across runs. You can also set `grain_cache_dir` (default: `grain_cache`) and `grain_cache_include_rng_seed` if you want separate cached grains per RNG seed.

## Astrophysical Campaign Setup

For ISM-like production runs (gas-kinetic adsorption mode), use:

- `config_astro_pilot.yaml` for short pilot checks.
- `config_astro_full.yaml` for large UV/T/nH sweeps.

Run locally:

```bash
python run_sweep.py config_astro_pilot.yaml
```

For Anvil CPU arrays:

1. Generate shard configs:
   `python anvil/generate_combo_configs.py --base-config config_astro_full.yaml --shards 64`
2. Upload project directory to Anvil.
3. Submit array:
   `bash anvil/submit_array.sh anvil/generated/manifest.txt "$PWD" h2-kmc`
4. Merge results:
   `python anvil/merge_results.py --input-glob 'results/anvil/aggregated_shard_*.csv' --output results/astro_full_merged.csv`

Detailed Anvil instructions are in `anvil/README.md`.

## Results and Field Comparison

- Current Status: Initial simulations reveal the inherent challenges of H₂ formation under realistic interstellar conditions. Our results show adsorption rates of ~10⁻⁶ s⁻¹ and H residence times of ~500,000 years at 10 K, consistent with theoretical expectations but highlighting the need for longer simulation times and higher gas densities for meaningful H₂ production.

- Field Context: Our work represents the first comprehensive KMC study of H₂ formation on carbonaceous grains with 3D surface structure and stochastic UV effects. While previous studies have focused primarily on silicate grains or simplified 2D models, our approach captures the unique properties of carbonaceous surfaces that may enable H₂ formation at higher temperatures than traditionally expected.

- Technical Validation: The simulation framework successfully reproduces expected physical behaviors: (1) temperature-dependent adsorption/desorption kinetics, (2) realistic binding energy distributions, and (3) proper stochastic event handling. The early termination (1 second) with zero H₂ formation is physically consistent given the extremely low gas densities (100-10,000 cm⁻³) and short simulation times relative to astrophysical timescales.

- Novel Contributions: Our methodology addresses three key limitations in existing literature: (1) the lack of 3D surface structure in previous models, (2) the absence of stochastic UV effects in most KMC studies, and (3) the focus on silicate rather than carbonaceous grain chemistry. The framework is designed to eventually test recent experimental findings suggesting H₂ formation on carbonaceous surfaces up to 250 K, which would resolve the "high-temperature H₂ formation paradox" currently debated in the field.

- Next Steps: Parameter optimization is required to achieve meaningful H₂ formation rates within computational constraints, including increased gas densities, extended simulation times, and larger ensemble sizes. Once optimized, the model will provide the first theoretical framework for understanding the relative importance of carbonaceous vs. silicate grains in interstellar H₂ chemistry and make testable predictions for JWST observations of H₂ formation in diverse astrophysical environments.


Aug 8th update: finished todos (2D to 3D lattice, physisorption site map, n-fold way event selection, UV stuff)
