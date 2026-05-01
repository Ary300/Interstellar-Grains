# KMC Literature Positioning Notes

This note is a working paper-positioning summary for the H2-on-carbonaceous-grains KMC project.

## Core positioning

- The strongest novelty claim is not "we did KMC on carbonaceous grains".
- The strongest novelty claim is "we present a KMC model calibrated against the Grieco et al. (2023) coronene experiment across 10-250 K, then extend that calibrated microphysics to ISM conditions."
- That is a real gap relative to the published KMC/MC literature surveyed below.

## Key KMC / MC papers

### Chang, Cuppen & Herbst (2005), A&A 434, 599-611

- Method: continuous-time random-walk Monte Carlo.
- Surface class: olivine, amorphous carbon, mixed/inhomogeneous surfaces.
- Main use in our framing: foundational lattice/CTRW H2-grain Monte Carlo reference.
- Main result relevant here: inhomogeneous surfaces broaden the temperature range for efficient H2 formation relative to homogeneous surfaces.

### Cuppen & Herbst (2005), MNRAS 361, 565-576

- Method: CTRW Monte Carlo on rough surfaces.
- Surface class: olivine and amorphous carbon with varying roughness.
- Main use in our framing: classic demonstration that roughness broadens the efficient H2-formation window.
- Relevance to this project: supports the idea that a heterogeneous carbonaceous surface can sustain broader-temperature formation than a flat one.

### Cuppen, Morata & Herbst (2006), MNRAS 367, 1757-1765

- Method: Monte Carlo with stochastic grain heating.
- Main use in our framing: the main citation for stochastic heating as a known but usually omitted extension.
- Relevance to this project: supports discussing stochastic heating as future work rather than a missing standard ingredient.

### Cuppen & Hornekær (2008), JCP 128, 174707

- Method: kinetic Monte Carlo constrained by DFT-derived barriers.
- Surface class: graphite.
- Mechanism focus: Eley-Rideal abstraction on chemisorbed hydrogen.
- Main use in our framing: direct precedent for KMC treatment of ER chemistry on graphitic/carbonaceous surfaces.

### Cazaux et al. (2011), A&A 535, A27

- Method: kinetic Monte Carlo plus sticking/barrier physics.
- Surface class: graphite.
- Main use in our framing: chemisorption-aware high-temperature formation precedent.
- Relevant result: realistic chemisorption barriers materially affect high-temperature H2 formation efficiency.

### Iqbal, Acharyya & Herbst (2012), ApJ 751, 58

- Method: Monte Carlo with physisorption and chemisorption.
- Surface class: olivine and amorphous carbon.
- Temperature range: extends to high temperatures.
- Main use in our framing: the closest older Monte Carlo precedent for chemisorption-inclusive H2 formation.
- Limitation relative to this project: predates Grieco et al. (2023), so it cannot be experimentally calibrated to that dataset.

### Satonkin et al. (2025), MNRAS 543, 2567-2574

- Method: off-lattice microscopic Monte Carlo.
- Surface class: rough carbonaceous grain.
- Temperature range: 5-35 K.
- Mechanisms: LH + ER with thermal diffusion and tunnelling.
- Main result relevant here: rough carbonaceous surfaces broaden the low-temperature efficient window; thermal hopping dominates mobility in their off-lattice model; reported E_hop / E_bind is around 0.5-0.6.
- Main limitation relative to this project: physisorption-focused low-temperature model, not a 10-250 K chemisorption-calibrated reproduction of Grieco.

## Grieco et al. (2023) experiment

- Main paper role: experimental anchor for the whole project.
- Critical observations for model comparison:
  - high-temperature efficiency plateau around 20% above about 80-100 K
  - intermediate-temperature maximum near 20 K
  - low-temperature suppression attributed to blocking by adsorbed D2
  - steady-state exposure after a few 10^15 atoms cm^-2
  - direct astrophysical framing around warm carbonaceous grains and high-z star formation

## What this project can plausibly claim if the referee-check jobs pass

- First KMC model in this repo ecosystem calibrated against the Grieco et al. (2023) coronene efficiency curve.
- Mechanism decomposition across the LH-to-ER transition that the experiment itself cannot directly measure.
- A controlled comparison between the validated lab-inspired setup and the ISM prediction model assumptions.
- ISM-rate and timescale estimates tied back to a laboratory-calibrated carbonaceous-grain model.

## Referee-sensitive points

These should be handled explicitly in the paper and checked against the new support jobs.

- Validation/prediction consistency:
  - This is the most important methodological issue.
  - The prediction-model Grieco reruns and LH consistency check are the direct response.

- Empirical transport knobs:
  - `lh_diffusion_factor`
  - `diffusion_rate_cap_s`
  - These need sensitivity bounds, not silence.

- Tunnelling:
  - Omission is most relevant at the lowest temperatures.
  - It should be discussed as a limitation, not treated as a blocker for the 100-250 K story.

- Stochastic grain heating:
  - Known in the literature but not a universal ingredient of KMC H2 papers.
  - Best handled as future work / scope control.

- Single-grain versus MRN-distributed grains:
  - The MRN support run is important for explaining why single-grain rates can sit below the canonical observational benchmark.

## Paper-comparison takeaways

- Against Chang/Cuppen/Cuppen-Morata:
  - our project is newer in target system and experimental anchoring
  - their main value is methodological precedent and low-temperature roughness/stochasticity context

- Against Cuppen & Hornekær / Cazaux et al.:
  - our project is less atomistically detailed in graphitic microphysics
  - our strength is the direct link to the Grieco coronene experiment and the broader 10-250 K validation target

- Against Iqbal et al.:
  - our project is later, carbonaceous-focused in the Grieco sense, and experimentally anchored to a modern benchmark

- Against Satonkin et al.:
  - their model is off-lattice and low-temperature focused
  - our project is broader in temperature and explicitly centered on the chemisorption-enabled warm regime highlighted by Grieco

## Plotting plan once the new jobs finish cleanly

The highest-value figure set for the paper should be:

1. Grieco Figure-2-style overlay with fit statistics.
2. ISM epsilon versus temperature at several densities.
3. H2 release-rate versus temperature on a log scale.
4. Mechanism decomposition (LH / ER fractions versus temperature).
5. LH-mode consistency comparison.
6. Transition-region distribution boxplots.
7. Ensemble convergence plot (running mean plus SEM).
8. Sensitivity-envelope plot.
9. Grain-characterization figure (binding energies and E_hop / E_bind).
10. MRN correction-factor figure.

The plotting style should emphasize:

- publication-clean typography
- consistent temperature-color mapping
- light gridlines, no clutter
- larger labels than the current default scripts
- panel titles that state the physical point, not just the variable name
- SI/cgs-consistent axis labels

## Source list used for this note

- Grieco et al. (2023), Nature Astronomy:
  https://www.nature.com/articles/s41550-023-01902-4
- Chang, Cuppen & Herbst (2005), A&A 434:
  https://www.aanda.org/articles/aa/full/2005/17/aa1842/aa1842.right.html
- Cuppen & Herbst (2005), MNRAS 361:
  https://academic.oup.com/mnras/article/361/2/565/1058966
- Cuppen, Morata & Herbst (2006), MNRAS 367:
  https://academic.oup.com/mnras/article-abstract/367/4/1757/1747618
- Cuppen & Hornekær (2008), JCP 128:
  https://pubmed.ncbi.nlm.nih.gov/18465936/
- Cazaux et al. (2011), A&A 535:
  https://www.aanda.org/articles/aa/full_html/2011/11/aa17220-11/aa17220-11.html
- Satonkin et al. (2025), MNRAS 543 / arXiv:
  https://academic.oup.com/mnras/article-abstract/543/3/2567/8262827
  https://arxiv.org/abs/2509.04913
