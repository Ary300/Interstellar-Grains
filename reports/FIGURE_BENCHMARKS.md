# Figure Benchmarks for the H2 KMC Paper

This note summarizes the most useful figure patterns from nearby H2-on-grain KMC / Monte Carlo papers and maps them onto the figure set for this project.

## What nearby papers tend to show

### 1. Efficiency or rate versus temperature is always central

Across the classic H2-grain Monte Carlo papers, the most important figure is almost always a temperature trend:

- recombination efficiency versus grain/surface temperature
- H2 formation or release rate versus temperature
- sometimes comparisons across different surface classes or roughness assumptions

This is the correct anchor for our paper too. The Grieco validation overlay and the ISM `epsilon(T)` / rate-versus-temperature panels should remain the visual core.

## 2. Mechanism separation is especially valuable when multiple channels exist

Older papers often compare:

- smooth versus rough surfaces
- homogeneous versus heterogeneous binding energies
- low- versus high-barrier assumptions
- with- versus without-ER contributions

Our closest equivalent is mechanism decomposition:

- LH fraction versus temperature
- ER fraction versus temperature
- showing the crossover clearly

This is one of the most publishable strengths of the current project because the experiment itself cannot directly provide it.

## 3. Surface-energy distributions are common in the more microscopic papers

Papers with more explicit grain microphysics often include:

- binding-energy distributions
- diffusion-barrier distributions
- site-class fractions
- ratios like `E_hop / E_bind`

That makes the grain-characterization figure important, especially for defending the empirical transport choices.

## 4. Robustness/sensitivity figures are often simpler than people expect

The best sensitivity figures in this area are usually not complicated:

- one envelope band
- one comparison curve set
- one concise ablation or ratio plot

For this project, the right robustness figures are:

- sensitivity envelope for chemisorption fraction and ER probability
- transport sensitivity for `lh_diffusion_factor` and `diffusion_rate_cap_s`
- LH-mode consistency
- small quick-check plots for porosity, grain size, and sticking

## 5. Distribution plots matter only when the stochastic claim is central

Most older papers report means and temperature trends, but if we want to claim transition-region stochasticity, then the paper should show actual distributions.

That makes the transition-region boxplots and convergence plot more than "supplementary polish" — they are part of the argument.

## Figure style lessons from the literature

The strongest papers in this subfield usually have figures that are:

- simple in structure
- physically labeled
- sparse in color
- built around a small number of comparisons per panel

The weakest figures in this literature tend to be:

- overloaded with too many curves
- inconsistent in color meaning across panels
- too small to read when embedded in a paper
- descriptive rather than interpretive in the title

## Style choices we should keep

The updated plotting scripts now aim for:

- one color per density across panels
- one color per mechanism where mechanism panels appear
- larger type and cleaner labels
- `.png` and `.pdf` output for paper use
- titles that say what the panel means physically

## Recommended final figure set

### Core figures

1. Grieco experiment versus KMC validation overlay
2. ISM `epsilon(T)` at several densities
3. ISM H2 release rate versus temperature
4. Mechanism decomposition versus temperature
5. CT10 comparison

### Robustness / referee figures

6. Sensitivity envelope
7. Transition-region boxplots
8. Ensemble convergence
9. LH-mode consistency
10. Grain characterization
11. MRN correction
12. Quick robustness checks: porosity, grain size, sticking

## What to avoid in our paper figures

- Do not put UV in the main figure set.
- Do not mix too many storylines into one panel.
- Do not use inconsistent density colors from one plot to the next.
- Do not rely on only one giant heatmap; line-based interpretation is clearer for this paper.

## Useful references for figure benchmarking

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
- Iqbal, Acharyya & Herbst (2012), ApJ 751:
  https://doi.org/10.1088/0004-637X/751/1/58
- Satonkin et al. (2025), MNRAS 543:
  https://academic.oup.com/mnras/article-abstract/543/3/2567/8262827
- Grieco et al. (2023), Nature Astronomy:
  https://www.nature.com/articles/s41550-023-01902-4
