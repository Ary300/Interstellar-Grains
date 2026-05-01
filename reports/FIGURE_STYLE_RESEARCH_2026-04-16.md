# Figure Style Research Notes

Updated: 2026-04-16

This note captures the figure-design constraints and nearby-paper patterns used
for the current manuscript figure refresh.

## Journal-facing constraints

### AAS / ApJ

From the AAS Graphics Guide:

- do not trust plotting defaults
- use vector formats where possible
- minimum font size is 6 pt
- line widths should be at least 0.5 pt
- color should not be the only separator; line style and marker shape should also differ

Source:

- https://journals.aas.org/graphics-guide/

### MNRAS

From the MNRAS instructions:

- prepare figures to survive reduction to a single-column width of about 84 mm
- final line weight should not be below 0.3 pt
- labels should still read at about 8 pt final size
- overly complex dash patterns should be avoided

Source:

- https://academic.oup.com/mnras/pages/general_instructions

## Nearby-paper patterns

### Cuppen and Herbst / stochastic-heating Monte Carlo tradition

The classic MNRAS Monte Carlo H2-formation papers organize figures around:

- temperature histories or temperature-response curves
- formation efficiency versus temperature
- grain-size dependence
- a limited number of comparisons per panel

This is a strong match to the present paper's baseline overview plus MRN and
mechanism-crossover panels.

Source:

- https://academic.oup.com/mnras/article/367/4/1757/1747618

### Recent off-lattice MNRAS example

Recent microscopic grain papers still combine:

- morphology or grain snapshots
- occupancy versus temperature
- energy-distribution figures
- mechanism-importance plots

That supports keeping the grain-characterization figure in the main or near-main
set rather than hiding it as a disposable diagnostic.

Source:

- https://academic.oup.com/mnras/article/543/3/2567/8262827

### Le Bourlot / broader ISM surface-chemistry framing

The larger ISM surface-chemistry papers emphasize compact physical stories:

- mechanism competition
- rate normalization against canonical ISM values
- parameter trends rather than decorative plot complexity

That argues for warm-regime annotations and simple regime bands instead of adding
more crowded overlays.

Source:

- https://arxiv.org/abs/1202.0374

## Applied design decisions

The refreshed manuscript plots therefore use:

- a white background instead of tinted axes fills
- reduction-safe line weights and type
- one color + one marker + one line style per density
- shaded physical regimes only where they add interpretation
- titles that state the conclusion rather than restate the axes
- a no-UV main manuscript build; UV remains supplementary if needed
