# Manuscript Figure Review — 2026-04-16

This note records the reference pass and production review used to push the manuscript figures toward journal-ready quality.

## External References

- AAS Graphics Guide: https://journals.aas.org/graphics-guide/
- MNRAS Instructions to Authors: https://academic.oup.com/mnras/pages/general_instructions
- IOP figure guidance: https://publishingsupport.iopscience.iop.org/questions/figures-journal-articles/
- NeurIPS author guidance on readability/accessibility: https://neurips.cc/Conferences/2021/PaperInformation/Author-Guidelines
- Cuppen & Herbst (2005), MNRAS: https://academic.oup.com/mnras/article/361/2/565/1058966
- Cuppen, Morata & Herbst (2006), MNRAS: https://academic.oup.com/mnras/article/367/4/1757/1747618

## Production Criteria Applied

- Single-column figures built at `3.33 in` width; double-column figures built at `7.0 in` width.
- Serif body-matched figure typography with embedded TrueType fonts in exported PDFs.
- Inward ticks on all four sides with visible minors on scientific axes.
- Distinguishable marker+line combinations so grayscale reading does not depend on color alone.
- Minimal in-panel titling; captions are expected to carry interpretation.
- Confidence regions kept light enough not to overpower the mean curve.
- Vector PDF export used for line-art figures; grayscale review copies generated for all manuscript figures.

## Review Outcome

### Ready or close to ready on visual grounds

- `fig04_grieco_validation`
- `fig07_efficiency_density`
- `fig08_mechanism_decomposition`
- `fig09_surface_inventory`
- `fig10_transition_zoom`
- `fig11_phase_map`
- `fig12_release_rate`
- `fig13_ct10_comparison`
- `fig16_ensemble_convergence`
- `fig17_sensitivity_envelope`
- `fig18_transition_boxplots`
- `fig19_lh_mode_consistency`
- `fig20_grain_size_dependence`
- `fig21_porosity_sticking`
- `fig22_timescale_phase_map`

### Improved substantially, but still more stylistic than journal-native

- `fig01_grain_structure`
  - Better framing and depth now, but still reads more like a careful methods render than a classic polished plate.
- `fig03_mechanism_schematic`
  - Much clearer and more publication-like than the earlier draft, but still a custom schematic rather than something drawn in a dedicated vector editor.

### Science caveat rather than style caveat

- `fig14_mrn_integration`
- `mrn_correction_factors`

These are visually much stronger than before, but the warm-regime MRN interpretation still needs provenance/science confirmation before they should be treated as final manuscript figures.

## Concrete Fixes Completed In This Pass

- Forced PDF exports away from Type 3 fonts and into embedded TrueType fonts.
- Reworked the grain render framing and perspective.
- Refined the mechanism schematic into a cleaner three-panel manuscript figure.
- Smoothed the phase-map presentation for `fig11` and `fig22` while keeping the underlying sampled structure intact.
- Added or adjusted marker labels in the astrophysical timescale map for readability at column width.
- Re-ran grayscale QC on the manuscript figure directory.

## Remaining Recommendations

- Use the numbered `figXX_*` outputs as the canonical manuscript figures.
- Treat legacy outputs such as `fig01_grieco_validation` and `fig02_ism_overview` as non-canonical helper products.
- If the paper is moving to submission soon, the next worthwhile step is a caption-and-LaTeX pass, not more styling churn.
