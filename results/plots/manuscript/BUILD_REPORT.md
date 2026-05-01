# Canonical MNRAS Figure Build

Generated with `mnras_style.mplstyle` + `mnras_figures.py` helpers and saved under `results/plots/manuscript/`.

## Generated
- `fig01_grain_lattice.pdf/.png`
- `fig02_binding_energies.pdf/.png`
- `fig03_mechanism_schematic.tex`
- `fig04_grieco_validation.pdf/.png`
- `fig05_physisorption_only.pdf/.png`
- `fig06_tau_sensitivity.pdf/.png`
- `fig07_epsilon_all_densities.pdf/.png`
- `fig08_mechanism_decomp.pdf/.png`
- `fig09_surface_h.pdf/.png`
- `fig10_transition_zoom.pdf/.png`
- `fig11_phase_map.pdf/.png`
- `fig12_release_rate.pdf/.png`
- `fig13_ct10_comparison.pdf/.png`
- `fig14_mrn_integration.pdf/.png`
- `fig16_ensemble_convergence.pdf/.png`
- `fig17_sensitivity_envelope.pdf/.png`
- `fig18_transition_distributions.pdf/.png`
- `fig19_lh_mode.pdf/.png`
- `fig20_grain_size.pdf/.png`
- `fig21_porosity_sticking.pdf/.png`
- `fig22_timescale_map.pdf/.png`
- `fig23_kinetic_trajectory.pdf/.png`

## Notes
- Figure 1 is rendered as an illustrative porous 3D grain consistent with the manuscript geometry parameters (radius, porosity, chemisorption fraction, defect fraction), rather than the shallow cached simulation slab used internally for kinetics.
- Figure 3 is emitted as TikZ source for manual Inkscape/LaTeX finishing rather than as a Matplotlib plot.
- Figure 6 uses the available 20 K and 200 K tau-sensitivity table; the local repo does not currently contain a matching 10/20/50 K set.
- MRN figure styling is final-form, but the local warm-regime MRN provenance still deserves a science-level recheck before submission.
- Figure 15 (UV suppression) was intentionally not rebuilt in this pass because the current manuscript direction is non-UV.
