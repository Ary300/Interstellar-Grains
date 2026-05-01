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
- `fig24_surface_energy_map.pdf/.png`
- `fig25_layer_gallery.pdf/.png`
- `fig26_pore_architecture.pdf/.png`

## Notes
- Figure 1 is a hybrid setup figure: a schematic shell cutaway consistent with the manuscript geometry parameters plus a discrete slice from the actual cached grain lattice.
- Figure 3 is emitted as TikZ source for manual Inkscape/LaTeX finishing rather than as a Matplotlib plot.
- Figure 6 uses the available 20 K and 200 K tau-sensitivity table; the local repo does not currently contain a matching 10/20/50 K set.
- MRN figure styling is final-form, but the local warm-regime MRN provenance still deserves a science-level recheck before submission.
- Figure 24 is an additional support-style figure showing the top-layer site classes and binding-energy map from the actual cached grain, useful for methods or supplementary placement.
- Figure 25 shows how site classes are distributed across surface, mid-depth, and deep cached grain layers; it is a non-line support figure tied directly to the cached lattice.
- Figure 26 shows pore architecture from the actual cached grain using occupied-layer count and a meridional material/void slice.
- Figure 15 (UV suppression) was intentionally not rebuilt in this pass because the current manuscript direction is non-UV.
