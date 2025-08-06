to do:
-Upgrade from 2D square grid model to 3D lattice or multi-layer surface representation for better amorphous structure, porosity, and defect modeling

-Introduce a separate site type map for chemisorption vs. physisorption
Adjust diffusion, reaction rates, and binding energies accordingly (1.5–2.0 eV for physi)

-Implement n-fold way event selection and time advancement (Δt = -ln(r)/k_total) algorithm for KMC

-No UV-driven processes (photodesorption, defect creation) in simulation logic, though parameters exist in scientific_data.py; add stochastic UV pulses with:
Photodesorption rate k_UV = G0 × 10⁻³ s⁻¹
UV-induced site modification (defect creation)
Couple uv_flux_factor from config.yaml to actual UV event generation

-Parallelization & performance??
