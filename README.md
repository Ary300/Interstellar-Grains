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


08/06 results:

✅ Temperature-dependent sticking: S(T) ∝ exp(-T/100, K)
  - 10 K: factor 0.905, rate 2.40e-04 s^-1
  - 100 K: factor 0.368, rate 9.78e-05 s^-1
  - ✓ Sticking decreases with temperature
  - ✓ Temperature dependence correct

✅ Surface diffusion: k_diff = ν exp(-E_diff / k_B T)
  - ν = 10^12 s^-1
  - E_diff = 0.3 × E_bind for all site types
  - ✓ Diffusion barriers correctly implemented

✅ H2 formation mechanisms:
  - Langmuir-Hinshelwood: Two adsorbed H atoms meet via diffusion
  - Eley-Rideal: Gas-phase H atom reacts with adsorbed H
  - ✓ Both mechanisms fully implemented
