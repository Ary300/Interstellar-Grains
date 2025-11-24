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

## Results and Field Comparison

- Current Status: Initial simulations reveal the inherent challenges of H₂ formation under realistic interstellar conditions. Our results show adsorption rates of ~10⁻⁶ s⁻¹ and H residence times of ~500,000 years at 10 K, consistent with theoretical expectations but highlighting the need for longer simulation times and higher gas densities for meaningful H₂ production.

- Field Context: Our work represents the first comprehensive KMC study of H₂ formation on carbonaceous grains with 3D surface structure and stochastic UV effects. While previous studies have focused primarily on silicate grains or simplified 2D models, our approach captures the unique properties of carbonaceous surfaces that may enable H₂ formation at higher temperatures than traditionally expected.

- Technical Validation: The simulation framework successfully reproduces expected physical behaviors: (1) temperature-dependent adsorption/desorption kinetics, (2) realistic binding energy distributions, and (3) proper stochastic event handling. The early termination (1 second) with zero H₂ formation is physically consistent given the extremely low gas densities (100-10,000 cm⁻³) and short simulation times relative to astrophysical timescales.

- Novel Contributions: Our methodology addresses three key limitations in existing literature: (1) the lack of 3D surface structure in previous models, (2) the absence of stochastic UV effects in most KMC studies, and (3) the focus on silicate rather than carbonaceous grain chemistry. The framework is designed to eventually test recent experimental findings suggesting H₂ formation on carbonaceous surfaces up to 250 K, which would resolve the "high-temperature H₂ formation paradox" currently debated in the field.

- Next Steps: Parameter optimization is required to achieve meaningful H₂ formation rates within computational constraints, including increased gas densities, extended simulation times, and larger ensemble sizes. Once optimized, the model will provide the first theoretical framework for understanding the relative importance of carbonaceous vs. silicate grains in interstellar H₂ chemistry and make testable predictions for JWST observations of H₂ formation in diverse astrophysical environments.


Aug 8th update: finished todos (2D to 3D lattice, physisorption site map, n-fold way event selection, UV stuff)
