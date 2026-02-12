import numpy as np

# Physical constants
EV_TO_KELVIN = 11604.525
EV_TO_J = 1.602176634e-19
EV_TO_ERG = EV_TO_J * 1e7
K_B = 8.617e-5  # eV/K
K_B_ERG = 1.380649e-16  # erg/K
H_BAR = 6.582e-16  # eV·s
H_BAR_J_S = H_BAR * EV_TO_J
M_H = 1.673e-24  # g
M_H2 = 3.348e-24  # g

# Literature-derived energetics for H atoms on carbonaceous surfaces
# Based on: Zecho et al. 2002, Sha et al. 2002, Cuppen et al. 2013, Morisset et al. 2005

# H diffusion barriers (eV)
H_DIFFUSION_BARRIERS = {
    "graphite": 0.03,      # Zecho et al. 2002
    "amorphous_carbon": 0.025,  # Estimated from graphite
    "defect": 0.015,       # Enhanced diffusion at defects
}

# H binding energies (eV)
H_BINDING_ENERGIES = {
    "physisorption": 0.045,    # 45 meV = ~520 K equivalent
    # Chemisorbed H on carbonaceous materials (C–H bond) is much stronger than physisorption
    # and enables long-lived surface H up to high temperatures (e.g., PAH/coronene films).
    "chemisorption": 1.75,
    "defect": 0.35,            # Enhanced physisorption-like binding at defects
}

# H2 formation barriers (eV)
H2_FORMATION_BARRIERS = {
    "LH_barrierless": 0.0,     # Cuppen et al. 2013 - often barrierless
    "LH_with_barrier": 0.02,   # Small barrier if present
    "ER_barrier": 0.15,        # Morisset et al. 2005
}

# Prefactors (attempt frequencies) from TST (s^-1)
TST_PREFACTORS = {
    "diffusion": 1e12,         # Surface phonon frequency
    "desorption": 1e13,        # Molecular vibration frequency
    "reaction": 1e12,          # Surface reaction attempt frequency
}

# Quantum tunneling parameters
TUNNELING_PARAMS = {
    # Effective barrier width for surface diffusion tunneling.
    # Values ~2–3 Å keep tunneling significant at low T without making diffusion unrealistically instantaneous.
    "barrier_width_angstroms": 3.0,
    "effective_mass": M_H,            # H atom mass
}

# UV photochemistry parameters
UV_PARAMS = {
    "absorption_cross_section_cm2": 1e-17,  # Typical for H atoms on surfaces
    "photodesorption_yield": 0.1,           # 10% yield per absorbed photon
    "photodissociation_yield": 1e-6,        # Very low for H2
}

# Surface chemistry data (existing)
surface_chemistry_data = {
    "er_cross_section_cm2": 1e-15,  # Eley-Rideal cross-section
    "uv_h2_formation_yield_per_pair": 1e-4,  # UV-assisted H2 formation yield
}

# UV photon flux (existing)
uv_photon_flux = {
    "integrated_fuv_photon_flux_photons_cm2_s": 1e8,  # Standard ISRF
}
