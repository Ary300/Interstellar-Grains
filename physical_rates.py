"""
Physically correct rate calculations using Transition State Theory (TST)
with quantum tunneling corrections for H atoms on interstellar grain surfaces.

Based on literature energetics from:
- Zecho et al. 2002 (diffusion barriers)
- Sha et al. 2002 (binding energies) 
- Cuppen et al. 2013 (H2 formation)
- Morisset et al. 2005 (Eley-Rideal)
"""

import numpy as np
from scientific_data import (
    K_B, H_BAR, M_H, EV_TO_KELVIN,
    H_DIFFUSION_BARRIERS, H_BINDING_ENERGIES, H2_FORMATION_BARRIERS,
    TST_PREFACTORS, TUNNELING_PARAMS, UV_PARAMS
)


def thermal_rate(prefactor, barrier_eV, temperature_k):
    """
    Calculate thermal rate using Transition State Theory.
    
    Args:
        prefactor: Attempt frequency (s^-1)
        barrier_eV: Energy barrier (eV)
        temperature_k: Temperature (K)
    
    Returns:
        Rate (s^-1)
    """
    if temperature_k <= 0:
        return 0.0
    
    kBT = K_B * temperature_k  # eV
    if barrier_eV <= 0:
        return prefactor
    
    return prefactor * np.exp(-barrier_eV / kBT)


def quantum_tunneling_rate(barrier_eV, barrier_width_angstroms, effective_mass, temperature_k):
    """
    Calculate quantum tunneling rate using WKB approximation.
    Based on Karssemeijer & Cuppen 2014, ApJ 781, 16
    
    Args:
        barrier_eV: Energy barrier (eV)
        barrier_width_angstroms: Barrier width (Å)
        effective_mass: Particle mass (g)
        temperature_k: Temperature (K)
    
    Returns:
        Tunneling rate (s^-1)
    """
    if barrier_eV <= 0:
        return 1e12
    
    barrier_width_m = barrier_width_angstroms * 1e-10
    barrier_joules = barrier_eV * 1.602e-19
    mass_kg = effective_mass * 1e-3
    
    hbar_si = 1.055e-34
    
    kappa = np.sqrt(2 * mass_kg * barrier_joules) / hbar_si
    
    transmission = np.exp(-2 * kappa * barrier_width_m)
    
    attempt_frequency = np.sqrt(barrier_joules / (2 * np.pi * mass_kg)) / barrier_width_m
    
    return attempt_frequency * transmission


def combined_rate(thermal_rate, tunneling_rate):
    """
    Combine thermal and tunneling rates using quantum statistics.
    
    Args:
        thermal_rate: Thermal rate (s^-1)
        tunneling_rate: Tunneling rate (s^-1)
    
    Returns:
        Combined rate (s^-1)
    """
    return thermal_rate + tunneling_rate


def h_diffusion_rate(site_type, temperature_k):
    """
    Calculate H diffusion rate using TST + quantum tunneling.
    
    Args:
        site_type: 0=void, 1=physisorption, 2=chemisorption, 3=defect
        temperature_k: Temperature (K)
    
    Returns:
        Diffusion rate (s^-1)
    """
    if site_type == 0:  # void
        return 0.0
    
    # Get appropriate barrier and prefactor
    if site_type == 1:  # physisorption
        barrier = H_DIFFUSION_BARRIERS["amorphous_carbon"]
        prefactor = TST_PREFACTORS["diffusion"]
    elif site_type == 2:  # chemisorption
        barrier = H_DIFFUSION_BARRIERS["amorphous_carbon"]
        prefactor = TST_PREFACTORS["diffusion"]
    elif site_type == 3:  # defect
        barrier = H_DIFFUSION_BARRIERS["defect"]
        prefactor = TST_PREFACTORS["diffusion"]
    else:
        return 0.0
    
    # Calculate thermal and tunneling rates
    thermal = thermal_rate(prefactor, barrier, temperature_k)
    tunneling = quantum_tunneling_rate(
        barrier, 
        TUNNELING_PARAMS["barrier_width_angstroms"],
        TUNNELING_PARAMS["effective_mass"],
        temperature_k
    )
    
    return combined_rate(thermal, tunneling)


def h_desorption_rate(binding_energy_eV, temperature_k):
    """
    Calculate H desorption rate using TST.
    
    Args:
        binding_energy_eV: Binding energy (eV)
        temperature_k: Temperature (K)
    
    Returns:
        Desorption rate (s^-1)
    """
    prefactor = TST_PREFACTORS["desorption"]
    return thermal_rate(prefactor, binding_energy_eV, temperature_k)


def h2_formation_lh_rate(temperature_k, adjacent_pairs_count):
    """
    Calculate H2 formation rate via Langmuir-Hinshelwood mechanism.
    
    Args:
        temperature_k: Temperature (K)
        adjacent_pairs_count: Number of adjacent H atom pairs
    
    Returns:
        H2 formation rate (s^-1)
    """
    if adjacent_pairs_count <= 0:
        return 0.0
    
    # LH formation is often barrierless (Cuppen et al. 2013)
    barrier = H2_FORMATION_BARRIERS["LH_barrierless"]
    prefactor = TST_PREFACTORS["reaction"]
    
    # Rate per pair
    rate_per_pair = thermal_rate(prefactor, barrier, temperature_k)
    
    return rate_per_pair * adjacent_pairs_count


def h2_formation_er_rate(gas_flux_cm2_s, h_atoms_count):
    """
    Calculate H2 formation rate via Eley-Rideal mechanism.
    
    Args:
        gas_flux_cm2_s: Gas flux (atoms cm^-2 s^-1)
        h_atoms_count: Number of H atoms on surface
    
    Returns:
        H2 formation rate (s^-1)
    """
    if h_atoms_count <= 0:
        return 0.0
    
    # ER cross-section from literature
    cross_section_cm2 = 1e-15  # cm^2
    
    # Reaction probability (can be temperature dependent)
    reaction_probability = 0.1  # 10% probability
    
    return gas_flux_cm2_s * cross_section_cm2 * reaction_probability * h_atoms_count


def uv_photodesorption_rate(uv_flux_photons_cm2_s, h_atoms_count):
    """
    Calculate UV photodesorption rate.
    
    Args:
        uv_flux_photons_cm2_s: UV photon flux (photons cm^-2 s^-1)
        h_atoms_count: Number of H atoms on surface
    
    Returns:
        Photodesorption rate (s^-1)
    """
    if h_atoms_count <= 0:
        return 0.0
    
    cross_section = UV_PARAMS["absorption_cross_section_cm2"]
    yield_param = UV_PARAMS["photodesorption_yield"]
    
    return uv_flux_photons_cm2_s * cross_section * yield_param * h_atoms_count


def uv_h2_formation_rate(uv_flux_photons_cm2_s, adjacent_pairs_count, site_area_cm2):
    """
    Calculate UV-assisted H2 formation rate.
    
    Args:
        uv_flux_photons_cm2_s: UV photon flux (photons cm^-2 s^-1)
        adjacent_pairs_count: Number of adjacent H atom pairs
        site_area_cm2: Site area (cm^2)
    
    Returns:
        UV H2 formation rate (s^-1)
    """
    if adjacent_pairs_count <= 0:
        return 0.0
    
    cross_section = UV_PARAMS["absorption_cross_section_cm2"]
    yield_param = UV_PARAMS["photodissociation_yield"]  # Very low for H2
    
    return uv_flux_photons_cm2_s * cross_section * yield_param * adjacent_pairs_count * site_area_cm2


def adsorption_rate(gas_density_cm3, temperature_k, sticking_probability, accessible_area_cm2):
    """
    Calculate adsorption rate from gas phase.
    Based on Hollenbach & McKee 1979 for H sticking on carbonaceous grains.
    
    Args:
        gas_density_cm3: Gas density (cm^-3)
        temperature_k: Gas temperature (K)
        sticking_probability: Sticking coefficient
        accessible_area_cm2: Accessible surface area (cm^2)
    
    Returns:
        Adsorption rate (s^-1)
    """
    v_thermal = np.sqrt(8 * K_B * temperature_k * 1.602e-19 / (np.pi * M_H))
    
    gas_flux = 0.25 * gas_density_cm3 * v_thermal
    
    if temperature_k > 250:
        E_ads_barrier = 0.01
    else:
        E_ads_barrier = 0.0
    temp_factor = np.exp(-E_ads_barrier / (K_B * temperature_k))
    effective_sticking = sticking_probability * temp_factor
    
    return gas_flux * accessible_area_cm2 * effective_sticking
