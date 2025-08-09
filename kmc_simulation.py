import numpy as np
import random
from scientific_data import surface_chemistry_data, uv_photon_flux, EV_TO_KELVIN, K_B, M_H
from physical_rates import (
    h_diffusion_rate, h_desorption_rate, h2_formation_lh_rate, h2_formation_er_rate,
    uv_photodesorption_rate, uv_h2_formation_rate, adsorption_rate
)

class KineticMonteCarlo:
    def __init__(self, simulation_parameters):
        self.simulation_parameters = simulation_parameters
        self.time = 0.0
        self.h2_molecules_formed = 0
        self.h2_molecules_formed_LH = 0
        self.h2_molecules_formed_ER = 0
        self.h2_molecules_formed_UV = 0
        self.total_adsorbed_h_atoms = 0
        self.total_desorbed_h_atoms = 0
        
        self.use_3d_lattice = self.simulation_parameters.get("use_3d_lattice", True)
        self.porosity_fraction = self.simulation_parameters.get("porosity_fraction", 0.2)
        
        self.uv_pulse_enabled = self.simulation_parameters.get("uv_pulse_enabled", True)
        self.uv_defect_creation_rate = self.simulation_parameters.get("uv_defect_creation_rate", 0.5)
        self.uv_pulse_duration = self.simulation_parameters.get("uv_pulse_duration", 1e-6)
        self.last_uv_pulse_time = 0.0
        self.uv_pulse_active = False
        
        self.E_bind_mean_meV = float(self.simulation_parameters.get("E_phys_mean_meV", 45.0))
        self.E_bind_sigma_meV = float(self.simulation_parameters.get("heterogeneity_E_bind_sigma_meV", 5.0))
        
        self.initialize_3d_lattice()

    def initialize_3d_lattice(self):
        grain_radius_um = self.simulation_parameters.get("grain_radius_um", 0.1)
        grain_radius_cm = grain_radius_um * 1e-4
        site_area_angstroms_sq = self.simulation_parameters.get("site_area_angstroms_sq", 9)
        site_area_cm2 = site_area_angstroms_sq * 1e-16
        grain_surface_area_cm2 = 4 * np.pi * (grain_radius_cm**2)
        
        calculated_sites = int(grain_surface_area_cm2 / site_area_cm2)
        self.surface_dimension = int(np.sqrt(calculated_sites))
        self.depth_layers = max(3, self.surface_dimension // 10)
        
        self.lattice = np.full((self.depth_layers, self.surface_dimension, self.surface_dimension), "C", dtype=object)
        
        self.E_bind_eV_map = np.zeros((self.depth_layers, self.surface_dimension, self.surface_dimension))
        self.E_diff_eV_map = np.zeros((self.depth_layers, self.surface_dimension, self.surface_dimension))
        
        self._generate_amorphous_structure()
        self._generate_energy_maps()
        
        initial_h_coverage = self.simulation_parameters.get("initial_h_coverage", 0.0)
        self.h_atoms_on_surface = 0
        self.adjacent_h_pairs_count = 0
        
        if initial_h_coverage > 0:
            self._initialize_h_atoms(initial_h_coverage)

    def _generate_amorphous_structure(self):
        rng = np.random.default_rng()
        porosity_mask = rng.random(self.lattice.shape) < self.porosity_fraction
        self.lattice[porosity_mask] = None
        
        surface_sites = self.lattice[0, :, :]
        if not np.any(surface_sites != None):
            num_accessible = max(1, int(self.surface_dimension * self.surface_dimension * 0.1))
            accessible_indices = np.random.choice(self.surface_dimension * self.surface_dimension, 
                                                num_accessible, replace=False)
            for idx in accessible_indices:
                r, c = idx // self.surface_dimension, idx % self.surface_dimension
                self.lattice[0, r, c] = "C"

    def _generate_energy_maps(self):
        rng = np.random.default_rng()
        
        # Generate site types: 0=void, 1=physisorption, 2=chemisorption, 3=defect
        self.site_types = np.zeros(self.lattice.shape, dtype=int)
        self.site_types = np.where(self.lattice != None, 1, 0)  # Default to physisorption
        
        # Add chemisorption sites (10% of accessible sites)
        chemisorption_fraction = self.simulation_parameters.get("chemisorption_fraction", 0.1)
        accessible_mask = (self.lattice != None)
        num_accessible = np.sum(accessible_mask)
        num_chemisorption = int(num_accessible * chemisorption_fraction)
        
        if num_chemisorption > 0:
            accessible_indices = np.where(accessible_mask)
            chemisorption_indices = rng.choice(len(accessible_indices[0]), num_chemisorption, replace=False)
            for idx in chemisorption_indices:
                d, r, c = accessible_indices[0][idx], accessible_indices[1][idx], accessible_indices[2][idx]
                self.site_types[d, r, c] = 2  # Chemisorption
        
        # Add defect sites (15% of accessible sites)
        defect_fraction = self.simulation_parameters.get("surface_defect_fraction", 0.15)
        num_defects = int(num_accessible * defect_fraction)
        
        if num_defects > 0:
            defect_indices = rng.choice(len(accessible_indices[0]), num_defects, replace=False)
            for idx in defect_indices:
                d, r, c = accessible_indices[0][idx], accessible_indices[1][idx], accessible_indices[2][idx]
                if self.site_types[d, r, c] == 1:  # Only convert physisorption sites to defects
                    self.site_types[d, r, c] = 3  # Defect
        
        mean_eV = self.E_bind_mean_meV / 1000.0
        sigma_eV = self.E_bind_sigma_meV / 1000.0
        
        E_bind_distribution = rng.normal(mean_eV, sigma_eV, self.lattice.shape)
        
        self.E_bind_eV_map = np.where(self.lattice != None, E_bind_distribution, 0.0)
        self.E_bind_eV_map = np.clip(self.E_bind_eV_map, 0.0, None)
        
        self.E_diff_eV_map = 0.3 * self.E_bind_eV_map
        self.E_diff_eV_map = np.where(self.lattice != None, self.E_diff_eV_map, 0.0)

    def _initialize_h_atoms(self, initial_coverage):
        accessible_sites = np.where(self.lattice[0, :, :] != None)
        if len(accessible_sites[0]) > 0:
            num_initial_h = int(len(accessible_sites[0]) * initial_coverage)
            if num_initial_h > 0:
                indices = random.sample(range(len(accessible_sites[0])), min(num_initial_h, len(accessible_sites[0])))
                for idx in indices:
                    r, c = accessible_sites[0][idx], accessible_sites[1][idx]
                    self.lattice[0, r, c] = "H"
                    self.h_atoms_on_surface += 1
                self._update_adjacent_h_pairs_count()

    def get_neighbors_3d(self, d, r, c):
        neighbors = []
        moves_2d = [(-1, 0), (1, 0), (0, -1), (0, 1)]
        for dr, dc in moves_2d:
            nr, nc = (r + dr) % self.surface_dimension, (c + dc) % self.surface_dimension
            if self.lattice[d, nr, nc] is not None:
                neighbors.append((d, nr, nc))
        
        if d < self.depth_layers - 1 and self.lattice[d + 1, r, c] is not None:
            neighbors.append((d + 1, r, c))
        
        if d > 0 and self.lattice[d - 1, r, c] is not None:
            neighbors.append((d - 1, r, c))
        
        return neighbors

    def _update_adjacent_h_pairs_count(self):
        self.adjacent_h_pairs_count = 0
        for d, r, c in self.get_occupied_sites():
            neighbors = self.get_neighbors_3d(d, r, c)
            for nd, nr, nc in neighbors:
                if self.lattice[nd, nr, nc] == "H" and (d, r, c) < (nd, nr, nc):
                    self.adjacent_h_pairs_count += 1

    def update_adjacent_h_pairs_count(self, d, r, c, add_atom):
        change = 1 if add_atom else -1
        neighbors = self.get_neighbors_3d(d, r, c)
        for nd, nr, nc in neighbors:
            if self.lattice[nd, nr, nc] == "H":
                self.adjacent_h_pairs_count += change

    def get_accessible_surface_sites(self):
        return np.where((self.lattice[0, :, :] != "H") & (self.lattice[0, :, :] != None))

    def get_occupied_sites(self):
        occupied_mask = (self.lattice == "H")
        if not np.any(occupied_mask):
            return []
        return list(zip(*np.where(occupied_mask)))

    def _calculate_diffusion_rate_with_tunneling(self, d, r, c, surface_temp_k):
        hbar = 1.055e-27
        mH = 1.674e-24
        kB = 1.381e-16
        eV_to_erg = 1.602e-12
        
        E_diff_eV = self.E_diff_eV_map[d, r, c]
        if E_diff_eV <= 0:
            return 0.0
        E_diff_erg = E_diff_eV * eV_to_erg
        
        a = 3.0e-8
        nu0 = 1e12
        
        classical_rate = nu0 * np.exp(-E_diff_erg / (kB * surface_temp_k))
        
        tunneling_prob = np.exp(-2 * np.sqrt(2 * mH * E_diff_erg) * a / hbar)
        tunneling_rate = nu0 * tunneling_prob
        
        return classical_rate + tunneling_rate

    def calculate_rates(self):
        """
        Calculate all rates using physically correct Transition State Theory
        with quantum tunneling corrections.
        """
        rates = {}
        surface_temp_k = self.simulation_parameters.get("surface_temperature_k")
        gas_temp_k = self.simulation_parameters.get("gas_temperature_k")
        h_gas_density = self.simulation_parameters.get("h_gas_density_cm3")
        uv_flux_factor = self.simulation_parameters.get("uv_flux_factor", 1.0)
        site_area_cm2 = self.simulation_parameters.get("site_area_angstroms_sq", 9) * 1e-16
        sticking_probability = self.simulation_parameters.get("sticking_probability", 0.3)

        accessible_sites = self.get_accessible_surface_sites()
        num_accessible_sites = len(accessible_sites[0])
        accessible_area_cm2 = num_accessible_sites * site_area_cm2

        # 1. Adsorption rate - physically derived from gas kinetics
        rates["adsorption"] = adsorption_rate(
            h_gas_density, gas_temp_k, sticking_probability, accessible_area_cm2
        )

        # 2. Eley-Rideal formation - gas impingement × cross-section × reaction probability
        if self.h_atoms_on_surface > 0:
            # Calculate gas flux for ER
            v_thermal = np.sqrt(8 * K_B * gas_temp_k / (np.pi * M_H * 1.602e-19))  # cm/s
            gas_flux_cm2_s = 0.25 * h_gas_density * v_thermal
            
            rates["h2_formation_ER"] = h2_formation_er_rate(
                gas_flux_cm2_s, self.h_atoms_on_surface
            )

        # 3. Desorption and diffusion rates - TST + quantum tunneling
        if self.h_atoms_on_surface > 0:
            occupied_sites = self.get_occupied_sites()
            if occupied_sites:
                total_desorption_rate = 0.0
                total_diffusion_rate = 0.0
                
                for d, r, c in occupied_sites:
                    # Get site-specific energetics
                    site_type = self.site_types[d, r, c]
                    binding_energy = self.E_bind_eV_map[d, r, c]
                    
                    # Desorption rate using TST
                    desorption_rate = h_desorption_rate(binding_energy, surface_temp_k)
                    total_desorption_rate += desorption_rate
                    
                    # Diffusion rate using TST + tunneling
                    diffusion_rate = h_diffusion_rate(site_type, surface_temp_k)
                    total_diffusion_rate += diffusion_rate
                
                rates["desorption"] = total_desorption_rate
                rates["diffusion"] = total_diffusion_rate

        # 4. Langmuir-Hinshelwood formation - TST for surface reactions
        if self.adjacent_h_pairs_count > 0:
            rates["h2_formation_LH"] = h2_formation_lh_rate(
                surface_temp_k, self.adjacent_h_pairs_count
            )

        # 5. UV processes - photon flux × cross-section × yield
        if uv_flux_factor > 0:
            uv_photon_flux_total = uv_photon_flux["integrated_fuv_photon_flux_photons_cm2_s"] * uv_flux_factor
            base_uv_rate = 5.0e-8
            
            # UV pulse rate (stochastic model)
            base_uv_rate = 5.0e-8  # 5 photons grain⁻¹ yr⁻¹
            if self.uv_pulse_enabled:
                rates["uv_pulse_start"] = base_uv_rate * uv_flux_factor
            
            if self.uv_pulse_active:
                if self.h_atoms_on_surface > 0:
                    # Photodesorption using physical cross-sections and yields
                    rates["uv_photodesorption"] = uv_photodesorption_rate(
                        uv_photon_flux_total, self.h_atoms_on_surface
                    )
                
                if self.adjacent_h_pairs_count > 0:
                    # UV-assisted H2 formation
                    rates["h2_formation_UV"] = uv_h2_formation_rate(
                        uv_photon_flux_total, self.adjacent_h_pairs_count, site_area_cm2
                    )
                
                # UV-induced defect creation (simplified model)
                accessible_surface_sites = self.get_accessible_surface_sites()
                if len(accessible_surface_sites[0]) > 0:
                    defect_creation_rate = self.uv_defect_creation_rate / (3.154e13)  # Convert Myr to s
                    rates["uv_defect_creation"] = defect_creation_rate * uv_flux_factor
                
                # UV-stimulated diffusion enhancement
                uv_diffusion_factor = float(self.simulation_parameters.get("uv_stimulated_diffusion_factor", 1.0))
                if uv_diffusion_factor > 1.0 and "diffusion" in rates:
                    rates["uv_stimulated_diffusion"] = (uv_diffusion_factor - 1.0) * rates["diffusion"]
        
        return {k: v for k, v in rates.items() if v > 0}

    def execute_event(self, event_type):
        if event_type == "adsorption":
            accessible_sites = self.get_accessible_surface_sites()
            if len(accessible_sites[0]) > 0:
                idx = random.randint(0, len(accessible_sites[0]) - 1)
                r, c = accessible_sites[0][idx], accessible_sites[1][idx]
                self.lattice[0, r, c] = "H"
                self.h_atoms_on_surface += 1
                self.total_adsorbed_h_atoms += 1
                self.update_adjacent_h_pairs_count(0, r, c, True)
                
        elif event_type in ["desorption", "uv_photodesorption"]:
            if self.h_atoms_on_surface > 0:
                occupied_sites = self.get_occupied_sites()
                if occupied_sites:
                    d, r, c = random.choice(occupied_sites)
                    self.lattice[d, r, c] = "C"
                    self.h_atoms_on_surface -= 1
                    self.total_desorbed_h_atoms += 1
                    self.update_adjacent_h_pairs_count(d, r, c, False)
                    
        elif event_type in ["diffusion", "uv_stimulated_diffusion"]:
            if self.h_atoms_on_surface > 0:
                occupied_sites = self.get_occupied_sites()
                if occupied_sites:
                    mobile_atoms = [
                        (d, r, c) for d, r, c in occupied_sites 
                        if any(self.lattice[nd, nr, nc] != "H" for nd, nr, nc in self.get_neighbors_3d(d, r, c))
                    ]
                    if mobile_atoms:
                        d_h, r_h, c_h = random.choice(mobile_atoms)
                        empty_neighbors = [
                            (nd, nr, nc) for nd, nr, nc in self.get_neighbors_3d(d_h, r_h, c_h) 
                            if self.lattice[nd, nr, nc] != "H"
                        ]
                        if empty_neighbors:
                            d_empty, r_empty, c_empty = random.choice(empty_neighbors)
                            self.lattice[d_h, r_h, c_h] = "C"
                            self.update_adjacent_h_pairs_count(d_h, r_h, c_h, False)
                            self.lattice[d_empty, r_empty, c_empty] = "H"
                            self.update_adjacent_h_pairs_count(d_empty, r_empty, c_empty, True)
                            
        elif event_type in ["h2_formation_LH", "h2_formation_UV"]:
            if self.adjacent_h_pairs_count > 0:
                pairs = []
                for d, r, c in self.get_occupied_sites():
                    for nd, nr, nc in self.get_neighbors_3d(d, r, c):
                        if self.lattice[nd, nr, nc] == "H" and (d, r, c) < (nd, nr, nc):
                            pairs.append(((d, r, c), (nd, nr, nc)))
                
                if pairs:
                    (d1, r1, c1), (d2, r2, c2) = random.choice(pairs)
                    self.lattice[d1, r1, c1] = "C"
                    self.update_adjacent_h_pairs_count(d1, r1, c1, False)
                    self.lattice[d2, r2, c2] = "C"
                    self.update_adjacent_h_pairs_count(d2, r2, c2, False)
                    self.h_atoms_on_surface -= 2
                    self.h2_molecules_formed += 1
                    if event_type == "h2_formation_LH":
                        self.h2_molecules_formed_LH += 1
                    else:
                        self.h2_molecules_formed_UV += 1
                                    
        elif event_type == "h2_formation_ER":
            if self.h_atoms_on_surface > 0:
                occupied_sites = self.get_occupied_sites()
                if occupied_sites:
                    d, r, c = random.choice(occupied_sites)
                    self.lattice[d, r, c] = "C"
                    self.h_atoms_on_surface -= 1
                    self.update_adjacent_h_pairs_count(d, r, c, False)
                    self.h2_molecules_formed_ER += 1
                    self.h2_molecules_formed += 1
                    
        elif event_type == "uv_pulse_start":
            self.uv_pulse_active = True
            self.last_uv_pulse_time = self.time

    def run_gillespie(self, max_time, max_steps=None):
        step_count = 0
        while self.time < max_time:
            if max_steps and step_count >= max_steps:
                break
            
            if self.uv_pulse_active and (self.time - self.last_uv_pulse_time) >= self.uv_pulse_duration:
                self.uv_pulse_active = False
            
            rates = self.calculate_rates()
            if not rates:
                break
            total_rate = sum(rates.values())
            if total_rate == 0:
                break
            
            delta_t = random.expovariate(total_rate)
            
            if self.time + delta_t > max_time:
                self.time = max_time
                break
            
            self.time += delta_t
            
            chosen_event = random.choices(list(rates.keys()), weights=list(rates.values()), k=1)[0]
            self.execute_event(chosen_event)
            step_count += 1
        return []

