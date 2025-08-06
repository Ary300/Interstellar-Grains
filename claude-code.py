#This could be useful code from Claude, just to help you Julia 
import numpy as np
import random
from scientific_data import surface_chemistry_data, uv_photon_flux, EV_TO_KELVIN, K_B, M_H

class EnhancedKineticMonteCarlo:
    """
    Enhanced KMC simulation for H2 formation on carbonaceous grains
    Includes high-temperature pathways and PAH surface chemistry
    """
    def __init__(self, simulation_parameters):
        self.simulation_parameters = simulation_parameters
        self.time = 0.0
        self.h2_molecules_formed = 0
        self.h2_molecules_formed_LH = 0
        self.h2_molecules_formed_ER = 0
        self.h2_molecules_formed_UV = 0
        self.h2_molecules_formed_chemi = 0  # NEW: Chemisorption pathway
        self.total_adsorbed_h_atoms = 0
        self.total_desorbed_h_atoms = 0
        
        # Enhanced surface chemistry parameters
        self.use_site_heterogeneity = self.simulation_parameters.get("use_site_heterogeneity", True)
        self.enable_high_temp_pathway = self.simulation_parameters.get("enable_high_temp_pathway", True)
        self.enable_pah_charging = self.simulation_parameters.get("enable_pah_charging", True)
        self.pah_charge_state = 0.0  # Average charge per grain
        
        # Temperature-dependent surface properties
        self.surface_temp = self.simulation_parameters.get("surface_temperature_k", 10.0)
        self.dynamic_surface = self.simulation_parameters.get("dynamic_surface", True)
        
        self.initialize_enhanced_surface()
        
    def initialize_enhanced_surface(self):
        """Initialize surface with enhanced carbonaceous grain properties"""
        grain_radius_um = self.simulation_parameters.get("grain_radius_um", 0.1)
        grain_radius_cm = grain_radius_um * 1e-4
        site_area_angstroms_sq = self.simulation_parameters.get("site_area_angstroms_sq", 9)
        site_area_cm2 = site_area_angstroms_sq * 1e-16
        grain_surface_area_cm2 = 4 * np.pi * (grain_radius_cm**2)
        calculated_sites = int(grain_surface_area_cm2 / site_area_cm2)
        
        self.surface_dimension = int(np.sqrt(calculated_sites))
        self.total_surface_sites = self.surface_dimension * self.surface_dimension
        self.surface_sites = np.full((self.surface_dimension, self.surface_dimension), None, dtype=object)
        
        # Enhanced site types for carbonaceous grains
        self.site_types = np.random.choice(['edge', 'bridge', 'hollow', 'defect'], 
                                          size=(self.surface_dimension, self.surface_dimension),
                                          p=[0.3, 0.4, 0.2, 0.1])  # PAH-like site distribution
        
        if self.use_site_heterogeneity:
            self._initialize_heterogeneous_sites()
        else:
            self._initialize_uniform_sites()
            
        # Initialize with some coverage if specified
        initial_h_coverage = self.simulation_parameters.get("initial_h_coverage", 0.0)
        num_initial_h = int(self.total_surface_sites * initial_h_coverage)
        self.h_atoms_on_surface = 0
        
        if num_initial_h > 0:
            available_indices = list(np.ndindex(self.surface_sites.shape))
            random.shuffle(available_indices)
            for i in range(num_initial_h):
                r, c = available_indices[i]
                self.surface_sites[r, c] = "H"
                self.h_atoms_on_surface += 1
                
        self.adjacent_h_pairs_count = 0
        self._update_all_adjacent_pairs()

    def _initialize_heterogeneous_sites(self):
        """Initialize heterogeneous binding energies for different site types"""
        rng = np.random.default_rng()
        
        # Site-dependent binding energies (carbonaceous grain specific)
        binding_energies = {
            'edge': 45,    # meV - weaker binding at edges
            'bridge': 60,  # meV - typical physisorption
            'hollow': 75,  # meV - stronger in hollow sites
            'defect': 120  # meV - very strong at defects (enables high-T formation)
        }
        
        # Diffusion barriers
        diffusion_barriers = {
            'edge': 0.035,    # eV
            'bridge': 0.045,  # eV  
            'hollow': 0.055,  # eV
            'defect': 0.025   # eV - easier diffusion from defects (paradoxical but real)
        }
        
        self.E_bind_eV_map = np.zeros((self.surface_dimension, self.surface_dimension))
        self.E_diff_eV_map = np.zeros((self.surface_dimension, self.surface_dimension))
        self.E_chemi_eV_map = np.zeros((self.surface_dimension, self.surface_dimension))  # NEW
        
        for site_type in ['edge', 'bridge', 'hollow', 'defect']:
            mask = (self.site_types == site_type)
            
            # Physisorption energies with heterogeneity
            base_bind = binding_energies[site_type] / 1000.0  # Convert to eV
            sigma = 0.005  # eV spread
            self.E_bind_eV_map[mask] = rng.normal(base_bind, sigma, size=np.sum(mask))
            
            # Diffusion barriers
            base_diff = diffusion_barriers[site_type]
            self.E_diff_eV_map[mask] = rng.normal(base_diff, 0.003, size=np.sum(mask))
            
            # Chemisorption energies (NEW - enables high-T pathway)
            if site_type == 'defect':
                self.E_chemi_eV_map[mask] = rng.normal(0.8, 0.1, size=np.sum(mask))  # Strong chemisorption
            elif site_type == 'edge':
                self.E_chemi_eV_map[mask] = rng.normal(0.5, 0.05, size=np.sum(mask))  # Moderate
            else:
                self.E_chemi_eV_map[mask] = rng.normal(0.3, 0.03, size=np.sum(mask))  # Weak
        
        # Ensure physical constraints
        self.E_bind_eV_map = np.clip(self.E_bind_eV_map, 0.01, 0.15)
        self.E_diff_eV_map = np.clip(self.E_diff_eV_map, 0.01, 0.08)
        self.E_chemi_eV_map = np.clip(self.E_chemi_eV_map, 0.2, 1.2)

    def _initialize_uniform_sites(self):
        """Fallback to uniform site properties"""
        base_bind = surface_chemistry_data["h_physisorption_binding_energy_mev_typical"] / 1000.0
        base_diff = surface_chemistry_data["h_physisorbed_amorphous_diffusion_barrier_ev_min"]
        
        self.E_bind_eV_map = np.full((self.surface_dimension, self.surface_dimension), base_bind)
        self.E_diff_eV_map = np.full((self.surface_dimension, self.surface_dimension), base_diff)
        self.E_chemi_eV_map = np.full((self.surface_dimension, self.surface_dimension), 0.4)

    def _update_all_adjacent_pairs(self):
        """Recalculate all adjacent H pairs count"""
        self.adjacent_h_pairs_count = 0
        H_positions = np.argwhere(self.surface_sites == "H")
        
        for r1, c1 in H_positions:
            for r2, c2 in self.get_neighbors(r1, c1):
                if self.surface_sites[r2, c2] == "H" and (r1, c1) < (r2, c2):
                    self.adjacent_h_pairs_count += 1

    def get_neighbors(self, r, c):
        """Get neighboring sites with periodic boundary conditions"""
        neighbors = []
        moves = [(-1, 0), (1, 0), (0, -1), (0, 1)]
        for dr, dc in moves:
            nr, nc = (r + dr) % self.surface_dimension, (c + dc) % self.surface_dimension
            neighbors.append((nr, nc))
        return neighbors

    def update_pah_charge_state(self):
        """Update PAH charge state based on UV flux and temperature"""
        if not self.enable_pah_charging:
            return
            
        uv_flux_factor = self.simulation_parameters.get("uv_flux_factor", 1.0)
        
        # Simplified charging model
        # Positive charging from UV photons, neutralization from electrons
        photoionization_rate = uv_flux_factor * 1e-3  # s^-1
        recombination_rate = 1e-7 * np.exp(-0.1 * 11604 / self.surface_temp)  # s^-1
        
        # Steady state charge (simplified)
        self.pah_charge_state = photoionization_rate / (photoionization_rate + recombination_rate)

    def calculate_enhanced_rates(self):
        """Calculate all reaction rates including high-temperature pathways"""
        rates = {}
        surface_temp_k = self.surface_temp
        gas_temp_k = self.simulation_parameters.get("gas_temperature_k")
        h_gas_density = self.simulation_parameters.get("h_gas_density_cm3")
        uv_flux_factor = self.simulation_parameters.get("uv_flux_factor", 1.0)
        site_area_cm2 = self.simulation_parameters.get("site_area_angstroms_sq", 9) * 1e-16
        pre_exp_frequency = 1e12

        # Update charge state
        self.update_pah_charge_state()

        # Enhanced adsorption (charge-dependent sticking)
        v_h_thermal = np.sqrt(8 * K_B * gas_temp_k / (np.pi * M_H))
        gas_flux_per_cm2_s = 0.25 * h_gas_density * v_h_thermal
        
        # Charge-enhanced sticking probability
        base_sticking = self.simulation_parameters.get("sticking_probability", 0.3)
        charge_enhancement = 1.0 + 0.5 * abs(self.pah_charge_state)  # Charged grains stick better
        effective_sticking = min(1.0, base_sticking * charge_enhancement)
        
        rates["adsorption"] = gas_flux_per_cm2_s * (self.total_surface_sites * site_area_cm2) * effective_sticking

        # Enhanced desorption (site-dependent)
        if self.h_atoms_on_surface > 0:
            occ_positions = np.argwhere(self.surface_sites == "H")
            total_desorption_rate = 0.0
            
            for r, c in occ_positions:
                # Temperature-dependent pathway selection
                E_phys = self.E_bind_eV_map[r, c]
                E_chemi = self.E_chemi_eV_map[r, c]
                
                # Physisorption desorption
                rate_phys = pre_exp_frequency * np.exp(-E_phys * EV_TO_KELVIN / surface_temp_k)
                
                # Chemisorption desorption (high-T pathway)
                if self.enable_high_temp_pathway and surface_temp_k > 30.0:
                    # Transition probability from physi to chemi
                    transition_prob = 1.0 / (1.0 + np.exp(-(surface_temp_k - 50.0) / 10.0))
                    rate_chemi = pre_exp_frequency * np.exp(-E_chemi * EV_TO_KELVIN / surface_temp_k)
                    
                    # Effective rate combines both pathways
                    effective_rate = (1 - transition_prob) * rate_phys + transition_prob * rate_chemi
                else:
                    effective_rate = rate_phys
                
                total_desorption_rate += effective_rate
            
            rates["desorption"] = float(total_desorption_rate)

        # Enhanced diffusion
        if self.h_atoms_on_surface > 0:
            occ_positions = np.argwhere(self.surface_sites == "H")
            total_diffusion_rate = 0.0
            
            for r, c in occ_positions:
                # Check if diffusion is possible
                empty_neighbors = [pos for pos in self.get_neighbors(r, c) 
                                 if self.surface_sites[pos] is None]
                if empty_neighbors:
                    E_diff = self.E_diff_eV_map[r, c]
                    diff_rate = pre_exp_frequency * np.exp(-E_diff * EV_TO_KELVIN / surface_temp_k)
                    total_diffusion_rate += diff_rate
            
            rates["diffusion"] = float(total_diffusion_rate)

        # Enhanced H2 formation mechanisms
        
        # 1. Langmuir-Hinshelwood (thermal)
        if self.adjacent_h_pairs_count > 0:
            pairs = self._get_adjacent_pairs()
            total_lh_rate = 0.0
            
            for (r1, c1), (r2, c2) in pairs:
                # Use average diffusion barrier as reaction barrier
                E_barrier = 0.5 * (self.E_diff_eV_map[r1, c1] + self.E_diff_eV_map[r2, c2])
                lh_rate = pre_exp_frequency * np.exp(-E_barrier * EV_TO_KELVIN / surface_temp_k)
                total_lh_rate += lh_rate
            
            rates["h2_formation_LH"] = float(total_lh_rate)

        # 2. Eley-Rideal (gas-surface)
        if self.h_atoms_on_surface > 0:
            er_cross_section = surface_chemistry_data["er_cross_section_cm2"]
            # Enhanced by charge state
            effective_cross_section = er_cross_section * (1.0 + 0.3 * abs(self.pah_charge_state))
            rates["h2_formation_ER"] = gas_flux_per_cm2_s * effective_cross_section * self.h_atoms_on_surface

        # 3. NEW: High-temperature chemisorption pathway
        if (self.enable_high_temp_pathway and self.adjacent_h_pairs_count > 0 and 
            surface_temp_k > 30.0):
            
            pairs = self._get_adjacent_pairs()
            total_chemi_rate = 0.0
            
            for (r1, c1), (r2, c2) in pairs:
                # Check if sites support chemisorption
                site1_type = self.site_types[r1, c1]
                site2_type = self.site_types[r2, c2]
                
                if site1_type in ['defect', 'edge'] or site2_type in ['defect', 'edge']:
                    # Lower barrier for chemisorbed pathway
                    E_chemi_barrier = 0.3  # eV, much lower than desorption
                    chemi_rate = pre_exp_frequency * np.exp(-E_chemi_barrier * EV_TO_KELVIN / surface_temp_k)
                    
                    # Temperature-dependent activation
                    activation = 1.0 / (1.0 + np.exp(-(surface_temp_k - 80.0) / 20.0))
                    total_chemi_rate += chemi_rate * activation
            
            rates["h2_formation_chemi"] = float(total_chemi_rate)

        # 4. UV-assisted pathways
        if uv_flux_factor > 0:
            uv_photon_flux_total = uv_photon_flux["integrated_fuv_photon_flux_photons_cm2_s"] * uv_flux_factor
            
            # UV photodesorption
            if self.h_atoms_on_surface > 0:
                uv_yield = surface_chemistry_data["uv_photodesorption_yield_h_atom"]
                # Charge-dependent enhancement
                effective_yield = uv_yield * (1.0 + 0.2 * abs(self.pah_charge_state))
                rates["uv_photodesorption"] = (uv_photon_flux_total * site_area_cm2 * 
                                             effective_yield * self.h_atoms_on_surface)
            
            # UV-assisted H2 formation
            if self.adjacent_h_pairs_count > 0:
                uv_h2_yield = surface_chemistry_data["uv_h2_formation_yield_per_pair"]
                effective_h2_yield = uv_h2_yield * (1.0 + 0.3 * abs(self.pah_charge_state))
                rates["h2_formation_UV"] = (uv_photon_flux_total * (2 * site_area_cm2) * 
                                          effective_h2_yield * self.adjacent_h_pairs_count)

        return {k: v for k, v in rates.items() if v > 0}

    def _get_adjacent_pairs(self):
        """Get all adjacent H atom pairs"""
        pairs = []
        H_positions = np.argwhere(self.surface_sites == "H")
        
        for r1, c1 in H_positions:
            for r2, c2 in self.get_neighbors(r1, c1):
                if self.surface_sites[r2, c2] == "H" and (r1, c1) < (r2, c2):
                    pairs.append(((r1, c1), (r2, c2)))
        return pairs

    def execute_enhanced_event(self, event_type):
        """Execute events including new high-temperature pathways"""
        
        if event_type == "adsorption":
            if self.h_atoms_on_surface < self.total_surface_sites:
                empty_sites_coords = np.argwhere(self.surface_sites == None)
                if len(empty_sites_coords) > 0:
                    r, c = random.choice(empty_sites_coords.tolist())
                    self.surface_sites[r, c] = "H"
                    self.h_atoms_on_surface += 1
                    self.total_adsorbed_h_atoms += 1
                    self._update_adjacent_h_pairs_local(r, c, True)

        elif event_type in ["desorption", "uv_photodesorption"]:
            if self.h_atoms_on_surface > 0:
                occupied_sites_coords = np.argwhere(self.surface_sites == "H").tolist()
                r, c = random.choice(occupied_sites_coords)
                self.surface_sites[r, c] = None
                self.h_atoms_on_surface -= 1
                self.total_desorbed_h_atoms += 1
                self._update_adjacent_h_pairs_local(r, c, False)

        elif event_type == "diffusion":
            if self.h_atoms_on_surface > 0:
                occupied_sites_coords = np.argwhere(self.surface_sites == "H").tolist()
                mobile_atoms = []
                
                for pos in occupied_sites_coords:
                    empty_neighbors = [n_pos for n_pos in self.get_neighbors(*pos) 
                                     if self.surface_sites[n_pos] is None]
                    if empty_neighbors:
                        mobile_atoms.append(pos)
                
                if mobile_atoms:
                    r_h, c_h = random.choice(mobile_atoms)
                    empty_neighbors = [n_pos for n_pos in self.get_neighbors(r_h, c_h) 
                                     if self.surface_sites[n_pos] is None]
                    if empty_neighbors:
                        r_empty, c_empty = random.choice(empty_neighbors)
                        
                        self.surface_sites[r_h, c_h] = None
                        self._update_adjacent_h_pairs_local(r_h, c_h, False)
                        self.surface_sites[r_empty, c_empty] = "H"
                        self._update_adjacent_h_pairs_local(r_empty, c_empty, True)

        elif event_type == "h2_formation_LH":
            self._execute_h2_formation("LH")
        elif event_type == "h2_formation_ER":
            self._execute_h2_formation("ER")
        elif event_type == "h2_formation_UV":
            self._execute_h2_formation("UV")
        elif event_type == "h2_formation_chemi":  # NEW
            self._execute_h2_formation("chemi")

    def _execute_h2_formation(self, mechanism):
        """Execute H2 formation by specified mechanism"""
        if mechanism == "ER":
            if self.h_atoms_on_surface > 0:
                occupied_sites_coords = np.argwhere(self.surface_sites == "H").tolist()
                r, c = random.choice(occupied_sites_coords)
                self.surface_sites[r, c] = None
                self.h_atoms_on_surface -= 1
                self._update_adjacent_h_pairs_local(r, c, False)
                self.h2_molecules_formed_ER += 1
                self.h2_molecules_formed += 1
        
        else:  # LH, UV, or chemi - all require adjacent pairs
            if self.adjacent_h_pairs_count > 0:
                pairs = self._get_adjacent_pairs()
                if pairs:
                    (r1, c1), (r2, c2) = random.choice(pairs)
                    
                    # Remove both H atoms
                    self.surface_sites[r1, c1] = None
                    self._update_adjacent_h_pairs_local(r1, c1, False)
                    self.surface_sites[r2, c2] = None
                    self._update_adjacent_h_pairs_local(r2, c2, False)
                    
                    self.h_atoms_on_surface -= 2
                    self.h2_molecules_formed += 1
                    
                    # Track mechanism
                    if mechanism == "LH":
                        self.h2_molecules_formed_LH += 1
                    elif mechanism == "UV":
                        self.h2_molecules_formed_UV += 1
                    elif mechanism == "chemi":
                        self.h2_molecules_formed_chemi += 1

    def _update_adjacent_h_pairs_local(self, r, c, add_atom):
        """Update adjacent H pairs count for local changes"""
        change = 1 if add_atom else -1
        for nr, nc in self.get_neighbors(r, c):
            if self.surface_sites[nr, nc] == "H":
                self.adjacent_h_pairs_count += change

    def run_enhanced_gillespie(self, max_time, max_steps=None):
        """Run Gillespie algorithm with enhanced rate calculations"""
        step_count = 0
        
        while self.time < max_time:
            if max_steps and step_count >= max_steps:
                break
                
            rates = self.calculate_enhanced_rates()
            
            if not rates:
                # No possible events - simulation stalled
                break
                
            total_rate = sum(rates.values())
            if total_rate == 0:
                break
            
            # Time advancement
            try:
                delta_t = random.expovariate(total_rate)
            except (ValueError, OverflowError):
                # Handle numerical issues
                delta_t = 1e-12
            
            if self.time + delta_t > max_time:
                delta_t = max_time - self.time
            
            self.time += delta_t
            
            # Event selection
            rand_val = random.random() * total_rate
            cumulative = 0.0
            chosen_event = None
            
            for event, rate in rates.items():
                cumulative += rate
                if rand_val <= cumulative:
                    chosen_event = event
                    break
            
            if chosen_event:
                self.execute_enhanced_event(chosen_event)
            
            step_count += 1
            
            # Optional: Dynamic temperature updates
            if self.dynamic_surface and step_count % 1000 == 0:
                self._update_surface_temperature()
        
        return []

    def _update_surface_temperature(self):
        """Update surface temperature due to UV heating (optional)"""
        if not self.simulation_parameters.get("enable_dynamic_heating", False):
            return
            
        uv_flux_factor = self.simulation_parameters.get("uv_flux_factor", 1.0)
        base_temp = self.simulation_parameters.get("surface_temperature_k", 10.0)
        
        # Simple heating model
        heating_increment = uv_flux_factor * 0.1  # K per UV flux unit
        self.surface_temp = min(base_temp + heating_increment, 300.0)  # Cap at 300K

    def get_results_summary(self):
        """Return comprehensive results including new mechanisms"""
        return {
            "final_time": self.time,
            "final_h_atoms_on_surface": self.h_atoms_on_surface,
            "h2_formed_LH": self.h2_molecules_formed_LH,
            "h2_formed_ER": self.h2_molecules_formed_ER,
            "h2_formed_UV": self.h2_molecules_formed_UV,
            "h2_formed_chemi": self.h2_molecules_formed_chemi,  # NEW
            "total_h2_formed": self.h2_molecules_formed,
            "total_adsorbed_h": self.total_adsorbed_h_atoms,
            "total_desorbed_h": self.total_desorbed_h_atoms,
            "h2_formation_efficiency": (self.h2_molecules_formed / max(self.total_adsorbed_h_atoms, 1)),
            "final_pah_charge": self.pah_charge_state,
            "final_surface_temp": self.surface_temp
        }

# Usage example:
if __name__ == "__main__":
    # Test the enhanced simulation
    test_params = {
        "grain_radius_um": 0.1,
        "site_area_angstroms_sq": 9,
        "surface_temperature_k": 100.0,  # High temperature test
        "gas_temperature_k": 100.0,
        "h_gas_density_cm3": 1000.0,
        "sticking_probability": 0.3,
        "initial_h_coverage": 0.1,
        "uv_flux_factor": 1.0,
        "enable_high_temp_pathway": True,
        "enable_pah_charging": True,
        "dynamic_surface": False
    }
    
    enhanced_kmc = EnhancedKineticMonteCarlo(test_params)
    enhanced_kmc.run_enhanced_gillespie(max_time=1e5)  # 100,000 seconds
    
    results = enhanced_kmc.get_results_summary()
    print("Enhanced KMC Results:")
    for key, value in results.items():
        print(f"{key}: {value}")
