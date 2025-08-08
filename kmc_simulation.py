import numpy as np
import random
from scientific_data import surface_chemistry_data, uv_photon_flux, EV_TO_KELVIN, K_B, M_H

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
        self.use_site_heterogeneity = self.simulation_parameters.get("use_site_heterogeneity", True)
        
        # 3D lattice parameters
        self.use_3d_lattice = self.simulation_parameters.get("use_3d_lattice", True)
        self.porosity_fraction = self.simulation_parameters.get("porosity_fraction", 0.2)
        self.surface_defect_fraction = self.simulation_parameters.get("surface_defect_fraction", 0.15)
        self.chemisorption_fraction = self.simulation_parameters.get("chemisorption_fraction", 0.1)
        
        # UV stochastic pulse parameters
        self.uv_pulse_enabled = self.simulation_parameters.get("uv_pulse_enabled", True)
        self.uv_defect_creation_rate = self.simulation_parameters.get("uv_defect_creation_rate", 0.5)  # sites/Myr
        self.uv_pulse_duration = self.simulation_parameters.get("uv_pulse_duration", 1e-6)  # s
        self.last_uv_pulse_time = 0.0
        self.uv_pulse_active = False
        
        # Energy parameters for different site types
        self.E_phys_mean_meV = float(self.simulation_parameters.get(
            "E_phys_mean_meV", 50.0  # ~500 K = 50 meV
        ))
        self.E_chem_mean_eV = float(self.simulation_parameters.get("E_chem_mean_eV", 1.75))  # 1.5-2.0 eV range
        self.E_bind_sigma_meV = float(self.simulation_parameters.get("heterogeneity_E_bind_sigma_meV", 5.0))
        self.E_diff_sigma_eV = float(self.simulation_parameters.get("heterogeneity_E_diff_sigma_eV", 0.005))
        
        self.initialize_3d_lattice()

    def initialize_3d_lattice(self):
        """Initialize 3D amorphous carbon lattice with porosity and defects"""
        grain_radius_um = self.simulation_parameters.get("grain_radius_um", 0.1)
        grain_radius_cm = grain_radius_um * 1e-4
        site_area_angstroms_sq = self.simulation_parameters.get("site_area_angstroms_sq", 9)
        site_area_cm2 = site_area_angstroms_sq * 1e-16
        grain_surface_area_cm2 = 4 * np.pi * (grain_radius_cm**2)
        
        # Calculate 3D lattice dimensions
        calculated_sites = int(grain_surface_area_cm2 / site_area_cm2)
        self.surface_dimension = int(np.cbrt(calculated_sites * 2))  # 3D cube root
        self.depth_layers = max(3, self.surface_dimension // 4)  # Multiple layers for 3D structure
        
        # Initialize 3D lattice: (depth, height, width)
        # Site types: None=void, "C"=carbon, "H"=hydrogen, "defect"=surface defect
        self.lattice = np.full((self.depth_layers, self.surface_dimension, self.surface_dimension), None, dtype=object)
        
        # Site type mapping: 0=void, 1=physisorption, 2=chemisorption, 3=defect
        self.site_types = np.zeros((self.depth_layers, self.surface_dimension, self.surface_dimension), dtype=int)
        
        # Energy maps for 3D lattice
        self.E_bind_eV_map = np.zeros((self.depth_layers, self.surface_dimension, self.surface_dimension))
        self.E_diff_eV_map = np.zeros((self.depth_layers, self.surface_dimension, self.surface_dimension))
        
        self._generate_amorphous_structure()
        self._assign_site_types()
        self._generate_energy_maps()
        
        # Initialize H atoms
        initial_h_coverage = self.simulation_parameters.get("initial_h_coverage", 0.0)
        self.h_atoms_on_surface = 0
        self.adjacent_h_pairs_count = 0
        
        if initial_h_coverage > 0:
            self._initialize_h_atoms(initial_h_coverage)

    def _generate_amorphous_structure(self):
        """Generate amorphous carbon structure with porosity and defects"""
        rng = np.random.default_rng()
        
        # Start with solid carbon structure
        self.lattice.fill("C")
        
        # Add porosity (voids) throughout the lattice
        porosity_mask = rng.random(self.lattice.shape) < self.porosity_fraction
        self.lattice[porosity_mask] = None
        
        # Add surface defects (stronger binding sites) on the top layer
        surface_defect_mask = rng.random((self.surface_dimension, self.surface_dimension)) < self.surface_defect_fraction
        self.lattice[0, surface_defect_mask] = "defect"
        
        # Ensure surface layer has accessible sites for adsorption
        surface_sites = self.lattice[0, :, :]
        accessible_mask = (surface_sites == "C") | (surface_sites == "defect")
        if not np.any(accessible_mask):
            # If no accessible sites, create some
            num_accessible = max(1, int(self.surface_dimension * self.surface_dimension * 0.1))
            accessible_indices = np.random.choice(self.surface_dimension * self.surface_dimension, 
                                                num_accessible, replace=False)
            for idx in accessible_indices:
                r, c = idx // self.surface_dimension, idx % self.surface_dimension
                self.lattice[0, r, c] = "C"

    def _assign_site_types(self):
        """Assign site types: 0=void, 1=physisorption, 2=chemisorption, 3=defect"""
        rng = np.random.default_rng()
        
        # Surface layer (layer 0) - accessible for adsorption
        surface_mask = self.lattice[0, :, :] != None
        accessible_sites = np.where(surface_mask)
        
        if len(accessible_sites[0]) > 0:
            # Assign chemisorption sites (stronger binding)
            num_chemisorption = int(len(accessible_sites[0]) * self.chemisorption_fraction)
            if num_chemisorption > 0:
                chem_indices = rng.choice(len(accessible_sites[0]), num_chemisorption, replace=False)
                for idx in chem_indices:
                    r, c = accessible_sites[0][idx], accessible_sites[1][idx]
                    if self.lattice[0, r, c] == "defect":
                        self.site_types[0, r, c] = 3  # defect site
                    else:
                        self.site_types[0, r, c] = 2  # chemisorption site
            
            # Remaining accessible sites are physisorption
            for r, c in zip(accessible_sites[0], accessible_sites[1]):
                if self.site_types[0, r, c] == 0:  # Not yet assigned
                    self.site_types[0, r, c] = 1  # physisorption site

    def _generate_energy_maps(self):
        """Generate binding and diffusion energy maps for 3D lattice"""
        rng = np.random.default_rng()
        
        for d in range(self.depth_layers):
            for r in range(self.surface_dimension):
                for c in range(self.surface_dimension):
                    site_type = self.site_types[d, r, c]
                    
                    if site_type == 0:  # void
                        self.E_bind_eV_map[d, r, c] = 0.0
                        self.E_diff_eV_map[d, r, c] = 0.0
                    elif site_type == 1:  # physisorption (450-550 K range)
                        # Convert K to eV: 500 K ≈ 43 meV
                        E_bind_kelvin = rng.normal(500.0, 50.0)  # 450-550 K range
                        E_bind_kelvin = np.clip(E_bind_kelvin, 450.0, 550.0)
                        E_bind_eV = E_bind_kelvin * 8.617e-5  # Convert K to eV
                        self.E_bind_eV_map[d, r, c] = E_bind_eV
                        # Diffusion barrier is 0.3 times binding energy
                        self.E_diff_eV_map[d, r, c] = 0.3 * E_bind_eV
                    elif site_type == 2:  # chemisorption (1.5-2.0 eV range)
                        E_bind_eV = rng.normal(self.E_chem_mean_eV, 0.25)  # 1.5-2.0 eV range
                        E_bind_eV = np.clip(E_bind_eV, 1.5, 2.0)
                        self.E_bind_eV_map[d, r, c] = E_bind_eV
                        # Diffusion barrier is 0.3 times binding energy
                        self.E_diff_eV_map[d, r, c] = 0.3 * E_bind_eV
                    elif site_type == 3:  # defect (enhanced chemisorption)
                        E_bind_eV = rng.normal(self.E_chem_mean_eV * 1.1, 0.25)  # Slightly stronger
                        E_bind_eV = np.clip(E_bind_eV, 1.6, 2.2)
                        self.E_bind_eV_map[d, r, c] = E_bind_eV
                        # Diffusion barrier is 0.3 times binding energy
                        self.E_diff_eV_map[d, r, c] = 0.3 * E_bind_eV
                    
                    # Ensure diffusion barrier doesn't exceed binding energy
                    self.E_diff_eV_map[d, r, c] = min(self.E_diff_eV_map[d, r, c], 
                                                      self.E_bind_eV_map[d, r, c] * 0.8)

    def _initialize_h_atoms(self, initial_coverage):
        """Initialize H atoms on accessible surface sites"""
        accessible_sites = np.where((self.site_types[0, :, :] > 0) & (self.lattice[0, :, :] != "H"))
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
        """Get 3D neighbors including surface and subsurface connections"""
        neighbors = []
        
        # Surface neighbors (same layer)
        moves_2d = [(-1, 0), (1, 0), (0, -1), (0, 1)]
        for dr, dc in moves_2d:
            nr, nc = (r + dr) % self.surface_dimension, (c + dc) % self.surface_dimension
            if self.site_types[d, nr, nc] > 0:  # Accessible site
                neighbors.append((d, nr, nc))
        
        # Subsurface connections (if not at bottom layer)
        if d < self.depth_layers - 1:
            if self.site_types[d + 1, r, c] > 0:
                neighbors.append((d + 1, r, c))
        
        # Surface connections (if not at top layer)
        if d > 0:
            if self.site_types[d - 1, r, c] > 0:
                neighbors.append((d - 1, r, c))
        
        return neighbors

    def _update_adjacent_h_pairs_count(self):
        """Count adjacent H atom pairs in 3D lattice"""
        self.adjacent_h_pairs_count = 0
        
        for d in range(self.depth_layers):
            for r in range(self.surface_dimension):
                for c in range(self.surface_dimension):
                    if self.lattice[d, r, c] == "H":
                        neighbors = self.get_neighbors_3d(d, r, c)
                        for nd, nr, nc in neighbors:
                            if self.lattice[nd, nr, nc] == "H":
                                # Count each pair only once
                                if (d, r, c) < (nd, nr, nc):
                                    self.adjacent_h_pairs_count += 1

    def update_adjacent_h_pairs_count(self, d, r, c, add_atom):
        """Update adjacent H pairs count when adding/removing an H atom"""
        change = 1 if add_atom else -1
        neighbors = self.get_neighbors_3d(d, r, c)
        for nd, nr, nc in neighbors:
            if self.lattice[nd, nr, nc] == "H":
                self.adjacent_h_pairs_count += change

    def get_accessible_surface_sites(self):
        """Get coordinates of all accessible surface sites (layer 0)"""
        accessible_mask = (self.site_types[0, :, :] > 0) & (self.lattice[0, :, :] != "H")
        return np.where(accessible_mask)

    def get_occupied_sites(self):
        """Get coordinates of all occupied sites (with H atoms) - optimized version"""
        # Use numpy operations for better performance
        occupied_mask = (self.lattice == "H")
        if not np.any(occupied_mask):
            return []
        
        # Get indices where H atoms are present
        occupied_indices = np.where(occupied_mask)
        occupied_sites = list(zip(occupied_indices[0], occupied_indices[1], occupied_indices[2]))
        return occupied_sites

    def calculate_rates(self):
        rates = {}
        surface_temp_k = self.simulation_parameters.get("surface_temperature_k")
        gas_temp_k = self.simulation_parameters.get("gas_temperature_k")
        h_gas_density = self.simulation_parameters.get("h_gas_density_cm3")
        uv_flux_factor = self.simulation_parameters.get("uv_flux_factor", 1.0)
        site_area_cm2 = self.simulation_parameters.get("site_area_angstroms_sq", 9) * 1e-16
        pre_exp_frequency = 1e12

        # Calculate accessible surface area (only surface layer)
        accessible_sites = self.get_accessible_surface_sites()
        num_accessible_sites = len(accessible_sites[0])
        accessible_area_cm2 = num_accessible_sites * site_area_cm2

        # Adsorption rate
        v_h_thermal = np.sqrt(8 * K_B * gas_temp_k / (np.pi * M_H))
        gas_flux_per_cm2_s = 0.25 * h_gas_density * v_h_thermal
        
        # Temperature-dependent sticking probability: S(T) ∝ exp(-T/100, K)
        base_sticking_prob = self.simulation_parameters.get("sticking_probability", 0.3)
        surface_temp_k = self.simulation_parameters.get("surface_temperature_k")
        temperature_factor = np.exp(-surface_temp_k / 100.0)
        sticking_prob = base_sticking_prob * temperature_factor
        
        rates["adsorption"] = gas_flux_per_cm2_s * accessible_area_cm2 * sticking_prob

        # Eley-Rideal formation
        if self.h_atoms_on_surface > 0:
            rates["h2_formation_ER"] = gas_flux_per_cm2_s * surface_chemistry_data["er_cross_section_cm2"] * self.h_atoms_on_surface

        # Desorption and diffusion rates
        if self.h_atoms_on_surface > 0:
            occupied_sites = self.get_occupied_sites()
            if occupied_sites:
                desorption_rates = []
                diffusion_rates = []
                
                for d, r, c in occupied_sites:
                    E_bind = self.E_bind_eV_map[d, r, c]
                    E_diff = self.E_diff_eV_map[d, r, c]
                    
                    desorption_rate = pre_exp_frequency * np.exp(-E_bind * EV_TO_KELVIN / surface_temp_k)
                    diffusion_rate = pre_exp_frequency * np.exp(-E_diff * EV_TO_KELVIN / surface_temp_k)
                    
                    desorption_rates.append(desorption_rate)
                    diffusion_rates.append(diffusion_rate)
                
                rates["desorption"] = float(np.sum(desorption_rates))
                rates["diffusion"] = float(np.sum(diffusion_rates))

        # Langmuir-Hinshelwood formation
        if self.adjacent_h_pairs_count > 0:
            occupied_sites = self.get_occupied_sites()
            if occupied_sites:
                E_diff_values = [self.E_diff_eV_map[d, r, c] for d, r, c in occupied_sites]
                E_diff_mean = float(np.mean(E_diff_values))
            else:
                E_diff_mean = float(np.mean(self.E_diff_eV_map[self.E_diff_eV_map > 0]))
            
            h2_formation_barrier_ev = E_diff_mean
            rates["h2_formation_LH"] = pre_exp_frequency * np.exp(-h2_formation_barrier_ev * EV_TO_KELVIN / surface_temp_k) * self.adjacent_h_pairs_count

        # UV processes
        if uv_flux_factor > 0:
            uv_photon_flux_total = uv_photon_flux["integrated_fuv_photon_flux_photons_cm2_s"] * uv_flux_factor
            
            # UV pulse rate based on simulation workflow: 1-10 photons grain⁻¹ yr⁻¹
            # Convert to per-second rate: 1-10 photons / (3.154e7 s) = 3.17e-8 to 3.17e-7 s⁻¹
            base_uv_rate = 5.0e-8  # 5 photons grain⁻¹ yr⁻¹ (middle of 1-10 range)
            
            # Stochastic UV pulse generation
            if self.uv_pulse_enabled:
                uv_pulse_rate = base_uv_rate * uv_flux_factor  # Scale with UV field strength
                rates["uv_pulse_start"] = uv_pulse_rate
            
            # UV processes during active pulse
            if self.uv_pulse_active:
                if self.h_atoms_on_surface > 0:
                    # Photodesorption rate: proportional to UV pulse rate and H atoms
                    # Use same base rate as UV pulses for consistency
                    photodesorption_rate = base_uv_rate * uv_flux_factor * self.h_atoms_on_surface
                    rates["uv_photodesorption"] = photodesorption_rate
                
                if self.adjacent_h_pairs_count > 0:
                    rates["h2_formation_UV"] = uv_photon_flux_total * (2 * site_area_cm2) * surface_chemistry_data["uv_h2_formation_yield_per_pair"] * self.adjacent_h_pairs_count
                    if self.simulation_parameters.get("enable_uv_photodissociation", True):
                        yield_pair = float(self.simulation_parameters.get("uv_photodissociation_yield_per_pair", 1e-6))
                        rates["uv_photodissociation_pair"] = uv_photon_flux_total * (2 * site_area_cm2) * yield_pair * self.adjacent_h_pairs_count
                
                # UV-induced defect creation
                accessible_surface_sites = self.get_accessible_surface_sites()
                if len(accessible_surface_sites[0]) > 0:
                    # Convert rate from sites/Myr to sites/s
                    defect_creation_rate = self.uv_defect_creation_rate / (3.154e13)  # Convert Myr to s
                    rates["uv_defect_creation"] = defect_creation_rate * uv_flux_factor
                
                uv_diffusion_factor = float(self.simulation_parameters.get("uv_stimulated_diffusion_factor", 1.0))
                if uv_diffusion_factor != 1.0 and self.h_atoms_on_surface > 0:
                    rates["uv_stimulated_diffusion"] = max(0.0, (uv_diffusion_factor - 1.0)) * rates.get("diffusion", 0.0)

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
                    self.lattice[d, r, c] = "C"  # Return to carbon site
                    self.h_atoms_on_surface -= 1
                    self.total_desorbed_h_atoms += 1
                    self.update_adjacent_h_pairs_count(d, r, c, False)
                    
        elif event_type == "diffusion":
            if self.h_atoms_on_surface > 0:
                occupied_sites = self.get_occupied_sites()
                if occupied_sites:
                    # Find mobile atoms (with accessible neighbors)
                    mobile_atoms = []
                    for d, r, c in occupied_sites:
                        neighbors = self.get_neighbors_3d(d, r, c)
                        empty_neighbors = [(nd, nr, nc) for nd, nr, nc in neighbors if self.lattice[nd, nr, nc] != "H"]
                        if empty_neighbors:
                            mobile_atoms.append((d, r, c))
                    
                    if mobile_atoms:
                        d_h, r_h, c_h = random.choice(mobile_atoms)
                        neighbors = self.get_neighbors_3d(d_h, r_h, c_h)
                        empty_neighbors = [(nd, nr, nc) for nd, nr, nc in neighbors if self.lattice[nd, nr, nc] != "H"]
                        
                        if empty_neighbors:
                            d_empty, r_empty, c_empty = random.choice(empty_neighbors)
                            self.lattice[d_h, r_h, c_h] = "C"
                            self.update_adjacent_h_pairs_count(d_h, r_h, c_h, False)
                            self.lattice[d_empty, r_empty, c_empty] = "H"
                            self.update_adjacent_h_pairs_count(d_empty, r_empty, c_empty, True)
                            
        elif event_type == "h2_formation_LH":
            if self.adjacent_h_pairs_count > 0:
                # Find all adjacent H pairs
                pairs = []
                for d in range(self.depth_layers):
                    for r in range(self.surface_dimension):
                        for c in range(self.surface_dimension):
                            if self.lattice[d, r, c] == "H":
                                neighbors = self.get_neighbors_3d(d, r, c)
                                for nd, nr, nc in neighbors:
                                    if self.lattice[nd, nr, nc] == "H":
                                        if (d, r, c) < (nd, nr, nc):
                                            pairs.append(((d, r, c), (nd, nr, nc)))
                
                if pairs:
                    (d1, r1, c1), (d2, r2, c2) = random.choice(pairs)
                    self.lattice[d1, r1, c1] = "C"
                    self.update_adjacent_h_pairs_count(d1, r1, c1, False)
                    self.lattice[d2, r2, c2] = "C"
                    self.update_adjacent_h_pairs_count(d2, r2, c2, False)
                    self.h_atoms_on_surface -= 2
                    self.h2_molecules_formed += 1
                    self.h2_molecules_formed_LH += 1
                    
        elif event_type == "h2_formation_UV":
            if self.adjacent_h_pairs_count > 0:
                pairs = []
                for d in range(self.depth_layers):
                    for r in range(self.surface_dimension):
                        for c in range(self.surface_dimension):
                            if self.lattice[d, r, c] == "H":
                                neighbors = self.get_neighbors_3d(d, r, c)
                                for nd, nr, nc in neighbors:
                                    if self.lattice[nd, nr, nc] == "H":
                                        if (d, r, c) < (nd, nr, nc):
                                            pairs.append(((d, r, c), (nd, nr, nc)))
                
                if pairs:
                    (d1, r1, c1), (d2, r2, c2) = random.choice(pairs)
                    self.lattice[d1, r1, c1] = "C"
                    self.update_adjacent_h_pairs_count(d1, r1, c1, False)
                    self.lattice[d2, r2, c2] = "C"
                    self.update_adjacent_h_pairs_count(d2, r2, c2, False)
                    self.h_atoms_on_surface -= 2
                    self.h2_molecules_formed += 1
                    self.h2_molecules_formed_UV += 1
                    
        elif event_type == "uv_photodissociation_pair":
            if self.adjacent_h_pairs_count > 0:
                pairs = []
                for d in range(self.depth_layers):
                    for r in range(self.surface_dimension):
                        for c in range(self.surface_dimension):
                            if self.lattice[d, r, c] == "H":
                                neighbors = self.get_neighbors_3d(d, r, c)
                                for nd, nr, nc in neighbors:
                                    if self.lattice[nd, nr, nc] == "H":
                                        if (d, r, c) < (nd, nr, nc):
                                            pairs.append(((d, r, c), (nd, nr, nc)))
                
                if pairs:
                    (d1, r1, c1), (d2, r2, c2) = random.choice(pairs)
                    self.lattice[d1, r1, c1] = "C"
                    self.update_adjacent_h_pairs_count(d1, r1, c1, False)
                    self.lattice[d2, r2, c2] = "C"
                    self.update_adjacent_h_pairs_count(d2, r2, c2, False)
                    self.h_atoms_on_surface -= 2
                    
        elif event_type == "uv_stimulated_diffusion":
            if self.h_atoms_on_surface > 0:
                occupied_sites = self.get_occupied_sites()
                if occupied_sites:
                    mobile_atoms = []
                    for d, r, c in occupied_sites:
                        neighbors = self.get_neighbors_3d(d, r, c)
                        empty_neighbors = [(nd, nr, nc) for nd, nr, nc in neighbors if self.lattice[nd, nr, nc] != "H"]
                        if empty_neighbors:
                            mobile_atoms.append((d, r, c))
                    
                    if mobile_atoms:
                        d_h, r_h, c_h = random.choice(mobile_atoms)
                        neighbors = self.get_neighbors_3d(d_h, r_h, c_h)
                        empty_neighbors = [(nd, nr, nc) for nd, nr, nc in neighbors if self.lattice[nd, nr, nc] != "H"]
                        
                        if empty_neighbors:
                            d_empty, r_empty, c_empty = random.choice(empty_neighbors)
                            self.lattice[d_h, r_h, c_h] = "C"
                            self.update_adjacent_h_pairs_count(d_h, r_h, c_h, False)
                            self.lattice[d_empty, r_empty, c_empty] = "H"
                            self.update_adjacent_h_pairs_count(d_empty, r_empty, c_empty, True)
                            
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
            # Start UV pulse
            self.uv_pulse_active = True
            self.last_uv_pulse_time = self.time
            # Schedule pulse end
            pulse_end_time = self.time + self.uv_pulse_duration
            
        elif event_type == "uv_pulse_end":
            # End UV pulse
            self.uv_pulse_active = False
            
        elif event_type == "uv_defect_creation":
            # Create new defect site on surface
            accessible_sites = self.get_accessible_surface_sites()
            if len(accessible_sites[0]) > 0:
                # Choose random accessible site
                idx = random.randint(0, len(accessible_sites[0]) - 1)
                r, c = accessible_sites[0][idx], accessible_sites[1][idx]
                
                # Convert to defect site
                if self.lattice[0, r, c] == "C":
                    self.lattice[0, r, c] = "defect"
                    self.site_types[0, r, c] = 3  # defect site type
                    
                    # Update energy map for the new defect site
                    rng = np.random.default_rng()
                    E_bind_eV = rng.normal(self.E_chem_mean_eV * 1.1, 0.25)
                    E_bind_eV = np.clip(E_bind_eV, 1.6, 2.2)
                    self.E_bind_eV_map[0, r, c] = E_bind_eV
                    self.E_diff_eV_map[0, r, c] = 0.3 * E_bind_eV

    def run_gillespie(self, max_time, max_steps=None):
        """
        Run Kinetic Monte Carlo simulation using n-fold way algorithm.
        
        Algorithm:
        - Precompute all event rates (adsorption, desorption, diffusion, reaction)
        - Calculate total rate k_total = sum of all rates
        - Time advancement: Δt = -ln(r)/k_total where r is random number in (0,1)
        - Event selection: Choose event i with probability k_i/k_total
        - Time resolution: 10^-6 seconds as specified in methodology
        """
        step_count = 0
        while self.time < max_time:
            if max_steps and step_count >= max_steps: break
            
            # Check if UV pulse should end
            if self.uv_pulse_active and (self.time - self.last_uv_pulse_time) >= self.uv_pulse_duration:
                self.uv_pulse_active = False
            
            # Precompute all event rates (n-fold way requirement)
            rates = self.calculate_rates()
            if not rates: break
            total_rate = sum(rates.values())
            if total_rate == 0: break
            
            # Time advancement: Δt = -ln(r)/k_total (n-fold way formula)
            # random.expovariate(total_rate) is equivalent to -ln(r)/k_total
            delta_t = random.expovariate(total_rate)
            
            # Time resolution: 10^-6 seconds as specified in methodology
            if delta_t < 1e-6:
                delta_t = 1e-6
            
            if self.time + delta_t > max_time:
                delta_t = max_time - self.time
            
            self.time += delta_t
            
            # Event selection: Choose event with probability proportional to its rate
            # This implements the n-fold way event selection
            chosen_event = random.choices(list(rates.keys()), weights=list(rates.values()), k=1)[0]
            self.execute_event(chosen_event)
            step_count += 1
        return []

