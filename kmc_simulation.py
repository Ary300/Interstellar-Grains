# import numpy as np
# import random
# from scientific_data import surface_chemistry_data, uv_photon_flux, EV_TO_KELVIN, K_B, M_H
# from physical_rates import (
#     h_diffusion_rate, h_desorption_rate, h2_formation_lh_rate, h2_formation_er_rate,
#     uv_photodesorption_rate, uv_h2_formation_rate, adsorption_rate
# )
# import pickle
# import hashlib
# import os
# from pathlib import Path

# class KineticMonteCarlo:
#     def __init__(self, simulation_parameters):
#         self.simulation_parameters = simulation_parameters
#         self.time = 0.0
#         self.h2_molecules_formed = 0
#         self.h2_molecules_formed_LH = 0
#         self.h2_molecules_formed_ER = 0
#         self.h2_molecules_formed_UV = 0
#         self.total_adsorbed_h_atoms = 0
#         self.total_desorbed_h_atoms = 0
        
#         self.use_3d_lattice = self.simulation_parameters.get("use_3d_lattice", True)
#         self.porosity_fraction = self.simulation_parameters.get("porosity_fraction", 0.2)
        
#         self.uv_pulse_enabled = self.simulation_parameters.get("uv_pulse_enabled", True)
#         self.uv_defect_creation_rate = self.simulation_parameters.get("uv_defect_creation_rate", 0.5)
#         self.uv_pulse_duration = self.simulation_parameters.get("uv_pulse_duration", 1e-6)
#         self.last_uv_pulse_time = 0.0
#         self.uv_pulse_active = False
        
#         self.E_bind_mean_meV = float(self.simulation_parameters.get("E_phys_mean_meV", 45.0))
#         self.E_bind_sigma_meV = float(self.simulation_parameters.get("heterogeneity_E_bind_sigma_meV", 5.0))
        
#         # CRITICAL: Initialize with caching
#         self.grain_cache = GrainCache()
#         self.initialize_off_lattice_structure_cached()

#     def initialize_off_lattice_structure_cached(self):
#         """
#         Initialize grain structure with caching support.
#         Checks cache first, generates only if needed.
#         """
#         grain_radius_um = self.simulation_parameters.get("grain_radius_um", 0.1)
#         site_area_angstroms_sq = self.simulation_parameters.get("site_area_angstroms_sq", 9)
        
#         # Generate cache key
#         cache_key = self.grain_cache._get_cache_key(
#             grain_radius_um, site_area_angstroms_sq,
#             self.E_bind_mean_meV, self.E_bind_sigma_meV
#         )
        
#         # Try to load from cache
#         cached_data = self.grain_cache.load_grain(cache_key)
        
#         if cached_data is not None:
#             print(f"Loading grain from cache (key: {cache_key[:8]}...)")
#             self._load_grain_from_cache(cached_data)
#         else:
#             print(f"Generating new grain (will be cached with key: {cache_key[:8]}...)")
#             self._generate_new_grain()
            
#             # Save to cache
#             grain_data = self._extract_grain_data()
#             self.grain_cache.save_grain(grain_data, cache_key)
        
#         # Initialize H atoms (always fresh, not cached)
#         initial_h_coverage = self.simulation_parameters.get("initial_h_coverage", 0.0)
#         self.h_atoms_on_surface = 0
#         self.adjacent_h_pairs_count = 0
        
#         if initial_h_coverage > 0:
#             self._initialize_h_atoms(initial_h_coverage)
    
#     def _generate_new_grain(self):
#         """Generate grain from scratch (original initialization)"""
#         grain_radius_um = self.simulation_parameters.get("grain_radius_um", 0.1)
#         grain_radius_cm = grain_radius_um * 1e-4
#         site_area_angstroms_sq = self.simulation_parameters.get("site_area_angstroms_sq", 9)
#         site_area_cm2 = site_area_angstroms_sq * 1e-16
#         grain_surface_area_cm2 = 4 * np.pi * (grain_radius_cm**2)
#         total_sites_needed = int(grain_surface_area_cm2 / site_area_cm2)
        
#         # CRITICAL: Initialize ALL attributes FIRST
#         self.molecular_coordinates = []
#         self.site_types = {}
#         self.site_occupancy = {}
#         self.site_energetics = {}
#         self.neighbor_cache = {}  # ADD THIS LINE
        
#         # Calculate grain radius in Ångstroms (needed by _generate_off_lattice_sites)
#         self.grain_radius_angstroms = grain_radius_um * 1e4  # ADD THIS LINE
        
#         # Now generate
#         self._generate_off_lattice_sites(total_sites_needed)
#         self._generate_realistic_surface_topology()
#         self._generate_energy_maps()
#         self._build_neighbor_cache()
    
#     def _extract_grain_data(self):
#         """Extract cacheable grain data (structure + energetics)"""
#         return {
#             'molecular_coordinates': self.molecular_coordinates.copy(),
#             'site_types': self.site_types.copy(),
#             'site_energetics': self.site_energetics.copy(),
#             'neighbor_cache': self.neighbor_cache.copy(),
#             'grain_radius_angstroms': self.grain_radius_angstroms
#         }
    
#     def _load_grain_from_cache(self, cached_data):
#         """Load grain structure from cached data"""
#         self.molecular_coordinates = cached_data['molecular_coordinates']
#         self.site_types = cached_data['site_types']
#         self.site_energetics = cached_data['site_energetics']
#         self.neighbor_cache = cached_data['neighbor_cache']
#         self.grain_radius_angstroms = cached_data['grain_radius_angstroms']
        
#         # Initialize fresh occupancy dictionary
#         self.site_occupancy = {i: None for i in range(len(self.molecular_coordinates))}
        
#         print(f"Loaded grain: {len(self.molecular_coordinates)} sites from cache")

#     def initialize_off_lattice_structure(self):
#         """
#         Initialize true off-lattice amorphous carbon structure using continuous molecular coordinates,
#         realistic surface topology, and dynamic nearest-neighbor algorithms.
#         Based on:
#         - Draine & Li 2007, ApJ 657, 810 (amorphous carbon properties)
#         - Mennella et al. 2003, ApJ 587, 727 (surface topology)
#         - Jones & Williams 1987, MNRAS 224, 473 (interstellar carbon grains)
#         """
#         grain_radius_um = self.simulation_parameters.get("grain_radius_um", 0.1)
#         grain_radius_cm = grain_radius_um * 1e-4
#         site_area_angstroms_sq = self.simulation_parameters.get("site_area_angstroms_sq", 9)
#         site_area_cm2 = site_area_angstroms_sq * 1e-16
#         grain_surface_area_cm2 = 4 * np.pi * (grain_radius_cm**2)
        
#         # Calculate total number of sites needed
#         total_sites_needed = int(grain_surface_area_cm2 / site_area_cm2)
        
#         # Create off-lattice structure with continuous molecular coordinates
#         self._generate_off_lattice_sites(total_sites_needed)
#         self._generate_realistic_surface_topology()
#         self._generate_energy_maps()
        
#         # CRITICAL: Build neighbor cache ONCE during initialization
#         self._build_neighbor_cache()
        
#         initial_h_coverage = self.simulation_parameters.get("initial_h_coverage", 0.0)
#         self.h_atoms_on_surface = 0
#         self.adjacent_h_pairs_count = 0
        
#         if initial_h_coverage > 0:
#             self._initialize_h_atoms(initial_h_coverage)

#     def _generate_off_lattice_sites(self, total_sites_needed):
#         """
#         Generate true off-lattice sites using continuous molecular coordinates.
#         Uses realistic amorphous carbon structure with irregular bond networks.
#         """
#         rng = np.random.default_rng()
        
#         # Calculate grain radius in Å
#         grain_radius_um = self.simulation_parameters.get("grain_radius_um", 0.1)
#         self.grain_radius_angstroms = grain_radius_um * 1e4
        
#         # Initialize off-lattice storage
#         self.molecular_coordinates = []  # List of (x, y, z) coordinates in Å
#         self.site_types = {}  # Dictionary mapping site index to type
#         self.site_occupancy = {}  # Dictionary mapping site index to occupant
#         self.site_energetics = {}  # Dictionary mapping site index to (E_bind, E_diff)
        
#         # Generate realistic amorphous carbon structure
#         # Based on Draine & Li 2007: amorphous carbon has irregular bond networks
#         sites_placed = 0
#         attempts = 0
#         max_attempts = total_sites_needed * 1000
        
#         while sites_placed < total_sites_needed and attempts < max_attempts:
#             # Generate random position in 3D space within grain volume
#             # Use spherical coordinates for proper grain geometry
#             r = self.grain_radius_angstroms * (rng.uniform(0.8, 1.0))**(1/3)  # Concentrate near surface
#             theta = rng.uniform(0, np.pi)
#             phi = rng.uniform(0, 2 * np.pi)
            
#             # Convert to Cartesian coordinates
#             x = r * np.sin(theta) * np.cos(phi)
#             y = r * np.sin(theta) * np.sin(phi)
#             z = r * np.cos(theta)
            
#             # Check minimum distance from existing sites
#             too_close = False
#             for existing_coords in self.molecular_coordinates:
#                 distance = np.sqrt((x - existing_coords[0])**2 + 
#                                  (y - existing_coords[1])**2 + 
#                                  (z - existing_coords[2])**2)
#                 if distance < 2.0:  # Minimum 2 Å separation
#                     too_close = True
#                     break
            
#             if not too_close:
#                 site_index = len(self.molecular_coordinates)
#                 self.molecular_coordinates.append((x, y, z))
#                 self.site_occupancy[site_index] = None
#                 sites_placed += 1
            
#             attempts += 1
        
#         # Ensure minimum surface coverage
#         surface_sites = [i for i, coords in enumerate(self.molecular_coordinates) 
#                         if coords[0]**2 + coords[1]**2 + coords[2]**2 > (0.95 * self.grain_radius_angstroms)**2]
        
#         if len(surface_sites) < 10:
#             # Add more surface sites if needed
#             additional_sites = 10 - len(surface_sites)
#             for _ in range(additional_sites):
#                 r = self.grain_radius_angstroms * 0.98  # Near surface
#                 theta = rng.uniform(0, np.pi)
#                 phi = rng.uniform(0, 2 * np.pi)
                
#                 x = r * np.sin(theta) * np.cos(phi)
#                 y = r * np.sin(theta) * np.sin(phi)
#                 z = r * np.cos(theta)
                
#                 site_index = len(self.molecular_coordinates)
#                 self.molecular_coordinates.append((x, y, z))
#                 self.site_occupancy[site_index] = None

#     def _generate_realistic_surface_topology(self):
#         """
#         Generate realistic amorphous carbon surface topology with irregular bond networks,
#         surface roughness, and proper nearest-neighbor algorithms.
#         """
#         rng = np.random.default_rng()
        
#         # Calculate surface roughness based on amorphous carbon properties
#         # Draine & Li 2007: amorphous carbon has significant surface roughness
#         surface_roughness = 0.1 * self.grain_radius_angstroms  # 10% of grain radius
        
#         # Generate surface topology with irregular features
#         for i, coords in enumerate(self.molecular_coordinates):
#             x, y, z = coords
            
#             # Calculate distance from center
#             distance_from_center = np.sqrt(x**2 + y**2 + z**2)
            
#             # Determine if site is on surface (within roughness layer)
#             is_surface_site = distance_from_center > (self.grain_radius_angstroms - surface_roughness)
            
#             if is_surface_site:
#                 # Surface sites: determine type based on local environment
#                 neighbors = self._find_nearest_neighbors(i, cutoff_distance=3.0)
#                 coordination_number = len(neighbors)
                
#                 # Physics-based site type assignment
#                 if coordination_number <= 2:
#                     # Low coordination: dangling bonds → defect sites
#                     self.site_types[i] = 3
#                 elif coordination_number <= 4:
#                     # Medium coordination: unsaturated bonds → chemisorption
#                     self.site_types[i] = 2
#                 else:
#                     # High coordination: saturated bonds → physisorption
#                     self.site_types[i] = 1
#             else:
#                 # Subsurface sites: bulk-like properties
#                 self.site_types[i] = 1

#     def _find_nearest_neighbors(self, site_index, cutoff_distance):
#         """
#         Dynamic nearest-neighbor algorithm using continuous molecular coordinates.
#         Finds all sites within cutoff_distance of the target site.
#         NOTE: This is O(N) and should only be called during initialization!
#         """
#         target_coords = self.molecular_coordinates[site_index]
#         neighbors = []
        
#         for i, coords in enumerate(self.molecular_coordinates):
#             if i == site_index:
#                 continue
            
#             distance = np.sqrt((target_coords[0] - coords[0])**2 + 
#                              (target_coords[1] - coords[1])**2 + 
#                              (target_coords[2] - coords[2])**2)
            
#             if distance <= cutoff_distance:
#                 neighbors.append(i)
        
#         return neighbors

#     def _build_neighbor_cache(self):
#         """
#         Build neighbor cache ONCE during initialization.
#         This is the KEY OPTIMIZATION that eliminates O(N²) lookups during simulation.
#         """
#         print("Building neighbor cache...")
#         self.neighbor_cache = {}
#         for i in range(len(self.molecular_coordinates)):
#             self.neighbor_cache[i] = self._find_nearest_neighbors(i, cutoff_distance=3.0)
#         print(f"Neighbor cache built: {len(self.neighbor_cache)} sites")

#     def _generate_energy_maps(self):
#         """
#         Generate energy maps for off-lattice structure with site-specific energetics.
#         """
#         rng = np.random.default_rng()
        
#         mean_eV = self.E_bind_mean_meV * 1e-3
#         sigma_eV = self.E_bind_sigma_meV * 1e-3
        
#         for i, coords in enumerate(self.molecular_coordinates):
#             site_type = self.site_types.get(i, 1)
#             neighbors = self._find_nearest_neighbors(i, cutoff_distance=3.0)
#             num_neighbors = len(neighbors)
            
#             # Calculate local environment effects
#             local_density = self._calculate_local_density(i)
#             surface_distance = self._calculate_surface_distance(coords)
            
#             # Assign binding energies based on site type and local environment
#             if site_type == 1:  # Physisorption
#                 base_energy = mean_eV
#                 # Local density effect: higher density → stronger binding
#                 density_factor = 1.0 + 0.2 * local_density
#                 # Surface distance effect: closer to surface → stronger binding
#                 surface_factor = 1.0 + 0.1 * (1.0 - surface_distance / self.grain_radius_angstroms)
                
#                 binding_energy = base_energy * density_factor * surface_factor + rng.normal(0, sigma_eV)
#                 diffusion_barrier = 0.025 + rng.normal(0, 0.005)
                
#             elif site_type == 2:  # Chemisorption
#                 binding_energy = 0.25 + rng.normal(0, 0.05)
#                 diffusion_barrier = 0.03 + rng.normal(0, 0.005)
                
#             elif site_type == 3:  # Defect
#                 binding_energy = 0.35 + rng.normal(0, 0.05)
#                 diffusion_barrier = 0.015 + rng.normal(0, 0.003)
            
#             # Ensure positive values
#             binding_energy = max(binding_energy, 0.01)
#             diffusion_barrier = max(diffusion_barrier, 0.01)
            
#             self.site_energetics[i] = (binding_energy, diffusion_barrier)

#     def _calculate_local_density(self, site_index):
#         """
#         Calculate local density around a site using continuous coordinates.
#         """
#         target_coords = self.molecular_coordinates[site_index]
#         local_volume = (4.0/3) * np.pi * (3.0)**3  # 3 Å radius sphere
#         neighbors_in_volume = 0
        
#         for coords in self.molecular_coordinates:
#             distance = np.sqrt((target_coords[0] - coords[0])**2 + 
#                              (target_coords[1] - coords[1])**2 + 
#                              (target_coords[2] - coords[2])**2)
#             if distance <= 3.0:
#                 neighbors_in_volume += 1
        
#         # Normalize to density (sites per Å³)
#         density = (neighbors_in_volume - 1) / local_volume  # -1 to exclude self
#         return density

#     def _calculate_surface_distance(self, coords):
#         """
#         Calculate distance from surface for a given coordinate.
#         """
#         distance_from_center = np.sqrt(coords[0]**2 + coords[1]**2 + coords[2]**2)
#         surface_distance = abs(distance_from_center - self.grain_radius_angstroms)
#         return surface_distance

#     def _initialize_h_atoms(self, initial_coverage):
#         """
#         Initialize H atoms on off-lattice surface sites.
#         """
#         surface_sites = [i for i, coords in enumerate(self.molecular_coordinates) 
#                         if self._is_surface_site(coords)]
        
#         if surface_sites:
#             num_initial_h = int(len(surface_sites) * initial_coverage)
#             if num_initial_h > 0:
#                 selected_sites = random.sample(surface_sites, min(num_initial_h, len(surface_sites)))
#                 for site_index in selected_sites:
#                     self.site_occupancy[site_index] = "H"
#                     self.h_atoms_on_surface += 1
#                 self._update_adjacent_h_pairs_count()

#     def _is_surface_site(self, coords):
#         """
#         Determine if a site is on the surface using continuous coordinates.
#         """
#         distance_from_center = np.sqrt(coords[0]**2 + coords[1]**2 + coords[2]**2)
#         surface_thickness = 0.05 * self.grain_radius_angstroms  # 5% of radius
#         return distance_from_center > (self.grain_radius_angstroms - surface_thickness)

#     def get_neighbors_3d(self, site_index):
#         """
#         Get neighbors for off-lattice structure using CACHED neighbor list.
#         This is now O(1) instead of O(N)!
#         """
#         return self.neighbor_cache.get(site_index, [])

#     def get_accessible_surface_sites(self):
#         """
#         Get accessible surface sites (not occupied by H atoms).
#         """
#         accessible_sites = []
#         for i, coords in enumerate(self.molecular_coordinates):
#             if self._is_surface_site(coords) and self.site_occupancy.get(i) != "H":
#                 accessible_sites.append(i)
#         return accessible_sites

#     def get_occupied_sites(self):
#         """
#         Get all sites occupied by H atoms in off-lattice structure.
#         """
#         occupied_sites = []
#         for i, occupancy in self.site_occupancy.items():
#             if occupancy == "H":
#                 occupied_sites.append(i)
#         return occupied_sites

#     def get_site_occupant(self, site_index):
#         """
#         Get the occupant of a site in the off-lattice structure.
#         """
#         return self.site_occupancy.get(site_index, None)

#     def occupy_site(self, site_index, species):
#         """
#         Occupy a site in the off-lattice structure.
#         """
#         if site_index < len(self.molecular_coordinates):
#             self.site_occupancy[site_index] = species

#     def vacate_site(self, site_index):
#         """
#         Vacate a site in the off-lattice structure.
#         """
#         if site_index < len(self.molecular_coordinates):
#             self.site_occupancy[site_index] = None

#     def _update_adjacent_h_pairs_count(self):
#         """
#         Initial full calculation of adjacent H pairs.
#         Only called during initialization!
#         """
#         self.adjacent_h_pairs_count = 0
#         for site_index in self.get_occupied_sites():
#             neighbors = self.get_neighbors_3d(site_index)
#             for neighbor_index in neighbors:
#                 if self.get_site_occupant(neighbor_index) == "H" and site_index < neighbor_index:
#                     self.adjacent_h_pairs_count += 1

#     def update_adjacent_h_pairs_count(self, site_index):
#         """
#         OPTIMIZED: Incremental update - only check neighbors of changed site.
#         This is now O(k) where k is number of neighbors, instead of O(N²)!
#         """
#         neighbors = self.get_neighbors_3d(site_index)
        
#         # Check if site is now occupied or was just vacated
#         is_occupied = (self.get_site_occupant(site_index) == "H")
        
#         # Update pair count based on neighbors
#         for neighbor_index in neighbors:
#             if self.get_site_occupant(neighbor_index) == "H":
#                 # Ensure we don't double count (only count each pair once)
#                 if site_index < neighbor_index:
#                     if is_occupied:
#                         self.adjacent_h_pairs_count += 1
#                     else:
#                         self.adjacent_h_pairs_count -= 1

#     def _calculate_realistic_lateral_interactions(self, site_index, surface_temp_k):
#         """
#         DISABLED FOR PERFORMANCE: Lateral interactions are O(N²) and have minimal impact.
#         Returns 0.0 to maintain interface compatibility.
#         """
#         # PERFORMANCE: This function was causing massive slowdowns (O(N²) per rate calculation)
#         # Lateral interactions have minimal physical impact on rates
#         # Can be re-enabled for detailed studies if needed
#         return 0.0

#     def _calculate_lh_coverage_enhancement(self, current_coverage, surface_temp_k):
#         """
#         Calculate coverage-dependent enhancement factor for LH H2 formation.
#         At low coverage: reduced formation (few neighbors)
#         At medium coverage: enhanced formation (optimal neighbor density)
#         At high coverage: reduced formation (site blocking)
#         """
#         if current_coverage <= 0.0:
#             return 0.0
        
#         # Optimal coverage around 0.1-0.3 ML (monolayer)
#         if current_coverage < 0.1:
#             # Low coverage: linearly increase with coverage
#             enhancement = current_coverage / 0.1
#         elif current_coverage < 0.3:
#             # Optimal range: full enhancement
#             enhancement = 1.0
#         else:
#             # High coverage: site blocking reduces formation
#             enhancement = 1.0 - 0.5 * (current_coverage - 0.3) / 0.7
        
#         # Ensure non-negative
#         enhancement = max(0.0, enhancement)
        
#         # Temperature-dependent boost (higher T = more mobility = better LH)
#         temp_factor = 1.0 + 0.1 * (surface_temp_k - 10.0) / 40.0
    
#         return enhancement * temp_factor
    
#     def _calculate_diffusion_rate_with_tunneling(self, site_index, surface_temp_k):
#         hbar = 1.055e-27
#         mH = 1.674e-24
#         kB = 1.381e-16
#         eV_to_erg = 1.602e-12
        
#         E_diff_eV = self.site_energetics[site_index][1]
#         if E_diff_eV <= 0:
#             return 0.0
#         E_diff_erg = E_diff_eV * eV_to_erg
        
#         a = 3.0e-8
#         nu0 = 1e12
        
#         classical_rate = nu0 * np.exp(-E_diff_erg / (kB * surface_temp_k))
        
#         tunneling_prob = np.exp(-2 * np.sqrt(2 * mH * E_diff_erg) * a / hbar)
#         tunneling_rate = nu0 * tunneling_prob
        
#         return classical_rate + tunneling_rate

#     def calculate_rates(self):
#         """
#         Calculate all rates using physically correct Transition State Theory
#         with quantum tunneling corrections.
#         OPTIMIZED: Removed O(N²) lateral interaction calculations.
#         """
#         rates = {}
#         surface_temp_k = self.simulation_parameters.get("surface_temperature_k")
#         gas_temp_k = self.simulation_parameters.get("gas_temperature_k")
#         h_gas_density = self.simulation_parameters.get("h_gas_density_cm3")
#         uv_flux_factor = self.simulation_parameters.get("uv_flux_factor", 1.0)
#         site_area_cm2 = self.simulation_parameters.get("site_area_angstroms_sq", 9) * 1e-16
#         sticking_probability = self.simulation_parameters.get("sticking_probability", 0.3)

#         accessible_sites = self.get_accessible_surface_sites()
#         num_accessible_sites = len(accessible_sites)
#         accessible_area_cm2 = num_accessible_sites * site_area_cm2
        
#         # Coverage effects: reduce available sites
#         current_coverage = self.h_atoms_on_surface / max(1, num_accessible_sites)
#         available_sites_fraction = max(0.0, 1.0 - current_coverage)
#         effective_area_cm2 = accessible_area_cm2 * available_sites_fraction

#         # 1. Adsorption rate - physically derived from gas kinetics with site blocking
#         rates["adsorption"] = adsorption_rate(
#             h_gas_density, gas_temp_k, sticking_probability, effective_area_cm2
#         )

#         # 2. Eley-Rideal formation - gas impingement × cross-section × reaction probability
#         if self.h_atoms_on_surface > 0:
#             # Calculate gas flux for ER
#             v_thermal = np.sqrt(8 * K_B * gas_temp_k * 1.602e-19 / (np.pi * M_H))
#             gas_flux_cm2_s = 0.25 * h_gas_density * v_thermal
            
#             rates["h2_formation_ER"] = h2_formation_er_rate(
#                 gas_flux_cm2_s, self.h_atoms_on_surface
#             )

#         # 3. Desorption and diffusion rates - TST + quantum tunneling
#         if self.h_atoms_on_surface > 0:
#             occupied_sites = self.get_occupied_sites()
#             if occupied_sites:
#                 total_desorption_rate = 0.0
#                 total_diffusion_rate = 0.0
                
#                 for site_index in occupied_sites:
#                     site_type = self.site_types.get(site_index, 1)
#                     base_binding_energy = self.site_energetics[site_index][0]
                    
#                     # PERFORMANCE: Lateral interactions disabled (was O(N²))
#                     lateral_energy = 0.0
#                     binding_energy = base_binding_energy + lateral_energy
                    
#                     # Desorption rate using TST with coverage effects
#                     desorption_rate = h_desorption_rate(binding_energy, surface_temp_k)
#                     total_desorption_rate += desorption_rate
                    
#                     # Diffusion rate using TST + tunneling
#                     diffusion_rate = h_diffusion_rate(site_type, surface_temp_k)
#                     total_diffusion_rate += diffusion_rate
                
#                 rates["desorption"] = total_desorption_rate
#                 rates["diffusion"] = total_diffusion_rate

#         # 4. Langmuir-Hinshelwood formation - TST for surface reactions with realistic coverage effects
#         if self.adjacent_h_pairs_count > 0:
#             # Physics-based coverage effects: electronic structure changes and site availability
#             coverage_enhancement = self._calculate_lh_coverage_enhancement(current_coverage, surface_temp_k)
#             enhanced_lh_rate = h2_formation_lh_rate(
#                 surface_temp_k, self.adjacent_h_pairs_count
#             ) * coverage_enhancement
#             rates["h2_formation_LH"] = enhanced_lh_rate

#         # 5. UV processes - photon flux × cross-section × yield
#         if uv_flux_factor > 0:
#             uv_photon_flux_total = uv_photon_flux["integrated_fuv_photon_flux_photons_cm2_s"] * uv_flux_factor
            
#             # UV pulse rate from Andersson et al. 2006, A&A 453, 459 (interstellar UV field)
#             base_uv_rate = 5.0e-8  # 5 photons grain⁻¹ yr⁻¹ in diffuse ISM
#             if self.uv_pulse_enabled:
#                 rates["uv_pulse_start"] = base_uv_rate * uv_flux_factor
            
#             if self.uv_pulse_active:
#                 pulse_end_rate = 1.0 / self.uv_pulse_duration
#                 rates["uv_pulse_end"] = pulse_end_rate
            
#             if self.uv_pulse_active:
#                 if self.h_atoms_on_surface > 0:
#                     # Photodesorption using physical cross-sections and yields
#                     rates["uv_photodesorption"] = uv_photodesorption_rate(
#                         uv_photon_flux_total, self.h_atoms_on_surface
#                     )
                
#                 if self.adjacent_h_pairs_count > 0:
#                     # UV-assisted H2 formation
#                     rates["h2_formation_UV"] = uv_h2_formation_rate(
#                         uv_photon_flux_total, self.adjacent_h_pairs_count, site_area_cm2
#                     )
                
#                 # UV-induced defect creation (simplified model)
#                 accessible_surface_sites = self.get_accessible_surface_sites()
#                 if len(accessible_surface_sites) > 0:
#                     defect_creation_rate = self.uv_defect_creation_rate / (3.154e13)  # Convert Myr to s
#                     rates["uv_defect_creation"] = defect_creation_rate * uv_flux_factor
                
#                 # UV-stimulated diffusion enhancement
#                 uv_diffusion_factor = float(self.simulation_parameters.get("uv_stimulated_diffusion_factor", 1.0))
#                 if uv_diffusion_factor > 1.0 and "diffusion" in rates:
#                     rates["uv_stimulated_diffusion"] = (uv_diffusion_factor - 1.0) * rates["diffusion"]
        
#         return {k: v for k, v in rates.items() if v > 0}

#     def execute_event(self, event_type):
#         if event_type == "adsorption":
#             accessible_sites = self.get_accessible_surface_sites()
#             if len(accessible_sites) > 0:
#                 site_index = random.choice(accessible_sites)
#                 self.occupy_site(site_index, "H")
#                 self.h_atoms_on_surface += 1
#                 self.total_adsorbed_h_atoms += 1
#                 self.update_adjacent_h_pairs_count(site_index)
                
#         elif event_type in ["desorption", "uv_photodesorption"]:
#             if self.h_atoms_on_surface > 0:
#                 occupied_sites = self.get_occupied_sites()
#                 if occupied_sites:
#                     site_index = random.choice(occupied_sites)
#                     self.vacate_site(site_index)
#                     self.h_atoms_on_surface -= 1
#                     self.total_desorbed_h_atoms += 1
#                     self.update_adjacent_h_pairs_count(site_index)
                    
#         elif event_type in ["diffusion", "uv_stimulated_diffusion"]:
#             if self.h_atoms_on_surface > 0:
#                 occupied_sites = self.get_occupied_sites()
#                 if occupied_sites:
#                     mobile_atoms = [
#                         site_index for site_index in occupied_sites 
#                         if any(self.get_site_occupant(nd) != "H" for nd in self.get_neighbors_3d(site_index))
#                     ]
#                     if mobile_atoms:
#                         source_site_index = random.choice(mobile_atoms)
#                         empty_neighbors = [
#                             neighbor_index for neighbor_index in self.get_neighbors_3d(source_site_index) 
#                             if self.get_site_occupant(neighbor_index) != "H"
#                         ]
#                         if empty_neighbors:
#                             target_site_index = random.choice(empty_neighbors)
#                             self.vacate_site(source_site_index)
#                             self.update_adjacent_h_pairs_count(source_site_index)
#                             self.occupy_site(target_site_index, "H")
#                             self.update_adjacent_h_pairs_count(target_site_index)
                            
#         elif event_type in ["h2_formation_LH", "h2_formation_UV"]:
#             if self.adjacent_h_pairs_count > 0:
#                 pairs = []
#                 pair_rates = []
#                 for site_index in self.get_occupied_sites():
#                     for neighbor_index in self.get_neighbors_3d(site_index):
#                         if self.get_site_occupant(neighbor_index) == "H" and site_index < neighbor_index:
#                             avg_binding = (self.site_energetics[site_index][0] + self.site_energetics[neighbor_index][0]) / 2
#                             site_rate = avg_binding * np.exp(-0.02 / (K_B * self.simulation_parameters.get("surface_temperature_k", 10)))
#                             pairs.append((site_index, neighbor_index))
#                             pair_rates.append(site_rate)
                
#                 if pairs:
#                     if len(pair_rates) > 1:
#                         total_rate = sum(pair_rates)
#                         probs = [r/total_rate for r in pair_rates]
#                         chosen_idx = random.choices(range(len(pairs)), weights=probs, k=1)[0]
#                         site1, site2 = pairs[chosen_idx]
#                     else:
#                         site1, site2 = pairs[0]
                        
#                     self.vacate_site(site1)
#                     self.update_adjacent_h_pairs_count(site1)
#                     self.vacate_site(site2)
#                     self.update_adjacent_h_pairs_count(site2)
#                     self.h_atoms_on_surface -= 2
#                     self.h2_molecules_formed += 1
#                     if event_type == "h2_formation_LH":
#                         self.h2_molecules_formed_LH += 1
#                     else:
#                         self.h2_molecules_formed_UV += 1
                                    
#         elif event_type == "h2_formation_ER":
#             if self.h_atoms_on_surface > 0:
#                 occupied_sites = self.get_occupied_sites()
#                 if occupied_sites:
#                     site_index = random.choice(occupied_sites)
#                     self.vacate_site(site_index)
#                     self.h_atoms_on_surface -= 1
#                     self.update_adjacent_h_pairs_count(site_index)
#                     self.h2_molecules_formed_ER += 1
#                     self.h2_molecules_formed += 1
                    
#         elif event_type == "uv_pulse_start":
#             self.uv_pulse_active = True
#             self.last_uv_pulse_time = self.time
            
#         elif event_type == "uv_pulse_end":
#             self.uv_pulse_active = False

#     def run_gillespie(self, max_time, max_steps=None):
#         step_count = 0
#         steady_state_check_interval = 1000
#         h2_formation_history = []
        
#         print(f"Starting Gillespie: max_time={max_time}, num_sites={len(self.molecular_coordinates)}")
#         print(f"Initial H atoms: {self.h_atoms_on_surface}, Initial H2: {self.h2_molecules_formed}")
        
#         while self.time < max_time:
#             if max_steps and step_count >= max_steps:
#                 print(f"Reached max_steps={max_steps}, breaking")
#                 break
            
#             # Diagnostic every 1000 steps
#             if step_count % 1000 == 0:
#                 print(f"Step {step_count}, time={self.time:.6f}, H_atoms={self.h_atoms_on_surface}, H2={self.h2_molecules_formed}")
            
#             rates = self.calculate_rates()
#             if not rates:
#                 print(f"No rates available at step {step_count}, breaking")
#                 break
#             total_rate = sum(rates.values())
#             if total_rate == 0:
#                 print(f"Total rate is zero at step {step_count}, breaking")
#                 break
            
#             delta_t = random.expovariate(total_rate)
#             self.time += delta_t
            
#             if self.time > max_time:
#                 self.time = max_time
#                 print(f"Reached max_time at step {step_count}")
#                 break
            
#             chosen_event = random.choices(list(rates.keys()), weights=list(rates.values()), k=1)[0]
#             self.execute_event(chosen_event)
#             step_count += 1
            
#             if step_count % steady_state_check_interval == 0:
#                 h2_formation_history.append(self.h2_molecules_formed)
#                 if len(h2_formation_history) >= 5:
#                     recent_rates = []
#                     for i in range(1, len(h2_formation_history)):
#                         rate = (h2_formation_history[i] - h2_formation_history[i-1]) / (steady_state_check_interval * delta_t)
#                         recent_rates.append(rate)
                    
#                     if len(recent_rates) >= 4:
#                         avg_rate = np.mean(recent_rates[-4:])
#                         if avg_rate > 0 and all(abs(r - avg_rate) / avg_rate < 0.1 for r in recent_rates[-4:]):
#                             print(f"Steady state reached at step {step_count}")
#                             break
#                     h2_formation_history = h2_formation_history[-5:]
        
#         print(f"Gillespie finished: {step_count} steps, final_time={self.time:.6f}, H2_formed={self.h2_molecules_formed}")
#         return []





# class GrainCache:
#     """
#     Cache system for pre-generated grain structures.
#     Avoids regenerating identical grains across multiple simulation runs.
#     """
    
#     def __init__(self, cache_dir="grain_cache"):
#         self.cache_dir = Path(cache_dir)
#         self.cache_dir.mkdir(exist_ok=True)
        
#     def _get_cache_key(self, grain_radius_um, site_area_angstroms_sq, 
#                        E_phys_mean_meV, heterogeneity_E_bind_sigma_meV):
#         """Generate unique cache key from grain parameters"""
#         params_str = f"{grain_radius_um}_{site_area_angstroms_sq}_{E_phys_mean_meV}_{heterogeneity_E_bind_sigma_meV}"
#         return hashlib.md5(params_str.encode()).hexdigest()
    
#     def _get_cache_path(self, cache_key):
#         """Get file path for cached grain"""
#         return self.cache_dir / f"grain_{cache_key}.pkl"
    
#     def save_grain(self, grain_data, cache_key):
#         """Save grain structure to disk"""
#         cache_path = self._get_cache_path(cache_key)
#         with open(cache_path, 'wb') as f:
#             pickle.dump(grain_data, f, protocol=pickle.HIGHEST_PROTOCOL)
#         print(f"Grain cached to {cache_path}")
    
#     def load_grain(self, cache_key):
#         """Load grain structure from disk"""
#         cache_path = self._get_cache_path(cache_key)
#         if cache_path.exists():
#             with open(cache_path, 'rb') as f:
#                 return pickle.load(f)
#         return None
    
#     def exists(self, cache_key):
#         """Check if grain is cached"""
#         return self._get_cache_path(cache_key).exists()


#     # Usage instructions:
#     """
#     1. First run with new grain parameters will generate and cache the grain
#     2. Subsequent runs with SAME parameters load from cache instantly
#     3. Different parameters (grain size, site area, etc.) create new cache entries
    
#     Expected speedup:
#     - First run: ~217s initialization (same as before)
#     - Subsequent runs: ~2-5s initialization (40-100× faster!)
    
#     Cache location: ./grain_cache/ directory
#     Cache persistence: Survives across Python sessions, can be reused indefinitely
#     """

#     def generate_grain_library():
#         """Pre-generate common grain configurations"""
#         cache = GrainCache()
        
#         # Common ISM grain sizes
#         grain_sizes = [0.005, 0.01, 0.02, 0.05, 0.1]  # μm
#         site_areas = [50, 60, 80]  # Å²
        
#         print("Generating grain library...")
#         for grain_radius_um in grain_sizes:
#             for site_area_angstroms_sq in site_areas:
#                 params = {
#                     'grain_radius_um': grain_radius_um,
#                     'site_area_angstroms_sq': site_area_angstroms_sq,
#                     'E_phys_mean_meV': 45.0,
#                     'heterogeneity_E_bind_sigma_meV': 5.0,
#                     'initial_h_coverage': 0.0,
#                     'use_3d_lattice': False,
#                     'porosity_fraction': 0.2,
#                 }
                
#                 cache_key = cache._get_cache_key(
#                     grain_radius_um, site_area_angstroms_sq, 45.0, 5.0
#                 )
                
#                 if not cache.exists(cache_key):
#                     print(f"Generating grain: r={grain_radius_um}μm, area={site_area_angstroms_sq}Ų")
#                     kmc = KineticMonteCarlo(params)
#                     # Grain is automatically cached during __init__
#                 else:
#                     print(f"Grain already cached: r={grain_radius_um}μm, area={site_area_angstroms_sq}Ų")
        
#         print(f"Grain library complete! Generated {len(grain_sizes) * len(site_areas)} configurations")

# if __name__ == "__main__":
#     # Pre-generate library (run once)
#     generate_grain_library()

import numpy as np
import random
from scientific_data import surface_chemistry_data, uv_photon_flux, EV_TO_KELVIN, K_B, M_H
from physical_rates import (
    h_diffusion_rate, h_desorption_rate, h2_formation_lh_rate, h2_formation_er_rate,
    uv_photodesorption_rate, uv_h2_formation_rate, adsorption_rate
)
import pickle
import hashlib
import os
from pathlib import Path


class GrainCache:
    """
    Cache system for pre-generated grain structures.
    Avoids regenerating identical grains across multiple simulation runs.
    """
    
    def __init__(self, cache_dir="grain_cache"):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)
        
    def _get_cache_key(self, grain_radius_um, site_area_angstroms_sq, 
                       E_phys_mean_meV, heterogeneity_E_bind_sigma_meV):
        """Generate unique cache key from grain parameters"""
        params_str = f"{grain_radius_um}_{site_area_angstroms_sq}_{E_phys_mean_meV}_{heterogeneity_E_bind_sigma_meV}"
        return hashlib.md5(params_str.encode()).hexdigest()
    
    def _get_cache_path(self, cache_key):
        """Get file path for cached grain"""
        return self.cache_dir / f"grain_{cache_key}.pkl"
    
    def save_grain(self, grain_data, cache_key):
        """Save grain structure to disk"""
        cache_path = self._get_cache_path(cache_key)
        with open(cache_path, 'wb') as f:
            pickle.dump(grain_data, f, protocol=pickle.HIGHEST_PROTOCOL)
        print(f"Grain cached to {cache_path}")
    
    def load_grain(self, cache_key):
        """Load grain structure from disk"""
        cache_path = self._get_cache_path(cache_key)
        if cache_path.exists():
            with open(cache_path, 'rb') as f:
                return pickle.load(f)
        return None
    
    def exists(self, cache_key):
        """Check if grain is cached"""
        return self._get_cache_path(cache_key).exists()


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
        
        # CRITICAL: Initialize with caching
        self.grain_cache = GrainCache()
        self.initialize_off_lattice_structure_cached()

    def _initialize_h_atoms(self, initial_coverage):
        """
        Initialize H atoms on off-lattice surface sites.
        
        Parameters:
        -----------
        initial_coverage : float
            Fraction of surface sites to populate with H atoms (0.0 to 1.0)
        """
        if initial_coverage <= 0:
            return
        
        # Get surface sites
        surface_sites = [i for i, coords in enumerate(self.molecular_coordinates)
                         if self._is_surface_site(coords)]
        
        if surface_sites:
            num_initial_h = int(len(surface_sites) * initial_coverage)
            if num_initial_h > 0:
                selected_sites = random.sample(surface_sites, min(num_initial_h, len(surface_sites)))
                for site_index in selected_sites:
                    self.site_occupancy[site_index] = "H"
                    self.h_atoms_on_surface += 1
                self._update_adjacent_h_pairs_count()

    def initialize_off_lattice_structure_cached(self):
        """
        Initialize grain structure with caching support.
        Checks cache first, generates only if needed.
        """
        grain_radius_um = self.simulation_parameters.get("grain_radius_um", 0.1)
        site_area_angstroms_sq = self.simulation_parameters.get("site_area_angstroms_sq", 9)
        
        # Generate cache key
        cache_key = self.grain_cache._get_cache_key(
            grain_radius_um, site_area_angstroms_sq,
            self.E_bind_mean_meV, self.E_bind_sigma_meV
        )
        
        # Try to load from cache
        cached_data = self.grain_cache.load_grain(cache_key)
        
        if cached_data is not None:
            print(f"Loading grain from cache (key: {cache_key[:8]}...)")
            self._load_grain_from_cache(cached_data)
        else:
            print(f"Generating new grain (will be cached with key: {cache_key[:8]}...)")
            self._generate_new_grain()
            
            # Save to cache
            grain_data = self._extract_grain_data()
            self.grain_cache.save_grain(grain_data, cache_key)
        
        # Initialize H atoms (always fresh, not cached)
        initial_h_coverage = self.simulation_parameters.get("initial_h_coverage", 0.0)
        self.h_atoms_on_surface = 0
        self.adjacent_h_pairs_count = 0
        
        if initial_h_coverage > 0:
            self._initialize_h_atoms(initial_h_coverage)

    def _generate_new_grain(self):
        """Generate grain from scratch (original initialization)"""
        grain_radius_um = self.simulation_parameters.get("grain_radius_um", 0.1)
        grain_radius_cm = grain_radius_um * 1e-4
        site_area_angstroms_sq = self.simulation_parameters.get("site_area_angstroms_sq", 9)
        site_area_cm2 = site_area_angstroms_sq * 1e-16
        grain_surface_area_cm2 = 4 * np.pi * (grain_radius_cm**2)
        total_sites_needed = int(grain_surface_area_cm2 / site_area_cm2)
        
        # CRITICAL: Initialize ALL attributes FIRST
        self.molecular_coordinates = []
        self.site_types = {}
        self.site_occupancy = {}
        self.site_energetics = {}
        self.neighbor_cache = {}
        
        # Calculate grain radius in Ångstroms (needed by _generate_off_lattice_sites)
        self.grain_radius_angstroms = grain_radius_um * 1e4
        
        # Now generate
        self._generate_off_lattice_sites(total_sites_needed)
        self._generate_realistic_surface_topology()
        self._generate_energy_maps()
        self._build_neighbor_cache()

    def _extract_grain_data(self):
        """Extract cacheable grain data (structure + energetics)"""
        return {
            'molecular_coordinates': self.molecular_coordinates.copy(),
            'site_types': self.site_types.copy(),
            'site_energetics': self.site_energetics.copy(),
            'neighbor_cache': self.neighbor_cache.copy(),
            'grain_radius_angstroms': self.grain_radius_angstroms
        }

    def _load_grain_from_cache(self, cached_data):
        """Load grain structure from cached data"""
        self.molecular_coordinates = cached_data['molecular_coordinates']
        self.site_types = cached_data['site_types']
        self.site_energetics = cached_data['site_energetics']
        self.neighbor_cache = cached_data['neighbor_cache']
        self.grain_radius_angstroms = cached_data['grain_radius_angstroms']
        
        # Initialize fresh occupancy dictionary
        self.site_occupancy = {i: None for i in range(len(self.molecular_coordinates))}
        
        print(f"Loaded grain: {len(self.molecular_coordinates)} sites from cache")

    def _generate_off_lattice_sites(self, total_sites_needed):
        """
        Generate true off-lattice sites using continuous molecular coordinates.
        Uses realistic amorphous carbon structure with irregular bond networks.
        """
        rng = np.random.default_rng()
        
        # Note: attributes already initialized in _generate_new_grain, don't re-initialize here
        
        # Generate realistic amorphous carbon structure
        # Based on Draine & Li 2007: amorphous carbon has irregular bond networks
        sites_placed = 0
        attempts = 0
        max_attempts = total_sites_needed * 1000
        
        while sites_placed < total_sites_needed and attempts < max_attempts:
            # Generate random position in 3D space within grain volume
            # Use spherical coordinates for proper grain geometry
            r = self.grain_radius_angstroms * (rng.uniform(0.8, 1.0))**(1/3)  # Concentrate near surface
            theta = rng.uniform(0, np.pi)
            phi = rng.uniform(0, 2 * np.pi)
            
            # Convert to Cartesian coordinates
            x = r * np.sin(theta) * np.cos(phi)
            y = r * np.sin(theta) * np.sin(phi)
            z = r * np.cos(theta)
            
            # Check minimum distance from existing sites
            too_close = False
            for existing_coords in self.molecular_coordinates:
                distance = np.sqrt((x - existing_coords[0])**2 + 
                                 (y - existing_coords[1])**2 + 
                                 (z - existing_coords[2])**2)
                if distance < 2.0:  # Minimum 2 Å separation
                    too_close = True
                    break
            
            if not too_close:
                site_index = len(self.molecular_coordinates)
                self.molecular_coordinates.append((x, y, z))
                self.site_occupancy[site_index] = None
                sites_placed += 1
            
            attempts += 1
        
        # Ensure minimum surface coverage
        surface_sites = [i for i, coords in enumerate(self.molecular_coordinates) 
                        if coords[0]**2 + coords[1]**2 + coords[2]**2 > (0.95 * self.grain_radius_angstroms)**2]
        
        if len(surface_sites) < 10:
            # Add more surface sites if needed
            additional_sites = 10 - len(surface_sites)
            for _ in range(additional_sites):
                r = self.grain_radius_angstroms * 0.98  # Near surface
                theta = rng.uniform(0, np.pi)
                phi = rng.uniform(0, 2 * np.pi)
                
                x = r * np.sin(theta) * np.cos(phi)
                y = r * np.sin(theta) * np.sin(phi)
                z = r * np.cos(theta)
                
                site_index = len(self.molecular_coordinates)
                self.molecular_coordinates.append((x, y, z))
                self.site_occupancy[site_index] = None

    def _generate_realistic_surface_topology(self):
        """
        Generate realistic amorphous carbon surface topology with irregular bond networks,
        surface roughness, and proper nearest-neighbor algorithms.
        """
        rng = np.random.default_rng()
        
        # Calculate surface roughness based on amorphous carbon properties
        # Draine & Li 2007: amorphous carbon has significant surface roughness
        surface_roughness = 0.1 * self.grain_radius_angstroms  # 10% of grain radius
        
        # Generate surface topology with irregular features
        for i, coords in enumerate(self.molecular_coordinates):
            x, y, z = coords
            
            # Calculate distance from center
            distance_from_center = np.sqrt(x**2 + y**2 + z**2)
            
            # Determine if site is on surface (within roughness layer)
            is_surface_site = distance_from_center > (self.grain_radius_angstroms - surface_roughness)
            
            if is_surface_site:
                # Surface sites: determine type based on local environment
                neighbors = self._find_nearest_neighbors(i, cutoff_distance=3.0)
                coordination_number = len(neighbors)
                
                # Physics-based site type assignment
                if coordination_number <= 2:
                    # Low coordination: dangling bonds → defect sites
                    self.site_types[i] = 3
                elif coordination_number <= 4:
                    # Medium coordination: unsaturated bonds → chemisorption
                    self.site_types[i] = 2
                else:
                    # High coordination: saturated bonds → physisorption
                    self.site_types[i] = 1
            else:
                # Subsurface sites: bulk-like properties
                self.site_types[i] = 1

    def _find_nearest_neighbors(self, site_index, cutoff_distance):
        """
        Dynamic nearest-neighbor algorithm using continuous molecular coordinates.
        Finds all sites within cutoff_distance of the target site.
        NOTE: This is O(N) and should only be called during initialization!
        """
        target_coords = self.molecular_coordinates[site_index]
        neighbors = []
        
        for i, coords in enumerate(self.molecular_coordinates):
            if i == site_index:
                continue
            
            distance = np.sqrt((target_coords[0] - coords[0])**2 + 
                             (target_coords[1] - coords[1])**2 + 
                             (target_coords[2] - coords[2])**2)
            
            if distance <= cutoff_distance:
                neighbors.append(i)
        
        return neighbors

    def _build_neighbor_cache(self):
        """
        Build neighbor cache ONCE during initialization.
        This is the KEY OPTIMIZATION that eliminates O(N²) lookups during simulation.
        """
        print("Building neighbor cache...")
        self.neighbor_cache = {}
        for i in range(len(self.molecular_coordinates)):
            self.neighbor_cache[i] = self._find_nearest_neighbors(i, cutoff_distance=3.0)
        print(f"Neighbor cache built: {len(self.neighbor_cache)} sites")

    def _generate_energy_maps(self):
        """
        Generate energy maps for off-lattice structure with site-specific energetics.
        """
        rng = np.random.default_rng()
        
        mean_eV = self.E_bind_mean_meV * 1e-3
        sigma_eV = self.E_bind_sigma_meV * 1e-3
        
        for i, coords in enumerate(self.molecular_coordinates):
            site_type = self.site_types.get(i, 1)
            neighbors = self._find_nearest_neighbors(i, cutoff_distance=3.0)
            num_neighbors = len(neighbors)
            
            # Calculate local environment effects
            local_density = self._calculate_local_density(i)
            surface_distance = self._calculate_surface_distance(coords)
            
            # Assign binding energies based on site type and local environment
            if site_type == 1:  # Physisorption
                base_energy = mean_eV
                # Local density effect: higher density → stronger binding
                density_factor = 1.0 + 0.2 * local_density
                # Surface distance effect: closer to surface → stronger binding
                surface_factor = 1.0 + 0.1 * (1.0 - surface_distance / self.grain_radius_angstroms)
                
                binding_energy = base_energy * density_factor * surface_factor + rng.normal(0, sigma_eV)
                diffusion_barrier = 0.025 + rng.normal(0, 0.005)
                
            elif site_type == 2:  # Chemisorption
                binding_energy = 0.25 + rng.normal(0, 0.05)
                diffusion_barrier = 0.03 + rng.normal(0, 0.005)
                
            elif site_type == 3:  # Defect
                binding_energy = 0.35 + rng.normal(0, 0.05)
                diffusion_barrier = 0.015 + rng.normal(0, 0.003)
            
            # Ensure positive values
            binding_energy = max(binding_energy, 0.01)
            diffusion_barrier = max(diffusion_barrier, 0.01)
            
            self.site_energetics[i] = (binding_energy, diffusion_barrier)

    def _calculate_local_density(self, site_index):
        """
        Calculate local density around a site using continuous coordinates.
        """
        target_coords = self.molecular_coordinates[site_index]
        local_volume = (4.0/3) * np.pi * (3.0)**3  # 3 Å radius sphere
        neighbors_in_volume = 0
        
        for coords in self.molecular_coordinates:
            distance = np.sqrt((target_coords[0] - coords[0])**2 + 
                             (target_coords[1] - coords[1])**2 + 
                             (target_coords[2] - coords[2])**2)
            if distance <= 3.0:
                neighbors_in_volume += 1
        
        # Normalize to density (sites per Å³)
        density = (neighbors_in_volume - 1) / local_volume  # -1 to exclude self
        return density

    def _calculate_surface_distance(self, coords):
        """
        Calculate distance from surface for a given coordinate.
        """
        distance_from_center = np.sqrt(coords[0]**2 + coords[1]**2 + coords[2]**2)
        surface_distance = abs(distance_from_center - self.grain_radius_angstroms)
        return surface_distance

    def _is_surface_site(self, coords):
        """
        Determine if a site is on the surface using continuous coordinates.
        """
        distance_from_center = np.sqrt(coords[0]**2 + coords[1]**2 + coords[2]**2)
        surface_thickness = 0.05 * self.grain_radius_angstroms  # 5% of radius
        return distance_from_center > (self.grain_radius_angstroms - surface_thickness)

    def get_neighbors_3d(self, site_index):
        """
        Get neighbors for off-lattice structure using CACHED neighbor list.
        This is now O(1) instead of O(N)!
        """
        return self.neighbor_cache.get(site_index, [])

    def get_accessible_surface_sites(self):
        """
        Get accessible surface sites (not occupied by H atoms).
        """
        accessible_sites = []
        for i, coords in enumerate(self.molecular_coordinates):
            if self._is_surface_site(coords) and self.site_occupancy.get(i) != "H":
                accessible_sites.append(i)
        return accessible_sites

    def get_occupied_sites(self):
        """
        Get all sites occupied by H atoms in off-lattice structure.
        """
        occupied_sites = []
        for i, occupancy in self.site_occupancy.items():
            if occupancy == "H":
                occupied_sites.append(i)
        return occupied_sites

    def get_site_occupant(self, site_index):
        """
        Get the occupant of a site in the off-lattice structure.
        """
        return self.site_occupancy.get(site_index, None)

    def occupy_site(self, site_index, species):
        """
        Occupy a site in the off-lattice structure.
        """
        if site_index < len(self.molecular_coordinates):
            self.site_occupancy[site_index] = species

    def vacate_site(self, site_index):
        """
        Vacate a site in the off-lattice structure.
        """
        if site_index < len(self.molecular_coordinates):
            self.site_occupancy[site_index] = None

    def _update_adjacent_h_pairs_count(self):
        """
        Initial full calculation of adjacent H pairs.
        Only called during initialization!
        """
        self.adjacent_h_pairs_count = 0
        for site_index in self.get_occupied_sites():
            neighbors = self.get_neighbors_3d(site_index)
            for neighbor_index in neighbors:
                if self.get_site_occupant(neighbor_index) == "H" and site_index < neighbor_index:
                    self.adjacent_h_pairs_count += 1

    def update_adjacent_h_pairs_count(self, site_index):
        """
        OPTIMIZED: Incremental update - only check neighbors of changed site.
        This is now O(k) where k is number of neighbors, instead of O(N²)!
        """
        neighbors = self.get_neighbors_3d(site_index)
        
        # Check if site is now occupied or was just vacated
        is_occupied = (self.get_site_occupant(site_index) == "H")
        
        # Update pair count based on neighbors
        for neighbor_index in neighbors:
            if self.get_site_occupant(neighbor_index) == "H":
                # Ensure we don't double count (only count each pair once)
                if site_index < neighbor_index:
                    if is_occupied:
                        self.adjacent_h_pairs_count += 1
                    else:
                        self.adjacent_h_pairs_count -= 1

    def _calculate_realistic_lateral_interactions(self, site_index, surface_temp_k):
        """
        DISABLED FOR PERFORMANCE: Lateral interactions are O(N²) and have minimal impact.
        Returns 0.0 to maintain interface compatibility.
        """
        # PERFORMANCE: This function was causing massive slowdowns (O(N²) per rate calculation)
        # Lateral interactions have minimal physical impact on rates
        # Can be re-enabled for detailed studies if needed
        return 0.0

    def _calculate_lh_coverage_enhancement(self, current_coverage, surface_temp_k):
        """
        Calculate coverage-dependent enhancement factor for LH H2 formation.
        At low coverage: reduced formation (few neighbors)
        At medium coverage: enhanced formation (optimal neighbor density)
        At high coverage: reduced formation (site blocking)
        """
        if current_coverage <= 0.0:
            return 0.0
        
        # Optimal coverage around 0.1-0.3 ML (monolayer)
        if current_coverage < 0.1:
            # Low coverage: linearly increase with coverage
            enhancement = current_coverage / 0.1
        elif current_coverage < 0.3:
            # Optimal range: full enhancement
            enhancement = 1.0
        else:
            # High coverage: site blocking reduces formation
            enhancement = 1.0 - 0.5 * (current_coverage - 0.3) / 0.7
        
        # Ensure non-negative
        enhancement = max(0.0, enhancement)
        
        # Temperature-dependent boost (higher T = more mobility = better LH)
        temp_factor = 1.0 + 0.1 * (surface_temp_k - 10.0) / 40.0
    
        return enhancement * temp_factor
    
    def _calculate_diffusion_rate_with_tunneling(self, site_index, surface_temp_k):
        hbar = 1.055e-27
        mH = 1.674e-24
        kB = 1.381e-16
        eV_to_erg = 1.602e-12
        
        E_diff_eV = self.site_energetics[site_index][1]
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
        OPTIMIZED: Removed O(N²) lateral interaction calculations.
        """
        rates = {}
        surface_temp_k = self.simulation_parameters.get("surface_temperature_k")
        gas_temp_k = self.simulation_parameters.get("gas_temperature_k")
        h_gas_density = self.simulation_parameters.get("h_gas_density_cm3")
        uv_flux_factor = self.simulation_parameters.get("uv_flux_factor", 1.0)
        site_area_cm2 = self.simulation_parameters.get("site_area_angstroms_sq", 9) * 1e-16
        sticking_probability = self.simulation_parameters.get("sticking_probability", 0.3)

        accessible_sites = self.get_accessible_surface_sites()
        num_accessible_sites = len(accessible_sites)
        accessible_area_cm2 = num_accessible_sites * site_area_cm2
        
        # Coverage effects: reduce available sites
        current_coverage = self.h_atoms_on_surface / max(1, num_accessible_sites)
        available_sites_fraction = max(0.0, 1.0 - current_coverage)
        effective_area_cm2 = accessible_area_cm2 * available_sites_fraction

        # 1. Adsorption rate - physically derived from gas kinetics with site blocking
        rates["adsorption"] = adsorption_rate(
            h_gas_density, gas_temp_k, sticking_probability, effective_area_cm2
        )

        # 2. Eley-Rideal formation - gas impingement × cross-section × reaction probability
        if self.h_atoms_on_surface > 0:
            # Calculate gas flux for ER
            v_thermal = np.sqrt(8 * K_B * gas_temp_k * 1.602e-19 / (np.pi * M_H))
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
                
                for site_index in occupied_sites:
                    site_type = self.site_types.get(site_index, 1)
                    base_binding_energy = self.site_energetics[site_index][0]
                    
                    # PERFORMANCE: Lateral interactions disabled (was O(N²))
                    lateral_energy = 0.0
                    binding_energy = base_binding_energy + lateral_energy
                    
                    # Desorption rate using TST with coverage effects
                    desorption_rate = h_desorption_rate(binding_energy, surface_temp_k)
                    total_desorption_rate += desorption_rate
                    
                    # Diffusion rate using TST + tunneling
                    diffusion_rate = h_diffusion_rate(site_type, surface_temp_k)
                    total_diffusion_rate += diffusion_rate
                
                rates["desorption"] = total_desorption_rate
                rates["diffusion"] = total_diffusion_rate

        # 4. Langmuir-Hinshelwood formation - TST for surface reactions with realistic coverage effects
        if self.adjacent_h_pairs_count > 0:
            # Physics-based coverage effects: electronic structure changes and site availability
            coverage_enhancement = self._calculate_lh_coverage_enhancement(current_coverage, surface_temp_k)
            enhanced_lh_rate = h2_formation_lh_rate(
                surface_temp_k, self.adjacent_h_pairs_count
            ) * coverage_enhancement
            rates["h2_formation_LH"] = enhanced_lh_rate

        # 5. UV processes - photon flux × cross-section × yield
        if uv_flux_factor > 0:
            uv_photon_flux_total = uv_photon_flux["integrated_fuv_photon_flux_photons_cm2_s"] * uv_flux_factor
            
            # UV pulse rate from Andersson et al. 2006, A&A 453, 459 (interstellar UV field)
            base_uv_rate = 5.0e-8  # 5 photons grain⁻¹ yr⁻¹ in diffuse ISM
            if self.uv_pulse_enabled:
                rates["uv_pulse_start"] = base_uv_rate * uv_flux_factor
            
            if self.uv_pulse_active:
                pulse_end_rate = 1.0 / self.uv_pulse_duration
                rates["uv_pulse_end"] = pulse_end_rate
            
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
                if len(accessible_surface_sites) > 0:
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
            if len(accessible_sites) > 0:
                site_index = random.choice(accessible_sites)
                self.occupy_site(site_index, "H")
                self.h_atoms_on_surface += 1
                self.total_adsorbed_h_atoms += 1
                self.update_adjacent_h_pairs_count(site_index)
                
        elif event_type in ["desorption", "uv_photodesorption"]:
            if self.h_atoms_on_surface > 0:
                occupied_sites = self.get_occupied_sites()
                if occupied_sites:
                    site_index = random.choice(occupied_sites)
                    self.vacate_site(site_index)
                    self.h_atoms_on_surface -= 1
                    self.total_desorbed_h_atoms += 1
                    self.update_adjacent_h_pairs_count(site_index)
                    
        elif event_type in ["diffusion", "uv_stimulated_diffusion"]:
            if self.h_atoms_on_surface > 0:
                occupied_sites = self.get_occupied_sites()
                if occupied_sites:
                    mobile_atoms = [
                        site_index for site_index in occupied_sites 
                        if any(self.get_site_occupant(nd) != "H" for nd in self.get_neighbors_3d(site_index))
                    ]
                    if mobile_atoms:
                        source_site_index = random.choice(mobile_atoms)
                        empty_neighbors = [
                            neighbor_index for neighbor_index in self.get_neighbors_3d(source_site_index) 
                            if self.get_site_occupant(neighbor_index) != "H"
                        ]
                        if empty_neighbors:
                            target_site_index = random.choice(empty_neighbors)
                            self.vacate_site(source_site_index)
                            self.update_adjacent_h_pairs_count(source_site_index)
                            self.occupy_site(target_site_index, "H")
                            self.update_adjacent_h_pairs_count(target_site_index)
                            
        elif event_type in ["h2_formation_LH", "h2_formation_UV"]:
            if self.adjacent_h_pairs_count > 0:
                pairs = []
                pair_rates = []
                for site_index in self.get_occupied_sites():
                    for neighbor_index in self.get_neighbors_3d(site_index):
                        if self.get_site_occupant(neighbor_index) == "H" and site_index < neighbor_index:
                            avg_binding = (self.site_energetics[site_index][0] + self.site_energetics[neighbor_index][0]) / 2
                            site_rate = avg_binding * np.exp(-0.02 / (K_B * self.simulation_parameters.get("surface_temperature_k", 10)))
                            pairs.append((site_index, neighbor_index))
                            pair_rates.append(site_rate)
                
                if pairs:
                    if len(pair_rates) > 1:
                        total_rate = sum(pair_rates)
                        probs = [r/total_rate for r in pair_rates]
                        chosen_idx = random.choices(range(len(pairs)), weights=probs, k=1)[0]
                        site1, site2 = pairs[chosen_idx]
                    else:
                        site1, site2 = pairs[0]
                        
                    self.vacate_site(site1)
                    self.update_adjacent_h_pairs_count(site1)
                    self.vacate_site(site2)
                    self.update_adjacent_h_pairs_count(site2)
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
                    site_index = random.choice(occupied_sites)
                    self.vacate_site(site_index)
                    self.h_atoms_on_surface -= 1
                    self.update_adjacent_h_pairs_count(site_index)
                    self.h2_molecules_formed_ER += 1
                    self.h2_molecules_formed += 1
                    
        elif event_type == "uv_pulse_start":
            self.uv_pulse_active = True
            self.last_uv_pulse_time = self.time
            
        elif event_type == "uv_pulse_end":
            self.uv_pulse_active = False

    def run_gillespie(self, max_time, max_steps=None):
        step_count = 0
        steady_state_check_interval = 1000
        h2_formation_history = []
        
        print(f"Starting Gillespie: max_time={max_time}, num_sites={len(self.molecular_coordinates)}")
        print(f"Initial H atoms: {self.h_atoms_on_surface}, Initial H2: {self.h2_molecules_formed}")
        
        while self.time < max_time:
            if max_steps and step_count >= max_steps:
                print(f"Reached max_steps={max_steps}, breaking")
                break
            
            # Diagnostic every 1000 steps
            if step_count % 1000 == 0:
                print(f"Step {step_count}, time={self.time:.6f}, H_atoms={self.h_atoms_on_surface}, H2={self.h2_molecules_formed}")
            
            rates = self.calculate_rates()
            if not rates:
                print(f"No rates available at step {step_count}, breaking")
                break
            total_rate = sum(rates.values())
            if total_rate == 0:
                print(f"Total rate is zero at step {step_count}, breaking")
                break
            
            delta_t = random.expovariate(total_rate)
            self.time += delta_t
            
            if self.time > max_time:
                self.time = max_time
                print(f"Reached max_time at step {step_count}")
                break
            
            chosen_event = random.choices(list(rates.keys()), weights=list(rates.values()), k=1)[0]
            self.execute_event(chosen_event)
            step_count += 1
            
            if step_count % steady_state_check_interval == 0:
                h2_formation_history.append(self.h2_molecules_formed)
                if len(h2_formation_history) >= 5:
                    recent_rates = []
                    for i in range(1, len(h2_formation_history)):
                        rate = (h2_formation_history[i] - h2_formation_history[i-1]) / (steady_state_check_interval * delta_t)
                        recent_rates.append(rate)
                    
                    if len(recent_rates) >= 4:
                        avg_rate = np.mean(recent_rates[-4:])
                        if avg_rate > 0 and all(abs(r - avg_rate) / avg_rate < 0.1 for r in recent_rates[-4:]):
                            print(f"Steady state reached at step {step_count}")
                            break
                    h2_formation_history = h2_formation_history[-5:]
        
        print(f"Gillespie finished: {step_count} steps, final_time={self.time:.6f}, H2_formed={self.h2_molecules_formed}")
        return []


def generate_grain_library():
    """Pre-generate common grain configurations"""
    cache = GrainCache()
    
    # Common ISM grain sizes
    grain_sizes = [0.005, 0.01, 0.02, 0.05, 0.1]  # μm
    site_areas = [50, 60, 80]  # Å²
    
    print("Generating grain library...")
    for grain_radius_um in grain_sizes:
        for site_area_angstroms_sq in site_areas:
            params = {
                'grain_radius_um': grain_radius_um,
                'site_area_angstroms_sq': site_area_angstroms_sq,
                'E_phys_mean_meV': 45.0,
                'heterogeneity_E_bind_sigma_meV': 5.0,
                'initial_h_coverage': 0.0,
                'use_3d_lattice': False,
                'porosity_fraction': 0.2,
                'surface_temperature_k': 10.0,
                'gas_temperature_k': 100.0,
                'h_gas_density_cm3': 1.0,
                'uv_flux_factor': 1.0,
                'sticking_probability': 0.3,
            }
            
            cache_key = cache._get_cache_key(
                grain_radius_um, site_area_angstroms_sq, 45.0, 5.0
            )
            
            if not cache.exists(cache_key):
                print(f"  Generating: radius={grain_radius_um}μm, area={site_area_angstroms_sq}Å²")
                kmc = KineticMonteCarlo(params)
                # Grain is automatically cached during __init__
            else:
                print(f"  Already cached: radius={grain_radius_um}μm, area={site_area_angstroms_sq}Å²")
    
    print(f"Grain library complete! Generated {len(grain_sizes) * len(site_areas)} configurations")


if __name__ == "__main__":
    # Pre-generate library (run once)
    generate_grain_library()
