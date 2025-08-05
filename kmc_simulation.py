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
        self.E_bind_mean_meV = float(self.simulation_parameters.get(
            "heterogeneity_E_bind_mean_meV",
            surface_chemistry_data["h_physisorption_binding_energy_mev_typical"]
        ))
        self.E_bind_sigma_meV = float(self.simulation_parameters.get("heterogeneity_E_bind_sigma_meV", 5.0))
        self.E_diff_mean_eV = float(self.simulation_parameters.get(
            "heterogeneity_E_diff_mean_eV",
            surface_chemistry_data["h_physisorbed_amorphous_diffusion_barrier_ev_min"]
        ))
        self.E_diff_sigma_eV = float(self.simulation_parameters.get("heterogeneity_E_diff_sigma_eV", 0.005))
        self.initialize_surface()

    def initialize_surface(self):
        grain_radius_um = self.simulation_parameters.get("grain_radius_um", 0.1)
        grain_radius_cm = grain_radius_um * 1e-4
        site_area_angstroms_sq = self.simulation_parameters.get("site_area_angstroms_sq", 9)
        site_area_cm2 = site_area_angstroms_sq * 1e-16
        grain_surface_area_cm2 = 4 * np.pi * (grain_radius_cm**2)
        calculated_sites = int(grain_surface_area_cm2 / site_area_cm2)
        self.surface_dimension = int(np.sqrt(calculated_sites))
        self.total_surface_sites = self.surface_dimension * self.surface_dimension
        self.surface_sites = np.full((self.surface_dimension, self.surface_dimension), None, dtype=object)

        if self.use_site_heterogeneity:
            rng = np.random.default_rng()
            E_bind_eV = (rng.normal(self.E_bind_mean_meV, self.E_bind_sigma_meV, size=(self.surface_dimension, self.surface_dimension))
                         .clip(min=1.0) / 1000.0)
            E_diff_eV = rng.normal(self.E_diff_mean_eV, self.E_diff_sigma_eV, size=(self.surface_dimension, self.surface_dimension))
            E_diff_eV = np.clip(E_diff_eV, a_min=0.0, a_max=np.maximum(0.0, E_bind_eV))
            self.E_bind_eV_map = E_bind_eV
            self.E_diff_eV_map = E_diff_eV
        else:
            self.E_bind_eV_map = np.full((self.surface_dimension, self.surface_dimension),
                                         surface_chemistry_data["h_physisorption_binding_energy_mev_typical"]/1000.0)
            self.E_diff_eV_map = np.full((self.surface_dimension, self.surface_dimension),
                                         surface_chemistry_data["h_physisorbed_amorphous_diffusion_barrier_ev_min"])
        
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

    def get_neighbors(self, r, c):
        neighbors = []
        moves = [(-1, 0), (1, 0), (0, -1), (0, 1)]
        for dr, dc in moves:
            nr, nc = (r + dr) % self.surface_dimension, (c + dc) % self.surface_dimension
            neighbors.append((nr, nc))
        return neighbors

    def update_adjacent_h_pairs_count(self, r, c, add_atom):
        change = 1 if add_atom else -1
        for nr, nc in self.get_neighbors(r, c):
            if self.surface_sites[nr, nc] == "H":
                self.adjacent_h_pairs_count += change

    def calculate_rates(self):
        rates = {}
        surface_temp_k = self.simulation_parameters.get("surface_temperature_k")
        gas_temp_k = self.simulation_parameters.get("gas_temperature_k")
        h_gas_density = self.simulation_parameters.get("h_gas_density_cm3")
        uv_flux_factor = self.simulation_parameters.get("uv_flux_factor", 1.0)
        site_area_cm2 = self.simulation_parameters.get("site_area_angstroms_sq", 9) * 1e-16
        pre_exp_frequency = 1e12

        v_h_thermal = np.sqrt(8 * K_B * gas_temp_k / (np.pi * M_H))
        gas_flux_per_cm2_s = 0.25 * h_gas_density * v_h_thermal
        sticking_prob = self.simulation_parameters.get("sticking_probability", 0.3)
        rates["adsorption"] = gas_flux_per_cm2_s * (self.total_surface_sites * site_area_cm2) * sticking_prob
        if self.h_atoms_on_surface > 0:
            rates["h2_formation_ER"] = gas_flux_per_cm2_s * surface_chemistry_data["er_cross_section_cm2"] * self.h_atoms_on_surface

        if self.h_atoms_on_surface > 0:
            occ = (self.surface_sites == "H")
            if np.any(occ):
                E_bind_occ = self.E_bind_eV_map[occ]
                desorption_per_atom = pre_exp_frequency * np.exp(-E_bind_occ * EV_TO_KELVIN / surface_temp_k)
                rates["desorption"] = float(np.sum(desorption_per_atom))
                E_diff_occ = self.E_diff_eV_map[occ]
                diffusion_per_atom = pre_exp_frequency * np.exp(-E_diff_occ * EV_TO_KELVIN / surface_temp_k)
                rates["diffusion"] = float(np.sum(diffusion_per_atom))
        if self.adjacent_h_pairs_count > 0:
            occ = (self.surface_sites == "H")
            if np.any(occ):
                E_diff_mean = float(np.mean(self.E_diff_eV_map[occ]))
            else:
                E_diff_mean = float(np.mean(self.E_diff_eV_map))
            h2_formation_barrier_ev = E_diff_mean
            rates["h2_formation_LH"] = pre_exp_frequency * np.exp(-h2_formation_barrier_ev * EV_TO_KELVIN / surface_temp_k) * self.adjacent_h_pairs_count
        
        if uv_flux_factor > 0:
            uv_photon_flux_total = uv_photon_flux["integrated_fuv_photon_flux_photons_cm2_s"] * uv_flux_factor
            if self.h_atoms_on_surface > 0:
                rates["uv_photodesorption"] = uv_photon_flux_total * site_area_cm2 * surface_chemistry_data["uv_photodesorption_yield_h_atom"] * self.h_atoms_on_surface
            if self.adjacent_h_pairs_count > 0:
                rates["h2_formation_UV"] = uv_photon_flux_total * (2 * site_area_cm2) * surface_chemistry_data["uv_h2_formation_yield_per_pair"] * self.adjacent_h_pairs_count
                if self.simulation_parameters.get("enable_uv_photodissociation", True):
                    yield_pair = float(self.simulation_parameters.get("uv_photodissociation_yield_per_pair", 1e-6))
                    rates["uv_photodissociation_pair"] = uv_photon_flux_total * (2 * site_area_cm2) * yield_pair * self.adjacent_h_pairs_count
            uv_diffusion_factor = float(self.simulation_parameters.get("uv_stimulated_diffusion_factor", 1.0))
            if uv_diffusion_factor != 1.0 and self.h_atoms_on_surface > 0:
                rates["uv_stimulated_diffusion"] = max(0.0, (uv_diffusion_factor - 1.0)) * rates.get("diffusion", 0.0)
        
        return {k: v for k, v in rates.items() if v > 0}

    def execute_event(self, event_type):
        if event_type == "adsorption":
            if self.h_atoms_on_surface < self.total_surface_sites:
                empty_sites_coords = np.argwhere(self.surface_sites == None)
                if len(empty_sites_coords) > 0:
                    r, c = random.choice(empty_sites_coords.tolist())
                    self.surface_sites[r, c] = "H"
                    self.h_atoms_on_surface += 1
                    self.total_adsorbed_h_atoms += 1
                    self.update_adjacent_h_pairs_count(r, c, True)
        elif event_type in ["desorption", "uv_photodesorption"]:
            if self.h_atoms_on_surface > 0:
                occupied_sites_coords = np.argwhere(self.surface_sites == "H").tolist()
                r, c = random.choice(occupied_sites_coords)
                self.surface_sites[r, c] = None
                self.h_atoms_on_surface -= 1
                self.total_desorbed_h_atoms += 1
                self.update_adjacent_h_pairs_count(r, c, False)
        elif event_type == "diffusion":
            if self.h_atoms_on_surface > 0:
                occupied_sites_coords = np.argwhere(self.surface_sites == "H").tolist()
                mobile_atoms = [pos for pos in occupied_sites_coords if any(self.surface_sites[n_pos] is None for n_pos in self.get_neighbors(*pos))]
                if mobile_atoms:
                    r_h, c_h = random.choice(mobile_atoms)
                    empty_neighbors = [n_pos for n_pos in self.get_neighbors(r_h, c_h) if self.surface_sites[n_pos] is None]
                    if empty_neighbors: 
                        r_empty, c_empty = random.choice(empty_neighbors)
                        self.surface_sites[r_h, c_h] = None
                        self.update_adjacent_h_pairs_count(r_h, c_h, False)
                        self.surface_sites[r_empty, c_empty] = "H"
                        self.update_adjacent_h_pairs_count(r_empty, c_empty, True)
        elif event_type == "h2_formation_LH":
            if self.adjacent_h_pairs_count > 0:
                pairs = []
                H = np.argwhere(self.surface_sites == "H")
                for r1, c1 in H:
                    for r2, c2 in self.get_neighbors(r1, c1):
                        if self.surface_sites[r2, c2] == "H":
                            if (r1, c1) < (r2, c2):
                                pairs.append(((r1, c1), (r2, c2)))
                if pairs:
                    (r1, c1), (r2, c2) = random.choice(pairs)
                    self.surface_sites[r1, c1] = None
                    self.update_adjacent_h_pairs_count(r1, c1, False)
                    self.surface_sites[r2, c2] = None
                    self.update_adjacent_h_pairs_count(r2, c2, False)
                    self.h_atoms_on_surface -= 2
                    self.h2_molecules_formed += 1
                    self.h2_molecules_formed_LH += 1
        elif event_type == "h2_formation_UV":
            if self.adjacent_h_pairs_count > 0:
                pairs = []
                H = np.argwhere(self.surface_sites == "H")
                for r1, c1 in H:
                    for r2, c2 in self.get_neighbors(r1, c1):
                        if self.surface_sites[r2, c2] == "H":
                            if (r1, c1) < (r2, c2):
                                pairs.append(((r1, c1), (r2, c2)))
                if pairs:
                    (r1, c1), (r2, c2) = random.choice(pairs)
                    self.surface_sites[r1, c1] = None
                    self.update_adjacent_h_pairs_count(r1, c1, False)
                    self.surface_sites[r2, c2] = None
                    self.update_adjacent_h_pairs_count(r2, c2, False)
                    self.h_atoms_on_surface -= 2
                    self.h2_molecules_formed += 1
                    self.h2_molecules_formed_UV += 1
        elif event_type == "uv_photodissociation_pair":
            if self.adjacent_h_pairs_count > 0:
                pairs = []
                H = np.argwhere(self.surface_sites == "H")
                for r1, c1 in H:
                    for r2, c2 in self.get_neighbors(r1, c1):
                        if self.surface_sites[r2, c2] == "H":
                            if (r1, c1) < (r2, c2):
                                pairs.append(((r1, c1), (r2, c2)))
                if pairs:
                    (r1, c1), (r2, c2) = random.choice(pairs)
                    self.surface_sites[r1, c1] = None
                    self.update_adjacent_h_pairs_count(r1, c1, False)
                    self.surface_sites[r2, c2] = None
                    self.update_adjacent_h_pairs_count(r2, c2, False)
                    self.h_atoms_on_surface -= 2
        elif event_type == "uv_stimulated_diffusion":
            if self.h_atoms_on_surface > 0:
                occupied_sites_coords = np.argwhere(self.surface_sites == "H").tolist()
                mobile_atoms = [pos for pos in occupied_sites_coords if any(self.surface_sites[n_pos] is None for n_pos in self.get_neighbors(*pos))]
                if mobile_atoms:
                    r_h, c_h = random.choice(mobile_atoms)
                    empty_neighbors = [n_pos for n_pos in self.get_neighbors(r_h, c_h) if self.surface_sites[n_pos] is None]
                    if empty_neighbors:
                        r_empty, c_empty = random.choice(empty_neighbors)
                        self.surface_sites[r_h, c_h] = None
                        self.update_adjacent_h_pairs_count(r_h, c_h, False)
                        self.surface_sites[r_empty, c_empty] = "H"
                        self.update_adjacent_h_pairs_count(r_empty, c_empty, True)
        elif event_type == "h2_formation_ER":
            if self.h_atoms_on_surface > 0:
                occupied_sites_coords = np.argwhere(self.surface_sites == "H").tolist()
                r, c = random.choice(occupied_sites_coords)
                self.surface_sites[r, c] = None; self.h_atoms_on_surface -= 1; self.update_adjacent_h_pairs_count(r, c, False); self.h2_molecules_formed_ER += 1; self.h2_molecules_formed += 1

    def run_gillespie(self, max_time, max_steps=None):
        step_count = 0
        while self.time < max_time:
            if max_steps and step_count >= max_steps: break
            rates = self.calculate_rates()
            if not rates: break
            total_rate = sum(rates.values())
            if total_rate == 0: break
            
            delta_t = random.expovariate(total_rate)
            if self.time + delta_t > max_time:
                delta_t = max_time - self.time
            
            self.time += delta_t
            chosen_event = random.choices(list(rates.keys()), weights=list(rates.values()), k=1)[0]
            self.execute_event(chosen_event)
            step_count += 1
        return []

