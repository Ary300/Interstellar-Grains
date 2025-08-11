# FIXED: Updated tests for amorphous structure and coverage effects
import unittest
import numpy as np
from kmc_simulation import KineticMonteCarlo
from scientific_data import surface_chemistry_data, uv_photon_flux, EV_TO_KELVIN, K_B, M_H

class TestKineticMonteCarlo(unittest.TestCase):

    def setUp(self):
        self.sim_params = {
            "grain_radius_um": 0.1,
            "site_area_angstroms_sq": 9,
            "surface_temperature_k": 10.0,
            "gas_temperature_k": 100.0,
            "h_gas_density_cm3": 1e2,
            "sticking_probability": 0.3,
            "initial_h_coverage": 0.0,
            "uv_flux_factor": 1.0,
            "surface_dimension": 20,  # Smaller for tests
            "coverage_effects_enabled": True,
            "lateral_interaction_eV": 0.005,
            "porosity_fraction": 0.2,
            "E_phys_mean_meV": 45.0,
            "E_phys_sigma_meV": 5.0,
        }
        self.kmc = KineticMonteCarlo(self.sim_params)

    def test_initialization(self):
        """Test that amorphous structure initializes correctly"""
        self.assertGreater(len(self.kmc.sites), 0)
        self.assertEqual(self.kmc.h_atoms_on_surface, 0)
        self.assertEqual(self.kmc.adjacent_h_pairs_count, 0)
        
        # Check amorphous structure components
        self.assertTrue(hasattr(self.kmc, 'voronoi_points'))
        self.assertTrue(hasattr(self.kmc, 'site_positions'))
        self.assertTrue(hasattr(self.kmc, 'site_neighbors'))
        self.assertTrue(hasattr(self.kmc, 'coordination_numbers'))
        
        # Check that sites have proper types
        for site_id in self.kmc.sites.keys():
            self.assertIn(self.kmc.site_types[site_id], [1, 2, 3])  # physisorption, chemisorption, defect

    def test_site_coordination(self):
        """Test coordination number calculation"""
        for site_id in self.kmc.sites.keys():
            coord = self.kmc.coordination_numbers[site_id]
            neighbors = self.kmc.site_neighbors.get(site_id, [])
            self.assertEqual(coord, len(neighbors))
            self.assertGreaterEqual(coord, 0)

    def test_coverage_effects(self):
        """Test coverage-dependent binding energies"""
        site_id = list(self.kmc.sites.keys())[0]
        
        # Test with no H neighbors
        base_energy = self.kmc.get_coverage_modified_binding_energy(site_id)
        
        # Add H to neighboring site if possible
        neighbors = self.kmc.site_neighbors.get(site_id, [])
        if neighbors:
            neighbor_id = neighbors[0]
            self.kmc.sites[neighbor_id] = "H"
            self.kmc.h_atoms_on_surface += 1
            
            # Energy should be different due to lateral interactions
            modified_energy = self.kmc.get_coverage_modified_binding_energy(site_id)
            self.assertNotEqual(base_energy, modified_energy)

    def test_adsorption_event(self):
        """Test H adsorption on amorphous surface"""
        initial_h = self.kmc.h_atoms_on_surface
        available_before = len(self.kmc.get_accessible_surface_sites())
        
        if available_before > 0:
            self.kmc.execute_event("adsorption")
            self.assertEqual(self.kmc.h_atoms_on_surface, initial_h + 1)
            
            # Check that a site changed from C to H
            h_sites = [sid for sid, occ in self.kmc.sites.items() if occ == "H"]
            self.assertEqual(len(h_sites), 1)

    def test_desorption_event(self):
        """Test H desorption from amorphous surface"""
        # Add an H atom first
        available_sites = self.kmc.get_accessible_surface_sites()
        if available_sites:
            site_id = available_sites[0]
            self.kmc.sites[site_id] = "H"
            self.kmc.h_atoms_on_surface = 1
            self.kmc._update_adjacent_h_pairs_count()
            
            initial_h = self.kmc.h_atoms_on_surface
            self.kmc.execute_event("desorption")
            self.assertEqual(self.kmc.h_atoms_on_surface, initial_h - 1)

    def test_diffusion_event(self):
        """Test H diffusion between amorphous sites"""
        # Add H atom to a site with empty neighbors
        available_sites = self.kmc.get_accessible_surface_sites()
        mobile_site = None
        
        for site_id in available_sites:
            neighbors = self.kmc.site_neighbors.get(site_id, [])
            empty_neighbors = [nid for nid in neighbors if self.kmc.sites[nid] == "C"]
            if empty_neighbors:
                mobile_site = site_id
                break
        
        if mobile_site:
            self.kmc.sites[mobile_site] = "H"
            self.kmc.h_atoms_on_surface = 1
            self.kmc._update_adjacent_h_pairs_count()
            
            # Record original position
            original_occupied = self.kmc.get_occupied_sites()
            
            # Execute diffusion
            self.kmc.execute_event("diffusion")
            
            # Check that H moved
            new_occupied = self.kmc.get_occupied_sites()
            self.assertEqual(len(new_occupied), 1)  # Still one H atom
            self.assertEqual(self.kmc.h_atoms_on_surface, 1)

    def test_h2_formation_lh(self):
        """Test Langmuir-Hinshelwood H2 formation"""
        # Find two neighboring sites
        for site_id in self.kmc.sites.keys():
            neighbors = self.kmc.site_neighbors.get(site_id, [])
            if neighbors:
                neighbor_id = neighbors[0]
                # Place H atoms on both sites
                self.kmc.sites[site_id] = "H"
                self.kmc.sites[neighbor_id] = "H"
                self.kmc.h_atoms_on_surface = 2
                self.kmc._update_adjacent_h_pairs_count()
                
                if self.kmc.adjacent_h_pairs_count > 0:
                    initial_h2 = self.kmc.h2_molecules_formed_LH
                    self.kmc.execute_event("h2_formation_LH")
                    
                    # Check that H2 formed and H atoms removed
                    self.assertEqual(self.kmc.h2_molecules_formed_LH, initial_h2 + 1)
                    self.assertEqual(self.kmc.h_atoms_on_surface, 0)
                    break

    def test_h2_formation_er(self):
        """Test Eley-Rideal H2 formation"""
        # Add one H atom
        available_sites = self.kmc.get_accessible_surface_sites()
        if available_sites:
            site_id = available_sites[0]
            self.kmc.sites[site_id] = "H"
            self.kmc.h_atoms_on_surface = 1
            
            initial_h2 = self.kmc.h2_molecules_formed_ER
            self.kmc.execute_event("h2_formation_ER")
            
            self.assertEqual(self.kmc.h2_molecules_formed_ER, initial_h2 + 1)
            self.assertEqual(self.kmc.h_atoms_on_surface, 0)

    def test_rate_calculation(self):
        """Test that rates are calculated correctly"""
        rates = self.kmc.calculate_rates()
        
        # Should always have adsorption rate
        self.assertIn("adsorption", rates)
        self.assertGreater(rates["adsorption"], 0)
        
        # Add H atoms and check for other rates
        available_sites = self.kmc.get_accessible_surface_sites()
        if len(available_sites) >= 2:
            # Add some H atoms
            self.kmc.sites[available_sites[0]] = "H"
            self.kmc.sites[available_sites[1]] = "H"
            self.kmc.h_atoms_on_surface = 2
            self.kmc._update_adjacent_h_pairs_count()
            
            rates = self.kmc.calculate_rates()
            
            # Should have desorption and diffusion
            self.assertIn("desorption", rates)
            self.assertIn("diffusion", rates)
            self.assertIn("h2_formation_ER", rates)

    def test_gillespie_algorithm(self):
        """Test basic Gillespie simulation"""
        max_time = 1000.0  # Short simulation
        
        initial_time = self.kmc.time
        self.kmc.run_gillespie(max_time, max_steps=100)
        
        # Time should advance
        self.assertGreaterEqual(self.kmc.time, initial_time)
        
        # Should have some statistics
        stats = self.kmc.get_statistics()
        self.assertIn("time", stats)
        self.assertIn("h_coverage", stats)
        self.assertIn("h2_formed_total", stats)

    def test_uv_pulse_events(self):
        """Test UV pulse start/end events"""
        self.assertFalse(self.kmc.uv_pulse_active)
        
        self.kmc.execute_event("uv_pulse_start")
        self.assertTrue(self.kmc.uv_pulse_active)
        
        self.kmc.execute_event("uv_pulse_end")
        self.assertFalse(self.kmc.uv_pulse_active)

    def test_adjacent_pairs_counting(self):
        """Test accurate counting of adjacent H pairs"""
        # Find neighboring sites
        for site_id in self.kmc.sites.keys():
            neighbors = self.kmc.site_neighbors.get(site_id, [])
            if len(neighbors) >= 2:
                # Place H on original site and two neighbors
                self.kmc.sites[site_id] = "H"
                self.kmc.sites[neighbors[0]] = "H"
                self.kmc.sites[neighbors[1]] = "H"
                self.kmc.h_atoms_on_surface = 3
                self.kmc._update_adjacent_h_pairs_count()
                
                # Should count pairs correctly
                self.assertGreater(self.kmc.adjacent_h_pairs_count, 0)
                self.assertLessEqual(self.kmc.adjacent_h_pairs_count, 3)  # Maximum possible pairs
                break

    def test_statistics_output(self):
        """Test statistics collection"""
        stats = self.kmc.get_statistics()
        
        expected_keys = [
            "time", "h_coverage", "h2_formed_total", 
            "h2_formed_LH", "h2_formed_ER", "h2_formed_UV",
            "h_atoms_surface", "h_pairs_adjacent",
            "total_adsorbed", "total_desorbed"
        ]
        
        for key in expected_keys:
            self.assertIn(key, stats)
            self.assertIsInstance(stats[key], (int, float))


if __name__ == '__main__':
    unittest.main()
