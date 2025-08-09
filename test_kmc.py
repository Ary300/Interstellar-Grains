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
            "use_3d_lattice": True,
            "porosity_fraction": 0.2,
            "E_phys_mean_meV": 45.0,
            "heterogeneity_E_bind_sigma_meV": 5.0,
        }
        self.kmc = KineticMonteCarlo(self.sim_params)

    def test_initialization(self):
        self.assertGreater(self.kmc.surface_dimension, 0)
        self.assertGreater(self.kmc.depth_layers, 0)
        self.assertEqual(self.kmc.h_atoms_on_surface, 0)
        self.assertEqual(self.kmc.adjacent_h_pairs_count, 0)
        self.assertEqual(self.kmc.lattice.shape, (self.kmc.depth_layers, self.kmc.surface_dimension, self.kmc.surface_dimension))
        self.assertTrue(hasattr(self.kmc, 'E_bind_eV_map'))

    def test_adsorption_event(self):
        initial_h = self.kmc.h_atoms_on_surface
        self.kmc.execute_event("adsorption")
        self.assertEqual(self.kmc.h_atoms_on_surface, initial_h + 1)

    def test_desorption_event(self):
        accessible_sites = self.kmc.get_accessible_surface_sites()
        if len(accessible_sites[0]) > 0:
            r, c = accessible_sites[0][0], accessible_sites[1][0]
            self.kmc.lattice[0, r, c] = "H"
            self.kmc.h_atoms_on_surface = 1
            self.kmc._update_adjacent_h_pairs_count()
            initial_h = self.kmc.h_atoms_on_surface
            self.kmc.execute_event("desorption")
            self.assertEqual(self.kmc.h_atoms_on_surface, initial_h - 1)

    def test_diffusion_event(self):
        accessible_sites = self.kmc.get_accessible_surface_sites()
        if len(accessible_sites[0]) > 1:
            r1, c1 = accessible_sites[0][0], accessible_sites[1][0]
            self.kmc.lattice[0, r1, c1] = "H"
            self.kmc.h_atoms_on_surface = 1
            self.kmc._update_adjacent_h_pairs_count()
            initial_h = self.kmc.h_atoms_on_surface
            self.kmc.execute_event("diffusion")
            self.assertEqual(self.kmc.h_atoms_on_surface, initial_h)
            self.assertIn("H", self.kmc.lattice)

    def test_h2_formation_LH_event(self):
        accessible_sites = np.argwhere(self.kmc.lattice[0] != None)
        if len(accessible_sites) > 1:
            r1, c1 = accessible_sites[0]
            r2, c2 = accessible_sites[1]
            
            self.kmc.lattice[0, r1, c1] = "H"
            self.kmc.lattice[0, r2, c2] = "H"
            self.kmc.h_atoms_on_surface = 2
            self.kmc._update_adjacent_h_pairs_count()
            
            if self.kmc.adjacent_h_pairs_count > 0:
                initial_h = self.kmc.h_atoms_on_surface
                initial_h2_lh = self.kmc.h2_molecules_formed_LH
                self.kmc.execute_event("h2_formation_LH")
                self.assertEqual(self.kmc.h_atoms_on_surface, initial_h - 2)
                self.assertEqual(self.kmc.h2_molecules_formed_LH, initial_h2_lh + 1)

    def test_h2_formation_ER_event(self):
        accessible_sites = self.kmc.get_accessible_surface_sites()
        if len(accessible_sites[0]) > 0:
            r, c = accessible_sites[0][0], accessible_sites[1][0]
            self.kmc.lattice[0, r, c] = "H"
            self.kmc.h_atoms_on_surface = 1
            self.kmc._update_adjacent_h_pairs_count()
            initial_h = self.kmc.h_atoms_on_surface
            initial_h2_er = self.kmc.h2_molecules_formed_ER
            self.kmc.execute_event("h2_formation_ER")
            self.assertEqual(self.kmc.h_atoms_on_surface, initial_h - 1)
            self.assertEqual(self.kmc.h2_molecules_formed_ER, initial_h2_er + 1)

    def test_h2_formation_UV_event(self):
        accessible_sites = np.argwhere(self.kmc.lattice[0] != None)
        if len(accessible_sites) > 1:
            r1, c1 = accessible_sites[0]
            r2, c2 = accessible_sites[1]
            
            self.kmc.lattice[0, r1, c1] = "H"
            self.kmc.lattice[0, r2, c2] = "H"
            self.kmc.h_atoms_on_surface = 2
            self.kmc._update_adjacent_h_pairs_count()

            if self.kmc.adjacent_h_pairs_count > 0:
                initial_h = self.kmc.h_atoms_on_surface
                initial_h2_uv = self.kmc.h2_molecules_formed_UV
                self.kmc.execute_event("h2_formation_UV")
                self.assertEqual(self.kmc.h_atoms_on_surface, initial_h - 2)
                self.assertEqual(self.kmc.h2_molecules_formed_UV, initial_h2_uv + 1)

    def test_calculate_rates_values(self):
        self.kmc.h_atoms_on_surface = 3
        self.kmc.adjacent_h_pairs_count = 1
        
        rates = self.kmc.calculate_rates()
        T_g = self.sim_params["gas_temperature_k"]
        n_H = self.sim_params["h_gas_density_cm3"]
        site_area_cm2 = self.sim_params["site_area_angstroms_sq"] * 1e-16
        
        accessible_sites = self.kmc.get_accessible_surface_sites()
        num_accessible_sites = len(accessible_sites[0])
        accessible_area_cm2 = num_accessible_sites * site_area_cm2
        v_th = np.sqrt(8 * K_B * T_g / (np.pi * M_H))
        flux = 0.25 * n_H * v_th
        sticking = self.sim_params["sticking_probability"]
        expected_adsorption = flux * accessible_area_cm2 * sticking
        self.assertAlmostEqual(rates["adsorption"], expected_adsorption, delta=expected_adsorption*1e-9)
        
        if self.kmc.h_atoms_on_surface > 0:
            self.assertIn("desorption", rates)
            self.assertIn("diffusion", rates)
            if self.kmc.adjacent_h_pairs_count > 0:
                self.assertIn("h2_formation_LH", rates)

    def test_run_gillespie_time_limit(self):
        max_time = 1e-5
        self.kmc.run_gillespie(max_time=max_time)
        self.assertLessEqual(self.kmc.time, max_time)

    def test_3d_structure_properties(self):
        void_sites = np.sum(self.kmc.lattice == None)
        total_sites = self.kmc.lattice.size
        actual_porosity = void_sites / total_sites
        self.assertAlmostEqual(actual_porosity, self.sim_params["porosity_fraction"], delta=0.05)
        
        self.assertTrue(np.all(self.kmc.E_bind_eV_map[self.kmc.lattice == None] == 0))
        self.assertTrue(np.all(self.kmc.E_diff_eV_map[self.kmc.lattice == None] == 0))
        self.assertTrue(np.all(self.kmc.E_diff_eV_map >= 0))
        self.assertTrue(np.all(self.kmc.E_bind_eV_map >= 0))

if __name__ == '__main__':
    unittest.main()
