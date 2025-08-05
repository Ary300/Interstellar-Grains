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
            "uv_flux_factor": 1.0
        }
        self.kmc = KineticMonteCarlo(self.sim_params)

    def test_initialization(self):
        self.assertGreater(self.kmc.total_surface_sites, 0)
        self.assertEqual(self.kmc.h_atoms_on_surface, 0)
        self.assertEqual(self.kmc.adjacent_h_pairs_count, 0)
        self.assertEqual(self.kmc.surface_sites.shape, (self.kmc.surface_dimension, self.kmc.surface_dimension))

    def test_adsorption_event(self):
        initial_h = self.kmc.h_atoms_on_surface
        self.kmc.execute_event("adsorption")
        self.assertEqual(self.kmc.h_atoms_on_surface, initial_h + 1)

    def test_desorption_event(self):
        self.kmc.surface_sites[0, 0] = "H"
        self.kmc.h_atoms_on_surface = 1
        self.kmc.update_adjacent_h_pairs_count(0, 0, True)
        initial_h = self.kmc.h_atoms_on_surface
        self.kmc.execute_event("desorption")
        self.assertEqual(self.kmc.h_atoms_on_surface, initial_h - 1)

    def test_diffusion_event(self):
        self.kmc.surface_sites[0, 0] = "H"
        self.kmc.h_atoms_on_surface = 1
        self.kmc.surface_sites[0, 1] = None
        self.kmc.update_adjacent_h_pairs_count(0, 0, True)
        initial_h = self.kmc.h_atoms_on_surface
        self.kmc.execute_event("diffusion")
        self.assertEqual(self.kmc.h_atoms_on_surface, initial_h)
        self.assertFalse(self.kmc.surface_sites[0,0] == "H" and self.kmc.surface_sites[0,1] == None)

    def test_h2_formation_LH_event(self):
        self.kmc.surface_sites[0, 0] = "H"
        self.kmc.surface_sites[0, 1] = "H"
        self.kmc.h_atoms_on_surface = 2
        self.kmc.update_adjacent_h_pairs_count(0, 0, True)
        self.kmc.update_adjacent_h_pairs_count(0, 1, True)
        initial_h = self.kmc.h_atoms_on_surface
        initial_h2_lh = self.kmc.h2_molecules_formed_LH
        self.kmc.execute_event("h2_formation_LH")
        self.assertEqual(self.kmc.h_atoms_on_surface, initial_h - 2)
        self.assertEqual(self.kmc.h2_molecules_formed_LH, initial_h2_lh + 1)

    def test_h2_formation_ER_event(self):
        self.kmc.surface_sites[0, 0] = "H"
        self.kmc.h_atoms_on_surface = 1
        self.kmc.update_adjacent_h_pairs_count(0, 0, True)
        initial_h = self.kmc.h_atoms_on_surface
        initial_h2_er = self.kmc.h2_molecules_formed_ER
        self.kmc.execute_event("h2_formation_ER")
        self.assertEqual(self.kmc.h_atoms_on_surface, initial_h - 1)
        self.assertEqual(self.kmc.h2_molecules_formed_ER, initial_h2_er + 1)

    def test_h2_formation_UV_event(self):
        self.kmc.surface_sites[0, 0] = "H"
        self.kmc.surface_sites[0, 1] = "H"
        self.kmc.h_atoms_on_surface = 2
        self.kmc.update_adjacent_h_pairs_count(0, 0, True)
        self.kmc.update_adjacent_h_pairs_count(0, 1, True)
        initial_h = self.kmc.h_atoms_on_surface
        initial_h2_uv = self.kmc.h2_molecules_formed_UV
        self.kmc.execute_event("h2_formation_UV")
        self.assertEqual(self.kmc.h_atoms_on_surface, initial_h - 2)
        self.assertEqual(self.kmc.h2_molecules_formed_UV, initial_h2_uv + 1)

    def test_calculate_rates_values(self):
        placed = 0
        coords = [(0,0),(0,1),(2,2)]
        for r,c in coords:
            self.kmc.surface_sites[r, c] = "H"
            self.kmc.h_atoms_on_surface += 1
            self.kmc.update_adjacent_h_pairs_count(r, c, True)
        rates = self.kmc.calculate_rates()
        T_s = self.sim_params["surface_temperature_k"]
        T_g = self.sim_params["gas_temperature_k"]
        n_H = self.sim_params["h_gas_density_cm3"]
        site_area_cm2 = self.sim_params["site_area_angstroms_sq"] * 1e-16
        v_th = np.sqrt(8 * K_B * T_g / (np.pi * M_H))
        flux = 0.25 * n_H * v_th
        sticking = self.sim_params["sticking_probability"]
        expected_adsorption = flux * (self.kmc.total_surface_sites * site_area_cm2) * sticking
        self.assertAlmostEqual(rates["adsorption"], expected_adsorption, delta=expected_adsorption*1e-12)
        nu0 = 1e12
        E_bind_ev = surface_chemistry_data["h_physisorption_binding_energy_mev_typical"]/1000.0
        E_diff_ev = surface_chemistry_data["h_physisorbed_amorphous_diffusion_barrier_ev_min"]
        expected_desorption = nu0 * np.exp(-E_bind_ev * EV_TO_KELVIN / T_s) * self.kmc.h_atoms_on_surface
        expected_diffusion = nu0 * np.exp(-E_diff_ev * EV_TO_KELVIN / T_s) * self.kmc.h_atoms_on_surface
        self.assertAlmostEqual(rates["desorption"], expected_desorption, delta=max(1e-30, expected_desorption*1e-12))
        self.assertAlmostEqual(rates["diffusion"], expected_diffusion, delta=max(1e-30, expected_diffusion*1e-12))
        self.assertIn("h2_formation_LH", rates)
        self.assertGreater(rates["h2_formation_LH"], 0.0)

    def test_run_gillespie_time_limit(self):
        max_time = 1e-5
        self.kmc.run_gillespie(max_time=max_time)
        self.assertLessEqual(self.kmc.time, max_time)

if __name__ == '__main__':
    unittest.main()
