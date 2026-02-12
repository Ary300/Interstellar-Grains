import unittest
import numpy as np
import tempfile
from pathlib import Path
from kmc_simulation import KineticMonteCarlo
from scientific_data import K_B_ERG, M_H

class TestKineticMonteCarlo(unittest.TestCase):

    def setUp(self):
        self.sim_params = {
            # Keep the lattice small so unit tests run quickly.
            "rng_seed": 123,
            "grain_radius_um": 0.01,
            "site_area_angstroms_sq": 25,
            "surface_temperature_k": 10.0,
            "gas_temperature_k": 100.0,
            "h_gas_density_cm3": 1e2,
            "sticking_probability": 0.3,
            "initial_h_coverage": 0.0,
            "uv_flux_factor": 0.0,
            "uv_pulse_enabled": False,
            "use_3d_lattice": True,
            "porosity_fraction": 0.2,
            "E_phys_mean_meV": 45.0,
            "heterogeneity_E_bind_sigma_meV": 5.0,
        }
        self.kmc = KineticMonteCarlo(self.sim_params)

    def _place_h(self, d: int, r: int, c: int):
        """Place an H atom while keeping KMC bookkeeping consistent."""
        self.assertNotEqual(self.kmc.lattice[d, r, c], None)
        self.assertNotEqual(self.kmc.lattice[d, r, c], "H")
        self.kmc.lattice[d, r, c] = "H"
        self.kmc.h_atoms_on_surface += 1
        self.kmc.occupied.add((d, r, c))
        if d == 0:
            self.kmc.empty_surface.discard((r, c))
        self.kmc.update_adjacent_h_pairs_count(d, r, c, True)

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
            r, c = int(accessible_sites[0][0]), int(accessible_sites[1][0])
            self._place_h(0, r, c)
            initial_h = self.kmc.h_atoms_on_surface
            self.kmc.execute_event("desorption")
            self.assertEqual(self.kmc.h_atoms_on_surface, initial_h - 1)

    def test_diffusion_event(self):
        accessible_sites = self.kmc.get_accessible_surface_sites()
        if len(accessible_sites[0]) > 1:
            r1, c1 = int(accessible_sites[0][0]), int(accessible_sites[1][0])
            self._place_h(0, r1, c1)
            initial_h = self.kmc.h_atoms_on_surface
            self.kmc.execute_event("diffusion")
            self.assertEqual(self.kmc.h_atoms_on_surface, initial_h)
            self.assertIn("H", self.kmc.lattice)
            self.assertEqual(len(self.kmc.occupied), 1)

    def test_h2_formation_LH_event(self):
        # Find any accessible top-layer site and an accessible neighbor.
        accessible_sites = np.argwhere(self.kmc.lattice[0] != None)
        for r1, c1 in accessible_sites:
            neighbors = [(r1 - 1, c1), (r1 + 1, c1), (r1, c1 - 1), (r1, c1 + 1)]
            for rr, cc in neighbors:
                r2, c2 = int(rr) % self.kmc.surface_dimension, int(cc) % self.kmc.surface_dimension
                if self.kmc.lattice[0, r2, c2] is None:
                    continue
                self._place_h(0, int(r1), int(c1))
                self._place_h(0, r2, c2)
                self.assertGreater(self.kmc.adjacent_h_pairs_count, 0)
                initial_h = self.kmc.h_atoms_on_surface
                initial_h2_lh = self.kmc.h2_molecules_formed_LH
                self.kmc.execute_event("h2_formation_LH")
                self.assertEqual(self.kmc.h_atoms_on_surface, initial_h - 2)
                self.assertEqual(self.kmc.h2_molecules_formed_LH, initial_h2_lh + 1)
                return

    def test_h2_formation_ER_event(self):
        accessible_sites = self.kmc.get_accessible_surface_sites()
        if len(accessible_sites[0]) > 0:
            r, c = int(accessible_sites[0][0]), int(accessible_sites[1][0])
            self._place_h(0, r, c)
            initial_h = self.kmc.h_atoms_on_surface
            initial_h2_er = self.kmc.h2_molecules_formed_ER
            self.kmc.execute_event("h2_formation_ER")
            self.assertEqual(self.kmc.h_atoms_on_surface, initial_h - 1)
            self.assertEqual(self.kmc.h2_molecules_formed_ER, initial_h2_er + 1)

    def test_h2_formation_UV_event(self):
        accessible_sites = np.argwhere(self.kmc.lattice[0] != None)
        for r1, c1 in accessible_sites:
            neighbors = [(r1 - 1, c1), (r1 + 1, c1), (r1, c1 - 1), (r1, c1 + 1)]
            for rr, cc in neighbors:
                r2, c2 = int(rr) % self.kmc.surface_dimension, int(cc) % self.kmc.surface_dimension
                if self.kmc.lattice[0, r2, c2] is None:
                    continue
                self._place_h(0, int(r1), int(c1))
                self._place_h(0, r2, c2)
                self.assertGreater(self.kmc.adjacent_h_pairs_count, 0)
                initial_h = self.kmc.h_atoms_on_surface
                initial_h2_uv = self.kmc.h2_molecules_formed_UV
                self.kmc.execute_event("h2_formation_UV")
                self.assertEqual(self.kmc.h_atoms_on_surface, initial_h - 2)
                self.assertEqual(self.kmc.h2_molecules_formed_UV, initial_h2_uv + 1)
                return

    def test_calculate_rates_values(self):
        accessible_sites = self.kmc.get_accessible_surface_sites()
        if len(accessible_sites[0]) < 3:
            self.skipTest("Not enough accessible sites on surface")

        # Place 2 adjacent H atoms + 1 isolated H atom.
        r1, c1 = int(accessible_sites[0][0]), int(accessible_sites[1][0])
        self._place_h(0, r1, c1)
        r2, c2 = (r1 + 1) % self.kmc.surface_dimension, c1
        if self.kmc.lattice[0, r2, c2] is None:
            r2, c2 = r1, (c1 + 1) % self.kmc.surface_dimension
        if self.kmc.lattice[0, r2, c2] is None:
            self.skipTest("Could not find adjacent accessible site")
        self._place_h(0, r2, c2)

        # Third atom somewhere else.
        r3, c3 = int(accessible_sites[0][2]), int(accessible_sites[1][2])
        if (r3, c3) in [(r1, c1), (r2, c2)]:
            r3, c3 = int(accessible_sites[0][-1]), int(accessible_sites[1][-1])
        if self.kmc.lattice[0, r3, c3] is None or (r3, c3) in [(r1, c1), (r2, c2)]:
            self.skipTest("Could not find a distinct third accessible site")
        self._place_h(0, r3, c3)

        self.assertEqual(self.kmc.h_atoms_on_surface, 3)
        self.assertGreater(self.kmc.adjacent_h_pairs_count, 0)
        
        rates = self.kmc.calculate_rates()
        T_g = self.sim_params["gas_temperature_k"]
        n_H = self.sim_params["h_gas_density_cm3"]
        site_area_cm2 = self.sim_params["site_area_angstroms_sq"] * 1e-16
        
        num_accessible_sites = self.kmc.get_num_accessible_surface_sites()
        accessible_area_cm2 = num_accessible_sites * site_area_cm2

        v_th = np.sqrt(8 * K_B_ERG * T_g / (np.pi * M_H))
        flux = 0.25 * n_H * v_th
        sticking = self.sim_params["sticking_probability"]
        expected_adsorption = flux * accessible_area_cm2 * sticking * np.exp(-T_g / 100.0)
        self.assertAlmostEqual(rates["adsorption"], expected_adsorption, delta=expected_adsorption*1e-9)
        
        if self.kmc.h_atoms_on_surface > 0:
            self.assertIn("h2_formation_ER", rates)
            self.assertIn("desorption", rates)
            self.assertIn("diffusion", rates)
            if self.kmc.adjacent_h_pairs_count > 0:
                self.assertIn("h2_formation_LH", rates)

    def test_arrival_rate_per_site_values(self):
        params = dict(self.sim_params)
        params["arrival_rate_s"] = None
        params["arrival_rate_per_site_s"] = 0.01
        kmc = KineticMonteCarlo(params)
        rates = kmc.calculate_rates()
        self.assertIn("arrival", rates)
        expected = 0.01 * kmc.total_accessible_surface_sites
        self.assertAlmostEqual(rates["arrival"], expected, delta=expected * 1e-12)
        self.assertNotIn("adsorption", rates)

    def test_h2_blocking_adsorption_and_desorption(self):
        params = dict(self.sim_params)
        params.update(
            {
                "enable_h2_blocking": True,
                "h2_stick_transition_K": 100.0,
                "h2_stick_prob_lowT": 1.0,
                "uv_flux_factor": 0.0,
            }
        )
        self.kmc = KineticMonteCarlo(params)
        kmc = self.kmc

        # Place 2 adjacent H atoms (guarantee a single pair exists).
        accessible_sites = kmc.get_accessible_surface_sites()
        if len(accessible_sites[0]) < 2:
            self.skipTest("Not enough accessible sites on surface")
        r1, c1 = int(accessible_sites[0][0]), int(accessible_sites[1][0])
        r2, c2 = (r1 + 1) % kmc.surface_dimension, c1
        if kmc.lattice[0, r2, c2] is None:
            r2, c2 = r1, (c1 + 1) % kmc.surface_dimension
        if kmc.lattice[0, r2, c2] is None:
            self.skipTest("Could not find adjacent accessible site")

        self._place_h(0, r1, c1)
        self._place_h(0, r2, c2)
        self.assertGreater(kmc.adjacent_h_pairs_count, 0)

        kmc.execute_event("h2_formation_LH")
        self.assertEqual(kmc.h2_molecules_formed, 1)
        self.assertEqual(kmc.h2_molecules_on_surface, 1)
        self.assertEqual(kmc.h2_molecules_desorbed, 0)
        self.assertIn("H2", kmc.lattice)

        kmc.execute_event("h2_desorption")
        self.assertEqual(kmc.h2_molecules_on_surface, 0)
        self.assertEqual(kmc.h2_molecules_desorbed, 0)
        self.assertEqual(getattr(kmc, "h2_molecules_released_formed", 0), 1)
        self.assertEqual(getattr(kmc, "h2_molecules_desorbed_beam", 0), 0)

    def test_temperature_ramp_updates_surface_temperature(self):
        params = dict(self.sim_params)
        params.update(
            {
                "arrival_rate_s": 10.0,
                "h_gas_density_cm3": 0.0,
                "temp_ramp": {"enabled": True, "T_start_K": 10.0, "T_end_K": 20.0, "rate_K_per_s": 1.0},
            }
        )
        kmc = KineticMonteCarlo(params)
        kmc.run_gillespie(max_time=10.0, max_steps=5000)
        self.assertGreaterEqual(float(kmc.simulation_parameters.get("surface_temperature_k", 0.0)), 19.0)

    def test_run_gillespie_time_limit(self):
        max_time = 1e-5
        # Guard against pathological parameterizations by providing a step cap.
        self.kmc.run_gillespie(max_time=max_time, max_steps=100000)
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

    def test_grain_cache_reuse_structure(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            params = dict(self.sim_params)
            params.update(
                {
                    "enable_grain_cache": True,
                    "grain_cache_dir": tmpdir,
                    "rng_seed": 123,
                }
            )
            kmc1 = KineticMonteCarlo(dict(params))
            cache_files = list(Path(tmpdir).glob("grain_*.pkl"))
            self.assertGreater(len(cache_files), 0)

            # Different RNG seed but same structural parameters should reuse cached grain.
            params["rng_seed"] = 999
            kmc2 = KineticMonteCarlo(dict(params))
            self.assertEqual(kmc1.lattice.shape, kmc2.lattice.shape)
            self.assertTrue(np.array_equal(kmc1.lattice, kmc2.lattice))
            self.assertTrue(np.array_equal(kmc1.site_types, kmc2.site_types))
            self.assertTrue(np.array_equal(kmc1.E_bind_eV_map, kmc2.E_bind_eV_map))
            self.assertTrue(np.array_equal(kmc1.E_diff_eV_map, kmc2.E_diff_eV_map))

if __name__ == '__main__':
    unittest.main()
