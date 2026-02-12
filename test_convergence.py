import unittest
import numpy as np
from kmc_simulation import KineticMonteCarlo

def run_many(params, max_time=1e-3, steps=1500, n=4):
    vals = []
    for i in range(n):
        p = params.copy()
        p["rng_seed"] = 1337 + i
        kmc = KineticMonteCarlo(p)
        initial_h = kmc.h_atoms_on_surface
        kmc.run_gillespie(max_time=max_time, max_steps=steps)
        total_h_in = initial_h + kmc.total_adsorbed_h_atoms
        # Fraction of H atoms converted to H2 (2 H atoms per H2).
        eff = (2.0 * kmc.h2_molecules_formed) / max(float(total_h_in), 1.0)
        vals.append(eff)
    return float(np.mean(vals))

class TestConvergence(unittest.TestCase):
    def setUp(self):
        self.base = {
            # Keep lattice sizes modest so convergence tests finish quickly.
            "grain_radius_um": 0.01,
            "site_area_angstroms_sq": 9,
            "surface_temperature_k": 15.0,
            "gas_temperature_k": 100.0,
            "h_gas_density_cm3": 300.0,
            "sticking_probability": 0.5,
            "initial_h_coverage": 0.02,
            "uv_flux_factor": 0.0,
            "use_3d_lattice": True,
            "porosity_fraction": 0.2,
            "E_phys_mean_meV": 45.0,
            "heterogeneity_E_bind_sigma_meV": 5.0,
        }

    def test_lattice_size_convergence(self):
        p1 = self.base.copy()
        p2 = self.base.copy()
        p3 = self.base.copy()
        p1["site_area_angstroms_sq"] = 9
        p2["site_area_angstroms_sq"] = 6
        p3["site_area_angstroms_sq"] = 12
        m1 = run_many(p1, max_time=1e-3, steps=1500, n=4)
        m2 = run_many(p2, max_time=1e-3, steps=1500, n=4)
        m3 = run_many(p3, max_time=1e-3, steps=1500, n=4)
        for m in (m1, m2, m3):
            self.assertTrue(np.isfinite(m))
            self.assertGreaterEqual(m, 0.0)
            self.assertLessEqual(m, 1.0)
        self.assertLessEqual(max(m1, m2, m3) - min(m1, m2, m3), 0.6)

    def test_time_horizon_convergence(self):
        p = self.base.copy()
        m_short = run_many(p, max_time=5e-4, steps=800, n=4)
        m_long = run_many(p, max_time=2e-3, steps=3000, n=4)
        self.assertTrue(np.isfinite(m_short))
        self.assertTrue(np.isfinite(m_long))
        self.assertGreaterEqual(m_short, 0.0)
        self.assertGreaterEqual(m_long, 0.0)
        self.assertLessEqual(m_short, 1.0)
        self.assertLessEqual(m_long, 1.0)
        # Longer horizons should not reduce conversion (allow tiny numerical slack).
        self.assertGreaterEqual(m_long + 1e-12, m_short)

    def test_ensemble_runs_convergence(self):
        p = self.base.copy()
        m_few = run_many(p, max_time=1e-3, steps=1500, n=3)
        m_more = run_many(p, max_time=1e-3, steps=1500, n=8)
        self.assertTrue(np.isfinite(m_few))
        self.assertTrue(np.isfinite(m_more))
        self.assertGreaterEqual(m_few, 0.0)
        self.assertGreaterEqual(m_more, 0.0)
        self.assertLessEqual(m_few, 1.0)
        self.assertLessEqual(m_more, 1.0)
        self.assertLessEqual(abs(m_few - m_more), 0.6)

    def test_3d_structure_convergence(self):
        p1 = self.base.copy()
        p2 = self.base.copy()
        p3 = self.base.copy()
        
        p1["porosity_fraction"] = 0.1
        p2["porosity_fraction"] = 0.2
        p3["porosity_fraction"] = 0.3
        
        m1 = run_many(p1, max_time=1e-3, steps=1500, n=3)
        m2 = run_many(p2, max_time=1e-3, steps=1500, n=3)
        m3 = run_many(p3, max_time=1e-3, steps=1500, n=3)
        
        for m in (m1, m2, m3):
            self.assertTrue(np.isfinite(m))
            self.assertGreaterEqual(m, 0.0)
            self.assertLessEqual(m, 1.0)
        self.assertLessEqual(max(m1, m2, m3) - min(m1, m2, m3), 0.7)

    def test_energy_model_convergence(self):
        p1 = self.base.copy()
        p2 = self.base.copy()
        
        p1["heterogeneity_E_bind_sigma_meV"] = 2.0
        p2["heterogeneity_E_bind_sigma_meV"] = 10.0
        
        m1 = run_many(p1, max_time=1e-3, steps=1500, n=3)
        m2 = run_many(p2, max_time=1e-3, steps=1500, n=3)
        
        for m in (m1, m2):
            self.assertTrue(np.isfinite(m))
            self.assertGreaterEqual(m, 0.0)
            self.assertLessEqual(m, 1.0)
        self.assertLessEqual(abs(m2 - m1), 0.7)

if __name__ == "__main__":
    unittest.main()
