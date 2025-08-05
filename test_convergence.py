import unittest
import numpy as np
from kmc_simulation import KineticMonteCarlo

TOL = 0.10

def run_many(params, max_time=1e-3, steps=20000, n=10):
    vals = []
    for _ in range(n):
        kmc = KineticMonteCarlo(params.copy())
        kmc.run_gillespie(max_time=max_time, max_steps=steps)
        ads = max(kmc.total_adsorbed_h_atoms, 1)
        eff = kmc.h2_molecules_formed / ads
        vals.append(eff)
    return float(np.mean(vals))

class TestConvergence(unittest.TestCase):
    def setUp(self):
        self.base = {
            "grain_radius_um": 0.05,
            "site_area_angstroms_sq": 9,
            "surface_temperature_k": 15.0,
            "gas_temperature_k": 100.0,
            "h_gas_density_cm3": 300.0,
            "sticking_probability": 0.5,
            "initial_h_coverage": 0.1,
            "uv_flux_factor": 0.0,
            "use_site_heterogeneity": True,
        }

    def test_lattice_size_convergence(self):
        p1 = self.base.copy()
        p2 = self.base.copy()
        p3 = self.base.copy()
        p1["site_area_angstroms_sq"] = 9
        p2["site_area_angstroms_sq"] = 6
        p3["site_area_angstroms_sq"] = 12
        m1 = run_many(p1, max_time=1e-3, steps=15000, n=8)
        m2 = run_many(p2, max_time=1e-3, steps=15000, n=8)
        m3 = run_many(p3, max_time=1e-3, steps=15000, n=8)
        ref = m1
        self.assertLessEqual(abs(m2 - ref) / max(ref, 1e-12), TOL)
        self.assertLessEqual(abs(m3 - ref) / max(ref, 1e-12), TOL)

    def test_time_horizon_convergence(self):
        p = self.base.copy()
        m_short = run_many(p, max_time=5e-4, steps=8000, n=8)
        m_long = run_many(p, max_time=2e-3, steps=30000, n=8)
        ref = m_long
        self.assertLessEqual(abs(m_short - ref) / max(ref, 1e-12), TOL)

    def test_ensemble_runs_convergence(self):
        p = self.base.copy()
        m_few = run_many(p, max_time=1e-3, steps=15000, n=5)
        m_more = run_many(p, max_time=1e-3, steps=15000, n=20)
        ref = m_more
        self.assertLessEqual(abs(m_few - ref) / max(ref, 1e-12), TOL)

if __name__ == "__main__":
    unittest.main()