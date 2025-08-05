import unittest
import numpy as np
from kmc_simulation import KineticMonteCarlo
from scientific_data import surface_chemistry_data

class TestQuantitativeValidation(unittest.TestCase):
    def test_efficiency_band_carbonaceous_baseline(self):
        params = {
            "grain_radius_um": 0.1,
            "site_area_angstroms_sq": 9,
            "surface_temperature_k": 15.0,
            "gas_temperature_k": 100.0,
            "h_gas_density_cm3": 300.0,
            "sticking_probability": 0.5,
            "initial_h_coverage": 0.2,
            "uv_flux_factor": 0.0,
            "use_site_heterogeneity": True
        }
        original_binding = surface_chemistry_data["h_physisorption_binding_energy_mev_typical"]
        surface_chemistry_data["h_physisorption_binding_energy_mev_typical"] = 60.0
        try:
            effs = []
            for _ in range(10):
                kmc = KineticMonteCarlo(params)
                kmc.run_gillespie(max_time=1e-3, max_steps=20000)
                ads = max(kmc.total_adsorbed_h_atoms, 1)
                eff = kmc.h2_molecules_formed / ads
                effs.append(eff)
            m = float(np.mean(effs))
            self.assertGreaterEqual(m, 1e-4)
            self.assertLessEqual(m, 5e-1)
        finally:
            surface_chemistry_data["h_physisorption_binding_energy_mev_typical"] = original_binding

if __name__ == "__main__":
    unittest.main()