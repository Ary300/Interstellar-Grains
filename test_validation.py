import unittest
import numpy as np
from kmc_simulation import KineticMonteCarlo
from scientific_data import surface_chemistry_data

class TestKMCValidation(unittest.TestCase):
    def setUp(self):
        self.benchmark_params = {
            "grain_radius_um": 0.1,
            "site_area_angstroms_sq": 9,
            "surface_temperature_k": 10.0,
            "gas_temperature_k": 100.0,
            "h_gas_density_cm3": 100.0,
            "sticking_probability": 1.0,
            "initial_h_coverage": 0.5,
            "uv_flux_factor": 0.0
        }
        self.temperatures = [10.0, 25.118864315095795, 50.11872336272722, 100.0, 316.22776601683796]

    def test_h2_efficiency_trend(self):
        efficiencies = []
        original_binding_energy = surface_chemistry_data["h_physisorption_binding_energy_mev_typical"]
        surface_chemistry_data["h_physisorption_binding_energy_mev_typical"] = 85.0

        try:
            for T in self.temperatures:
                params = self.benchmark_params.copy()
                params["surface_temperature_k"] = T
                kmc = KineticMonteCarlo(params)
                kmc.run_gillespie(max_time=1e-3, max_steps=20000)
                total_ads = max(kmc.total_adsorbed_h_atoms, 1)
                eff = kmc.h2_molecules_formed / total_ads
                efficiencies.append((T, eff))
        finally:
            surface_chemistry_data["h_physisorption_binding_energy_mev_typical"] = original_binding_energy

        lowT_eff = [e for (_, e) in efficiencies[:3]]
        self.assertGreaterEqual(lowT_eff[0], lowT_eff[1] - 1e-12, "Efficiency should not increase from 10K to ~25K")
        self.assertGreaterEqual(lowT_eff[1], lowT_eff[2] - 1e-12, "Efficiency should not increase from ~25K to ~50K")
        self.assertLess(efficiencies[-1][1], 1e-2, "Efficiency at ~300 K should be very low")

if __name__ == '__main__':
    unittest.main()


