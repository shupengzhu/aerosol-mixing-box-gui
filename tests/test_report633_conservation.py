"""Conservation metadata distinguishes a closed box from an open scenario."""
import json
from pathlib import Path
import unittest

import numpy as np

from scrambox import case as casemod
from scrambox.run import conservation, conservation_diagnostic


ROOT = Path(__file__).resolve().parents[1]


class ConservationDiagnosticTest(unittest.TestCase):
    def test_exact_first_last_definition(self):
        grid = {"species": {"aerosol": ["PSO4"], "gases": ["SULF"],
                            "gas_of": {"PSO4": "SULF"}}}
        result = {"mass_ugm3": np.array([[[[1.0]]], [[[2.0]]]]),
                  "gas_ugm3": np.array([[2.0], [2.0]])}
        self.assertAlmostEqual(conservation(result, grid)["PSO4+SULF"], 1.0 / 3.0)

    def test_source_free_default_is_eligible_without_pass_threshold(self):
        diagnostic = conservation_diagnostic(casemod.resolve({}))
        self.assertTrue(diagnostic["applicable"])
        self.assertEqual(diagnostic["reasons"], [])
        self.assertIsNone(diagnostic["pass_threshold"])
        self.assertIn("first", diagnostic["definition"])
        self.assertIn("gas", diagnostic["definition"])

    def test_urban_plume_open_system_is_not_conservation_error(self):
        scenario = json.loads((ROOT / "scenarios/urban_plume_aging.json").read_text())
        diagnostic = conservation_diagnostic(casemod.resolve(scenario["case"]))
        self.assertFalse(diagnostic["applicable"])
        self.assertIn("particle_emissions", diagnostic["reasons"])
        self.assertIn("gas_sources", diagnostic["reasons"])
        self.assertIn("dilution", diagnostic["reasons"])

    def test_time_varying_dilution_is_detected(self):
        case = casemod.resolve({"environment": {"time_series": {
            "time_s": [0, 60], "dilution_s1": [0, 1e-5]}}})
        diagnostic = conservation_diagnostic(case)
        self.assertFalse(diagnostic["applicable"])
        self.assertIn("dilution", diagnostic["reasons"])


if __name__ == "__main__":
    unittest.main()
