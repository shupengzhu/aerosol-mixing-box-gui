"""PartMC rejects unsupported particle input before a run is queued."""
import copy
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

from scrambox.case import DEFAULTS
from scrambox.server import Workspace, partmc_options, start_run


class PartMCInputPreflightTest(unittest.TestCase):
    def test_supported_lognormal_case_passes(self):
        with patch("scrambox.partmc.available", return_value={"v23": True, "v1": True}):
            self.assertEqual(partmc_options({}, {"n_part": 200})["n_part"], 200)

    def test_unsupported_inputs_fail_before_job_creation(self):
        with patch("scrambox.partmc.available", return_value={"v23": True, "v1": True}):
            for section in ("initial", "emissions"):
                with self.subTest(section=section):
                    case = copy.deepcopy(DEFAULTS)
                    case[section]["particles"] = {"mode": "bin_mass", "mass_ugm3": []}
                    with self.assertRaisesRegex(ValueError, "PartMC.*bin_mass.*lognormal"):
                        partmc_options(case, None)

            case = copy.deepcopy(DEFAULTS)
            case["environment"]["dilution"]["background"] = {
                "particles": {"mode": "bin_mass", "mass_ugm3": []}, "gas_ugm3": {}}
            with self.assertRaisesRegex(ValueError, "PartMC.*bin_mass.*lognormal"):
                partmc_options(case, None)

    def test_unsupported_run_does_not_create_or_submit_job(self):
        case = copy.deepcopy(DEFAULTS)
        case["initial"]["particles"] = {"mode": "bin_mass", "mass_ugm3": []}
        with tempfile.TemporaryDirectory() as tmp, \
             patch("scrambox.partmc.available", return_value={"v23": True, "v1": True}), \
             patch("scrambox.server.POOL.submit") as submit:
            ws = Workspace(Path(tmp))
            with self.assertRaisesRegex(ValueError, "PartMC.*bin_mass.*lognormal"):
                start_run(ws, case, model="partmc")
            self.assertEqual(list(ws.runs.iterdir()), [])
            submit.assert_not_called()


if __name__ == "__main__":
    unittest.main()
