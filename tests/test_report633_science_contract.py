"""Report 633 scientific input-contract regressions for the pinned kernel."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scrambox import case, grid, server
from scrambox.compile import compile_case


class ClassPartitionContractTest(unittest.TestCase):
    def test_full_singleton_and_one_group_defaults_remain_valid(self):
        full = grid.composition_classes(3, {"preset": "explicit", "classes": [[[0, 1], [0, 1]]]})
        self.assertEqual(full, [[(0.0, 1.0)] * 3])
        self.assertEqual(grid.composition_classes(1, {"preset": "NC6"}), [[(0.0, 1.0)]])
        self.assertEqual(len(grid.composition_classes(3, {"preset": "NC6"})), 3)

    def test_partial_singleton_cannot_leave_compositions_unclassified(self):
        with self.assertRaisesRegex(ValueError, "not closed|no class"):
            grid.composition_classes(3, {"preset": "explicit", "classes": [[[0, .5], [0, 1]]]})

    def test_simplex_sliver_is_checked_even_when_rectangle_midpoint_is_outside(self):
        with self.assertRaisesRegex(ValueError, "not closed|no class"):
            grid.composition_classes(3, {"preset": "explicit", "classes": [[[0, .7], [0, .7]]]})

    def test_singleton_invalid_ranges_are_rejected(self):
        for bad in ([-.1, 1], [0, 1.1], [.7, .2], [0, 0]):
            with self.subTest(bad=bad), self.assertRaisesRegex(ValueError, "0 <= lo < hi <= 1"):
                grid.composition_classes(3, {"preset": "explicit", "classes": [[bad, [0, 1]]]})

    def test_grid_and_start_run_reject_gap_before_creating_job(self):
        invalid = {"mixing": {"preset": "explicit", "classes": [[[0, .5]] + [[0, 1]] * 4]}}
        with self.assertRaises(ValueError):
            server.grid_info(invalid)
        with tempfile.TemporaryDirectory() as td:
            ws = server.Workspace(Path(td))
            with patch.object(server.POOL, "submit", side_effect=AssertionError("worker submitted")):
                with self.assertRaises(ValueError):
                    server.start_run(ws, invalid)
            self.assertEqual(list(ws.runs.iterdir()), [])


class SplitKernelContractTest(unittest.TestCase):
    invalid = {"processes": {"class_redistribution": "split"}}

    def test_split_is_rejected_with_pinned_kernel_guidance(self):
        with self.assertRaisesRegex(ValueError, "split.*unsupported.*kernel.*whole"):
            case.resolve(self.invalid)
        self.assertEqual(case.resolve({})["processes"]["class_redistribution"], "whole")

    def test_compile_rejects_split_before_writing_inputs(self):
        with tempfile.TemporaryDirectory() as td:
            target = Path(td) / "split-run"
            with self.assertRaisesRegex(ValueError, "split.*unsupported.*kernel.*whole"):
                compile_case(self.invalid, target)
            self.assertFalse(target.exists())

    def test_start_run_rejects_split_without_job_or_run_data(self):
        with tempfile.TemporaryDirectory() as td:
            ws = server.Workspace(Path(td))
            with patch.object(server.POOL, "submit", side_effect=AssertionError("worker submitted")):
                with self.assertRaisesRegex(ValueError, "split.*unsupported.*kernel.*whole"):
                    server.start_run(ws, self.invalid)
            self.assertEqual(list(ws.runs.iterdir()), [])

    def test_whole_compilation_uses_existing_namelist_only(self):
        with tempfile.TemporaryDirectory() as td:
            run = Path(td) / "whole-run"
            compile_case({}, run)
            nml = (run / "run.nml").read_text()
            self.assertIn("&box_proc", nml)
            self.assertNotIn("class_split", nml)
        self.assertEqual(server.meta()["defaults"]["processes"]["class_redistribution"], "whole")
        self.assertNotIn("class_split_alpha", server.meta()["defaults"]["processes"])


if __name__ == "__main__":
    unittest.main()
