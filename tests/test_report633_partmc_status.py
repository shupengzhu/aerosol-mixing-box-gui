"""PartMC timeout cleanup and observable progress without a fabricated percentage."""
import json
import os
from pathlib import Path
import signal
import tempfile
import time
import unittest
from unittest.mock import patch

from scrambox import partmc
from scrambox.server import Workspace, _status, list_runs


def _live(pid: int) -> bool:
    stat = Path(f"/proc/{pid}/stat")
    return stat.exists() and stat.read_text().split()[2] not in {"Z", "X"}


class PartMCStatusTest(unittest.TestCase):
    def test_sparse_multi_gigabyte_log_tail_is_bounded(self):
        from scrambox.partmc import _log_tail
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "huge.log"
            with path.open("wb") as log:
                log.seek(2 * 1024 * 1024 * 1024)
                log.write(b"TAIL_ERROR_MARKER")
            self.assertIn("TAIL_ERROR_MARKER", _log_tail(path, 64))

    def test_large_debug_stdout_is_filtered_and_log_is_bounded(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            fake = directory / "fake_partmc"
            fake.write_text("#!/usr/bin/env python3\n"
                            "import sys\n"
                            "row=b'1.000000E+00 2.000000E+00 3.000000E+00 4.000000E+00 5.000000E+00\\n'\n"
                            "sys.stdout.buffer.write(row*120000)\n"
                            "sys.stdout.buffer.write((b'non-debug-noise'*4000+b'\\n')*110)\n"
                            "sys.stdout.buffer.write(b'TAIL_ERROR_MARKER\\n')\n"
                            "sys.stdout.buffer.flush()\n"
                            "sys.exit(7)\n")
            fake.chmod(0o755)
            with patch("scrambox.partmc.write_partmc", return_value={"kernel": {"isorropia": "v23"}}), \
                 patch.dict(partmc.PARTMC_BIN, {"v23": str(fake)}):
                with self.assertRaisesRegex(RuntimeError, "TAIL_ERROR_MARKER"):
                    partmc.run_partmc({}, directory, timeout_s=10)
            log = directory / "pmc.log"
            self.assertLess(log.stat().st_size, 5 * 1024 * 1024)
            self.assertIn("suppressed", log.read_text())

    def test_normal_completion_keeps_progress_log_and_runs_projection(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            fake = directory / "fake_partmc"
            fake.write_text("#!/usr/bin/env python3\nprint('t_sim 120s done', flush=True)\n")
            fake.chmod(0o755)
            with patch("scrambox.partmc.write_partmc", return_value={"kernel": {"isorropia": "v23"}}), \
                 patch.dict(partmc.PARTMC_BIN, {"v23": str(fake)}), \
                 patch("scrambox.partmc.project", return_value={"summary": {"n_output": 3}}) as project, \
                 patch("scrambox.partmc.truth") as truth:
                result = partmc.run_partmc({}, directory, n_part=200, timeout_s=5)
            project.assert_called_once_with(directory)
            truth.assert_called_once_with(directory)
            self.assertEqual(result["summary"]["n_part"], 200)
            self.assertIn("t_sim 120s done", (directory / "pmc.log").read_text())

    def test_progress_has_output_and_wall_time_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            ws = Workspace(Path(tmp))
            run = ws.runs / "synthetic"
            out = run / "out"
            out.mkdir(parents=True)
            (run / "case.resolved.json").write_text(json.dumps({"time": {
                "duration_s": 120, "output_every_s": 60}}))
            (run / "pmc.log").write_text("synthetic progress\n")
            start = time.time() - 15
            (run / "status.json").write_text(json.dumps({"state": "running", "model": "partmc",
                                                        "started": start}))
            for index in (1, 2):
                file = out / f"pmc_0001_{index:08d}.nc"
                file.write_bytes(b"synthetic")
                os.utime(file, (time.time() - 5, time.time() - 5))
            status = _status(run)
            self.assertAlmostEqual(status["progress"], 2 / 3)
            self.assertEqual(status["output_count"], 2)
            self.assertEqual(status["expected_output_count"], 3)
            self.assertGreaterEqual(status["elapsed_s"], 15)
            self.assertTrue(4 <= status["last_output_age_s"] <= 8)
            self.assertGreater(status["log_bytes"], 0)
            listed = list_runs(ws)[0]
            self.assertEqual(listed["output_count"], 2)
            self.assertIn("last_output_age_s", listed)

    def test_timeout_terminates_entire_owned_process_group(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            fake = directory / "fake_partmc"
            fake.write_text("#!/usr/bin/env python3\n"
                            "import subprocess,time\n"
                            "child=subprocess.Popen(['sleep','30'])\n"
                            "open('child.pid','w').write(str(child.pid))\n"
                            "print('synthetic child started',flush=True)\n"
                            "time.sleep(30)\n")
            fake.chmod(0o755)
            child_pid = None
            try:
                with patch("scrambox.partmc.write_partmc", return_value={"kernel": {"isorropia": "v23"}}), \
                     patch.dict(partmc.PARTMC_BIN, {"v23": str(fake)}):
                    with self.assertRaisesRegex(TimeoutError, "PartMC timed out"):
                        partmc.run_partmc({}, directory, timeout_s=0.5)
                child_pid = int((directory / "child.pid").read_text())
                for _ in range(30):
                    if not _live(child_pid):
                        break
                    time.sleep(0.1)
                self.assertFalse(_live(child_pid), "owned PartMC child survived timeout")
            finally:
                if child_pid is None and (directory / "child.pid").exists():
                    child_pid = int((directory / "child.pid").read_text())
                if child_pid is not None and _live(child_pid):
                    os.kill(child_pid, signal.SIGKILL)


if __name__ == "__main__":
    unittest.main()
