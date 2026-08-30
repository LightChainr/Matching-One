from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class P154LocalSingletPilotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        compiler = os.environ.get("CXX") or shutil.which("g++") or shutil.which("clang++")
        if compiler is None:
            raise unittest.SkipTest("no C++17 compiler")
        cls.temporary = tempfile.TemporaryDirectory()
        cls.binary = Path(cls.temporary.name) / "p154_local_singlet_pilot"
        command = [compiler, "-O2", "-std=c++17"]
        if sys.platform != "darwin" and "clang" not in Path(compiler).name:
            command.append("-fopenmp")
        command += [str(ROOT / "src" / "p154_local_singlet_pilot.cpp"), "-o", str(cls.binary)]
        subprocess.run(command, check=True, cwd=ROOT)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temporary.cleanup()

    def test_self_test(self) -> None:
        completed = subprocess.run([str(self.binary), "--self-test"], check=True,
                                   text=True, capture_output=True)
        self.assertIn("local-singlet self-test passed", completed.stdout)

    def test_tiny_same_stream_output(self) -> None:
        prefix = Path(self.temporary.name) / "tiny"
        subprocess.run([
            str(self.binary), "--n", "65", "--samples", "20", "--batches", "2",
            "--seed", "202615465", "--replica-offset", "15465000000", "--threads", "1",
            "--git-commit", "freeze-test", "--output-prefix", str(prefix),
        ], check=True)
        with Path(str(prefix) + ".batches.csv").open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 4)
        self.assertEqual({row["orientation"] for row in rows}, {"first", "second"})
        for row in rows:
            self.assertEqual(int(row["samples"]), 10)
            self.assertEqual(int(row["sum_i0"]) + int(row["sum_i1"]) + int(row["sum_i2"]), 10)
            self.assertLessEqual(int(row["sum_black_axis_pairs"]), 2 * 65 * 10)
            self.assertLessEqual(int(row["sum_white_matching_axis_pairs"]), 2 * 65 * 10)
            self.assertIn("sum_even_numerator_squared", row)
            self.assertIn("sum_i0_even_numerator", row)
            self.assertIn("sum_i2_even_numerator", row)
        metadata = json.loads(Path(str(prefix) + ".metadata.json").read_text())
        self.assertEqual(metadata["replica_counter_first"], 15465000000)
        self.assertEqual(metadata["replica_counter_last_exclusive"], 15465000020)
        self.assertEqual(metadata["git_commit"], "freeze-test")

    def test_locked_score_runs_and_reports_all_candidates(self) -> None:
        output = Path(self.temporary.name) / "score.json"
        report = Path(self.temporary.name) / "REPORT.md"
        subprocess.run([
            sys.executable, str(ROOT / "scripts" / "score_p154_local_singlet_pilot.py"),
            "--raw-dir", str(ROOT / "results" / "p154-phase-e-local-singlet-pilot" / "raw"),
            "--json", str(output), "--report", str(report),
        ], check=True)
        payload = json.loads(output.read_text())
        self.assertEqual(set(payload["common_ray_candidates"]), {"A/E/C", "A/E/J_even", "A/E/J_odd"})
        self.assertIn(payload["promotion_gate"]["decision"], {"extend_both_to_100k", "stop_at_20k"})
        self.assertIn("continuum-energy identification", payload["interpretation_boundary"])


if __name__ == "__main__":
    unittest.main()
