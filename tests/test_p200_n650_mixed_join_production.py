from __future__ import annotations

import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_p200_n650_mixed_join import render  # noqa: E402


class P200N650MixedJoinProductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.temp = tempfile.TemporaryDirectory()
        cls.directory = Path(cls.temp.name)
        cls.binary = cls.directory / "p200_n650_mixed_join_mc"
        compiler = shutil.which("c++")
        if compiler is None:
            raise unittest.SkipTest("C++ compiler unavailable")
        subprocess.run(
            [compiler, "-O2", "-std=c++17", str(ROOT / "src" / "p200_n650_mixed_join_mc.cpp"), "-o", str(cls.binary)],
            check=True,
        )

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temp.cleanup()

    def run_smoke(self, name: str, threads: int = 1, samples: int = 2000) -> Path:
        prefix = self.directory / name
        subprocess.run(
            [
                str(self.binary), "--samples", str(samples), "--batches", "20",
                "--seed", "2026102003", "--replica-offset", "18000000000",
                "--threads", str(threads), "--git-commit", "TEST", "--output-prefix", str(prefix),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        return prefix

    def test_exact_tiny_gate(self) -> None:
        result = subprocess.run([str(self.binary), "--self-test"], check=True, capture_output=True, text=True)
        self.assertIn("N10 exhaustive typed joins", result.stdout)
        self.assertIn("N650 HNF/lifts", result.stdout)

    def test_worker_count_does_not_change_batch_sums(self) -> None:
        first = self.run_smoke("worker1", threads=1, samples=200)
        second = self.run_smoke("worker2", threads=2, samples=200)
        self.assertEqual(
            (first.with_suffix(".batches.csv")).read_bytes(),
            (second.with_suffix(".batches.csv")).read_bytes(),
        )

    def test_scorer_reads_runner_contract_and_full_covariance(self) -> None:
        prefix = self.run_smoke("score")
        payload = render(
            prefix.with_suffix(".batches.csv"),
            prefix.with_suffix(".metadata.json"),
            ROOT / "predictions" / "p200_n650_mixed_join_phaseB_20260829.json",
        )
        self.assertEqual(payload["primary"]["state_order"], ["ES", "ED", "OS", "OD"])
        covariance = payload["primary"]["delete_one_covariance"]
        self.assertEqual(len(covariance), 4)
        for row in range(4):
            for column in range(4):
                self.assertAlmostEqual(covariance[row][column], covariance[column][row])
        self.assertEqual(payload["primary"]["joint_GLS"]["degrees_of_freedom"], 4)
        self.assertTrue(payload["decision"]["stop_recommended"])

    def test_large_run_requires_explicit_contract_gate(self) -> None:
        result = subprocess.run(
            [str(self.binary), "--samples", "20100", "--batches", "100", "--output-prefix", str(self.directory / "blocked")],
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("requires explicit --production", result.stderr)


if __name__ == "__main__":
    unittest.main()
