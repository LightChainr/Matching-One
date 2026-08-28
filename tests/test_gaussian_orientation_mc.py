#!/usr/bin/env python3
"""Integration tests for the bounded same-N Gaussian discovery engine."""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class GaussianOrientationMCTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        compiler = os.environ.get("CXX") or shutil.which("g++") or shutil.which("clang++")
        if compiler is None:
            raise unittest.SkipTest("no C++17 compiler found")
        cls.temp = tempfile.TemporaryDirectory()
        cls.binary = Path(cls.temp.name) / "gaussian_orientation_mc"
        command = [compiler, "-O2", "-std=c++17"]
        if sys.platform != "darwin" and "clang" not in Path(compiler).name:
            command.append("-fopenmp")
        command += [str(ROOT / "src" / "gaussian_orientation_mc.cpp"), "-o", str(cls.binary)]
        subprocess.run(command, check=True, cwd=ROOT)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temp.cleanup()

    def test_exact_reference_self_test(self) -> None:
        completed = subprocess.run(
            [str(self.binary), "--self-test"], check=True, text=True, capture_output=True
        )
        self.assertIn("exhaustive N=5,13", completed.stdout)

    def test_reproducible_batches_and_sector_analysis(self) -> None:
        first = Path(self.temp.name) / "first"
        second = Path(self.temp.name) / "second"
        common = [
            "--samples", "400", "--batches", "4", "--n", "65",
            "--p-ref", "0.59274605", "--seed", "17", "--threads", "1",
            "--git-commit", "test-sha",
        ]
        subprocess.run([str(self.binary), *common, "--output-prefix", str(first)], check=True)
        second_common = list(common)
        second_common[second_common.index("--threads") + 1] = "2"
        subprocess.run(
            [str(self.binary), *second_common, "--output-prefix", str(second)], check=True
        )
        self.assertEqual(
            Path(str(first) + ".batches.csv").read_bytes(),
            Path(str(second) + ".batches.csv").read_bytes(),
        )
        metadata = json.loads(Path(str(first) + ".metadata.json").read_text(encoding="utf-8"))
        self.assertEqual(metadata["git_commit"], "test-sha")
        self.assertEqual(metadata["designs"][0]["N"], 65)

        analysis_json = Path(self.temp.name) / "analysis.json"
        analysis_csv = Path(self.temp.name) / "analysis.csv"
        subprocess.run(
            [
                sys.executable, str(ROOT / "scripts" / "analyze_gaussian_orientation_mc.py"),
                "--batches", str(first) + ".batches.csv",
                "--metadata", str(first) + ".metadata.json",
                "--json", str(analysis_json), "--csv", str(analysis_csv),
            ],
            check=True,
        )
        payload = json.loads(analysis_json.read_text(encoding="utf-8"))
        summaries = payload["summaries"]
        self.assertEqual(len(summaries), 10)  # two channels x five sectors
        by_key = {(row["channel"], row["sector"]): row for row in summaries}
        for channel in ("either", "cross"):
            odd = by_key[(channel, "odd")]
            matching_function = by_key[(channel, "matching_function")]
            self.assertAlmostEqual(
                matching_function["difference_first_minus_second"],
                2 * odd["difference_first_minus_second"],
            )
            self.assertAlmostEqual(
                matching_function["difference_batch_se"],
                2 * odd["difference_batch_se"],
            )
        self.assertTrue(analysis_csv.exists())


if __name__ == "__main__":
    unittest.main()
