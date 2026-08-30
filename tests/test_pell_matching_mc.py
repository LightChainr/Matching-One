#!/usr/bin/env python3
"""Integration smoke tests for the discovery-stage Pell Monte Carlo engine."""

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


class PellMatchingMCTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        compiler = os.environ.get("CXX") or shutil.which("g++") or shutil.which("clang++")
        if compiler is None:
            raise unittest.SkipTest("no C++17 compiler found")
        cls.temp = tempfile.TemporaryDirectory()
        cls.binary = Path(cls.temp.name) / "pell_matching_mc"
        command = [compiler, "-O2", "-std=c++17"]
        if sys.platform != "darwin" and "clang" not in Path(compiler).name:
            command.append("-fopenmp")
        command += [str(ROOT / "src" / "pell_matching_mc.cpp"), "-o", str(cls.binary)]
        subprocess.run(command, check=True, cwd=ROOT)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temp.cleanup()

    def test_exact_self_test(self) -> None:
        completed = subprocess.run(
            [str(self.binary), "--self-test"], check=True, text=True, capture_output=True
        )
        self.assertIn("self-test passed", completed.stdout)

    def test_simulation_is_reproducible_and_analyzable(self) -> None:
        prefix_a = Path(self.temp.name) / "run_a"
        prefix_b = Path(self.temp.name) / "run_b"
        common = [
            "--axis", "3", "--diamond", "2", "--samples", "400",
            "--batches", "4", "--p-ref", "0.59", "--h", "0.02",
            "--seed", "17", "--threads", "1",
        ]
        subprocess.run([str(self.binary), *common, "--output-prefix", str(prefix_a)], check=True)
        subprocess.run([str(self.binary), *common, "--output-prefix", str(prefix_b)], check=True)
        self.assertEqual(
            Path(f"{prefix_a}.batches.csv").read_bytes(),
            Path(f"{prefix_b}.batches.csv").read_bytes(),
        )

        analysis_json = Path(self.temp.name) / "analysis.json"
        analysis_csv = Path(self.temp.name) / "analysis.csv"
        subprocess.run(
            [
                sys.executable, str(ROOT / "scripts" / "analyze_pell_mc.py"),
                "--batches", f"{prefix_a}.batches.csv",
                "--metadata", f"{prefix_a}.metadata.json",
                "--json", str(analysis_json), "--csv", str(analysis_csv),
            ],
            check=True,
        )
        result = json.loads(analysis_json.read_text(encoding="utf-8"))
        self.assertEqual(result["batch_count"], 4)
        self.assertIn("axis_linear_root", result["derived"])
        self.assertIn("orientation_root_gap_diamond_minus_axis", result["derived"])
        self.assertTrue(analysis_csv.exists())


if __name__ == "__main__":
    unittest.main()
