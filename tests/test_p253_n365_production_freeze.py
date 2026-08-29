#!/usr/bin/env python3

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
import yaml


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "experiments/p253_n365_annulus_recurrence_20260829.yaml"


def cos4(a: int, b: int) -> Fraction:
    n = a * a + b * b
    return Fraction(a**4 - 6 * a * a * b * b + b**4, n * n)


class P253N365ProductionFreezeTests(unittest.TestCase):
    def test_orientation_arithmetic_and_fresh_counter_domain(self) -> None:
        payload = yaml.safe_load(MANIFEST.read_text(encoding="utf-8"))
        first, second = cos4(14, 13), cos4(19, 2)
        self.assertEqual(first, Fraction(-131767, 133225))
        self.assertEqual(second, Fraction(121673, 133225))
        self.assertEqual(first - second, Fraction(-50688, 26645))
        self.assertEqual(payload["orientation_contract"]["Delta_cos4"], "-50688/26645")
        production = payload["production"]
        self.assertEqual(
            production["replica_counter_last_exclusive"]
            - production["replica_counter_first"],
            production["samples_per_design"],
        )
        self.assertNotEqual(production["seed"], 22550260829)

    def test_engine_accepts_exact_frozen_design_and_radii(self) -> None:
        compiler = shutil.which("c++") or shutil.which("g++")
        if compiler is None:
            self.skipTest("C++ compiler unavailable")
        with tempfile.TemporaryDirectory() as directory:
            binary = Path(directory) / "engine"
            subprocess.run(
                [compiler, "-O1", "-std=c++17",
                 str(ROOT / "src/matching_multiradius_pivotal_mc.cpp"),
                 "-o", str(binary)],
                check=True,
            )
            completed = subprocess.run(
                [str(binary), "--validate-only", "--radii", "2,4,7,8",
                 "--cutoff", "euclidean", "--design", "n365_first,14,13",
                 "--design", "n365_second,19,2"],
                check=True, capture_output=True, text=True,
            )
            self.assertIn("validated 2 designs at radii 2 4 7 8", completed.stdout)

    def test_production_entrypoint_is_frozen_but_not_executed(self) -> None:
        script = (ROOT / "scripts/run_p253_n365_annulus_20260829.sh").read_text(
            encoding="utf-8")
        for token in (
            "--samples 200000", "--batches 200", "--threads 16",
            "--seed 25336560829", "--replica-offset 25336500000",
            "--design n365_first,14,13", "--design n365_second,19,2",
        ):
            self.assertIn(token, script)


if __name__ == "__main__":
    unittest.main()
