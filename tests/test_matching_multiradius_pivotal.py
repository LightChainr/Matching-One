#!/usr/bin/env python3
"""Execution-contract tests for the Issue #225 norm-5 pilot."""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from analyze_norm5_multiradius_pivotal import analyze, read_rows  # noqa: E402


def cos4(a: int, b: int) -> Fraction:
    n = a * a + b * b
    return Fraction(a**4 - 6 * a * a * b * b + b**4, n * n)


class MatchingMultiRadiusPivotalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        compiler = shutil.which("c++") or shutil.which("g++")
        if compiler is None:
            raise unittest.SkipTest("C++ compiler unavailable")
        cls.temporary = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temporary.name)
        cls.binary = cls.root / "engine"
        subprocess.run(
            [compiler, "-std=c++17", "-O1",
             str(ROOT / "src" / "matching_multiradius_pivotal_mc.cpp"),
             "-o", str(cls.binary)],
            check=True,
        )

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temporary.cleanup()

    def test_norm5_radii_are_injective_only_under_frozen_disk_scheme(self) -> None:
        good = subprocess.run(
            [str(self.binary), "--validate-only", "--radii", "2,4,8",
             "--cutoff", "euclidean"], capture_output=True, text=True)
        self.assertEqual(good.returncode, 0, good.stderr)
        bad = subprocess.run(
            [str(self.binary), "--validate-only", "--radii", "2,4,8",
             "--cutoff", "chebyshev"], capture_output=True, text=True)
        self.assertEqual(bad.returncode, 2)
        self.assertIn("N=425,R=8", bad.stderr)

    def test_registry_order_has_no_cross_size_sign_reversal(self) -> None:
        delta325 = cos4(17, 6) - cos4(18, 1)
        delta425 = cos4(16, 13) - cos4(19, 8)
        self.assertEqual(delta325, Fraction(-16128, 21125))
        self.assertEqual(delta425, Fraction(-32256, 36125))
        self.assertLess(delta325, 0)
        self.assertLess(delta425, 0)

    def test_small_stream_and_frozen_scorer_shapes(self) -> None:
        prefix = self.root / "pilot"
        subprocess.run(
            [str(self.binary), "--samples", "4000", "--batches", "20",
             "--radii", "2,4,7,8", "--cutoff", "euclidean",
             "--seed", "22550260829", "--replica-offset", "15000000000",
             "--git-commit", "test", "--output-prefix", str(prefix)],
            check=True, capture_output=True, text=True)
        rows, batches = read_rows(Path(f"{prefix}.batches.csv"))
        self.assertEqual(len(rows), 4 * 4 * 20)
        self.assertEqual(batches, list(range(20)))
        result = analyze(
            Path(f"{prefix}.batches.csv"), Path(f"{prefix}.metadata.json"),
            require_production=False)
        self.assertEqual(len(result["same_R_UV"]), 6)
        self.assertEqual(len(result["dyadic_shells"]), 4)
        self.assertEqual(len(result["contrast_vector"]["order"]), 16)
        self.assertEqual(len(result["shell_vector"]["order"]), 8)
        self.assertLess(result["matched_delta"]["relative_delta_mismatch"], 0.001)

    def test_runtime_custom_designs(self) -> None:
        completed = subprocess.run(
            [str(self.binary), "--validate-only", "--radii", "1,2",
             "--design", "first,8,1", "--design", "second,7,4"],
            capture_output=True, text=True)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("validated 2 designs", completed.stdout)


if __name__ == "__main__":
    unittest.main()
