#!/usr/bin/env python3
"""Focused engine/analyzer tests for Issue #225."""

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
from analyze_c4_multiradius_pivotal import analyze, read_rows  # noqa: E402


class C4MultiRadiusPivotalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        compiler = shutil.which("c++") or shutil.which("g++")
        if compiler is None:
            raise unittest.SkipTest("C++ compiler unavailable")
        cls.temporary = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temporary.name)
        cls.multi = cls.root / "multi"
        cls.single = cls.root / "single"
        for source, output in (
            ("c4_multiradius_pivotal_mc.cpp", cls.multi),
            ("c4_local_odd_pivotal_mc.cpp", cls.single),
        ):
            subprocess.run(
                [compiler, "-std=c++17", "-O1", str(ROOT / "src" / source), "-o", str(output)],
                check=True,
            )

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temporary.cleanup()

    def run_engine(self, binary: Path, prefix: str, *extra: str) -> Path:
        output = self.root / prefix
        subprocess.run(
            [str(binary), "--samples", "400", "--batches", "10", "--seed", "225",
             "--output-prefix", str(output), *extra],
            check=True,
            capture_output=True,
            text=True,
        )
        return output

    def test_exact_oracle_and_multiradius_alignment(self) -> None:
        subprocess.run([str(self.multi), "--self-test"], check=True, capture_output=True, text=True)
        prefix = self.run_engine(self.multi, "dyadic", "--radii", "1,2,4")
        rows, sizes, radii, batches = read_rows(Path(f"{prefix}.batches.csv"))
        self.assertEqual(sizes, [130, 170])
        self.assertEqual(radii, [1, 2, 4])
        self.assertEqual(batches, list(range(10)))
        self.assertEqual(len(rows), 60)
        metadata = json.loads(Path(f"{prefix}.metadata.json").read_text(encoding="utf-8"))
        self.assertEqual(metadata["radii"], [1, 2, 4])
        self.assertEqual(metadata["cross_radius_coupling"], "same configuration and pivotal flag")

    def test_single_radius_reproduces_frozen_R3_stream(self) -> None:
        old = self.run_engine(self.single, "old", "--radius", "3")
        new = self.run_engine(self.multi, "new", "--radii", "3")
        with Path(f"{old}.batches.csv").open(newline="", encoding="utf-8") as handle:
            old_rows = {(int(row["n"]), int(row["batch"])): row for row in csv.DictReader(handle)}
        with Path(f"{new}.batches.csv").open(newline="", encoding="utf-8") as handle:
            new_rows = {(int(row["n"]), int(row["batch"])): row for row in csv.DictReader(handle)}
        self.assertEqual(set(old_rows), set(new_rows))
        mapping = {
            "sum_local_twice": "h4_minus",
            "local_twice_score_t": "h4_minus_score_t",
            "local_twice_score_lambda": "h4_minus_score_lambda",
            "sum_global_twice": "sum_global_twice",
            "global_twice_score_t": "global_twice_score_t",
            "global_twice_score_lambda": "global_twice_score_lambda",
        }
        for key in old_rows:
            for old_field, new_field in mapping.items():
                with self.subTest(key=key, field=old_field):
                    self.assertEqual(int(old_rows[key][old_field]), int(new_rows[key][new_field]))

    def test_analyzer_emits_full_aligned_covariance_and_shells(self) -> None:
        prefix = self.run_engine(self.multi, "analysis", "--radii", "1,2,4")
        result = analyze(Path(f"{prefix}.batches.csv"), Path(f"{prefix}.metadata.json"))
        self.assertEqual(len(result["points"]), 6)
        self.assertEqual(len(result["shell_increments"]), 4)
        base = result["base_vector"]
        shell = result["shell_vector"]
        self.assertEqual(len(base["order"]), 12)
        self.assertEqual([len(row) for row in base["covariance"]], [12] * 12)
        self.assertEqual(len(shell["order"]), 8)
        self.assertEqual([len(row) for row in shell["covariance"]], [8] * 8)

    def test_aliased_large_radius_fails_before_sampling(self) -> None:
        completed = subprocess.run(
            [str(self.multi), "--samples", "20", "--batches", "2", "--radii", "2,4,8",
             "--output-prefix", str(self.root / "aliased")],
            capture_output=True,
            text=True,
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIn("not injective", completed.stderr)


if __name__ == "__main__":
    unittest.main()
