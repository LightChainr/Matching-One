#!/usr/bin/env python3
"""Integration tests for the bounded same-N Gaussian discovery engine."""

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
        self.assertIn("matching channels", completed.stdout)

    def test_frozen_confirmation_designs_include_new_sizes(self) -> None:
        prefix = Path(self.temp.name) / "n170"
        subprocess.run(
            [
                str(self.binary), "--samples", "20", "--batches", "2",
                "--n", "170", "--seed", "19", "--threads", "1",
                "--git-commit", "test-sha", "--output-prefix", str(prefix),
            ],
            check=True,
        )
        metadata = json.loads(
            Path(str(prefix) + ".metadata.json").read_text(encoding="utf-8")
        )
        self.assertEqual(metadata["designs"][0]["N"], 170)
        self.assertEqual(metadata["designs"][0]["first"], [13, 1])
        self.assertEqual(metadata["designs"][0]["second"], [11, 7])

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
        covariance_json = Path(self.temp.name) / "covariance.json"
        delta_m_csv = Path(self.temp.name) / "delta_M_by_size_seed.csv"
        frozen_weights = Path(self.temp.name) / "frozen_weights.json"
        subprocess.run(
            [
                sys.executable, str(ROOT / "scripts" / "analyze_gaussian_orientation_mc.py"),
                "--batches", str(first) + ".batches.csv",
                "--metadata", str(first) + ".metadata.json",
                "--json", str(analysis_json), "--csv", str(analysis_csv),
                "--covariance-json", str(covariance_json),
                "--delta-m-csv", str(delta_m_csv),
                "--freeze-gls", str(frozen_weights),
            ],
            check=True,
        )
        payload = json.loads(analysis_json.read_text(encoding="utf-8"))
        summaries = payload["summaries"]
        channels = ("cross", "both", "either", "direction_0", "direction_1")
        self.assertEqual(len(summaries), 25)  # five channels x five sectors
        by_key = {(row["channel"], row["sector"]): row for row in summaries}
        for channel in channels:
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
        with Path(str(first) + ".batches.csv").open(newline="", encoding="utf-8") as handle:
            batch_rows = list(csv.DictReader(handle))
        self.assertEqual({row["channel"] for row in batch_rows}, set(channels))
        by_batch = {}
        for row in batch_rows:
            by_batch.setdefault(int(row["batch"]), []).append(row)
        for rows in by_batch.values():
            first_differences = {
                int(row["first_primal_sum"]) - int(row["first_matching_sum"])
                for row in rows
            }
            second_differences = {
                int(row["second_primal_sum"]) - int(row["second_matching_sum"])
                for row in rows
            }
            self.assertEqual(len(first_differences), 1)
            self.assertEqual(len(second_differences), 1)

        covariance = json.loads(covariance_json.read_text(encoding="utf-8"))["by_N"]["65"]
        with delta_m_csv.open(newline="", encoding="utf-8") as handle:
            delta_rows = list(csv.DictReader(handle))
        self.assertEqual(len(delta_rows), 1)
        self.assertEqual(delta_rows[0]["row_id"], "65:17")
        self.assertEqual(delta_rows[0]["N"], "65")
        self.assertEqual(delta_rows[0]["seed"], "17")
        self.assertAlmostEqual(float(delta_rows[0]["delta_M"]), 0.04)
        self.assertEqual(len(covariance["raw_orientation_channel_matrix"]["labels"]), 20)
        self.assertEqual(len(covariance["orientation_sector_effect_matrix"]["labels"]), 25)
        frozen = json.loads(frozen_weights.read_text(encoding="utf-8"))
        self.assertEqual(frozen["channel_names"], list(channels))
        self.assertAlmostEqual(sum(frozen["by_N"]["65"]["weights"]), 1.0)
        self.assertTrue(frozen["by_N"]["65"]["all_D_channels_identical_batchwise"])
        self.assertEqual(frozen["by_N"]["65"]["weights"], [0.2] * 5)

        overlap = subprocess.run(
            [
                sys.executable, str(ROOT / "scripts" / "analyze_gaussian_orientation_mc.py"),
                "--batches", str(first) + ".batches.csv",
                "--metadata", str(first) + ".metadata.json",
                "--json", str(Path(self.temp.name) / "overlap_analysis.json"),
                "--csv", str(Path(self.temp.name) / "overlap_sectors.csv"),
                "--frozen-gls", str(frozen_weights),
                "--evaluation-json", str(Path(self.temp.name) / "overlap_eval.json"),
                "--evaluation-csv", str(Path(self.temp.name) / "overlap_eval.csv"),
            ],
            text=True,
            capture_output=True,
        )
        self.assertNotEqual(overlap.returncode, 0)
        self.assertIn("overlap the pilot", overlap.stderr)

        evaluation = Path(self.temp.name) / "evaluation"
        evaluation_common = list(common)
        evaluation_common += ["--replica-offset", "400"]
        subprocess.run(
            [str(self.binary), *evaluation_common, "--output-prefix", str(evaluation)],
            check=True,
        )
        evaluation_json = Path(self.temp.name) / "evaluation.json"
        evaluation_csv = Path(self.temp.name) / "evaluation.csv"
        subprocess.run(
            [
                sys.executable, str(ROOT / "scripts" / "analyze_gaussian_orientation_mc.py"),
                "--batches", str(evaluation) + ".batches.csv",
                "--metadata", str(evaluation) + ".metadata.json",
                "--json", str(Path(self.temp.name) / "evaluation_analysis.json"),
                "--csv", str(Path(self.temp.name) / "evaluation_sectors.csv"),
                "--frozen-gls", str(frozen_weights),
                "--evaluation-json", str(evaluation_json),
                "--evaluation-csv", str(evaluation_csv),
            ],
            check=True,
        )
        evaluated = json.loads(evaluation_json.read_text(encoding="utf-8"))["by_N"]["65"]
        self.assertEqual(evaluated["target"], "orientation_difference")
        ratio = evaluated["estimators"]["equal_weight"]["variance_reduction_vs_optimized"]
        self.assertAlmostEqual(ratio, 1.0, places=10)
        self.assertTrue(evaluation_csv.exists())
        self.assertTrue(analysis_csv.exists())


if __name__ == "__main__":
    unittest.main()
