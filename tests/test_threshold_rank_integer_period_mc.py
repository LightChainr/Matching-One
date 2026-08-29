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

import yaml


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from integer_period_torus import integer_torus_geometry  # noqa: E402
from threshold_rank_nz import ThresholdRankEngine, counter_permutation  # noqa: E402


DESIGNS = {
    260: (
        ((16, -2), (2, 16)),
        ((14, -8), (8, 14)),
        (16, 2),
        (14, 8),
        (2, 130),
    ),
    340: (
        ((18, -4), (4, 18)),
        ((14, -12), (12, 14)),
        (18, 4),
        (14, 12),
        (2, 170),
    ),
    580: (
        ((24, -2), (2, 24)),
        ((18, -16), (16, 18)),
        (24, 2),
        (18, 16),
        (2, 290),
    ),
}


class ThresholdRankIntegerPeriodMCTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        compiler = os.environ.get("CXX") or shutil.which("g++") or shutil.which("clang++")
        if compiler is None:
            raise unittest.SkipTest("no C++17 compiler found")
        cls.temporary = tempfile.TemporaryDirectory()
        cls.binary = Path(cls.temporary.name) / "threshold_rank_integer_period_mc"
        command = [compiler, "-O2", "-std=c++17"]
        if sys.platform != "darwin" and "clang" not in Path(compiler).name:
            command.append("-fopenmp")
        command += [
            str(ROOT / "src" / "threshold_rank_integer_period_mc.cpp"),
            "-o", str(cls.binary),
        ]
        subprocess.run(command, check=True, cwd=ROOT)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temporary.cleanup()

    def run_predefined(self, n: int, *, samples: int = 8, threads: int = 1) -> Path:
        prefix = Path(self.temporary.name) / f"n{n}_{samples}_{threads}"
        subprocess.run(
            [
                str(self.binary), "--samples", str(samples), "--batches", "2",
                "--n", str(n), "--seed", "17", "--replica-offset", "5",
                "--threads", str(threads), "--git-commit", "test-sha",
                "--output-prefix", str(prefix),
            ],
            check=True,
        )
        return prefix

    def test_exact_self_test(self) -> None:
        completed = subprocess.run(
            [str(self.binary), "--self-test"],
            check=True, text=True, capture_output=True,
        )
        self.assertIn("arbitrary integer periods", completed.stdout)
        self.assertIn("basis invariance", completed.stdout)
        self.assertIn("Smith(2,130)/(2,170)", completed.stdout)
        self.assertIn("(2,290)", completed.stdout)

    def test_norm4_designs_and_smith_metadata(self) -> None:
        for n, (first_matrix, second_matrix, first_rep, second_rep, smith) in DESIGNS.items():
            prefix = self.run_predefined(n)
            metadata = json.loads(
                Path(str(prefix) + ".metadata.json").read_text(encoding="utf-8")
            )
            design = metadata["designs"][0]
            self.assertEqual(design["N"], n)
            self.assertEqual(design["first"], list(first_rep))
            self.assertEqual(design["second"], list(second_rep))
            self.assertEqual(design["first_period_matrix"], [list(row) for row in first_matrix])
            self.assertEqual(design["second_period_matrix"], [list(row) for row in second_matrix])
            self.assertEqual(design["first_smith_invariants"], list(smith))
            self.assertEqual(design["second_smith_invariants"], list(smith))
            self.assertEqual(metadata["channel"], "rank-2 cross wrapping")

    def test_cpp_n260_matches_general_period_python_oracle(self) -> None:
        samples = 24
        prefix = self.run_predefined(260, samples=samples)
        metadata = json.loads(
            Path(str(prefix) + ".metadata.json").read_text(encoding="utf-8")
        )
        design_metadata = metadata["designs"][0]

        cpp = {
            "first": {"minus": [0] * 261, "plus": [0] * 261},
            "second": {"minus": [0] * 261, "plus": [0] * 261},
        }
        with Path(str(prefix) + ".hist.csv").open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                cpp[row["orientation"]][row["kind"]][int(row["k"])] += int(row["count"])

        for orientation, matrix_key, hnf_key in (
            ("first", "first_period_matrix", "first_HNF"),
            ("second", "second_period_matrix", "second_HNF"),
        ):
            matrix = tuple(tuple(value for value in row) for row in design_metadata[matrix_key])
            geometry = integer_torus_geometry(matrix)
            engine = ThresholdRankEngine(geometry)
            h11 = design_metadata[hnf_key][0][0]
            expected = {"minus": [0] * 261, "plus": [0] * 261}
            for counter in range(5, 5 + samples):
                labels = counter_permutation(260, 17, counter)
                permutation = tuple(
                    geometry.vertex((label % h11, label // h11)) for label in labels
                )
                k_minus, k_plus = engine.threshold_ranks(permutation)
                expected["minus"][k_minus] += 1
                expected["plus"][k_plus] += 1
            self.assertEqual(cpp[orientation], expected)

    def test_n580_matches_python_oracle_and_is_thread_reproducible(self) -> None:
        samples = 8
        first_prefix = self.run_predefined(580, samples=samples, threads=1)
        second_prefix = self.run_predefined(580, samples=samples, threads=2)
        for suffix in (".hist.csv", ".moments.csv"):
            self.assertEqual(
                Path(str(first_prefix) + suffix).read_bytes(),
                Path(str(second_prefix) + suffix).read_bytes(),
            )

        metadata = json.loads(
            Path(str(first_prefix) + ".metadata.json").read_text(encoding="utf-8")
        )
        design = metadata["designs"][0]
        cpp = {
            "first": {"minus": [0] * 581, "plus": [0] * 581},
            "second": {"minus": [0] * 581, "plus": [0] * 581},
        }
        with Path(str(first_prefix) + ".hist.csv").open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                cpp[row["orientation"]][row["kind"]][int(row["k"])] += int(row["count"])

        for orientation, matrix_key, hnf_key in (
            ("first", "first_period_matrix", "first_HNF"),
            ("second", "second_period_matrix", "second_HNF"),
        ):
            matrix = tuple(tuple(value for value in row) for row in design[matrix_key])
            geometry = integer_torus_geometry(matrix)
            engine = ThresholdRankEngine(geometry)
            h11 = design[hnf_key][0][0]
            expected = {"minus": [0] * 581, "plus": [0] * 581}
            for counter in range(5, 5 + samples):
                labels = counter_permutation(580, 17, counter)
                permutation = tuple(
                    geometry.vertex((label % h11, label // h11)) for label in labels
                )
                k_minus, k_plus = engine.threshold_ranks(permutation)
                expected["minus"][k_minus] += 1
                expected["plus"][k_plus] += 1
            self.assertEqual(cpp[orientation], expected)

    def test_custom_arbitrary_matrices_and_thread_reproducibility(self) -> None:
        custom_prefix = Path(self.temporary.name) / "custom"
        subprocess.run(
            [
                str(self.binary), "--samples", "8", "--batches", "2",
                "--first-matrix", "3", "1", "1", "2",
                "--second-matrix", "3", "4", "1", "3",
                "--threads", "1", "--output-prefix", str(custom_prefix),
            ],
            check=True,
        )
        custom = json.loads(
            Path(str(custom_prefix) + ".metadata.json").read_text(encoding="utf-8")
        )["designs"][0]
        self.assertEqual(custom["N"], 5)
        self.assertEqual(custom["first_smith_invariants"], [1, 5])
        self.assertEqual(custom["second_smith_invariants"], [1, 5])

        first = self.run_predefined(260, samples=40, threads=1)
        second = self.run_predefined(260, samples=40, threads=2)
        for suffix in (".hist.csv", ".moments.csv"):
            self.assertEqual(
                Path(str(first) + suffix).read_bytes(),
                Path(str(second) + suffix).read_bytes(),
            )

    def test_negative_unsigned_cli_value_is_rejected(self) -> None:
        completed = subprocess.run(
            [
                str(self.binary), "--samples", "-2", "--batches", "2",
                "--n", "260", "--output-prefix",
                str(Path(self.temporary.name) / "must-not-run"),
            ],
            text=True, capture_output=True,
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIn("negative value for --samples", completed.stderr)

    def test_p200_production_contract_has_disjoint_counter_domain(self) -> None:
        manifest = yaml.safe_load(
            (ROOT / "experiments" / "p200_n580_phaseA_100m_20260829.yaml").read_text()
        )
        self.assertEqual(manifest["scope"], "N580_only_no_N650_path_flags")
        self.assertEqual(manifest["sampling"]["samples_per_pair"], 100_000_000)
        self.assertEqual(manifest["sampling"]["batches"], 100)
        self.assertEqual(manifest["sampling"]["threads"], 8)
        self.assertEqual(manifest["sampling"]["seed"], 2026102001)
        self.assertEqual(
            manifest["sampling"]["replica_counter_last_exclusive"]
            - manifest["sampling"]["replica_counter_first"],
            100_000_000,
        )
        self.assertIn("threshold_rank_integer_period_mc.cpp", manifest["build"]["command"])
        self.assertNotIn("path", manifest.get("geometry", {}))

    def test_p200_score_input_is_one_joint_four_state_view(self) -> None:
        prediction = json.loads(
            (ROOT / "predictions" / "p200_n580_q2_jordan_score_input_20260829.json").read_text()
        )
        self.assertEqual(prediction["state_order"], ["I_S", "I_Du", "T_D", "T_Su"])
        self.assertEqual(set(prediction["models"]), {"ordinary_q2", "rank2_Jordan"})
        for model in prediction["models"].values():
            covariance = model["N580_state_prediction_covariance"]
            self.assertEqual([len(row) for row in covariance], [4, 4, 4, 4])
            self.assertEqual(model["covariance_coordinate_order"], prediction["state_order"])
        self.assertIn("joint four-coordinate", prediction["score_semantics"]["primary"])


if __name__ == "__main__":
    unittest.main()
