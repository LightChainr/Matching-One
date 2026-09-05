
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


if __name__ == "__main__":
    unittest.main()
