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

    def run_predefined(
        self, n: int, *, samples: int = 8, threads: int = 1, marked: bool = False
    ) -> Path:
        prefix = Path(self.temporary.name) / f"n{n}_{samples}_{threads}_{int(marked)}"
        command = [
                str(self.binary), "--samples", str(samples), "--batches", "2",
                "--n", str(n), "--seed", "17", "--replica-offset", "5",
                "--threads", str(threads), "--git-commit", "test-sha",
                "--output-prefix", str(prefix),
        ]
        if marked:
            command.append("--marked-births")
        subprocess.run(command, check=True)
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

    def test_marked_path_full_source_frame_and_index_schema(self) -> None:
        prefix = self.run_predefined(65, samples=20, marked=True)
        metadata = json.loads(
            Path(str(prefix) + ".metadata.json").read_text(encoding="utf-8")
        )
        self.assertTrue(metadata["marked_birth_schema"])
        self.assertIn("active_S+inactive_S", metadata["full_source"])
        self.assertIn("primitive(P*ell)", metadata["chi4_frame"])
        self.assertIn("raw winding coefficients", metadata["saturation_index"])
        self.assertIn("C_black_NN-C_white_matching-q", metadata["external_observer"])
        self.assertIn("q*J_D4 retained only as contact control", metadata["external_products"])

        with Path(str(prefix) + ".path.csv").open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 2 * 2 * 65)
        required_external = {
            "sum_O_ext", "sum_O_ext2", "sum_O_ext_J_S_re", "sum_O_ext_J_S_im",
            "sum_O_ext_J_D_re", "sum_O_ext_J_D_im", "sum_J_D_conj_J_S_re",
            "sum_J_D_conj_J_S_im", "sum_abs_J_S2",
        }
        self.assertTrue(required_external.issubset(rows[0]))
        for row in rows:
            active_s = int(row["sum_active_S"])
            inactive_s = int(row["sum_inactive_S"])
            active_d = int(row["sum_active_D"])
            inactive_d = int(row["sum_inactive_D"])
            self.assertEqual(2 * int(row["sum_site_S"]), active_s + inactive_s)
            self.assertEqual(2 * int(row["sum_site_D"]), active_d - inactive_d)
            self.assertGreaterEqual(int(row["sum_O_ext2"]), 0)
            self.assertAlmostEqual(float(row["sum_J_D_conj_J_S_im"]), 0.0, places=12)

        with Path(str(prefix) + ".marked_births.csv").open(
            newline="", encoding="utf-8"
        ) as handle:
            sparse = list(csv.DictReader(handle))
        self.assertTrue(sparse)
        for row in sparse:
            if row["direct_0_to_2"] == "1":
                self.assertEqual((row["line_null"], row["iota01"], row["iota12"]), ("1", "0", "0"))
            else:
                self.assertEqual(row["line_null"], "0")
                self.assertGreaterEqual(max(int(row["iota01"]), int(row["iota12"])), 1)
                x, y = int(row["physical_x"]), int(row["physical_y"])
                self.assertNotEqual((x, y), (0, 0))

        with Path(str(prefix) + ".complement_audit.csv").open(
            newline="", encoding="utf-8"
        ) as handle:
            audits = list(csv.DictReader(handle))
        for row in audits:
            for name in (
                "endpoint_failures", "site_failures", "line_failures",
                "local_mark_failures", "index_mismatches",
            ):
                self.assertEqual(int(row[name]), 0)

        score_path = Path(self.temporary.name) / "external-score.json"
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts" / "score_marked_birth_path.py"),
                "--prefix", str(prefix),
                "--output", str(score_path),
            ],
            check=True,
        )
        score = json.loads(score_path.read_text(encoding="utf-8"))
        self.assertEqual(score["schema"], "matching-one/marked-birth-path-score/v2")
        for name in (
            "P4_connected_O_ext_J_D_re",
            "P4_connected_O_ext_J_S_re",
            "P4_Gram_J_D_conj_J_S_re",
            "P4_Gram_abs_J_S2",
        ):
            self.assertIn(name, score["P4_point"])
            self.assertIn(name, score["covariance_metric_order"])


if __name__ == "__main__":
    unittest.main()
