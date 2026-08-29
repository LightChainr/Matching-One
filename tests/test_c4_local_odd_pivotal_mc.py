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

from analyze_c4_local_odd_pivotal_mc import render  # noqa: E402


class C4LocalOddPivotalMCTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler = (shutil.which("g++-16") or shutil.which("g++") or
                    shutil.which("clang++"))
        if compiler is None:
            raise unittest.SkipTest("no C++17 compiler found")
        cls.temporary = tempfile.TemporaryDirectory()
        cls.binary = Path(cls.temporary.name) / "c4_local_odd_pivotal_mc"
        command = [compiler, "-O2", "-std=c++17"]
        if Path(compiler).name.startswith("g++-") or (
            sys.platform != "darwin" and "clang" not in Path(compiler).name
        ):
            command.append("-fopenmp")
        command += [str(ROOT / "src" / "c4_local_odd_pivotal_mc.cpp"),
                    "-o", str(cls.binary)]
        subprocess.run(command, check=True, cwd=ROOT)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_cpp_matches_exact_n10_oracle(self):
        completed = subprocess.run(
            [str(self.binary), "--self-test"], check=True,
            text=True, capture_output=True,
        )
        self.assertIn("exact N=10 local/global response matrix", completed.stdout)

    def make_synthetic(self, matrices):
        root = Path(self.temporary.name)
        batches = root / "synthetic.batches.csv"
        metadata = root / "synthetic.metadata.json"
        fields = [
            "n", "a", "b", "batch", "counter_first", "counter_last_exclusive",
            "samples", "sum_score_t", "sum_score_lambda", "sum_score_t2",
            "sum_score_lambda2", "sum_score_cross", "sum_global_twice",
            "sum_local_twice", "global_twice_score_t",
            "global_twice_score_lambda", "local_twice_score_t",
            "local_twice_score_lambda",
        ]
        with batches.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            for batch in range(100):
                for n, a, b in ((130, 11, 3), (170, 13, 1)):
                    matrix = matrices[n]
                    samples = 2000
                    writer.writerow({
                        "n": n, "a": a, "b": b, "batch": batch,
                        "counter_first": 2000 * batch,
                        "counter_last_exclusive": 2000 * (batch + 1),
                        "samples": samples, "sum_score_t": 0,
                        "sum_score_lambda": 0, "sum_score_t2": 4 * n * samples,
                        "sum_score_lambda2": 4 * n * samples,
                        "sum_score_cross": 0, "sum_global_twice": 0,
                        "sum_local_twice": 0,
                        "global_twice_score_t": int(2 * samples * matrix[0][0]),
                        "global_twice_score_lambda": int(2 * samples * matrix[0][1]),
                        "local_twice_score_t": int(2 * samples * matrix[1][0]),
                        "local_twice_score_lambda": int(2 * samples * matrix[1][1]),
                    })
        metadata.write_text(json.dumps({
            "schema": "matching-one/c4-local-odd-pivotal-score-stream/v1",
            "samples_per_size": 200000, "batches": 100, "radius": 3,
            "cross_size_coupling": "same seed/counter and prefix-coupled site bits",
            "designs": [
                {"N": 130, "a": 11, "b": 3},
                {"N": 170, "a": 13, "b": 1},
            ],
        }), encoding="utf-8")
        return batches, metadata

    def test_full_rank_gate_computes_matrix_pencil(self):
        batches, metadata = self.make_synthetic({
            130: [[2.0, 0.0], [0.0, 1.0]],
            170: [[3.0, 0.0], [0.0, 0.5]],
        })
        result = render(batches, metadata)
        self.assertTrue(result["sizes"]["130"]["gate"]["passes"])
        self.assertEqual(len(result["sizes"]["170"]["delete_one"]), 100)
        pencil = result["generalized_eigensystem"]
        self.assertTrue(pencil["computed"])
        self.assertAlmostEqual(pencil["branches_sorted_by_modulus"][0]["modulus"], 1.5)
        self.assertAlmostEqual(pencil["branches_sorted_by_modulus"][1]["modulus"], 0.5)

    def test_rank_one_gate_suppresses_second_eigenvalue(self):
        batches, metadata = self.make_synthetic({
            130: [[2.0, 1.0], [4.0, 2.0]],
            170: [[3.0, 1.5], [6.0, 3.0]],
        })
        result = render(batches, metadata)
        self.assertFalse(result["generalized_eigensystem"]["computed"])
        self.assertEqual(result["generalized_eigensystem"]["resolved_dimension"], 1)


if __name__ == "__main__":
    unittest.main()
