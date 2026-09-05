
from __future__ import annotations
import csv
from collections import Counter
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from integer_period_torus import axis_integer_torus  # noqa: E402
from threshold_rank_nz import ThresholdRankEngine, counter_permutation  # noqa: E402


class ThresholdRankAxisCppTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        compiler = shutil.which("g++") or shutil.which("c++")
        if compiler is None:
            raise unittest.SkipTest("no C++ compiler")
        cls.temporary = tempfile.TemporaryDirectory()
        cls.binary = Path(cls.temporary.name) / "threshold_rank_axis_mc"
        subprocess.run(
            [
                compiler,
                "-O2",
                "-std=c++17",
                str(ROOT / "src" / "threshold_rank_axis_mc.cpp"),
                "-o",
                str(cls.binary),
            ],
            check=True,
        )

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temporary.cleanup()

    def test_exact_self_test(self) -> None:
        completed = subprocess.run(
            [str(self.binary), "--self-test"],
            check=True,
            text=True,
            capture_output=True,
        )
        self.assertIn("axis L=2 all 24 permutations", completed.stdout)

    @staticmethod
    def read_histogram(path: Path):
        result = {}
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                key = (int(row["batch"]), row["kind"])
                result.setdefault(key, Counter())[int(row["k"])] = int(row["count"])
        return result

    def test_cpp_matches_python_axis_oracle_for_counter_stream(self) -> None:
        L = 3
        samples = 20
        batches = 2
        seed = 1701
        offset = 900
        prefix = Path(self.temporary.name) / "axis_l3"
        subprocess.run(
            [
                str(self.binary),
                "--L", str(L),
                "--samples", str(samples),
                "--batches", str(batches),
                "--seed", str(seed),
                "--replica-offset", str(offset),
                "--threads", "1",
                "--git-commit", "test",
                "--output-prefix", str(prefix),
            ],
            check=True,
            text=True,
            capture_output=True,
        )
        observed = self.read_histogram(Path(str(prefix) + ".hist.csv"))

        geometry = axis_integer_torus(L)
        engine = ThresholdRankEngine(geometry)
        per_batch = samples // batches
        expected = {}
        for batch in range(batches):
            minus = Counter()
            plus = Counter()
            begin = offset + batch * per_batch
            for counter in range(begin, begin + per_batch):
                permutation = counter_permutation(geometry.n, seed, counter)
                k_minus, k_plus = engine.threshold_ranks(permutation)
                minus[k_minus] += 1
                plus[k_plus] += 1
            expected[(batch, "minus")] = minus
            expected[(batch, "plus")] = plus

        self.assertEqual(observed, expected)

    def test_metadata_schema_marks_axis_geometry(self) -> None:
        prefix = Path(self.temporary.name) / "axis_l2_meta"
        subprocess.run(
            [
                str(self.binary),
                "--L", "2",
                "--samples", "4",
                "--batches", "2",
                "--seed", "19",
                "--threads", "1",
                "--git-commit", "abc123",
                "--output-prefix", str(prefix),
            ],
            check=True,
            text=True,
            capture_output=True,
        )
        import json
        payload = json.loads(Path(str(prefix) + ".metadata.json").read_text(encoding="utf-8"))
        self.assertEqual(payload["engine"], "axis threshold-rank Newman-Ziff")
        self.assertEqual(payload["L"], 2)
        self.assertEqual(payload["N"], 4)
        self.assertEqual(payload["git_commit"], "abc123")


if __name__ == "__main__":
    unittest.main()
