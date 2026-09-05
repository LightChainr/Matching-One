
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


class ThresholdRankAxisPairCppTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        compiler = shutil.which("g++") or shutil.which("c++")
        if compiler is None:
            raise unittest.SkipTest("no C++ compiler")
        cls.temporary = tempfile.TemporaryDirectory()
        cls.binary = Path(cls.temporary.name) / "threshold_rank_axis_pair_mc"
        subprocess.run(
            [
                compiler,
                "-O2",
                "-std=c++17",
                str(ROOT / "src" / "threshold_rank_axis_pair_mc.cpp"),
                "-o",
                str(cls.binary),
            ],
            check=True,
        )

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temporary.cleanup()

    def test_self_test(self) -> None:
        completed = subprocess.run(
            [str(self.binary), "--self-test"],
            check=True,
            text=True,
            capture_output=True,
        )
        self.assertIn("exact permutation restriction is uniform", completed.stdout)

    @staticmethod
    def restrict_to_lower(permutation, L):
        lower_L = L - 1
        result = []
        for vertex in permutation:
            x = vertex % L
            y = vertex // L
            if x < lower_L and y < lower_L:
                result.append(x + lower_L * y)
        return tuple(result)

    @staticmethod
    def read_histogram(path: Path):
        result = {}
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                key = (int(row["batch"]), row["role"], row["kind"])
                result.setdefault(key, Counter())[int(row["k"])] = int(row["count"])
        return result

    def test_cpp_matches_python_for_both_exact_marginals(self) -> None:
        L = 3
        samples = 20
        batches = 2
        seed = 311
        offset = 700
        prefix = Path(self.temporary.name) / "pair_l3"
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

        upper_geometry = axis_integer_torus(L)
        lower_geometry = axis_integer_torus(L - 1)
        upper_engine = ThresholdRankEngine(upper_geometry)
        lower_engine = ThresholdRankEngine(lower_geometry)
        per_batch = samples // batches
        expected = {}
        for batch in range(batches):
            counters = {
                ("upper", "minus"): Counter(),
                ("upper", "plus"): Counter(),
                ("lower", "minus"): Counter(),
                ("lower", "plus"): Counter(),
            }
            begin = offset + batch * per_batch
            for replica in range(begin, begin + per_batch):
                upper = counter_permutation(upper_geometry.n, seed, replica)
                lower = self.restrict_to_lower(upper, L)
                self.assertEqual(set(lower), set(range(lower_geometry.n)))
                ku_minus, ku_plus = upper_engine.threshold_ranks(upper)
                kl_minus, kl_plus = lower_engine.threshold_ranks(lower)
                counters[("upper", "minus")][ku_minus] += 1
                counters[("upper", "plus")][ku_plus] += 1
                counters[("lower", "minus")][kl_minus] += 1
                counters[("lower", "plus")][kl_plus] += 1
            for (role, kind), counts in counters.items():
                expected[(batch, role, kind)] = counts

        self.assertEqual(observed, expected)

    def test_restriction_of_uniform_permutation_has_uniform_small_order(self) -> None:
        # Exact combinatorial oracle independent of the Monte Carlo engine.
        import itertools
        orders = Counter()
        for permutation in itertools.permutations(range(4)):
            orders[tuple(value for value in permutation if value != 3)] += 1
        self.assertEqual(len(orders), 6)
        self.assertEqual(set(orders.values()), {4})


if __name__ == "__main__":
    unittest.main()
