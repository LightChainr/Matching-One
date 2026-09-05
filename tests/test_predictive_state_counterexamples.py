from __future__ import annotations

from collections import Counter
from fractions import Fraction
import json
from math import gcd
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from predictive_state_counterexamples import (  # noqa: E402
    N10BirthOracle,
    build_artifact,
    conditional_arm_arithmetic,
    direct_birth_priority_control,
    exact_determinant,
    exact_rank,
    generator_controls,
    hankel,
    hilbert_cauchy_determinant,
    one_clock_covariance,
    one_clock_controls,
)


def union_find_mark(mask: int):
    """Independent weighted-DSU winding check using only positive NN edges."""
    parent = list(range(10))
    weight = [(0, 0) for _ in range(10)]
    cycles = []

    def find(vertex):
        x = y = 0
        while parent[vertex] != vertex:
            dx, dy = weight[vertex]
            x, y = x+dx, y+dy
            vertex = parent[vertex]
        return vertex, (x, y)

    for u in range(10):
        if not (mask >> u) & 1:
            continue
        for offset, delta in ((3, (1, 0)), (1, (0, 1))):
            v = (u+offset) % 10
            if not (mask >> v) & 1:
                continue
            ru, du = find(u)
            rv, dv = find(v)
            difference = (du[0]+delta[0]-dv[0], du[1]+delta[1]-dv[1])
            if ru != rv:
                parent[rv] = ru
                weight[rv] = difference
            elif difference != (0, 0):
                x, y = difference
                numerator = (3*x+y, -x+3*y)
                if numerator[0] % 10 or numerator[1] % 10:
                    raise AssertionError("DSU cycle does not close on the declared torus")
                cycles.append((numerator[0]//10, numerator[1]//10))
    if not cycles:
        return 0, None
    x, y = cycles[0]
    if any(x*b-y*a for a, b in cycles[1:]):
        return 2, None
    scale = gcd(abs(x), abs(y))
    line = (x//scale, y//scale)
    if line[0] < 0 or (line[0] == 0 and line[1] < 0):
        line = (-line[0], -line[1])
    return 1, line


class GeneratorControlsTests(unittest.TestCase):
    def test_hilbert_determinants_match_independent_product(self):
        for size in range(1, 11):
            matrix = hankel(lambda n: Fraction(1, n+1), size)
            determinant = exact_determinant(matrix)
            self.assertEqual(determinant, hilbert_cauchy_determinant(size))
            self.assertGreater(determinant, 0)
            self.assertEqual(exact_rank(matrix), size)

    def test_other_rational_power_has_unbounded_finite_ranks(self):
        for size in range(1, 8):
            matrix = hankel(lambda n: Fraction(1, (n+2)**2), size)
            self.assertEqual(exact_rank(matrix), size)
            self.assertGreater(exact_determinant(matrix), 0)

    def test_geometric_power_and_log_pair(self):
        for size in range(1, 11):
            self.assertEqual(exact_rank(hankel(lambda n: Fraction(1, 2**n), size)), 1)
            self.assertEqual(exact_rank(hankel(lambda n: Fraction(n+1, 2**n), size)), min(size, 2))

    def test_log_pair_has_exact_repeated_root_recurrence(self):
        g = lambda n: Fraction(n+1, 2**n)
        for n in range(20):
            self.assertEqual(g(n+2)-g(n+1)+Fraction(1, 4)*g(n), 0)
        self.assertNotEqual(g(1)-Fraction(1, 2)*g(0), 0)

    def test_variable_not_constant_coefficient_closure(self):
        g = lambda n: Fraction(1, n+1)
        for n in range(20):
            self.assertEqual((n+2)*g(n+1)-(n+1)*g(n), 0)
        self.assertEqual(g(2)-(g(1)/g(0))*g(1), Fraction(1, 12))
        self.assertEqual(generator_controls(3)["false_constant_rank_one_first_heldout_residual"], "1/12")

    def test_exact_linear_algebra_edge_cases(self):
        self.assertEqual(exact_rank([]), 0)
        self.assertEqual(exact_rank([[], []]), 0)
        self.assertEqual(exact_rank([[1, 2, 3], [2, 4, 6]]), 1)
        self.assertEqual(exact_determinant([[0, 1], [2, 0]]), -2)
        self.assertEqual(exact_determinant([[1, 2], [2, 4]]), 0)
        self.assertEqual(exact_determinant([]), 1)
        with self.assertRaises(ValueError):
            exact_rank([[1], [2, 3]])
        with self.assertRaises(ValueError):
            exact_determinant([[1, 2]])

    def test_single_clock_covariance_is_full_rank(self):
        for row in one_clock_controls()["rows"]:
            self.assertEqual(row["size"], row["covariance_rank"])
            self.assertGreater(Fraction(row["determinant"]), 0)

    def test_single_clock_determinant_is_product_of_quantile_gaps(self):
        values = [Fraction(1, 10), Fraction(1, 3), Fraction(4, 5)]
        gaps = [values[0]] + [b-a for a, b in zip(values, values[1:])] + [1-values[-1]]
        product = Fraction(1)
        for gap in gaps:
            product *= gap
        self.assertEqual(exact_determinant(one_clock_covariance(values)), product)
        for invalid in ([], [0], [1], [Fraction(1, 2)]*2, [Fraction(2, 3), Fraction(1, 3)]):
            with self.assertRaises(ValueError):
                one_clock_covariance(invalid)

    def test_invalid_sizes_fail(self):
        for size in (0, -1, True, 1.5):
            with self.assertRaises(ValueError):
                hankel(lambda n: Fraction(1), size)
            with self.assertRaises(ValueError):
                generator_controls(size)


class BirthControlsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.oracle = N10BirthOracle()

    def test_bfs_matches_independent_dsu_on_every_configuration(self):
        for mask in range(1024):
            self.assertEqual(self.oracle.marks[mask], union_find_mark(mask), msg=str(mask))

    def test_exact_rank_census(self):
        self.assertEqual(Counter(rank for rank, _ in self.oracle.marks), {0: 548, 1: 310, 2: 166})
        self.assertEqual(self.oracle.marks[0], (0, None))
        self.assertEqual(self.oracle.marks[1023], (2, None))

    def test_equal_one_step_state_different_two_step_survival(self):
        for mask in (31, 47):
            self.assertEqual(self.oracle.marks[mask], (1, (0, 1)))
            self.assertEqual(len(self.oracle.vacant(mask)), 5)
            self.assertEqual(len(self.oracle.triggers(mask)[0]), 1)
            self.assertEqual(self.oracle.survival(mask, 1), Fraction(4, 5))
        self.assertEqual(self.oracle.survival(31, 2), Fraction(2, 5))
        self.assertEqual(self.oracle.survival(47, 2), Fraction(3, 10))

    def test_exact_trigger_pairs(self):
        self.assertEqual(self.oracle.triggers(31), ([7], [(5, 8), (6, 9)]))
        self.assertEqual(self.oracle.triggers(47), ([8], [(4, 7), (6, 7), (6, 9)]))

    def test_independent_future_permutation_counts(self):
        self.assertEqual(self.oracle.exit_counts(31), {1: 24, 2: 48, 3: 48})
        self.assertEqual(self.oracle.exit_counts(47), {1: 24, 2: 60, 3: 36})

    def test_subset_and_killed_kernel_for_all_rank_one_states(self):
        checked = 0
        for mask, (rank, _) in enumerate(self.oracle.marks):
            if rank == 1:
                for steps in range(len(self.oracle.vacant(mask))+1):
                    self.assertEqual(self.oracle.survival(mask, steps), self.oracle.subset_survival(mask, steps))
                    checked += 1
        self.assertEqual(checked, 1650)

    def test_direct_priority_integral_and_permutation_dp(self):
        result = direct_birth_priority_control(self.oracle)
        self.assertEqual(result["directed_edge_counts_by_k"], {"5": 40, "6": 40})
        self.assertEqual(result["directed_edge_count"], 80)
        self.assertEqual(result["permutation_paths_with_direct_birth"], 288000)
        self.assertEqual(result["total_permutation_paths"], 3628800)
        self.assertEqual(result["probability_from_beta_weights"], "5/63")
        self.assertEqual(result["probability_from_permutation_dp"], "5/63")

    def test_invalid_inputs_and_killing(self):
        self.oracle.survival(31, 1)  # Cache must not bypass bool/type validation.
        for mask in (-1, 1024, True, 31.0):
            with self.assertRaises(ValueError):
                self.oracle.survival(mask, 1)
        for steps in (-1, True, 1.5):
            with self.assertRaises(ValueError):
                self.oracle.survival(31, steps)
            with self.assertRaises(ValueError):
                self.oracle.subset_survival(31, steps)
        self.assertEqual(self.oracle.survival(1023, 0), 0)
        self.assertEqual(self.oracle.survival(31, 6), 0)
        with self.assertRaises(ValueError):
            self.oracle.triggers(0)


class ArtifactTests(unittest.TestCase):
    def test_checked_artifact_reproduces(self):
        checked = json.loads((ROOT / "results/exact-predictive-state-controls/20260830.json").read_text())
        self.assertEqual(build_artifact(), checked)
        self.assertEqual(checked["issues"], [400, 403, 405])
        self.assertIn("no Monte Carlo", checked["data_class"])

    def test_arm_arithmetic_remains_explicitly_conditional(self):
        result = conditional_arm_arithmetic()
        self.assertEqual(result["rows"], [
            {"arms": 6, "arm_exponent": "35/12", "area_decay_exponent": "5/6"},
            {"arms": 8, "arm_exponent": "21/4", "area_decay_exponent": "2"},
        ])
        self.assertEqual(result["relative_area_decay_exponent"], "7/6")
        self.assertIn("not an arm correspondence theorem", result["claim"])


if __name__ == "__main__":
    unittest.main()
