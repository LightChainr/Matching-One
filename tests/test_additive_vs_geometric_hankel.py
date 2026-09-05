import importlib.util
import json
from fractions import Fraction
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts" / "additive_vs_geometric_hankel.py"
ARTIFACT = ROOT / "analysis" / "additive_vs_geometric_hankel_certificate.json"
SPEC = importlib.util.spec_from_file_location("additive_vs_geometric_hankel", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class AdditiveVsGeometricHankelTests(unittest.TestCase):
    def test_committed_artifact_regenerates_exactly(self):
        self.assertEqual(json.loads(ARTIFACT.read_text()), MODULE.build_certificate())

    def test_hilbert_determinants_match_closed_formula(self):
        sequence = MODULE.additive_power_sequence(15)
        for size in range(1, 8):
            with self.subTest(size=size):
                self.assertEqual(MODULE.exact_determinant(MODULE.hankel(sequence, size)), MODULE.hilbert_determinant(size))

    def test_additive_hankel_ranks_are_full(self):
        self.assertEqual(MODULE.ranks_through(MODULE.additive_power_sequence(15), 8), list(range(1, 9)))

    def test_variable_coefficient_recurrence_is_exact(self):
        sequence = MODULE.additive_power_sequence(33)
        residuals = [(n + 2) * sequence[n + 1] - (n + 1) * sequence[n] for n in range(32)]
        self.assertEqual(residuals, [Fraction(0)] * 32)

    def test_geometric_power_has_rank_one(self):
        self.assertEqual(MODULE.ranks_through(MODULE.geometric_power_sequence(15), 8), [1] * 8)

    def test_logarithmic_partner_has_rank_two(self):
        self.assertEqual(MODULE.ranks_through(MODULE.logarithmic_partner_sequence(15), 8), [1] + [2] * 7)

    def test_three_distinct_modes_have_rank_three(self):
        sequence = MODULE.exponential_sum_sequence(15, ("1/2", "1/3", "1/5"))
        self.assertEqual(MODULE.ranks_through(sequence, 8), [1, 2] + [3] * 6)

    def test_deliberately_wrong_rank_one_additive_claim_fails(self):
        sequence = MODULE.additive_power_sequence(5)
        self.assertNotEqual(MODULE.exact_determinant(MODULE.hankel(sequence, 3)), 0)
        self.assertEqual(MODULE.exact_rank(MODULE.hankel(sequence, 3)), 3)


if __name__ == "__main__":
    unittest.main()
