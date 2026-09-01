import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import modular_homology_channel_oracle as oracle  # noqa: E402


class ModularHomologyChannelOracleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        config = json.loads(
            (ROOT / "analysis" / "modular_homology_channel_manifest.json").read_text()
        )
        cls.result = oracle.analyze(config)

    def test_rank_channels_are_exactly_scalar(self):
        scalar = self.result["theorem"]["scalar_channels"]
        self.assertEqual(set(scalar), {"rank", "either", "cross"})
        self.assertTrue(
            self.result["finite_regression"]["all_rank_channels_preserved"]
        )
        self.assertEqual(self.result["finite_regression"]["case_count"], 42)

    def test_sl2_action_preserves_rank_on_bounded_vectors(self):
        matrices = [
            oracle.parse_matrix(((0, -1), (1, 0))),
            oracle.parse_matrix(((1, 1), (0, 1))),
            oracle.parse_matrix(((1, 0), (-1, 1))),
        ]
        vectors = [
            (x, y)
            for x in range(-3, 4)
            for y in range(-3, 4)
            if (x, y) != (0, 0)
        ]
        for matrix in matrices:
            for first in vectors:
                for second in vectors:
                    basis = (first, second)
                    self.assertEqual(
                        oracle.subgroup_rank(basis),
                        oracle.subgroup_rank(oracle.transform_basis(matrix, basis)),
                    )

    def test_direction_and_both_have_exact_counterexamples(self):
        first = self.result["counterexamples"]["direction_0_and_both"]
        second = self.result["counterexamples"]["direction_1_and_both"]
        self.assertEqual(first["determinant"], 1)
        self.assertEqual(second["determinant"], 1)
        self.assertEqual(first["changed_channels"], ["both", "direction_0"])
        self.assertEqual(second["changed_channels"], ["both", "direction_1"])
        for example in (first, second):
            self.assertEqual(example["flags_before"]["rank"], 1)
            self.assertEqual(example["flags_after"]["rank"], 1)
            self.assertTrue(example["flags_before"]["either"])
            self.assertTrue(example["flags_after"]["either"])
            self.assertFalse(example["flags_before"]["cross"])
            self.assertFalse(example["flags_after"]["cross"])

    def test_scalar_classification_lifts_to_matching_combinations(self):
        lift = self.result["theorem"]["matching_combination_lift"]
        self.assertEqual(set(lift), {"cross", "either"})
        for combinations in lift.values():
            self.assertEqual(
                set(combinations), {"primal", "matching", "even", "odd"}
            )
            self.assertEqual(set(combinations.values()), {"modular_scalar"})

    def test_elliptic_filter_kills_h4_h8_not_h12(self):
        rows = {
            row["spin"]: row
            for row in self.result["elliptic_spin_filter"]["candidates"]
        }
        self.assertFalse(rows[4]["hexagonal_rho_allowed"])
        self.assertFalse(rows[8]["hexagonal_rho_allowed"])
        self.assertTrue(rows[12]["hexagonal_rho_allowed"])
        self.assertTrue(rows[24]["square_lattice_and_hexagonal_allowed"])
        self.assertEqual(
            self.result["elliptic_spin_filter"]["intersection_period"], 12
        )

    def test_invalid_matrix_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "not in SL"):
            oracle.parse_matrix(((2, 0), (0, 1)))


if __name__ == "__main__":
    unittest.main()
