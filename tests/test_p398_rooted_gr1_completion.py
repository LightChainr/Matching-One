import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from p398_rooted_gr1_completion import (  # noqa: E402
    ROOTED_SEEDS,
    build_result,
    orbit,
    selected_completion_families,
)


RESULT = ROOT / "results/p398-rooted-gr1-completion/latest.json"


class TestP398RootedGr1Completion(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = json.loads(RESULT.read_text(encoding="utf-8"))

    def test_frozen_rooted_orbits(self):
        self.assertEqual(
            {name: len(orbit(seed)) for name, seed in ROOTED_SEEDS.items()},
            {"AP": 4, "OP": 2, "DP": 2},
        )
        self.assertEqual(orbit(ROOTED_SEEDS["AP"])[0], (0, 0, 1, 2))
        self.assertEqual(orbit(ROOTED_SEEDS["OP"])[0], (0, 1, 0, 2))
        self.assertEqual(orbit(ROOTED_SEEDS["DP"])[0], (0, 0, 1, 1))

    def test_exactly_seven_projected_coordinates(self):
        families = selected_completion_families()
        self.assertEqual(
            {name: len(family["mark_action"]) for name, family in families.items()},
            {"rooted_trivial": 2, "rooted_charge1": 2, "rooted_charge2": 3},
        )
        self.assertEqual(self.result["new_coordinate_count"], 7)
        self.assertEqual(self.result["dimension_lower_bound"], 7)

    def test_sector_deficits_are_filled_exactly(self):
        expected = {
            "trivial": (3, 2, 5),
            "charge1_rational": (2, 2, 4),
            "charge2": (1, 3, 4),
        }
        for sector, (old_rank, increment, combined) in expected.items():
            row = self.result["sector_completion"][sector]
            self.assertEqual(row["old_rank"], old_rank)
            self.assertEqual(row["incremental_rank"], increment)
            self.assertEqual(row["combined_rank"], combined)
            self.assertTrue(row["complete"])
            self.assertIsNone(row["left_null_counterexample"])

    def test_combined_response_is_a_basis_of_grade_one(self):
        combined = self.result["combined"]
        self.assertEqual(combined["old_rank"], 6)
        self.assertEqual(combined["new_raw_rank"], 7)
        self.assertEqual(combined["old_plus_new_rank"], 13)
        self.assertEqual(combined["old_plus_new_B_coordinate_determinant"], 3072)
        self.assertTrue(combined["full_grade_one"])

    def test_all_new_families_are_translation_covariant(self):
        for projection in self.result["new_mark_projections"].values():
            self.assertEqual(projection["translation_covariance"]["raw_response_residual_rank"], 0)
            self.assertEqual(projection["translation_covariance"]["H_dual_state_residual_rank"], 0)

    def test_minimality_and_claim_boundary(self):
        self.assertEqual(self.result["decision"], "minimal_exact_completion")
        self.assertTrue(self.result["minimality"]["proved"])
        self.assertIn("not yet an affine", self.result["claim_boundary"][1])
        self.assertIn("Promote this seven-coordinate rooted registry", self.result["next_gate"])

    def test_full_certificate_recomputes(self):
        self.assertEqual(build_result(), self.result)


if __name__ == "__main__":
    unittest.main()
