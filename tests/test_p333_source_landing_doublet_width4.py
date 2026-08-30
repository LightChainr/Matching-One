import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from p333_source_landing_doublet_width4 import (  # noqa: E402
    build_result,
    landing_emission_vector,
    landing_reference_state,
    landing_rotation,
)


RESULT = ROOT / "results/p333-source-landing-doublet-width4/latest.json"


class TestP333SourceLandingDoubletWidth4(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = json.loads(RESULT.read_text(encoding="utf-8"))

    def test_landing_reference_states(self):
        self.assertEqual(
            [landing_reference_state(site) for site in range(4)],
            [(0, 1, 1, 1), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)],
        )

    def test_c4_emission_orbit(self):
        rotation = landing_rotation()
        vectors = [landing_emission_vector(site) for site in range(4)]
        for site, vector in enumerate(vectors):
            rotated = tuple(
                sum(rotation[row][pivot] * vector[pivot] for pivot in range(2))
                for row in range(2)
            )
            self.assertEqual(rotated, vectors[(site + 1) % 4])
        self.assertEqual(vectors, [(1, 0), (0, 1), (-1, 0), (0, -1)])

    def test_non_scalar_and_translation_gates(self):
        self.assertTrue(self.result["non_scalar_gate"]["passes"])
        self.assertEqual(self.result["non_scalar_gate"]["new_non_scalar_rank"], 2)
        self.assertEqual(self.result["translation_covariance"]["G0_residual_rank"], 0)
        self.assertEqual(self.result["translation_covariance"]["G1_residual_rank"], 0)

    def test_dimension_ladder_and_exact_failure(self):
        stages = self.result["stages"]
        self.assertEqual(stages["affine_q_jet"]["affine_tangent_dimension"], 4)
        self.assertEqual(stages["endpoint_radical_normalized"]["affine_tangent_dimension"], 2)
        self.assertFalse(stages["gram_self_adjoint"]["consistent"])
        witness = self.result["first_empty_restriction"]
        self.assertEqual(witness["added_stage"], "gram_self_adjoint")
        self.assertEqual(witness["coefficient_rank"], 0)
        self.assertEqual(witness["augmented_rank"], 1)
        self.assertEqual(witness["inconsistency_witness"]["left_times_rhs"], 1)

    def test_canonical_residual_and_typed_lower_bound(self):
        self.assertEqual(self.result["canonical_restricted_gram_skew_rank"], 2)
        bound = self.result["next_type_lower_bound_if_failure"]
        self.assertEqual(bound["minimum_non_scalar_mark_dimension"], 3)
        self.assertIn("charge-two", bound["missing_irrep"])
        self.assertFalse(bound["tested"])

    def test_width_four_stops(self):
        self.assertEqual(self.result["decision"], "width4_doublet_fails")
        self.assertIsNone(self.result["final_velocity"])

    def test_full_certificate_recomputes(self):
        self.assertEqual(build_result(), self.result)


if __name__ == "__main__":
    unittest.main()
