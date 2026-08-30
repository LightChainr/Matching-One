import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from p333_source_landing_doublet import (  # noqa: E402
    build_result,
    landing_emission_vector,
    landing_pair_state,
    landing_rotation,
    source_landing_detach_jet,
)


RESULT = ROOT / "results/p333-source-landing-doublet/latest.json"


class TestP333SourceLandingDoublet(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = json.loads(RESULT.read_text(encoding="utf-8"))

    def test_landing_pair_states(self):
        self.assertEqual(
            [landing_pair_state(site) for site in range(3)],
            [(0, 1, 1), (0, 1, 0), (0, 0, 1)],
        )

    def test_emission_vectors_form_c3_orbit(self):
        rotation = landing_rotation()
        vectors = [landing_emission_vector(site) for site in range(3)]
        for site in range(3):
            column = tuple((vectors[site][0], vectors[site][1]))
            rotated = tuple(
                sum(rotation[row][pivot] * column[pivot] for pivot in range(2))
                for row in range(2)
            )
            self.assertEqual(rotated, vectors[(site + 1) % 3])
        self.assertEqual(
            tuple(sum(vector[index] for vector in vectors) for index in range(2)),
            (0, 0),
        )

    def test_detach_mark_block_is_genuine(self):
        _, velocity = source_landing_detach_jet(0)
        mark_block = velocity[-2:]
        self.assertTrue(any(value for row in mark_block for value in row))
        self.assertEqual(mark_block[0][-2:], (0, 0))
        self.assertEqual(mark_block[1][-2:], (0, 0))

    def test_non_scalar_and_translation_gates(self):
        self.assertTrue(self.result["non_scalar_gate"]["passes"])
        self.assertEqual(self.result["non_scalar_gate"]["new_non_scalar_rank"], 2)
        self.assertEqual(self.result["translation_covariance"]["G0_residual_rank"], 0)
        self.assertEqual(self.result["translation_covariance"]["G1_residual_rank"], 0)

    def test_exact_dimension_ladder(self):
        stages = self.result["stages"]
        self.assertEqual(
            [
                stages[name]["affine_tangent_dimension"]
                for name in (
                    "affine_q_jet",
                    "endpoint_radical_normalized",
                    "gram_self_adjoint",
                    "source_landing_normalized",
                )
            ],
            [4, 2, 2, 0],
        )
        self.assertTrue(all(stage["consistent"] for stage in stages.values()))

    def test_unique_canonical_zero_velocity(self):
        self.assertEqual(self.result["decision"], "doublet_breaks_obstruction_unique")
        self.assertEqual(self.result["canonical_restricted_gram_skew_rank"], 0)
        self.assertTrue(self.result["final_velocity"]["particular_is_zero"])
        self.assertTrue(self.result["final_velocity"]["all_tangent_velocities_zero"])

    def test_full_certificate_recomputes(self):
        self.assertEqual(build_result(), self.result)


if __name__ == "__main__":
    unittest.main()
