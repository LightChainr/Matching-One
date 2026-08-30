import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from noncrossing_connectivity_codec import noncrossing_states  # noqa: E402
from p333_gram_source_intertwiner import join_block_count  # noqa: E402
from p398_qadic_jantzen import build_result, gram_coefficient  # noqa: E402


RESULT = ROOT / "results/p398-qadic-jantzen/latest.json"


class TestP398QAdicJantzen(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = json.loads(RESULT.read_text(encoding="utf-8"))
        cls.by_width = {row["width"]: row for row in cls.result["widths"]}

    def test_exact_polynomial_coefficients(self):
        for width in (3, 4):
            states = noncrossing_states(width)
            coefficients = [gram_coefficient(states, order) for order in range(width + 1)]
            for row, left in enumerate(states):
                for column, right in enumerate(states):
                    blocks = join_block_count(left, right)
                    self.assertEqual(
                        sum(matrix[row][column] for matrix in coefficients),
                        2 ** blocks,
                    )

    def test_invariant_valuations_and_filtration(self):
        self.assertEqual(self.by_width[3]["local_invariant_factor_valuations"], [0, 1, 1, 1, 1])
        self.assertEqual(self.by_width[4]["local_invariant_factor_valuations"], [0] + [1] * 13)
        self.assertEqual(self.by_width[3]["jantzen_filtration_dimensions"], {"J0": 5, "J1": 4, "J2": 0})
        self.assertEqual(self.by_width[4]["jantzen_filtration_dimensions"], {"J0": 14, "J1": 13, "J2": 0})

    def test_leading_radical_forms_are_unimodular(self):
        for width in (3, 4):
            leading = self.by_width[width]["leading_radical_form"]
            self.assertEqual(leading["rank"], self.by_width[width]["radical_dimension"])
            self.assertEqual(leading["determinant"], -1)
            self.assertTrue(leading["unimodular"])
            self.assertFalse(self.by_width[width]["higher_coefficients_needed"])

    def test_translation_sector_decomposition(self):
        self.assertEqual(
            self.by_width[3]["grade_one_sector_dimensions"],
            {"trivial": 2, "charge1_rational": 2},
        )
        self.assertEqual(
            self.by_width[4]["grade_one_sector_dimensions"],
            {"trivial": 5, "charge1_rational": 4, "charge2": 4},
        )

    def test_mark_projection_ranks_and_coverage(self):
        self.assertEqual(self.by_width[3]["combined_tested_projection_rank"], 4)
        self.assertEqual(self.by_width[3]["combined_uncovered_dimension"], 0)
        self.assertEqual(self.by_width[4]["combined_tested_projection_rank"], 6)
        self.assertEqual(self.by_width[4]["combined_uncovered_dimension"], 7)
        self.assertEqual(
            {name: row["tested_projection_rank"] for name, row in self.by_width[4]["sector_coverage"].items()},
            {"trivial": 3, "charge1_rational": 2, "charge2": 1},
        )

    def test_all_mark_projections_are_translation_covariant(self):
        for width in (3, 4):
            for projection in self.by_width[width]["mark_projections"].values():
                self.assertGreater(projection["grade_one_functional_rank"], 0)
                self.assertEqual(projection["translation_covariance"]["raw_response_residual_rank"], 0)
                self.assertEqual(projection["translation_covariance"]["H_dual_state_residual_rank"], 0)

    def test_nilpotent_types_are_not_identified(self):
        separation = self.result["nilpotent_separation"]
        self.assertTrue(separation["base_parameter_nilpotent"]["square_zero"])
        self.assertEqual(separation["base_parameter_nilpotent"]["rank_after_Q_equals_1_specialization"], 0)
        self.assertTrue(separation["fixed_Q_extension"]["not_implied_by_nonzero_Jantzen_projection"])

    def test_full_certificate_recomputes(self):
        self.assertEqual(build_result(), self.result)


if __name__ == "__main__":
    unittest.main()
