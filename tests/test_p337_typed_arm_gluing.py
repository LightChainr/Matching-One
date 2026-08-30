import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p337_typed_arm_gluing import (  # noqa: E402
    DEFAULT_OUTPUT,
    build_result,
    local_template,
    ordinary_six_arm_counterexample,
)


class TypedArmGluingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = build_result()

    def test_ordinary_six_arm_counterexample_is_contractible(self):
        row = ordinary_six_arm_counterexample()
        self.assertEqual(row["before_components"], 3)
        self.assertEqual(row["before_rank"], 0)
        self.assertEqual(row["after_adding_origin_components"], 1)
        self.assertEqual(row["after_adding_origin_rank"], 0)
        self.assertEqual(len(row["open_arms"]), 3)
        self.assertEqual(len(row["vacant_matching_arms"]), 3)

    def test_fixed_templates_have_scale_independent_cost(self):
        theta = local_template("theta")
        figure = local_template("figure_eight")
        self.assertEqual((theta["forced_open"], theta["forced_closed"]), (6, 18))
        self.assertEqual((figure["forced_open"], figure["forced_closed"]), (8, 16))
        self.assertEqual(theta["forced_total"], figure["forced_total"])
        self.assertEqual(theta["uniform_finite_energy_bound"], "eta^24 for p in [eta,1-eta]")

    def test_small_quotient_fixtures_realize_both_types(self):
        fixtures = self.result["small_quotient_fixtures"]
        self.assertEqual(fixtures["theta"]["descriptor"]["topology"], "one_carrier_theta")
        self.assertEqual(
            fixtures["figure_eight"]["descriptor"]["topology"],
            "two_carrier_figure_eight",
        )

    def test_minimal_global_fields_are_not_ordinary_arm_fields(self):
        fields = self.result["minimal_typed_fields"]["theta"]
        self.assertIn("ambient_rank_zero", fields)
        self.assertTrue(any("component_id" in field for field in fields))
        self.assertTrue(any("deck_address" in field for field in fields))

    def test_checked_in_result_reproduces(self):
        self.assertEqual(json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8")), self.result)


if __name__ == "__main__":
    unittest.main()
