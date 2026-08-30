import json
from pathlib import Path
import unittest

from p334_dual_hazard_ulc import build_result


ROOT = Path(__file__).resolve().parents[1]


class P334DualHazardULCTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = build_result()

    def test_complement_swaps_birth_and_exit_without_failure(self):
        audit = self.result["complement_duality"]
        self.assertEqual(audit["matrices_checked"], 83)
        self.assertEqual(audit["rank_one_states_checked"], 59922)
        self.assertEqual(audit["state_line_failures"], 0)
        self.assertEqual(audit["birth_exit_swap_failures"], 0)
        self.assertTrue(audit["all_pass"])

    def test_dual_hazard_criterion_passes_entire_bounded_atlas(self):
        for carrier in ("primal", "matching"):
            audit = self.result["bounded_exact_audit"][carrier]
            self.assertEqual(audit["fixed_line_sequences"], 240)
            self.assertEqual(audit["adjacent_layer_pairs"], 492)
            self.assertEqual(audit["ulc_comparisons"], 346)
            self.assertEqual(audit["exit_hazard_nondecreasing"], 240)
            self.assertEqual(audit["birth_hazard_nonincreasing"], 240)
            self.assertEqual(audit["q_ratio_strictly_decreasing"], 240)
            self.assertEqual(audit["pointwise_exit_nesting_failures"], 0)
            self.assertEqual(audit["pointwise_birth_nesting_failures"], 0)

    def test_degree_bias_is_the_real_unproved_term(self):
        primal = self.result["bounded_exact_audit"]["primal"]
        matching = self.result["bounded_exact_audit"]["matching"]
        self.assertEqual(primal["negative_exit_bias_corrections"], 208)
        self.assertEqual(matching["negative_exit_bias_corrections"], 276)
        worst = self.result["worst_negative_exit_degree_bias"]
        self.assertEqual(worst["correction"], "-1/12")
        self.assertEqual(worst["exit_edge_slack"], "5/12")
        self.assertEqual(worst["exit_uniform_hazard_delta"], "1/3")

    def test_current_and_ratio_identities_are_exercised(self):
        examples = self.result["representative_sequences"]
        self.assertTrue(examples)
        for example in examples:
            audit = example["audit"]
            self.assertTrue(audit["current_rows"])
            self.assertTrue(audit["ratio_rows"])
            self.assertTrue(
                all(
                    row["q_ratio"] == row["dual_hazard_ratio"]
                    for row in audit["ratio_rows"]
                )
            )

    def test_checked_in_artifact_reproduces(self):
        checked = json.loads(
            (ROOT / "results/p334-dual-hazard-ulc/latest.json").read_text()
        )
        self.assertEqual(checked, self.result)


if __name__ == "__main__":
    unittest.main()
