from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from q1_matching_intertwiner import build_oracle  # noqa: E402


class Q1MatchingIntertwinerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = build_oracle()

    def test_exchange_and_pullthrough_identity(self) -> None:
        site = self.result["site_matching_square_face"]
        self.assertTrue(site["bare_exchange"]["J_squared_is_identity"])
        identity = site["operator_identity"]
        self.assertTrue(identity["Ad_J_S_equals_S"])
        self.assertTrue(identity["Ad_J_D_equals_minus_D"])
        self.assertTrue(identity["identity_verified"])
        self.assertGreater(identity["nonzero_residual_rank"], 0)

    def test_planar_dual_has_exact_local_repair(self) -> None:
        control = self.result["edge_FK_positive_control"]
        self.assertTrue(control["adjusted_pullthrough_tangent_zero_configurationwise"])
        self.assertEqual(control["raw_defect_values_by_occupied_count"], {
            "0": [2], "1": [1], "2": [0], "3": [-1]
        })

    def test_site_matching_has_occupancy_counterterm_no_go(self) -> None:
        no_go = self.result["site_matching_square_face"]["occupancy_counterterm_no_go"]
        self.assertFalse(no_go["exists"])
        self.assertEqual(no_go["defect_values_by_occupied_count"]["2"], [0, 1])

    def test_equal_measure_score_witness_separates_derivatives(self) -> None:
        witnesses = self.result["site_matching_square_face"]["occupancy_counterterm_no_go"]["witnesses"]
        self.assertEqual(witnesses[0]["Bernoulli_measure_score_at_p_half"], 0)
        self.assertEqual(witnesses[1]["Bernoulli_measure_score_at_p_half"], 0)
        self.assertNotEqual(witnesses[0]["defect"], witnesses[1]["defect"])

    def test_even_odd_trace_selection_zero(self) -> None:
        selection = self.result["site_matching_square_face"]["trace_selection"]
        self.assertEqual(selection["uniform_doubled_trace_D"], "0")
        self.assertEqual(selection["uniform_doubled_trace_S_times_D"], "0")

    def test_committed_artifact_is_reproducible(self) -> None:
        committed = json.loads(
            (ROOT / "results" / "exact-q1-matching-intertwiner" / "latest.json").read_text()
        )
        self.assertEqual(committed, self.result)


if __name__ == "__main__":
    unittest.main()
