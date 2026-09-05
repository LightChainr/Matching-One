
from __future__ import annotations
import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import pc_block_event_exact as oracle  # noqa: E402


class PcBlockEventExactTest(unittest.TestCase):
    def test_open_boundary_edge_sets(self) -> None:
        self.assertEqual(len(oracle.graph_edges(1, "square")), 1)
        self.assertEqual(len(oracle.graph_edges(1, "matching")), 1)
        self.assertEqual(len(oracle.graph_edges(2, "square")), 10)
        self.assertEqual(len(oracle.graph_edges(2, "matching")), 16)

    def test_s1_requires_both_sites(self) -> None:
        for graph in oracle.GRAPHS:
            self.assertEqual(oracle.decide_event(0, 1, graph).reason, "left_empty")
            self.assertFalse(oracle.decide_event(1, 1, graph).success)
            self.assertFalse(oracle.decide_event(2, 1, graph).success)
            self.assertTrue(oracle.decide_event(3, 1, graph).success)

    def test_diagonal_tie_distinguishes_square_and_matching(self) -> None:
        # s=2 row-major indices: (1,0)=1, (2,0)=2, (0,1)=4.
        base = (1 << 1) | (1 << 2)
        enlarged = base | (1 << 4)
        self.assertTrue(oracle.decide_event(base, 2, "square").success)
        square = oracle.decide_event(enlarged, 2, "square")
        matching = oracle.decide_event(enlarged, 2, "matching")
        self.assertFalse(square.success)
        self.assertEqual(square.reason, "left_tie")
        self.assertTrue(matching.success)

    def test_disconnected_selected_clusters_fail(self) -> None:
        # Open opposite outer corners only.
        mask = (1 << 0) | (1 << 3)
        for graph in oracle.GRAPHS:
            decision = oracle.decide_event(mask, 2, graph)
            self.assertFalse(decision.success)
            self.assertEqual(decision.reason, "largest_clusters_disconnected")

    def test_square_s2_reliability_polynomial(self) -> None:
        result = oracle.enumerate_reliability(2, "square")
        self.assertEqual(result["success_by_occupied"], [0, 0, 2, 8, 19, 24, 20, 8, 1])
        self.assertEqual(result["success_count"], 82)
        self.assertEqual(result["probability_at_half"], "41/128")

    def test_matching_s2_reliability_polynomial(self) -> None:
        result = oracle.enumerate_reliability(2, "matching")
        self.assertEqual(result["success_by_occupied"], [0, 0, 4, 20, 41, 44, 26, 8, 1])
        self.assertEqual(result["success_count"], 144)
        self.assertEqual(result["probability_at_half"], "9/16")

    def test_unique_largest_event_is_not_assumed_monotone(self) -> None:
        square = oracle.single_site_addition_counterexamples(2, "square")
        matching = oracle.single_site_addition_counterexamples(2, "matching")
        self.assertEqual(len(square), 28)
        self.assertEqual(len(matching), 0)
        self.assertIn({"mask": 6, "opened_index": 4, "enlarged_mask": 22}, square)

    def test_all_enumerated_reasons_partition_configuration_space(self) -> None:
        for graph in oracle.GRAPHS:
            result = oracle.enumerate_reliability(2, graph)
            self.assertEqual(sum(result["reason_counts"].values()), 256)
            self.assertEqual(result["reason_counts"]["success"], result["success_count"])

    def test_checked_in_results_are_reproducible(self) -> None:
        artifact = oracle.build_artifact()
        checked_json = json.loads(
            (ROOT / "results/pc-block-event-exact/latest.json").read_text(encoding="utf-8")
        )
        checked_markdown = (
            ROOT / "results/pc-block-event-exact/latest.md"
        ).read_text(encoding="utf-8")
        self.assertEqual(checked_json, artifact)
        self.assertEqual(checked_markdown, oracle.render_markdown(artifact))


if __name__ == "__main__":
    unittest.main()
