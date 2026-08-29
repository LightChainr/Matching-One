from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import pc_block_event_exact as exact  # noqa: E402
import pc_block_event_linear as linear  # noqa: E402


class PcBlockEventLinearTests(unittest.TestCase):
    def test_all_tiny_configurations_match_exact_oracle(self) -> None:
        comparisons = 0
        for graph in linear.GRAPHS:
            for s in (1, 2):
                vertices = exact.ordered_vertices(s)
                lookup = {vertex: index for index, vertex in enumerate(vertices)}
                for mask in range(1 << len(vertices)):
                    expected = exact.decide_event(mask, s, graph)
                    observed = linear.evaluate_open_flags(
                        s, graph, linear.flags_from_mask(mask, len(vertices))
                    )
                    self.assertEqual(observed.success, expected.success)
                    self.assertEqual(observed.reason, expected.reason)
                    self.assertEqual(
                        observed.left_largest,
                        tuple(sorted(lookup[v] for v in expected.left_largest)),
                    )
                    self.assertEqual(
                        observed.right_largest,
                        tuple(sorted(lookup[v] for v in expected.right_largest)),
                    )
                    comparisons += 1
        self.assertEqual(comparisons, 520)

    def test_invalid_inputs_fail_closed(self) -> None:
        with self.assertRaises(ValueError):
            linear.evaluate_open_flags(0, "square", ())
        with self.assertRaises(ValueError):
            linear.evaluate_open_flags(1, "triangular", (False, False))
        with self.assertRaises(ValueError):
            linear.evaluate_open_flags(2, "square", (False,) * 7)
        with self.assertRaises(ValueError):
            linear.flags_from_mask(4, 2)

    def test_edge_check_counts_match_open_rectangle_formulas(self) -> None:
        for s in (1, 2, 7):
            square = linear.evaluate_open_flags(s, "square", (False,) * linear.site_count(s))
            matching = linear.evaluate_open_flags(s, "matching", (False,) * linear.site_count(s))
            self.assertEqual(square.full_edge_checks, 4 * s * s - 3 * s)
            self.assertEqual(square.half_edge_checks, 4 * s * (s - 1))
            self.assertEqual(matching.full_edge_checks, 8 * s * s - 9 * s + 2)
            self.assertEqual(matching.half_edge_checks, 8 * s * s - 12 * s + 4)

    def test_large_all_open_controls(self) -> None:
        for graph in linear.GRAPHS:
            decision = linear.evaluate_open_flags(64, graph, (True,) * linear.site_count(64))
            self.assertTrue(decision.success)
            self.assertEqual(len(decision.left_largest), 64 * 64)
            self.assertEqual(len(decision.right_largest), 64 * 64)

    def test_nonmonotone_square_witness_is_preserved(self) -> None:
        before = linear.evaluate_open_flags(2, "square", linear.flags_from_mask(6, 8))
        after = linear.evaluate_open_flags(2, "square", linear.flags_from_mask(22, 8))
        matching_after = linear.evaluate_open_flags(2, "matching", linear.flags_from_mask(22, 8))
        self.assertTrue(before.success)
        self.assertEqual(after.reason, "left_tie")
        self.assertTrue(matching_after.success)

    def test_checked_in_results_reproduce(self) -> None:
        artifact = linear.build_artifact()
        checked_json = json.loads(
            (ROOT / "results/pc-block-event-linear/latest.json").read_text(encoding="utf-8")
        )
        checked_markdown = (
            ROOT / "results/pc-block-event-linear/latest.md"
        ).read_text(encoding="utf-8")
        self.assertEqual(checked_json, artifact)
        self.assertEqual(checked_markdown, linear.render_markdown(artifact))


if __name__ == "__main__":
    unittest.main()
