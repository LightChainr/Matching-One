from __future__ import annotations

import copy
import json
import unittest

from scripts.gadget_complement_involution_certificate import (
    DEFAULT_OUTPUT,
    build_artifact,
    build_row,
    complement_commutes_with_relabeling,
    complement_graph,
    validate_artifact,
)
from scripts.gadget_graph_canonical import enumerate_graphs
from scripts.terminal_partition_canonical import full_symmetric_group


class GadgetComplementInvolutionCertificateTests(unittest.TestCase):
    def test_checked_artifact_reproduces(self) -> None:
        checked = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(checked, build_artifact())
        self.assertEqual(
            validate_artifact(checked)["status"],
            "valid_exact_complement_involution_certificate",
        )

    def test_labeled_complement_is_an_involution(self) -> None:
        for terminal_count in (3, 4):
            vertex_count = terminal_count + 1
            for graph in enumerate_graphs(vertex_count, terminal_count):
                self.assertEqual(
                    complement_graph(vertex_count, complement_graph(vertex_count, graph)),
                    graph,
                )

    def test_complement_commutes_with_every_audited_relabeling(self) -> None:
        for terminal_count in (3, 4):
            vertex_count = terminal_count + 1
            for graph in enumerate_graphs(vertex_count, terminal_count):
                for terminal_map in full_symmetric_group(terminal_count):
                    vertex_map = tuple(terminal_map) + (terminal_count,)
                    self.assertTrue(
                        complement_commutes_with_relabeling(
                            vertex_count, terminal_count, graph, vertex_map
                        )
                    )

    def test_exact_pair_counts(self) -> None:
        three = build_row(3)
        four = build_row(4)
        self.assertEqual((three["self_complementary_orbits"], three["complement_pairs"]), (0, 10))
        self.assertEqual((four["self_complementary_orbits"], four["complement_pairs"]), (2, 44))
        self.assertEqual(three["orbit_accounting"], 20)
        self.assertEqual(four["orbit_accounting"], 90)

    def test_histograms_are_palindromic(self) -> None:
        for terminal_count in (3, 4):
            histogram = build_row(terminal_count)["orbit_edge_count_histogram"]
            self.assertEqual(histogram, list(reversed(histogram)))

    def test_invalid_edges_and_scope_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "invalid edge"):
            complement_graph(3, ((0, 3),))
        with self.assertRaisesRegex(ValueError, "three or four"):
            build_row(2)

    def test_tampering_fails_closed(self) -> None:
        tampered = copy.deepcopy(build_artifact())
        tampered["rows"][1]["self_complementary_orbits"] = 0
        with self.assertRaisesRegex(ValueError, "does not exactly reproduce"):
            validate_artifact(tampered)


if __name__ == "__main__":
    unittest.main()
