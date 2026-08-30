from __future__ import annotations

import copy
import json
import unittest

from scripts.gadget_burnside_orbit_certificate import (
    DEFAULT_OUTPUT,
    build_artifact,
    build_row,
    cycle_type,
    fixed_graph_edge_polynomial,
    induced_edge_cycle_lengths,
    validate_artifact,
)


class GadgetBurnsideOrbitCertificateTests(unittest.TestCase):
    def test_checked_artifact_reproduces(self) -> None:
        checked = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(checked, build_artifact())
        self.assertEqual(validate_artifact(checked)["orbit_counts"], [20, 90])

    def test_terminal_cycle_type(self) -> None:
        self.assertEqual(cycle_type((0, 1, 2)), (1, 1, 1))
        self.assertEqual(cycle_type((1, 2, 0)), (3,))
        self.assertEqual(cycle_type((1, 0, 3, 2)), (2, 2))

    def test_s3_class_actions(self) -> None:
        self.assertEqual(induced_edge_cycle_lengths((0, 1, 2)), (1, 1, 1, 1, 1, 1))
        self.assertEqual(induced_edge_cycle_lengths((1, 0, 2)), (2, 2, 1, 1))
        self.assertEqual(induced_edge_cycle_lengths((1, 2, 0)), (3, 3))

    def test_fixed_graph_polynomial(self) -> None:
        self.assertEqual(fixed_graph_edge_polynomial((3, 3)), [1, 0, 0, 2, 0, 0, 1])
        self.assertEqual(sum(fixed_graph_edge_polynomial((2, 2, 1, 1))), 16)

    def test_burnside_rows_match_canonical_census(self) -> None:
        three = build_row(3)
        four = build_row(4)
        self.assertEqual((three["fixed_graph_sum"], three["group_order"]), (120, 6))
        self.assertEqual((four["fixed_graph_sum"], four["group_order"]), (2160, 24))
        self.assertEqual(three["canonical_orbits_by_burnside"], 20)
        self.assertEqual(four["canonical_orbits_by_burnside"], 90)

    def test_invalid_scope_fails_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "three or four"):
            build_row(2)

    def test_tampering_fails_closed(self) -> None:
        tampered = copy.deepcopy(build_artifact())
        tampered["rows"][0]["fixed_graph_sum"] += 1
        with self.assertRaisesRegex(ValueError, "does not exactly reproduce"):
            validate_artifact(tampered)


if __name__ == "__main__":
    unittest.main()
