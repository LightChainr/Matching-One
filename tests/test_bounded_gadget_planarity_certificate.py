
from __future__ import annotations
from itertools import combinations
import json
import unittest

from scripts.bounded_gadget_planarity_certificate import (
    DEFAULT_OUTPUT,
    build_artifact,
    build_row,
    cyclic_orders,
    minimum_orientable_genus,
    validate_artifact,
    validate_simple_graph,
)


class BoundedGadgetPlanarityCertificateTests(unittest.TestCase):
    def test_checked_artifact_reproduces(self) -> None:
        checked = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(checked, build_artifact())
        self.assertEqual(
            validate_artifact(checked)["status"],
            "valid_exact_bounded_planarity_certificate",
        )

    def test_cyclic_orders_remove_rotation_duplicates(self) -> None:
        self.assertEqual(len(cyclic_orders((0, 1, 2))), 2)
        self.assertEqual(len(cyclic_orders((0, 1, 2, 3))), 6)

    def test_k5_is_genus_one_after_exhaustion(self) -> None:
        k5 = tuple(combinations(range(5), 2))
        genus, checked, _ = minimum_orientable_genus(5, k5)
        self.assertEqual(genus, 1)
        self.assertEqual(checked, 7776)

    def test_k5_minus_edge_is_planar(self) -> None:
        k5_minus_edge = tuple(
            edge for edge in combinations(range(5), 2) if edge != (0, 1)
        )
        genus, checked, _ = minimum_orientable_genus(5, k5_minus_edge)
        self.assertEqual(genus, 0)
        self.assertGreater(checked, 0)

    def test_k33_control_is_genus_one(self) -> None:
        k33 = tuple((left, right) for left in range(3) for right in range(3, 6))
        genus, checked, _ = minimum_orientable_genus(6, k33)
        self.assertEqual((genus, checked), (1, 64))

    def test_exact_bounded_census(self) -> None:
        three = build_row(3)
        four = build_row(4)
        self.assertEqual((three["all_canonical_orbits"]["planar_orbits"], three["all_canonical_orbits"]["nonplanar_orbits"]), (20, 0))
        self.assertEqual((four["all_canonical_orbits"]["planar_orbits"], four["all_canonical_orbits"]["nonplanar_orbits"]), (89, 1))
        self.assertEqual(four["all_canonical_orbits"]["nonplanar_witnesses"][0]["encoding"], "4:5:1111111111")

    def test_filtered_counts(self) -> None:
        four = build_row(4)
        self.assertEqual((four["connected_carrier"]["planar_orbits"], four["connected_carrier"]["nonplanar_orbits"]), (57, 1))
        self.assertEqual((four["connected_internal_degree_at_least_3"]["planar_orbits"], four["connected_internal_degree_at_least_3"]["nonplanar_orbits"]), (26, 1))

    def test_invalid_graph_and_scope_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "self-loop"):
            validate_simple_graph(3, ((0, 0),))
        with self.assertRaisesRegex(ValueError, "duplicate"):
            validate_simple_graph(3, ((0, 1), (1, 0)))
        with self.assertRaisesRegex(ValueError, "three or four"):
            build_row(2)


if __name__ == "__main__":
    unittest.main()
