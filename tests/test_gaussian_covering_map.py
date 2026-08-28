from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from gaussian_covering_map import (  # noqa: E402
    CoveringMap,
    GaussianPair,
    canonical_cover,
    covering_units,
    det,
)


class GaussianCoveringMapTests(unittest.TestCase):
    def test_norm2_and_norm5_lineages_are_exact_covers(self) -> None:
        examples = (
            ((8, 1), (9, 7), 2),
            ((7, 4), (11, 3), 2),
            ((9, 2), (11, 7), 2),
            ((7, 6), (13, 1), 2),
            ((12, 1), (13, 11), 2),
            ((9, 8), (17, 1), 2),
            ((8, 1), (17, 6), 5),
            ((7, 4), (18, 1), 5),
            ((9, 2), (16, 13), 5),
            ((7, 6), (19, 8), 5),
        )
        for parent_pair, child_pair, degree in examples:
            with self.subTest(parent=parent_pair, child=child_pair):
                parent = GaussianPair(*parent_pair)
                child = GaussianPair(*child_pair)
                candidates = covering_units(parent, child)
                self.assertEqual(len(candidates), 4)
                cover = canonical_cover(parent, child)
                self.assertEqual(cover.degree, degree)
                self.assertEqual(abs(det(cover.homology_matrix())), degree)
                cover.verify_partition()
                cover.verify_edges()
                cover.verify_homology()

    def test_each_fiber_has_one_label_per_parent_residue_class(self) -> None:
        cover = canonical_cover(GaussianPair(8, 1), GaussianPair(17, 6))
        n = cover.parent.n
        for parent_label in range(n):
            fiber = cover.fiber(parent_label)
            self.assertEqual(len(fiber), 5)
            self.assertEqual(
                {cover.map_vertex(child_label) for child_label in fiber},
                {parent_label},
            )
            self.assertEqual(
                {child_label % n for child_label in fiber},
                {cover.inverse_t * parent_label % n},
            )

    def test_inverse_coordinates_are_a_bijection(self) -> None:
        cover = canonical_cover(GaussianPair(9, 2), GaussianPair(16, 13))
        coordinates = {
            cover.inverse_coordinates(child_label)
            for child_label in range(cover.child.n)
        }
        expected = {
            (parent_label, kernel_index)
            for parent_label in range(cover.parent.n)
            for kernel_index in range(cover.degree)
        }
        self.assertEqual(coordinates, expected)

    def test_all_d4_covering_units_preserve_edges_and_homology_degree(self) -> None:
        parent = GaussianPair(8, 1)
        child = GaussianPair(17, 6)
        for t in covering_units(parent, child):
            cover = CoveringMap(parent, child, t)
            self.assertEqual(len(set(cover.direction_map().values())), 2)
            self.assertEqual(abs(det(cover.direction_matrix())), 1)
            self.assertEqual(abs(det(cover.homology_matrix())), 5)
            cover.verify_partition()
            cover.verify_edges()
            cover.verify_homology()

    def test_global_winding_map_is_not_just_the_local_d4_direction_map(self) -> None:
        cover = canonical_cover(GaussianPair(8, 1), GaussianPair(9, 7))
        # For the canonical minimal t=28, local lattice directions are a D4
        # reflection, but child fundamental torus cycles map to combinations
        # of parent cycles. This guards against treating directional wrapping
        # classes as if they transformed only by the local D4 map.
        self.assertEqual(cover.t, 28)
        self.assertEqual(cover.direction_matrix(), ((-1, 0), (0, 1)))
        self.assertEqual(cover.homology_matrix(), ((-1, 1), (1, 1)))
        self.assertEqual(cover.map_winding((1, 0)), (-1, 1))
        self.assertEqual(cover.map_winding((0, 1)), (1, 1))
        self.assertEqual(abs(det(cover.homology_matrix())), 2)


if __name__ == "__main__":
    unittest.main()
