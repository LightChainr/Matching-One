from __future__ import annotations

import random
import sys
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from gaussian_circulant_geometry import GaussianTorus  # noqa: E402
from matched_torus_reference import (  # noqa: E402
    axis_geometry,
    cluster_stats,
    diamond_geometry,
    diamond_xy_geometry,
    gaussian_geometry,
    integer_period_geometry,
)
from torus_homology import (  # noqa: E402
    HomologyUnionFind,
    apply_matrix,
    classify_configuration,
    component_homologies,
    determinant,
    exhaustive_channel_counts,
    gaussian_period_matrix,
    matching_channel_differences,
    matmul2,
    qspan_contains,
    transport_basis,
    unimodular_inverse,
    winding_coefficients,
    wrapping_channels,
)


# Frozen PR #21 regressions. These counts must remain exact.
PR21_AXIS_L3 = {"rank0": 259, "rank1": 162, "rank2": 91, "d0": 175, "d1": 175}
PR21_DIAMOND_L2 = {"rank0": 143, "rank1": 68, "rank2": 45, "d0": 81, "d1": 81}
GAUSSIAN_2_1 = {"rank0": 16, "rank1": 10, "rank2": 6, "d0": 11, "d1": 11}


def _random_unimodular(rng: random.Random) -> tuple[tuple[int, int], tuple[int, int]]:
    """Return a random det +/-1 integer matrix generated from elementary shears."""

    matrix = ((1, 0), (0, 1))
    for _ in range(rng.randint(3, 8)):
        shift = rng.randint(-3, 3)
        if rng.randrange(2) == 0:
            shear = ((1, shift), (0, 1))
        else:
            shear = ((1, 0), (shift, 1))
        matrix = matmul2(matrix, shear)
    if rng.randrange(2) == 0:
        matrix = matmul2(matrix, ((-1, 0), (0, 1)))
    return matrix


class HomologyUnionFindTests(unittest.TestCase):
    def test_rank_and_primitive_spiral_direction(self) -> None:
        union_find = HomologyUnionFind(1, (4, 6))
        union_find.add_edge(0, 0, 8, 6)
        component = union_find.component(0)
        self.assertEqual(component.rank, 1)
        self.assertEqual(component.basis, ((2, 1),))
        self.assertTrue(component.direction_0)
        self.assertTrue(component.direction_1)
        self.assertTrue(component.both)
        self.assertFalse(component.cross)

        # A parallel cycle does not increase rank; an independent one does.
        union_find.add_edge(0, 0, -8, -6)
        self.assertEqual(union_find.component(0).rank, 1)
        union_find.add_edge(0, 0, 0, 6)
        self.assertEqual(union_find.component(0).rank, 2)
        self.assertTrue(union_find.component(0).cross)

    def test_rejects_cycle_outside_period_lattice(self) -> None:
        union_find = HomologyUnionFind(1, (4, 6))
        with self.assertRaisesRegex(ValueError, "period lattice"):
            union_find.add_edge(0, 0, 1, 0)

    def test_diagonal_tuple_matches_explicit_matrix(self) -> None:
        first = HomologyUnionFind(1, (4, 6))
        second = HomologyUnionFind(1, ((4, 0), (0, 6)))
        first.add_edge(0, 0, 8, 6)
        second.add_edge(0, 0, 8, 6)
        self.assertEqual(first.component(0).basis, second.component(0).basis)

    def test_adjugate_winding_is_exact_for_nondiagonal_matrix(self) -> None:
        period = gaussian_period_matrix(2, 1)
        self.assertEqual(winding_coefficients(2, 1, period), (1, 0))
        self.assertEqual(winding_coefficients(-1, 2, period), (0, 1))
        self.assertEqual(winding_coefficients(4, 2, period), (2, 0))
        self.assertEqual(winding_coefficients(1, 3, period), (1, 1))
        with self.assertRaisesRegex(ValueError, "period lattice"):
            winding_coefficients(1, 0, period)


class ExactTinyGeometryTests(unittest.TestCase):
    def _exhaustive_check(self, geometry: object) -> dict[str, int]:
        counts = {"rank0": 0, "rank1": 0, "rank2": 0, "d0": 0, "d1": 0}
        for mask in range(1 << geometry.n):
            black = [bool((mask >> vertex) & 1) for vertex in range(geometry.n)]
            white = [not value for value in black]
            black_channels, black_components = classify_configuration(geometry, black)
            white_channels, white_components = classify_configuration(
                geometry, white, matching=True
            )

            # Regression against the previous Boolean implementation.
            _, old_black_wrap = cluster_stats(black, geometry.primal_edges)
            _, old_white_wrap = cluster_stats(white, geometry.matching_edges)
            self.assertEqual(black_channels.either, old_black_wrap)
            self.assertEqual(white_channels.either, old_white_wrap)

            for channels, components in (
                (black_channels, black_components),
                (white_channels, white_components),
            ):
                self.assertEqual(channels.cross, any(c.rank == 2 for c in components))
                self.assertEqual(channels.either, channels.max_rank > 0)
                self.assertTrue(not channels.cross or channels.both)
                self.assertTrue(not channels.both or channels.either)

            # Configuration-level matching topology makes all five difference
            # channels identical on these symmetric periodic quotients.
            differences = matching_channel_differences(geometry, black)
            self.assertEqual(len(set(differences.values())), 1)

            counts[f"rank{black_channels.max_rank}"] += 1
            counts["d0"] += int(black_channels.direction_0)
            counts["d1"] += int(black_channels.direction_1)
        return counts

    def test_axis_l3_exhaustive_channels(self) -> None:
        counts = self._exhaustive_check(axis_geometry(3))
        self.assertEqual(counts, PR21_AXIS_L3)

    def test_diamond_l2_exhaustive_channels(self) -> None:
        counts = self._exhaustive_check(diamond_geometry(2))
        self.assertEqual(counts, PR21_DIAMOND_L2)

    def test_axis_via_general_integer_period_matrix(self) -> None:
        geometry = integer_period_geometry(((3, 0), (0, 3)), name="axis", L=3)
        self.assertEqual(exhaustive_channel_counts(geometry), PR21_AXIS_L3)
        self.assertEqual(self._exhaustive_check(geometry), PR21_AXIS_L3)

    def test_diamond_xy_matches_uv_embedding(self) -> None:
        counts = self._exhaustive_check(diamond_xy_geometry(2))
        self.assertEqual(counts, PR21_DIAMOND_L2)

    def test_gaussian_2_1_exhaustive_channels(self) -> None:
        geometry = gaussian_geometry(2, 1)
        self.assertEqual(geometry.n, 5)
        self.assertEqual(geometry.period_matrix, ((2, -1), (1, 2)))
        counts = self._exhaustive_check(geometry)
        self.assertEqual(counts, GAUSSIAN_2_1)
        self.assertEqual(exhaustive_channel_counts(geometry), GAUSSIAN_2_1)


class GaussianCirculantAgreementTests(unittest.TestCase):
    def test_cyclic_labels_and_winding_map(self) -> None:
        geometry = gaussian_geometry(2, 1)
        cyclic = GaussianTorus(2, 1)
        labels = [cyclic.label(x, y) for x, y in geometry.coordinates]
        self.assertEqual(sorted(labels), list(range(5)))
        for edge in geometry.primal_edges + geometry.matching_edges:
            start = geometry.coordinates[edge.i]
            dest = geometry.coordinates[edge.j]
            expected = (
                cyclic.label(*start) + cyclic.a * edge.dx + cyclic.b * edge.dy
            ) % cyclic.n
            self.assertEqual(cyclic.label(*dest), expected)
        self.assertEqual(cyclic.winding_coordinates(2, 1), (1, 0))
        self.assertEqual(
            winding_coefficients(2, 1, geometry.period_matrix),
            cyclic.winding_coordinates(2, 1),
        )
        self.assertEqual(
            winding_coefficients(-1, 2, geometry.period_matrix),
            cyclic.winding_coordinates(-1, 2),
        )


class UnimodularBasisChangeTests(unittest.TestCase):
    def _assert_invariant(
        self,
        geometry: object,
        period_matrix,
        change,
        active: list[bool],
        edges,
    ) -> None:
        changed = matmul2(period_matrix, change)
        inverse = unimodular_inverse(change)
        original = component_homologies(active, edges, period_matrix)
        transformed = component_homologies(active, edges, changed)
        original_channels = wrapping_channels(original)
        transformed_channels = wrapping_channels(transformed)
        self.assertEqual(original_channels.max_rank, transformed_channels.max_rank)
        self.assertEqual(original_channels.either, transformed_channels.either)
        self.assertEqual(original_channels.cross, transformed_channels.cross)
        self.assertEqual(len(original), len(transformed))
        for left, right in zip(original, transformed):
            self.assertEqual(left.root, right.root)
            self.assertEqual(left.size, right.size)
            self.assertEqual(left.rank, right.rank)
            transported = transport_basis(left.basis, inverse)
            self.assertEqual(len(transported), len(right.basis))
            for vector in transported:
                self.assertTrue(qspan_contains(vector, right.basis))
            for vector in right.basis:
                self.assertTrue(qspan_contains(apply_matrix(change, vector), left.basis))

    def _exhaustive_invariance(self, geometry: object, seed: int, n_changes: int) -> None:
        rng = random.Random(seed)
        period_matrix = geometry.period_matrix
        self.assertNotEqual(determinant(period_matrix), 0)
        changes = [_random_unimodular(rng) for _ in range(n_changes)]
        self.assertTrue(any(abs(determinant(change)) == 1 for change in changes))
        self.assertTrue(any(determinant(change) == 1 for change in changes))
        self.assertTrue(any(determinant(change) == -1 for change in changes))
        for change in changes:
            self.assertEqual(abs(determinant(change)), 1)
            for mask in range(1 << geometry.n):
                active = [bool((mask >> vertex) & 1) for vertex in range(geometry.n)]
                self._assert_invariant(
                    geometry, period_matrix, change, active, geometry.primal_edges
                )
                self._assert_invariant(
                    geometry,
                    period_matrix,
                    change,
                    [not value for value in active],
                    geometry.matching_edges,
                )

    def test_axis_l3_random_basis_change(self) -> None:
        self._exhaustive_invariance(axis_geometry(3), seed=20260828, n_changes=8)

    def test_diamond_l2_random_basis_change(self) -> None:
        self._exhaustive_invariance(diamond_geometry(2), seed=20260828, n_changes=8)

    def test_gaussian_2_1_random_basis_change(self) -> None:
        self._exhaustive_invariance(gaussian_geometry(2, 1), seed=20260828, n_changes=12)

    def test_shear_changes_generator_relative_channels_but_not_rank(self) -> None:
        union_find = HomologyUnionFind(1, ((3, 0), (0, 3)))
        union_find.add_edge(0, 0, 3, 0)
        self.assertEqual(union_find.component(0).basis, ((1, 0),))
        self.assertTrue(union_find.component(0).direction_0)
        self.assertFalse(union_find.component(0).direction_1)

        shear = ((1, 0), (1, 1))
        sheared = HomologyUnionFind(1, matmul2(((3, 0), (0, 3)), shear))
        sheared.add_edge(0, 0, 3, 0)
        self.assertEqual(sheared.component(0).rank, 1)
        self.assertTrue(sheared.component(0).direction_0)
        self.assertTrue(sheared.component(0).direction_1)


if __name__ == "__main__":
    unittest.main()
