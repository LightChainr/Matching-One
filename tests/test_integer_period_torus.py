from __future__ import annotations

import sys
import unittest
from pathlib import Path

import mpmath as mp


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from gaussian_circulant_geometry import GaussianTorus  # noqa: E402
from integer_period_torus import (  # noqa: E402
    IntegerPeriods,
    axis_integer_torus,
    classify_configuration,
    diamond_integer_torus,
    gaussian_integer_torus,
    integer_torus_geometry,
    matrix_vector,
    unimodular_inverse,
)
from matched_torus_reference import (  # noqa: E402
    axis_geometry,
    cluster_stats,
    diamond_geometry,
    exact_check,
)
from torus_homology import classify_configuration as old_classify  # noqa: E402


def remap_active(source, target, active):
    mapped = [False] * target.n
    for vertex, point in enumerate(source.coordinates):
        mapped[target.vertex(point)] = active[vertex]
    return mapped


class IntegerPeriodArithmeticTests(unittest.TestCase):
    def test_exact_adjugate_winding_and_rejection(self) -> None:
        periods = IntegerPeriods(((3, 1), (1, 2)))
        for winding in ((1, 0), (0, 1), (2, -3), (-4, 7)):
            displacement = periods.period_vector(winding)
            self.assertEqual(periods.winding(displacement), winding)
        with self.assertRaisesRegex(ValueError, "period lattice"):
            periods.winding((1, 0))

    def test_quotient_has_abs_determinant_vertices(self) -> None:
        for matrix in (
            ((3, 1), (1, 2)),
            ((2, -1), (1, 3)),
            ((-2, 1), (1, 2)),
            ((1, 0), (0, 1)),
        ):
            geometry = integer_torus_geometry(matrix)
            self.assertEqual(geometry.n, abs(geometry.periods.det))
            self.assertEqual(
                len({geometry.periods.quotient_key(p) for p in geometry.coordinates}),
                geometry.n,
            )

    def test_unimodular_basis_change_vector_law(self) -> None:
        periods = IntegerPeriods(((4, 1), (1, 3)))
        change = ((1, 1), (0, 1))
        changed = periods.change_basis(change)
        inverse = unimodular_inverse(change)
        for winding in ((1, 0), (0, 1), (2, -3)):
            displacement = periods.period_vector(winding)
            self.assertEqual(
                changed.winding(displacement), matrix_vector(inverse, winding)
            )
            self.assertTrue(changed.equivalent(displacement, (0, 0)))


class ReferenceRegressionTests(unittest.TestCase):
    def _compare_reference(self, old, new) -> None:
        self.assertEqual(old.n, new.n)
        for mask in range(1 << old.n):
            old_active = [bool((mask >> vertex) & 1) for vertex in range(old.n)]
            new_active = [False] * new.n
            for old_vertex, coordinate in enumerate(old.coordinates):
                if old.name == "diamond":
                    u, v = coordinate
                    # The old reference uses u=x+y, v=y-x.
                    x = (u - v) // 2
                    y = (u + v) // 2
                    point = (x, y)
                else:
                    point = coordinate
                new_active[new.vertex(point)] = old_active[old_vertex]

            for matching in (False, True):
                old_edges = old.matching_edges if matching else old.primal_edges
                _, old_wrap = cluster_stats(old_active, old_edges)
                old_channels, _ = old_classify(old, old_active, matching=matching)
                new_channels, _ = classify_configuration(
                    new, new_active, matching=matching
                )
                self.assertEqual(new_channels.either, old_wrap)
                self.assertEqual(new_channels.max_rank, old_channels.max_rank)
                self.assertEqual(new_channels.cross, old_channels.cross)

    def test_axis_l3_exhaustive_regression(self) -> None:
        self._compare_reference(axis_geometry(3), axis_integer_torus(3))

    def test_diamond_l2_exhaustive_regression(self) -> None:
        self._compare_reference(diamond_geometry(2), diamond_integer_torus(2))

    def test_gaussian_circulant_labels_and_edges(self) -> None:
        reference = GaussianTorus(2, 1)
        geometry = gaussian_integer_torus(2, 1)
        labels = {reference.label(x, y) for x, y in geometry.coordinates}
        self.assertEqual(labels, set(range(reference.n)))

        for matching in (False, True):
            edges = geometry.matching_edges if matching else geometry.primal_edges
            forward = {
                (reference.label(*geometry.coordinates[edge.j])
                 - reference.label(*geometry.coordinates[edge.i])) % reference.n
                for edge in edges
            }
            residues = forward | {(-residue) % reference.n for residue in forward}
            self.assertEqual(residues, reference.edge_residues(matching))

    def test_tiny_exact_matching_identities(self) -> None:
        mp.mp.dps = 80
        p = mp.mpf("0.317")
        for geometry in (
            axis_integer_torus(2),
            # L=1 collapses matching edges into degenerate self-loops in both
            # reference implementations, so use the first nondegenerate case.
            diamond_integer_torus(2),
            gaussian_integer_torus(2, 1),
        ):
            result = exact_check(geometry, p)
            self.assertLess(abs(result["difference"]), mp.mpf("1e-65"))


class UnimodularInvarianceTests(unittest.TestCase):
    def test_exhaustive_rank_invariance(self) -> None:
        matrix = ((3, 1), (1, 2))  # determinant 5
        change = ((1, 1), (0, 1))
        first = integer_torus_geometry(matrix)
        second = integer_torus_geometry(
            first.periods.change_basis(change).matrix
        )
        self.assertEqual(first.n, second.n)

        for mask in range(1 << first.n):
            active_first = [bool((mask >> vertex) & 1) for vertex in range(first.n)]
            active_second = remap_active(first, second, active_first)
            for matching in (False, True):
                channels_first, components_first = classify_configuration(
                    first, active_first, matching=matching
                )
                channels_second, components_second = classify_configuration(
                    second, active_second, matching=matching
                )
                self.assertEqual(channels_first.max_rank, channels_second.max_rank)
                self.assertEqual(channels_first.either, channels_second.either)
                self.assertEqual(channels_first.cross, channels_second.cross)
                self.assertEqual(
                    sorted((c.size, c.rank) for c in components_first),
                    sorted((c.size, c.rank) for c in components_second),
                )


if __name__ == "__main__":
    unittest.main()
