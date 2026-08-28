from __future__ import annotations

import sys
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from matched_torus_reference import (  # noqa: E402
    axis_geometry,
    cluster_stats,
    diamond_geometry,
)
from torus_homology import (  # noqa: E402
    HomologyUnionFind,
    classify_configuration,
    matching_channel_differences,
)


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
        self.assertEqual(counts, {"rank0": 259, "rank1": 162, "rank2": 91, "d0": 175, "d1": 175})

    def test_diamond_l2_exhaustive_channels(self) -> None:
        counts = self._exhaustive_check(diamond_geometry(2))
        self.assertEqual(counts, {"rank0": 143, "rank1": 68, "rank2": 45, "d0": 81, "d1": 81})


if __name__ == "__main__":
    unittest.main()
