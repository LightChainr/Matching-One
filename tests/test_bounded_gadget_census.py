from __future__ import annotations

import copy
import json
from pathlib import Path
import unittest

from scripts.bounded_gadget_census import (
    DEFAULT_OUTPUT,
    build_artifact,
    build_census_row,
    is_connected,
    validate_artifact,
    vertex_degrees,
)


class BoundedGadgetCensusTests(unittest.TestCase):
    def test_checked_artifact_reproduces(self) -> None:
        checked = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(checked, build_artifact())
        self.assertEqual(
            validate_artifact(checked)["status"],
            "valid_exact_bounded_candidate_space_census",
        )

    def test_three_terminal_census(self) -> None:
        row = build_census_row(3)
        self.assertEqual(row["labeled_simple_graphs"], 64)
        self.assertEqual(row["canonical_orbits"], 20)
        self.assertEqual(row["connected_carrier_orbits"], 11)
        self.assertEqual(row["connected_internal_degree_at_least_3_orbits"], 4)
        self.assertEqual(row["orbit_multiplicity_sum"], 64)

    def test_four_terminal_census(self) -> None:
        row = build_census_row(4)
        self.assertEqual(row["labeled_simple_graphs"], 1024)
        self.assertEqual(row["canonical_orbits"], 90)
        self.assertEqual(row["connected_carrier_orbits"], 58)
        self.assertEqual(row["connected_internal_degree_at_least_3_orbits"], 27)
        self.assertEqual(row["orbit_multiplicity_sum"], 1024)

    def test_edge_histograms_sum_to_orbit_counts(self) -> None:
        for terminal_count in (3, 4):
            row = build_census_row(terminal_count)
            for count_key, histogram_key in (
                ("canonical_orbits", "all_orbit_edge_count_histogram"),
                ("connected_carrier_orbits", "connected_carrier_edge_count_histogram"),
                (
                    "connected_internal_degree_at_least_3_orbits",
                    "connected_internal_degree_at_least_3_edge_count_histogram",
                ),
            ):
                self.assertEqual(sum(row[histogram_key].values()), row[count_key])

    def test_connectivity_and_degrees_are_structural(self) -> None:
        self.assertFalse(is_connected(4, ((0, 1), (1, 2))))
        self.assertTrue(is_connected(4, ((0, 3), (1, 3), (2, 3))))
        self.assertEqual(
            vertex_degrees(4, ((0, 3), (1, 3), (2, 3))),
            (1, 1, 1, 3),
        )

    def test_invalid_scope_fails_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "three or four"):
            build_census_row(2)

    def test_artifact_tampering_fails_closed(self) -> None:
        tampered = copy.deepcopy(build_artifact())
        tampered["rows"][1]["connected_carrier_orbits"] += 1
        with self.assertRaisesRegex(ValueError, "does not exactly reproduce"):
            validate_artifact(tampered)


if __name__ == "__main__":
    unittest.main()
