#!/usr/bin/env python3
"""Lock the P2 manuscript evidence assembly to the committed census artifacts."""

from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p2_manuscript_evidence_table import (  # noqa: E402
    DEFAULT_OUTPUT,
    build_result,
    decimal_prefix,
    decimal_string,
    _reflect,
    render_markdown,
    validate_result,
)

MANUSCRIPT = ROOT / "docs" / "manuscripts" / "p2-algebraic-exclusion"

from fractions import Fraction  # noqa: E402


class P2ManuscriptEvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = build_result()
        cls.committed = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))

    def test_committed_artifact_reproduces_exactly(self) -> None:
        self.assertEqual(self.committed, self.result)
        validate_result(self.committed)

    def test_decimal_rendering_is_exact_and_truncating(self) -> None:
        self.assertEqual(decimal_string(Fraction(1, 3), 6), "0.333333")
        self.assertEqual(decimal_string(Fraction(2, 3), 6), "0.666666")
        self.assertEqual(decimal_prefix(Fraction(1, 3), Fraction(34, 100), 6), "0.3")
        self.assertEqual(decimal_prefix(Fraction(1, 3), Fraction(334, 1000), 6), "0.33")

    def test_method_intervals_are_pairwise_disjoint(self) -> None:
        disjointness = self.result["interval_disjointness"]
        self.assertTrue(disjointness["pairwise_disjoint"])
        self.assertEqual(
            disjointness["ascending_order"],
            ["mertens-2022-p-cell", "mertens-2022-p-med", "yang-zhou-2024-corrected", "jacobsen-2015-eigenvalue"],
        )
        for gap in disjointness["adjacent_gaps"]:
            self.assertGreater(Fraction(gap["gap_text"]), 0, gap)

    def test_search_class_sizes_match_the_look_elsewhere_ledger(self) -> None:
        counts = self.result["search_class"]["primitive_polynomial_counts_by_degree"]
        self.assertEqual(counts["1"], 12175)
        self.assertEqual(counts["2"], 3355121)
        self.assertEqual(counts["3"], 749507743)
        self.assertEqual(counts["4"], 157309446881)
        self.assertEqual(sum(counts.values()), self.result["search_class"]["polynomials_per_interval"])

    def test_every_degree_and_interval_pair_is_reported(self) -> None:
        pairs = {(row["degree"], row["interval_id"]) for row in self.result["exclusion_table"]}
        intervals = {row["interval_id"] for row in self.result["intervals"]}
        self.assertEqual(len(intervals), 4)
        self.assertEqual(pairs, {(degree, interval) for degree in (1, 2, 3, 4) for interval in intervals})

    def test_degrees_one_to_three_exclude_every_interval(self) -> None:
        for row in self.result["exclusion_table"]:
            if row["degree"] <= 3:
                with self.subTest(degree=row["degree"], interval=row["interval_id"]):
                    self.assertTrue(row["excluded"])
                    self.assertEqual(row["root_containing_polynomials"], 0)

    def test_two_narrowest_intervals_are_excluded_at_every_degree(self) -> None:
        self.assertEqual(
            sorted(self.result["intervals_excluded_at_every_degree"]),
            ["jacobsen-2015-eigenvalue", "yang-zhou-2024-corrected"],
        )

    def test_approach_resolution_is_a_strictly_decreasing_staircase(self) -> None:
        approach = self.result["closest_approach_by_degree"]
        self.assertTrue(approach["approach_resolution_strictly_decreasing_in_degree"])
        resolution = approach["approach_resolution_by_degree"]
        self.assertEqual(sorted(resolution), ["1", "2", "3", "4"])
        floors = [Fraction(resolution[str(degree)]["closest_approach_text"]) for degree in (1, 2, 3, 4)]
        self.assertEqual(floors[3], 0)
        for lower_degree, higher_degree in zip(floors, floors[1:]):
            self.assertGreater(lower_degree, higher_degree)

    def test_degree_four_is_the_boundary_degree_for_the_frozen_intervals(self) -> None:
        approach = self.result["closest_approach_by_degree"]
        self.assertEqual(approach["degrees_whose_closest_polynomial_stays_further_than_one_interval_width"], [1, 2, 3])
        self.assertEqual(approach["boundary_degree"], 4)
        widths = {row["interval_id"]: Fraction(row["width_text"]) for row in self.result["intervals"]}
        for row in approach["rows"]:
            ratio = Fraction(row["floor_to_interval_width_ratio_text"])
            expected = Fraction(row["root_distance_lower_bound_text"]) / widths[row["interval_id"]]
            with self.subTest(degree=row["degree"], interval=row["interval_id"]):
                self.assertEqual(ratio, expected)
                if row["degree"] <= 3:
                    self.assertGreater(ratio, 1)
                self.assertEqual(row["root_inside_interval"], ratio == 0)

    def test_approach_floor_is_consistent_with_the_exclusion_table(self) -> None:
        excluded = {(row["degree"], row["interval_id"]): row["excluded"] for row in self.result["exclusion_table"]}
        for row in self.result["closest_approach_by_degree"]["rows"]:
            with self.subTest(degree=row["degree"], interval=row["interval_id"]):
                self.assertEqual(excluded[(row["degree"], row["interval_id"])], not row["root_inside_interval"])

    def test_historical_forms_are_low_complexity_and_certified(self) -> None:
        historical = self.result["historical_form_complexity"]
        self.assertEqual(historical["max_degree"], 6)
        self.assertEqual(historical["max_height"], 3)
        self.assertFalse(historical["all_inside_census_class"])
        self.assertEqual(historical["outside_census_class"], ["(3,12^2) site"])
        native = json.loads((ROOT / "results" / "pslq-lattice-native-candidates" / "latest.json").read_text(encoding="utf-8"))
        certified = {tuple(row["minimal_polynomial_coefficients_ascending"]) for row in native["candidates"]}
        for row in historical["rows"]:
            coefficients = row["minimal_polynomial_ascending"]
            with self.subTest(lattice=row["lattice"]):
                self.assertEqual(row["degree"], len(coefficients) - 1)
                self.assertEqual(row["height"], max(abs(value) for value in coefficients))
                self.assertGreater(coefficients[-1], 0)
                self.assertEqual(row["inside_census_class"], row["degree"] <= 4 and row["height"] <= 100)
        # every row is either certified by the native artifact, a matching reflection of one, or 1/2
        for row in historical["rows"]:
            coefficients = tuple(row["minimal_polynomial_ascending"])
            with self.subTest(lattice=row["lattice"]):
                self.assertTrue(
                    coefficients in certified
                    or tuple(_reflect(list(coefficients))) in certified
                    or coefficients == (-1, 2),
                    coefficients,
                )

    def test_matching_reflection_is_an_involution_and_matches_sykes_essam(self) -> None:
        kagome = [1, 0, -3, 1]
        triangular = _reflect(kagome)
        self.assertEqual(triangular, [1, -3, 0, 1])
        self.assertEqual(_reflect(triangular), kagome)

    def test_recommended_extension_is_costed_and_not_run(self) -> None:
        extension = self.result["historical_form_complexity"]["recommended_extension"]
        self.assertEqual(extension["status"], "not run here")
        self.assertEqual(extension["total_polynomials_per_interval"], sum(extension["primitive_counts_by_degree"].values()))
        self.assertEqual(sorted(extension["primitive_counts_by_degree"], key=int), ["1", "2", "3", "4", "5", "6"])
        self.assertLess(
            extension["total_polynomials_per_interval"],
            self.result["search_class"]["polynomials_per_interval"],
        )

    def test_sensitivity_control_is_digested_as_a_source(self) -> None:
        paths = {artifact["path"] for artifact in self.result["source_artifacts"]}
        self.assertIn("results/pslq-degree4-synthetic-boundary-control/latest.json", paths)
        control = json.loads(
            (ROOT / "results" / "pslq-degree4-synthetic-boundary-control" / "latest.json").read_text(encoding="utf-8")
        )
        self.assertTrue(control["conclusion"]["all_trials_passed"])

    def test_every_quartic_survivor_survives_exactly_one_interval(self) -> None:
        census = self.result["quartic_survivor_census"]
        self.assertEqual(census["distinct_surviving_quartics"], 16)
        self.assertEqual(census["max_intervals_per_survivor"], 1)
        self.assertEqual(census["survivors_per_interval"]["mertens-2022-p-cell"], 15)
        self.assertEqual(census["survivors_per_interval"]["mertens-2022-p-med"], 1)
        for survivor in census["survivors"]:
            with self.subTest(coefficients=survivor["coefficients_ascending"]):
                self.assertEqual(len(survivor["excluded_by_interval_ids"]), 3)
                self.assertEqual(len(survivor["separations"]), 3)
                self.assertLessEqual(survivor["height"], 100)
                self.assertGreater(survivor["coefficients_ascending"][-1], 0)

    def test_survivor_roots_are_certified_outside_the_other_intervals(self) -> None:
        bounds = {row["interval_id"]: (Fraction(row["lower"]), Fraction(row["upper"])) for row in self.result["intervals"]}
        for survivor in self.result["quartic_survivor_census"]["survivors"]:
            low, high = (Fraction(value) for value in survivor["root_bracket"])
            self.assertLess(low, high)
            home = bounds[survivor["surviving_interval_id"]]
            self.assertGreaterEqual(low, home[0])
            self.assertLessEqual(high, home[1])
            for separation in survivor["separations"]:
                gap = Fraction(separation["separation_lower_bound_text"])
                lower, upper = bounds[separation["interval_id"]]
                with self.subTest(coefficients=survivor["coefficients_ascending"], interval=separation["interval_id"]):
                    self.assertGreater(gap, 0)
                    self.assertTrue(high < lower or low > upper)

    def test_survivor_witnesses_agree_with_the_source_census_artifacts(self) -> None:
        census = self.result["quartic_survivor_census"]
        for interval_id, expected in census["survivors_per_interval"].items():
            source = json.loads((ROOT / "results" / f"pslq-degree4-{interval_id}" / "latest.json").read_text(encoding="utf-8"))
            witnesses = source["interval_result"]["root_witnesses"]
            self.assertEqual(len(witnesses), expected, interval_id)
            committed = {tuple(row["coefficients_ascending"]) for row in witnesses}
            derived = {
                tuple(row["coefficients_ascending"])
                for row in census["survivors"]
                if row["surviving_interval_id"] == interval_id
            }
            self.assertEqual(committed, derived, interval_id)

    def test_source_artifact_digests_are_recorded(self) -> None:
        import hashlib

        self.assertGreaterEqual(len(self.result["source_artifacts"]), 16)
        for artifact in self.result["source_artifacts"]:
            path = ROOT / artifact["path"]
            with self.subTest(path=artifact["path"]):
                self.assertTrue(path.is_file())
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), artifact["sha256"])

    def test_claim_boundary_refuses_transcendence_and_p_values(self) -> None:
        excluded = self.result["claim_boundary"]["excluded"]
        for forbidden in ("transcendence", "p-values", "closed forms", "degree/height expansion"):
            self.assertIn(forbidden, excluded)

    def test_rendered_tables_do_not_drift_from_the_artifact(self) -> None:
        committed = (MANUSCRIPT / "tables.md").read_text(encoding="utf-8")
        self.assertEqual(committed, render_markdown(self.result))
        self.assertIn("Do not edit by hand", committed.splitlines()[0])

    def test_manuscript_states_the_bounded_scope(self) -> None:
        draft = (MANUSCRIPT / "manuscript.md").read_text(encoding="utf-8")
        for required in ("!= transcendental", "!= non-algebraic", "bounded certified exclusion"):
            with self.subTest(required=required):
                self.assertIn(required, draft)
        for forbidden in ("proves that p_c is transcendental", "closed form for p_c is"):
            with self.subTest(forbidden=forbidden):
                self.assertNotIn(forbidden, draft)

    def test_manuscript_cites_only_generated_tables(self) -> None:
        draft = (MANUSCRIPT / "manuscript.md").read_text(encoding="utf-8")
        for table in ("tables.md#table-1", "tables.md#table-2", "tables.md#table-3", "tables.md#table-4", "tables.md#table-5"):
            with self.subTest(table=table):
                self.assertIn(table, draft)

    def test_tampering_fails(self) -> None:
        changed = copy.deepcopy(self.result)
        changed["quartic_survivor_census"]["distinct_surviving_quartics"] = 15
        with self.assertRaisesRegex(ValueError, "does not exactly reproduce"):
            validate_result(changed)


if __name__ == "__main__":
    unittest.main()
