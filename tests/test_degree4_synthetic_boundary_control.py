#!/usr/bin/env python3
"""Lock the quartic-census sensitivity control where the census returned a null."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from degree4_synthetic_boundary_control import (  # noqa: E402
    DEFAULT_OUTPUT,
    PLANTED,
    build_result,
    validate_result,
)


class Degree4BoundaryControlTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = build_result()
        cls.committed = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))

    def test_committed_artifact_reproduces_exactly(self) -> None:
        self.assertEqual(self.committed, self.result)
        validate_result(self.committed)

    def test_control_uses_the_unmodified_census_path(self) -> None:
        path = self.result["census_path"]
        self.assertFalse(path["modified_for_this_control"])
        self.assertEqual(path["driver"], "scripts/degree4_interval_exclusion.py::run_search")

    def test_planted_quartics_are_committed_census_witnesses(self) -> None:
        import hashlib

        self.assertEqual(len(self.result["planted_targets"]), len(PLANTED))
        for target in self.result["planted_targets"]:
            source = ROOT / target["committed_by"]
            payload = json.loads(source.read_text(encoding="utf-8"))
            committed = {tuple(row["coefficients_ascending"]) for row in payload["interval_result"]["root_witnesses"]}
            with self.subTest(coefficients=target["coefficients_ascending"]):
                self.assertIn(tuple(target["coefficients_ascending"]), committed)
                self.assertLessEqual(target["height"], 100)
                self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(), target["committed_by_sha256"])

    def test_control_covers_exactly_the_intervals_the_census_nulled(self) -> None:
        contract = json.loads((ROOT / "analysis" / "pslq_search_contract.json").read_text(encoding="utf-8"))
        widths = {row["id"]: Fraction(row["upper"]) - Fraction(row["lower"]) for row in contract["intervals"]}
        nulled, found = set(), set()
        for interval_id in widths:
            census = json.loads(
                (ROOT / "results" / f"pslq-degree4-{interval_id}" / "latest.json").read_text(encoding="utf-8")
            )
            target = nulled if census["interval_result"]["root_containing_polynomials"] == 0 else found
            target.add(interval_id)
        selection = self.result["width_selection"]
        self.assertEqual(set(selection["covered"]), nulled)
        self.assertEqual(set(selection["not_covered"]), found)
        self.assertEqual(set(self.result["widths_tested"]), nulled)
        for width_id, width_text in self.result["widths_tested"].items():
            self.assertEqual(Fraction(width_text), widths[width_id])
        covered = {(row["width_id"], row["polarity"], tuple(row["planted_coefficients_ascending"])) for row in self.result["trials"]}
        expected = {
            (width_id, polarity, tuple(planted["coefficients_ascending"]))
            for width_id in nulled
            for polarity in ("positive", "negative")
            for planted in PLANTED
        }
        self.assertEqual(covered, expected)

    def test_covered_widths_are_the_narrowest(self) -> None:
        contract = json.loads((ROOT / "analysis" / "pslq_search_contract.json").read_text(encoding="utf-8"))
        widths = {row["id"]: Fraction(row["upper"]) - Fraction(row["lower"]) for row in contract["intervals"]}
        covered = set(self.result["width_selection"]["covered"])
        self.assertEqual(covered, {"jacobsen-2015-eigenvalue", "yang-zhou-2024-corrected"})
        self.assertLess(max(widths[i] for i in covered), min(widths[i] for i in widths if i not in covered))

    def test_trial_intervals_have_the_declared_width_and_placement(self) -> None:
        for row in self.result["trials"]:
            lower, upper = Fraction(row["lower"]), Fraction(row["upper"])
            low, high = (Fraction(value) for value in row["planted_root_bracket"])
            with self.subTest(width=row["width_id"], polarity=row["polarity"], planted=row["planted_coefficients_ascending"]):
                self.assertEqual(upper - lower, Fraction(row["width_text"]))
                if row["polarity"] == "positive":
                    self.assertLess(lower, low)
                    self.assertLess(high, upper)
                else:
                    self.assertLess(high, lower)

    def test_planted_root_is_recovered_at_every_width(self) -> None:
        positive = [row for row in self.result["trials"] if row["polarity"] == "positive"]
        self.assertEqual(len(positive), 4)
        for row in positive:
            with self.subTest(width=row["width_id"], planted=row["planted_coefficients_ascending"]):
                self.assertTrue(row["planted_quartic_detected"])
                self.assertIn(row["planted_coefficients_ascending"], row["witness_coefficients"])
                self.assertGreaterEqual(row["root_containing_polynomials"], 1)

    def test_planted_root_is_not_reported_when_absent(self) -> None:
        negative = [row for row in self.result["trials"] if row["polarity"] == "negative"]
        self.assertEqual(len(negative), 4)
        for row in negative:
            with self.subTest(width=row["width_id"], planted=row["planted_coefficients_ascending"]):
                self.assertFalse(row["planted_quartic_detected"])
                self.assertNotIn(row["planted_coefficients_ascending"], row["witness_coefficients"])

    def test_narrowest_widths_are_sensitive_and_specific(self) -> None:
        """The two intervals the census excludes must be the ones proved sensitive."""
        for width_id in ("jacobsen-2015-eigenvalue", "yang-zhou-2024-corrected"):
            rows = {row["polarity"]: row for row in self.result["trials"]
                    if row["width_id"] == width_id and row["planted_coefficients_ascending"] == [-84, 99, -7, 99, 58]}
            with self.subTest(width=width_id):
                self.assertEqual(rows["positive"]["root_containing_polynomials"], 1)
                self.assertEqual(rows["negative"]["root_containing_polynomials"], 0)

    def test_conclusion_matches_the_trials(self) -> None:
        conclusion = self.result["conclusion"]
        self.assertTrue(conclusion["all_trials_passed"])
        self.assertTrue(conclusion["all_positive_trials_detected_the_planted_quartic"])
        self.assertTrue(conclusion["all_negative_trials_excluded_the_planted_quartic"])
        self.assertEqual(conclusion["positive_trials"], conclusion["negative_trials"])
        self.assertTrue(all(row["passed"] for row in self.result["trials"]))

    def test_claim_boundary_refuses_promotion_of_the_planted_quartics(self) -> None:
        excluded = self.result["claim_boundary"]["excluded"]
        for forbidden in ("candidate formulas", "p-value", "transcendence", "degree/height expansion"):
            self.assertIn(forbidden, excluded)


if __name__ == "__main__":
    unittest.main()
