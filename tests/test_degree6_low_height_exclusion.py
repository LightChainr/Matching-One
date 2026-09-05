#!/usr/bin/env python3
"""Lock the degree-1..6 exclusion at both heights, and its sensitivity control."""

from __future__ import annotations

from fractions import Fraction
import json
from math import gcd
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from degree6_low_height_exclusion import (  # noqa: E402
    DEGREE_MAX,
    HEIGHT,
    build_result,
    derivative_bound,
    enumerate_class,
    output_path,
    validate_result,
)

HEIGHT4 = 4
HEIGHT4_CLASS_SIZES = {1: 23, 2: 265, 3: 2639, 4: 24913, 5: 229703, 6: 2093785}
from degree6_low_height_control import (  # noqa: E402
    DEFAULT_OUTPUT as CONTROL_OUTPUT,
    build_result as build_control,
    planted_polynomial,
    validate_result as validate_control,
)
from pslq_look_elsewhere_ledger import primitive_polynomial_count  # noqa: E402

CONTRACT = json.loads((ROOT / "analysis" / "pslq_search_contract.json").read_text(encoding="utf-8"))
INTERVAL_IDS = [row["id"] for row in CONTRACT["intervals"]]


class Degree6LowHeightExclusionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.results = {interval_id: build_result(interval_id) for interval_id in INTERVAL_IDS}

    def test_committed_artifacts_reproduce_exactly(self) -> None:
        for interval_id, result in self.results.items():
            committed = json.loads(output_path(interval_id).read_text(encoding="utf-8"))
            with self.subTest(interval=interval_id):
                self.assertEqual(committed, result)
                validate_result(committed)

    def test_enumeration_is_primitive_sign_normalized_and_exact_degree(self) -> None:
        for degree in range(1, DEGREE_MAX + 1):
            members = list(enumerate_class(degree))
            with self.subTest(degree=degree):
                self.assertEqual(len(members), primitive_polynomial_count(degree, HEIGHT))
                self.assertEqual(len(set(members)), len(members))
                for coefficients in members:
                    self.assertEqual(len(coefficients) - 1, degree)
                    self.assertGreaterEqual(coefficients[-1], 1)
                    self.assertLessEqual(max(abs(v) for v in coefficients), HEIGHT)
                    common = 0
                    for value in coefficients:
                        common = gcd(common, value)
                    self.assertEqual(common, 1)

    def test_known_thresholds_are_inside_the_enumerated_class(self) -> None:
        """The class must actually contain the forms it claims to cover."""
        native = json.loads((ROOT / "results" / "pslq-lattice-native-candidates" / "latest.json").read_text(encoding="utf-8"))
        for row in native["candidates"]:
            coefficients = tuple(row["minimal_polynomial_coefficients_ascending"])
            degree = len(coefficients) - 1
            with self.subTest(candidate=row["candidate_id"]):
                self.assertLessEqual(degree, DEGREE_MAX)
                self.assertLessEqual(max(abs(v) for v in coefficients), HEIGHT)
                self.assertIn(coefficients, set(enumerate_class(degree)))

    def test_every_interval_is_excluded_at_every_degree(self) -> None:
        for interval_id, result in self.results.items():
            row = result["interval_result"]
            with self.subTest(interval=interval_id):
                self.assertTrue(row["excluded"])
                self.assertEqual(row["degrees_excluded"], list(range(1, DEGREE_MAX + 1)))
                self.assertEqual(row["degrees_with_roots"], [])
                self.assertEqual(row["polynomials_per_interval"], 409584)
                for entry in row["by_degree"]:
                    self.assertEqual(entry["root_containing_polynomials"], 0)
                    self.assertEqual(entry["root_witnesses"], [])

    def test_screen_bound_is_certified_and_consistent(self) -> None:
        """Every reported closest polynomial must clear the screen it was judged by."""
        for interval_id, result in self.results.items():
            row = result["interval_result"]
            lower, upper = Fraction(row["lower"]), Fraction(row["upper"])
            width = upper - lower
            for entry in row["by_degree"]:
                bound = entry["derivative_bound_on_unit_interval"]
                residual = Fraction(entry["closest_polynomial"]["minimum_absolute_endpoint_residual"])
                floor = Fraction(entry["root_distance_lower_bound_text"])
                with self.subTest(interval=interval_id, degree=entry["degree"]):
                    self.assertEqual(bound, derivative_bound(entry["degree"]))
                    self.assertEqual(floor, residual / bound)
                    self.assertEqual(Fraction(entry["floor_to_interval_width_ratio_text"]), floor / width)
                    # zero screened candidates <=> nothing came within D*(u-l) of the interval
                    if entry["screen_candidates_exactly_decided"] == 0:
                        self.assertGreater(residual, bound * width)

    def test_low_degree_rows_agree_with_the_merged_height_100_census(self) -> None:
        """Degrees 1-4 here are a strict subclass of the merged census; both must be null."""
        for interval_id, result in self.results.items():
            for entry in result["interval_result"]["by_degree"]:
                if entry["degree"] > 4:
                    continue
                with self.subTest(interval=interval_id, degree=entry["degree"]):
                    self.assertEqual(entry["root_containing_polynomials"], 0)
        # every committed height-100 quartic survivor is far outside the height bound tested here
        for interval_id in INTERVAL_IDS:
            census = json.loads((ROOT / "results" / f"pslq-degree4-{interval_id}" / "latest.json").read_text(encoding="utf-8"))
            for witness in census["interval_result"]["root_witnesses"]:
                height = max(abs(v) for v in witness["coefficients_ascending"])
                self.assertGreater(height, HEIGHT, witness["coefficients_ascending"])

    def test_claim_boundary_refuses_transcendence_and_expansion(self) -> None:
        for result in self.results.values():
            excluded = result["claim_boundary"]["excluded"]
            for forbidden in ("transcendence", "p-values", "closed forms", "higher degree or height"):
                self.assertIn(forbidden, excluded)


class Degree6LowHeightControlTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = build_control()

    def test_committed_artifact_reproduces_exactly(self) -> None:
        committed = json.loads(CONTROL_OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(committed, self.result)
        validate_control(committed)

    def test_control_uses_the_unmodified_scan_path(self) -> None:
        path = self.result["scan_path"]
        self.assertFalse(path["modified_for_this_control"])
        self.assertEqual(path["function"], "scripts/degree6_low_height_exclusion.py::_scan_degree")

    def test_planted_form_is_the_certified_three_twelve_site_polynomial(self) -> None:
        import hashlib

        planted = self.result["planted"]
        self.assertEqual(planted["coefficients_ascending"], [1, 0, 0, 0, -3, 0, 1])
        self.assertEqual(planted["coefficients_ascending"], list(planted_polynomial()))
        self.assertEqual(planted["degree"], 6)
        self.assertEqual(planted["height"], 3)
        source = ROOT / planted["source"]
        self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(), planted["source_sha256"])

    def test_planted_root_bracket_is_the_three_twelve_site_value(self) -> None:
        low, high = (Fraction(v) for v in self.result["planted"]["root_bracket"])
        self.assertLess(low, high)
        # sqrt(1 - 2 sin(pi/18)) = 0.80790...
        self.assertLess(Fraction(807900, 10**6), low)
        self.assertLess(high, Fraction(807901, 10**6))

    def test_every_frozen_width_is_covered_at_both_polarities(self) -> None:
        widths = {row["id"]: Fraction(row["upper"]) - Fraction(row["lower"]) for row in CONTRACT["intervals"]}
        self.assertEqual(set(self.result["widths_tested"]), set(widths))
        for width_id, width_text in self.result["widths_tested"].items():
            self.assertEqual(Fraction(width_text), widths[width_id])
        covered = {(row["width_id"], row["polarity"]) for row in self.result["trials"]}
        self.assertEqual(covered, {(w, p) for w in widths for p in ("positive", "negative")})

    def test_trial_geometry_matches_its_declared_polarity(self) -> None:
        low, high = (Fraction(v) for v in self.result["planted"]["root_bracket"])
        for row in self.result["trials"]:
            lower, upper = Fraction(row["lower"]), Fraction(row["upper"])
            with self.subTest(width=row["width_id"], polarity=row["polarity"]):
                self.assertEqual(upper - lower, Fraction(row["width_text"]))
                if row["polarity"] == "positive":
                    self.assertLess(lower, low)
                    self.assertLess(high, upper)
                else:
                    self.assertLess(high, lower)

    def test_the_scan_finds_a_planted_degree_six_root_and_only_then(self) -> None:
        conclusion = self.result["conclusion"]
        self.assertTrue(conclusion["all_trials_passed"])
        self.assertTrue(conclusion["all_positive_trials_detected_the_planted_polynomial"])
        self.assertTrue(conclusion["all_negative_trials_excluded_the_planted_polynomial"])
        self.assertEqual(conclusion["positive_trials"], 4)
        self.assertEqual(conclusion["negative_trials"], 4)
        for row in self.result["trials"]:
            with self.subTest(width=row["width_id"], polarity=row["polarity"]):
                self.assertGreaterEqual(row["screen_candidates_exactly_decided"], 1)
                if row["polarity"] == "positive":
                    self.assertEqual(row["root_containing_polynomials"], 1)
                    self.assertIn([1, 0, 0, 0, -3, 0, 1], row["witness_coefficients"])
                else:
                    self.assertEqual(row["root_containing_polynomials"], 0)

    def test_claim_boundary_refuses_promoting_the_planted_form(self) -> None:
        excluded = self.result["claim_boundary"]["excluded"]
        for forbidden in ("candidate for square-site p_c", "p-value", "transcendence", "degree/height expansion"):
            self.assertIn(forbidden, excluded)


class Degree6Height4ExclusionTests(unittest.TestCase):
    """The corrected historical range (manuscript section 6.6).

    The wrong number these tests exist to stop us believing is the exclusion
    verdict on a class that is 5.7 times larger than the height-3 one.  The
    A-lattice quintic has height 4, so this is the class that actually covers
    the exactly-known planar thresholds, and the paper's closing sentence now
    rests on it rather than on the height-3 result.
    """

    NARROWEST = "jacobsen-2015-eigenvalue"

    @classmethod
    def setUpClass(cls) -> None:
        # One interval is rebuilt from scratch; the other three are checked
        # against their own committed contents.  A full four-interval rebuild is
        # 20 s, and the Python 3.9 CI job has no room for it -- the same
        # scoping the replication test of #551 uses, and for the same reason.
        cls.rebuilt = build_result(cls.NARROWEST, height=HEIGHT4)
        cls.committed = {
            interval_id: json.loads(output_path(interval_id, HEIGHT4).read_text(encoding="utf-8"))
            for interval_id in INTERVAL_IDS
        }

    def test_the_rebuilt_interval_reproduces_its_committed_artifact(self) -> None:
        # No validate_result() here: it would rebuild the same interval a second
        # time for the same check.
        self.assertEqual(self.committed[self.NARROWEST], self.rebuilt)
        self.assertEqual(self.rebuilt["search"]["coefficient_height_max"], HEIGHT4)

    def test_every_interval_is_excluded_at_every_degree(self) -> None:
        for interval_id, artifact in self.committed.items():
            row = artifact["interval_result"]
            with self.subTest(interval=interval_id):
                self.assertTrue(row["excluded"])
                self.assertEqual(row["degrees_excluded"], list(range(1, DEGREE_MAX + 1)))
                self.assertEqual(row["degrees_with_roots"], [])
                self.assertEqual(row["polynomials_per_interval"], 2351328)
                for entry in row["by_degree"]:
                    self.assertEqual(entry["root_containing_polynomials"], 0)
                    self.assertEqual(entry["root_witnesses"], [])
                    self.assertEqual(entry["screen_candidates_exactly_decided"], 0)

    def test_class_sizes_are_the_independently_counted_ones(self) -> None:
        """A silently smaller class would exclude by not looking.

        The counts on the right are from a direct brute-force count of
        primitive sign-normalized tuples at height 4, done outside this code
        path; the enumerator has to agree with them degree by degree.
        """
        for degree, expected in HEIGHT4_CLASS_SIZES.items():
            with self.subTest(degree=degree):
                self.assertEqual(len(list(enumerate_class(degree, HEIGHT4))), expected)
        for interval_id, artifact in self.committed.items():
            rows = artifact["interval_result"]["by_degree"]
            with self.subTest(interval=interval_id):
                self.assertEqual({row["degree"]: row["polynomials_in_class"] for row in rows},
                                 HEIGHT4_CLASS_SIZES)

    def test_the_a_lattice_quintic_is_inside_this_class_and_outside_the_old_one(self) -> None:
        """The one thing height 4 was raised to cover.

        If this polynomial were not in the enumerated class, section 6.6 would
        be claiming to close the historical range while still missing the row
        that showed it was open.
        """
        a_lattice = tuple(json.loads(
            (ROOT / "results" / "ziff-a-lattice-complexity" / "latest.json").read_text(encoding="utf-8")
        )["polynomial_ascending"])
        self.assertEqual(max(abs(v) for v in a_lattice), HEIGHT4)
        self.assertIn(a_lattice, set(enumerate_class(len(a_lattice) - 1, HEIGHT4)))
        self.assertNotIn(a_lattice, set(enumerate_class(len(a_lattice) - 1, HEIGHT)))

    def test_the_height_3_artifacts_are_untouched_by_the_height_parameter(self) -> None:
        """Threading a parameter through must not move the committed result.

        The height-3 artifacts are cited by the manuscript and pinned by digest
        in the implementation-agreement comparison; a default-argument slip
        here would rewrite them under a green suite.
        """
        for interval_id in INTERVAL_IDS:
            self.assertEqual(output_path(interval_id, HEIGHT), output_path(interval_id))
            self.assertNotEqual(output_path(interval_id, HEIGHT4), output_path(interval_id))

    def test_the_screen_bound_scales_with_the_height(self) -> None:
        for degree in range(1, DEGREE_MAX + 1):
            with self.subTest(degree=degree):
                self.assertEqual(derivative_bound(degree, HEIGHT4), 2 * degree * (degree + 1))
                self.assertEqual(derivative_bound(degree, HEIGHT4),
                                 (HEIGHT4 * derivative_bound(degree, HEIGHT)) // HEIGHT)


if __name__ == "__main__":
    unittest.main()
