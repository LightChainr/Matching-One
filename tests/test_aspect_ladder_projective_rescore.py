#!/usr/bin/env python3
"""Tests for the denominator-free rescore of the N=580 ladder.

The wrong numbers here are the ones the ratio test produced, and one new way to
get them wrong: filling in a covariance entry nobody measured.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import aspect_ladder_projective_rescore as rescore  # noqa: E402


class CurvatureFunctionalTests(unittest.TestCase):
    def test_it_is_the_second_divided_difference_and_nothing_else(self) -> None:
        """Anchor: on r^2 it must return exactly 1, on any line exactly 0.

        The wrong number this stops us believing is a "curvature" that is
        really a rescaled first difference, which would be nonzero for the two
        families that are linear in r and would then appear to exclude them for
        free.
        """
        identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        square = rescore.curvature([1.0, 4.0, 16.0], identity)
        self.assertAlmostEqual(square["value"], 1.0, places=12)
        for intercept, slope in ((0.0, 1.0), (5.0, 0.0), (-2.0, 3.5)):
            line = [intercept + slope * r for r in (1.0, 2.0, 4.0)]
            self.assertAlmostEqual(
                rescore.curvature(line, identity)["value"], 0.0, places=12
            )

    def test_every_frozen_competitor_predicts_a_non_negative_curvature(self) -> None:
        """The claim that makes the measured sign interesting.

        If some competitor were concave, a negative measurement would be
        evidence *for* it rather than against the whole list, and the headline
        would be wrong. This checks the list, not our memory of it.
        """
        identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        for name, ray in rescore.load_competitors().items():
            predicted = rescore.curvature(list(ray), identity)["value"]
            self.assertGreaterEqual(predicted, -1e-12, f"{name} is concave in r")

    def test_the_measurement_is_negative_across_the_unknown_covariance(self) -> None:
        """Pins the finding, and pins how far it is from being decisive.

        The wrong number is a sign that turns out to depend on the covariance
        entry the first scoring run threw away. It does not; the significance
        does.
        """
        response = rescore.load_response()
        report = rescore.rescore(rescore.load_competitors(), response)
        low, high = report["curvature"]["z_over_admissible_correlations"]
        self.assertLess(high, 0.0, "the sign must not flip anywhere in the range")
        self.assertLess(report["curvature"]["value"], 0.0)
        # Reported honestly: not decisive at one end of the range.
        self.assertLess(abs(low), 4.0)
        self.assertGreater(abs(high), 2.0)


class MissingCovarianceTests(unittest.TestCase):
    def test_an_unmeasured_covariance_entry_is_declared_not_assumed(self) -> None:
        """Stops a zero being read as a measurement.

        The committed artifact has no cov(r2, r4). Filling it with zero is a
        choice, and an undeclared choice is the same failure as the ratio test:
        a number nobody checked doing load-bearing work.
        """
        response = rescore.load_response()
        self.assertFalse(response["covariance_is_complete"])
        self.assertEqual(response["covariance"][1][2], 0.0)
        report = rescore.rescore(rescore.load_competitors(), response)
        self.assertFalse(report["covariance_is_complete"])
        for name, row in report["competitors"].items():
            low, high = row["sigma_over_admissible_correlations"]
            self.assertLessEqual(low, row["equivalent_sigma"] + 1e-9, name)
            self.assertGreaterEqual(high, row["equivalent_sigma"] - 1e-9, name)

    def test_exactly_one_verdict_depends_on_the_missing_entry(self) -> None:
        """Says what a replay would buy, so nobody pays for it twice.

        Seven of the eight verdicts hold across every positive-definite value of
        cov(r2, r4). Only bare_aspect_ratio moves, so that entry is worth one
        deterministic replay and nothing more.
        """
        report = rescore.rescore(rescore.load_competitors(), rescore.load_response())
        unstable = sorted(
            name for name, row in report["competitors"].items()
            if not row["verdict_survives_the_missing_covariance"]
        )
        self.assertEqual(unstable, ["bare_aspect_ratio"])

    def test_the_admissible_range_really_is_positive_definite(self) -> None:
        """The sweep must not wander outside the set of real covariances.

        A correlation the data could never have produced would let the sweep
        manufacture either significance or its absence.
        """
        covariance = rescore.load_response()["covariance"]
        low, high = rescore.admissible_correlation_range(covariance)
        scale = math.sqrt(covariance[1][1] * covariance[2][2])
        for correlation in (low, 0.5 * (low + high), high):
            filled = rescore._with_corr(covariance, correlation, scale)
            self.assertGreater(rescore._determinant(filled), 0.0)
        for outside in (low - 0.05, high + 0.05):
            filled = rescore._with_corr(covariance, outside, scale)
            self.assertLess(rescore._determinant(filled), 0.0)


class Spin8ReconciliationTests(unittest.TestCase):
    def test_a_model_through_the_data_needs_no_spin8(self) -> None:
        """Anchor for the reconciliation arithmetic.

        The wrong number is a required |A8/A4| that is nonzero even for a model
        that already fits, which would make every entry in that column look
        alarming for a reason that has nothing to do with the data.
        """
        response = rescore.load_response()
        vector = response["vector"]
        exact = [1.0, vector[1] / vector[0], vector[2] / vector[0]]
        self.assertLess(
            rescore.required_spin8_ratio(vector, response["covariance"], exact), 1e-9
        )

    def test_every_survivor_of_the_clean_pair_needs_a_large_spin8(self) -> None:
        """The dichotomy the run actually produced.

        no_modulus_dependence is the only competitor the r=2 rung could accept
        within the assumed bound, and it is the one the clean r=1 / r=4 pair
        excludes outright. Everything else needs |A8/A4| far above 1 -- which
        is the assumption the whole two-orientation estimator rests on.
        """
        report = rescore.rescore(rescore.load_competitors(), rescore.load_response())
        cheap = {name for name, row in report["competitors"].items()
                 if row["required_abs_A8_over_A4_to_reach_r2"] < 1.0}
        self.assertEqual(cheap, {"no_modulus_dependence"})
        self.assertTrue(report["competitors"]["no_modulus_dependence"]["excluded_at_3_sigma"])
        for name, row in report["competitors"].items():
            if name != "no_modulus_dependence":
                self.assertGreater(row["required_abs_A8_over_A4_to_reach_r2"], 5.0, name)

    def test_the_bound_provenance_is_recorded_rather_than_asserted(self) -> None:
        """Keeps the weak link visible.

        The entire dichotomy rests on |A8/A4| being small, and that is not a
        measurement anywhere in this repository. If someone later writes it down
        as one, this is where the trade has to happen.
        """
        report = rescore.rescore(rescore.load_competitors(), rescore.load_response())
        provenance = report["spin8_bound_provenance"]
        self.assertIn("H4 0.4163/2", provenance["traces_to"])
        self.assertIn("not a measurement", provenance["what_that_actually_is"])
        self.assertIn("N=650", report["what_this_does_not_separate"])


class CommittedArtifactTests(unittest.TestCase):
    def test_the_committed_rescore_reproduces(self) -> None:
        committed = json.loads(
            (ROOT / "results" / "aspect-ladder-n580-projective" / "latest.json")
            .read_text(encoding="utf-8")
        )
        fresh = rescore.rescore(rescore.load_competitors(), rescore.load_response())
        self.assertEqual(json.loads(json.dumps(fresh)), committed)

    def test_the_competitors_come_from_the_freeze(self) -> None:
        """Stops the list drifting between the freeze and the rescore.

        The wrong result is one computed against a competitor set that quietly
        gained or lost a member after the run.
        """
        frozen = json.loads(json.dumps(rescore.load_competitors()))
        self.assertEqual(len(frozen), 8)
        self.assertEqual(frozen["q4_jordan_weight4"], [1.0, 2.75, 10.9908008589])
        self.assertEqual(frozen["bare_aspect_ratio"], [1.0, 2.0, 4.0])
        for ray in frozen.values():
            self.assertEqual(ray[0], 1.0)


if __name__ == "__main__":
    unittest.main()
