#!/usr/bin/env python3
"""Tests for the #612 corrections: the chart identity and the N725 score.

Two kinds of test, kept strictly apart:

* **ArtifactLock** — the committed result JSONs are pinned. These do not
  recompute anything; they stop a future refactor from silently moving a
  published number, and they pin the numbers #609/#612 pre-registered so that
  a re-run has to agree with them.
* **PureMachinery** — the estimator conventions themselves, on small synthetic
  problems with no histogram loading: the GLS covector really is the dual
  basis, the second-difference weights really annihilate constants, and the
  three-standard-error decision rule really sorts intervals the way #612
  specifies.
"""

from __future__ import annotations

import json
import math
import random
import unittest
from pathlib import Path

try:  # pragma: no cover
    from scripts.p612_chart_identity import gls_fit
    from scripts.p612_n725_score import decision
    from scripts.p582_amplitude_law import second_difference_weights
except ModuleNotFoundError:  # pragma: no cover
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    from p612_chart_identity import gls_fit
    from p612_n725_score import decision
    from p582_amplitude_law import second_difference_weights

ROOT = Path(__file__).resolve().parents[1]
CHART = ROOT / "results" / "p612-chart-identity" / "latest.json"
SCORE = ROOT / "results" / "p612-n725-score" / "latest.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text())


class ArtifactLockChartIdentity(unittest.TestCase):
    """The chart identity, and the ratios #609 published."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.report = _load(CHART)

    def test_the_identity_closes_inside_the_declared_tolerance(self) -> None:
        for name, block in self.report["weightings"].items():
            self.assertTrue(block["every_identity_closes"], name)
            self.assertLess(block["largest_absolute_identity_error"], 1e-9, name)

    def test_the_published_curvature_ratios_are_reproduced(self) -> None:
        """#609's 1.5368 / 1.5574 and the equal-weighting control."""
        expected = {"spin0": [1.5368, 1.5574], "equal": [1.5789, 1.4438]}
        for name, expected_ratios in expected.items():
            block = self.report["weightings"][name]
            measured = [entry["decomposition"]["ratio_measured_over_one_exponent"]
                        for entry in block["per_lineage"] if entry.get("usable")]
            self.assertEqual(len(measured), 2, name)
            for mine, theirs in zip(measured, expected_ratios):
                self.assertAlmostEqual(mine, theirs, places=3,
                                       msg=f"{name}: {mine} vs {theirs}")

    def test_the_chart_term_carries_most_of_the_excess(self) -> None:
        """#612 asked for roughly 95-98% on the primary weighting."""
        for name, low, high in (("spin0", 0.90, 1.02), ("equal", 0.85, 1.20)):
            for fraction in self.report["weightings"][name]["chart_fraction_of_excess"]:
                self.assertGreater(fraction, low, name)
                self.assertLess(fraction, high, name)

    def test_the_width_contraction_is_below_one_in_every_chart(self) -> None:
        """k = 1 + h0*beta0 < 1 is what makes the transport non-trivial."""
        for block in self.report["weightings"].values():
            for entry in block["per_lineage"]:
                if entry.get("usable"):
                    self.assertLess(entry["k"], 1.0)
                    self.assertGreater(entry["k"], 0.5)


class ArtifactLockN725Score(unittest.TestCase):
    """The one N725 block, and the numbers #612 stated in its own body."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.report = _load(SCORE)

    def test_the_declared_forecast_is_reproduced_by_the_frozen_fit(self) -> None:
        new = self.report["weightings"]["spin0"]["new_transition"]
        self.assertAlmostEqual(new["predicted_from_the_five_amplitude_fit"],
                               new["forecast_target"], delta=1e-9)

    def test_the_pre_registered_constants_are_recovered(self) -> None:
        """#612 states k=0.7767512 and the chart-corrected 1.03047842e-03."""
        curvature = self.report["weightings"]["spin0"]["p50_curvature"]
        self.assertAlmostEqual(curvature["k_145_to_290"], 0.7767512, places=6)
        self.assertAlmostEqual(curvature["predictions"]["chart_corrected"],
                               1.03047842e-03, delta=1e-8)
        self.assertAlmostEqual(curvature["predictions"]["naive_one_exponent"],
                               6.70567e-04, delta=1e-9)

    def test_the_chart_calculation_does_not_separate_the_two_targets(self) -> None:
        """#612: they differ by about 0.7%, not by 55%."""
        predictions = self.report["weightings"]["spin0"]["p50_curvature"]["predictions"]
        relative = abs(predictions["chart_corrected_minus_factor_1547_target"]) \
            / predictions["factor_1547_target_from_609"]
        self.assertLess(relative, 0.02)
        self.assertTrue(predictions["the_two_targets_are_separated_by_the_chart"])

    def test_the_pre_registered_second_difference_weights_are_reproduced(self) -> None:
        for block in self.report["weightings"].values():
            self.assertTrue(block["p50_curvature"]["weights_agree_with_pre_registration"])
        weights = self.report["weightings"]["spin0"]["p50_curvature"][
            "second_difference_weights"]
        self.assertAlmostEqual(sum(weights), 0.0, places=9)

    def test_the_primary_verdict_is_support_and_the_sensitivity_is_not(self) -> None:
        self.assertEqual(self.report["weightings"]["spin0"]["forecast_verdict"]["outcome"],
                         "supports_this_finite_forecast")
        self.assertEqual(self.report["weightings"]["equal"]["forecast_verdict"]["outcome"],
                         "stops_this_forecast")

    def test_the_point_forecast_miss_is_reported_and_is_more_than_three_se(self) -> None:
        """The pass is a 5% tolerance pass, not a nominal-error pass."""
        verdict = self.report["weightings"]["spin0"]["forecast_verdict"]
        self.assertGreater(abs(verdict["displacement_in_standard_errors"]), 3.0)
        self.assertLess(abs(verdict["relative_displacement"]), 0.05)

    def test_the_curvature_identity_also_closes_on_the_new_rung(self) -> None:
        for block in self.report["weightings"].values():
            self.assertLess(block["p50_curvature"]["identity"]["absolute_error"], 1e-9)

    def test_the_residual_transport_is_negligible_but_not_zero(self) -> None:
        identity = self.report["weightings"]["spin0"]["p50_curvature"]["identity"]
        curvature = self.report["weightings"]["spin0"]["p50_curvature"]["measured"]
        self.assertNotEqual(identity["residual_transport"], 0.0)
        self.assertLess(abs(identity["residual_transport"]) / abs(curvature), 0.01)

    def test_725_breaks_the_committed_label_identification(self) -> None:
        """#609's design claim: v_5 = 2 and an interpolating spin-0 combination."""
        crossing = self.report["weightings"]["spin0"]["label_crossing"]
        self.assertEqual(crossing["five_adic_valuation"], 2)
        self.assertTrue(crossing["spin0_is_an_interpolation"])
        self.assertTrue(crossing["breaks_the_committed_identification"])
        weights = crossing["spin0_weights"]
        self.assertAlmostEqual(sum(weights.values()), 1.0, places=12)

    def test_cross_size_coupling_was_measured_not_assumed(self) -> None:
        coupling = self.report["weightings"]["spin0"]["cross_size_coupling_290_725"]
        self.assertEqual(coupling["batches_compared"], 100)
        self.assertTrue(coupling["cross_term_droppable"])

    def test_the_sample_count_is_the_corrected_one(self) -> None:
        production = self.report["production"]
        self.assertEqual(production["batches"], 100)
        self.assertEqual(production["samples_per_batch"], 1_000_000)
        self.assertEqual(production["samples_per_orientation"], 100_000_000)

    def test_the_run_is_labelled_a_reanalysis(self) -> None:
        self.assertIn("reanalysis", self.report["design"])


class PureMachinery(unittest.TestCase):
    """The estimator conventions, on synthetic problems."""

    def test_the_covectors_are_the_dual_basis(self) -> None:
        rng = random.Random(20260907)
        size, width = 6, 3
        basis = [[rng.uniform(-1.0, 1.0) for _ in range(size)] for _ in range(width)]
        covariance = [[0.0] * size for _ in range(size)]
        for i in range(size):
            covariance[i][i] = 1.0 + 0.5 * i
            for j in range(i + 1, size):
                value = 0.1 * math.exp(-abs(i - j))
                covariance[i][j] = value
                covariance[j][i] = value
        fit = gls_fit([rng.uniform(-1.0, 1.0) for _ in range(size)],
                      covariance, basis)
        for i in range(width):
            for j in range(width):
                applied = math.fsum(c * v for c, v in zip(fit["covectors"][i], basis[j]))
                self.assertAlmostEqual(applied, 1.0 if i == j else 0.0, places=9)

    def test_a_vector_in_the_span_is_fit_exactly(self) -> None:
        rng = random.Random(11)
        size = 5
        basis = [[1.0] * size, [rng.uniform(0.0, 1.0) for _ in range(size)]]
        covariance = [[1.0 if i == j else 0.0 for j in range(size)] for i in range(size)]
        coefficients = [0.7, -1.3]
        observation = [math.fsum(c * b[level] for c, b in zip(coefficients, basis))
                       for level in range(size)]
        fit = gls_fit(observation, covariance, basis)
        for mine, theirs in zip(fit["amplitudes"], coefficients):
            self.assertAlmostEqual(mine, theirs, places=9)
        for value in fit["residual"]:
            self.assertAlmostEqual(value, 0.0, places=9)

    def test_the_second_difference_weights_annihilate_constants(self) -> None:
        weights = second_difference_weights([0.0, math.log(2.0), math.log(5.0)])
        self.assertAlmostEqual(sum(weights), 0.0, places=12)

    def test_the_decision_rule_sorts_intervals_as_specified(self) -> None:
        band = (-0.0004914884611524859, -0.0004446800362808205)
        se = 2.3e-06
        # well inside
        self.assertEqual(decision(-0.000468, se)["outcome"],
                         "supports_this_finite_forecast")
        # well outside, above the band
        self.assertEqual(decision(-0.000420, se)["outcome"], "stops_this_forecast")
        # straddling the lower boundary
        self.assertEqual(decision(band[0], se)["outcome"], "unresolved")
        # straddling the upper boundary
        self.assertEqual(decision(band[1], se)["outcome"], "unresolved")
        # an absurdly precise measurement displaced from the point forecast
        # still reads as "supports" while it is inside the band -- the rule is a
        # tolerance check, not a consistency test.
        verdict = decision(-0.000480, 1.0e-09)
        self.assertEqual(verdict["outcome"], "supports_this_finite_forecast")
        self.assertGreater(abs(verdict["displacement_in_standard_errors"]), 3.0)


if __name__ == "__main__":
    unittest.main()
