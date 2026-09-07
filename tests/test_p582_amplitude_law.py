"""Tests for the amplitude law of #582's frozen shape direction (#584).

Every test below names the wrong number it would stop us believing.  The
headline claims are an exponent measured on five amplitudes, a held-out
prediction, and an independent second-difference check that falsifies the
one-exponent model; each of those has a way of being wrong that arithmetic can
catch, and this file catches it.
"""

from __future__ import annotations

import json
import math
import unittest
from pathlib import Path

from scripts.p582_amplitude_law import (CANDIDATE_SIZES, DEFAULT_OUTPUT,
                                        OMEGA_WINDOW, candidate_sizes, decide,
                                        error_budget, fit_exponent,
                                        first_difference_image,
                                        five_adic_valuation,
                                        leave_one_transition_out,
                                        primitive_representatives,
                                        second_difference_image,
                                        second_difference_weights)

#: The real log geometry of #582's five transitions.  Synthetic amplitudes are
#: placed on it so that a recovery test exercises the same unequal steps and
#: midpoints the production fit sees.
GEOMETRY = [
    ("65->130", (math.log(65) + math.log(130)) / 2, math.log(2.0)),
    ("130->325", (math.log(130) + math.log(325)) / 2, math.log(2.5)),
    ("85->170", (math.log(85) + math.log(170)) / 2, math.log(2.0)),
    ("170->425", (math.log(170) + math.log(425)) / 2, math.log(2.5)),
    ("145->290", (math.log(145) + math.log(290)) / 2, math.log(2.0)),
]


def synthetic(scale: float, omega: float, sigma: float = 1e-9):
    return [{"label": label, "sbar": sbar, "h": step,
             "amplitude": first_difference_image(scale, omega, sbar, step),
             "standard_error": sigma}
            for label, sbar, step in GEOMETRY]


class FiniteDifferenceImage(unittest.TestCase):
    def test_the_image_is_the_difference_and_not_the_derivative(self):
        """Stops us believing an amplitude fitted with the derivative shortcut.

        The wrong number is a 1.9% error at ``h = log 2`` and 3.3% at
        ``h = log 2.5``.  Those are step-size dependent and the same size as the
        misfit being measured, so dropping ``sinh(x)/x`` would manufacture a
        difference between the m=2 and m=2.5 transitions -- exactly the
        `multiplier` label Gate 3 screened.
        """
        omega, sbar = 0.97, 5.0
        for step, expected in ((math.log(2.0), 0.018942), (math.log(2.5), 0.033242)):
            derivative = -omega * math.exp(-omega * sbar)
            image = first_difference_image(1.0, omega, sbar, step)
            self.assertAlmostEqual(image / derivative - 1.0, expected, places=5)

    def test_the_image_equals_the_exact_difference_of_the_exponential(self):
        """Stops us believing a model that is not the thing #582 measures.

        The wrong number is any image that disagrees with ``[f(s_t) - f(s_b)]/h``
        evaluated directly on ``lambda exp(-omega s)``.
        """
        scale, omega = 0.1783, 0.9702
        for _, sbar, step in GEOMETRY:
            base, target = sbar - step / 2, sbar + step / 2
            direct = (scale * math.exp(-omega * target)
                      - scale * math.exp(-omega * base)) / step
            self.assertAlmostEqual(first_difference_image(scale, omega, sbar, step),
                                   direct, places=15)

    def test_a_zero_step_is_refused(self):
        """Stops us believing an amplitude read off a transition with no lever arm."""
        with self.assertRaises(ValueError):
            first_difference_image(1.0, 1.0, 5.0, 0.0)


class SecondDifferenceWeights(unittest.TestCase):
    def test_the_weights_annihilate_constants(self):
        """Stops us believing a curvature that the limit shape leaked into.

        The wrong number is a nonzero weight sum, which would put the whole
        threshold law ``A(u)`` -- of order 0.5 -- into a curvature amplitude of
        order 2e-3, swamping it by three orders of magnitude.
        """
        nodes = [math.log(65), math.log(130), math.log(325)]
        self.assertAlmostEqual(sum(second_difference_weights(nodes)), 0.0, places=14)

    def test_the_weights_return_two_on_s_squared(self):
        """Stops us believing weights that are off by the unequal-step factor.

        ``2 f[s0,s1,s2]`` of ``s^2`` is exactly 2 on any grid.  The wrong number
        is the equal-step formula ``(f0 - 2f1 + f2)/h^2`` applied to steps that
        are not equal here (log 2 against log 2.5), which returns 1.87.
        """
        nodes = [math.log(65), math.log(130), math.log(325)]
        weights = second_difference_weights(nodes)
        self.assertAlmostEqual(sum(w * s * s for w, s in zip(weights, nodes)),
                               2.0, places=12)
        equal_step = (nodes[0] - 2 * nodes[1] + nodes[2]) / ((nodes[1] - nodes[0]) ** 2)
        self.assertNotAlmostEqual(equal_step, 2.0, places=1)

    def test_the_second_difference_image_matches_direct_evaluation(self):
        """Stops us believing a prediction that is not the model's own image."""
        nodes = [math.log(85), math.log(170), math.log(425)]
        scale, omega = 0.1783, 0.9702
        weights = second_difference_weights(nodes)
        direct = sum(w * scale * math.exp(-omega * s) for w, s in zip(weights, nodes))
        self.assertAlmostEqual(second_difference_image(scale, omega, nodes),
                               direct, places=18)

    def test_increasing_nodes_are_required(self):
        """Stops us believing a curvature built from an out-of-order lineage."""
        with self.assertRaises(ValueError):
            second_difference_weights([math.log(130), math.log(65), math.log(325)])


class ExponentFit(unittest.TestCase):
    def test_a_noiseless_one_exponent_series_is_recovered(self):
        """Stops us believing an exponent biased by the fitter itself.

        The wrong number is any ``omega`` that a noiseless series generated at
        0.9702 does not return, which would mean the reported 0.970 is a
        property of the golden-section search rather than of the data.
        """
        fitted = fit_exponent(synthetic(0.1783476, 0.9702))
        self.assertAlmostEqual(fitted["omega"], 0.9702, places=6)
        self.assertAlmostEqual(fitted["scale_amplitude"], 0.1783476, places=8)
        self.assertLess(fitted["statistic"], 1e-6)
        self.assertEqual(fitted["degrees_of_freedom"], 3)

    def test_the_search_window_does_not_bind_at_the_reported_exponent(self):
        """Stops us believing an exponent pinned by the search window's edge."""
        fitted = fit_exponent(synthetic(0.1783476, 0.9702))
        self.assertFalse(fitted["window_binds"])
        self.assertGreater(fitted["distance_to_the_nearest_window_edge"], 0.5)

    def test_the_window_binds_flag_fires_when_it_should(self):
        """Stops us believing a silent boundary hit reported as a measurement."""
        fitted = fit_exponent(synthetic(1.0, OMEGA_WINDOW[1] + 1.0))
        self.assertTrue(fitted["window_binds"])


class HeldOutPrediction(unittest.TestCase):
    def test_a_noiseless_series_is_predicted_exactly_out_of_sample(self):
        """Stops us believing a held-out check that secretly saw the held-out point.

        The wrong number is a nonzero prediction error on data that obeys the
        model exactly; if the fold leaked, the four-point exponents would still
        agree but the machinery would not be testing anything.
        """
        folds = leave_one_transition_out(synthetic(0.1783476, 0.9702))
        self.assertLess(folds["largest_relative_prediction_error"], 1e-6)
        self.assertLess(folds["omega_half_spread"], 1e-6)
        self.assertEqual(len(folds["folds"]), 5)
        for fold in folds["folds"]:
            self.assertEqual(len(fold["held_out"].split("->")), 2)


class ErrorBudget(unittest.TestCase):
    """The published exponent must not exclude unity on a statistics-only bar."""

    NUMBERS = {"fit": {"omega": 0.970179, "statistic": 277.4},
               "folds": {"omega_half_spread": 0.0154},
               "weighting": {"omega_shift": 0.0267}}

    def test_unity_survives_once_the_weighting_systematic_is_included(self):
        """Stops us believing ``omega != 1``.

        The wrong number is 0.970 quoted against the leave-one-out spread alone.
        The orientation weighting moves the exponent by 0.027 on its own -- #582
        warns that the equal-weight residue alternates in sign along a lineage --
        so the combined bar is 0.031 and unity is 0.96 bars away, not 2.
        """
        budget = error_budget(self.NUMBERS["fit"], self.NUMBERS["folds"],
                              self.NUMBERS["weighting"])
        self.assertFalse(budget["unit_exponent_is_excluded"])
        self.assertLess(budget["distance_from_unity_in_combined_uncertainties"], 2.0)

    def test_dropping_the_systematic_would_have_excluded_unity(self):
        """Names the exact mistake the test above prevents.

        With the weighting shift set to zero the combined bar is 0.0154 and
        unity sits 1.94 bars away -- still short of the two-bar rule, but the
        margin more than doubles, and that is the number a statistics-only
        budget would have published.
        """
        budget = error_budget(self.NUMBERS["fit"], self.NUMBERS["folds"],
                              {"omega_shift": 0.0})
        self.assertGreater(budget["distance_from_unity_in_combined_uncertainties"],
                           1.9)

    def test_the_fit_is_recorded_as_rejected(self):
        """Stops us believing a chi-square of 277 on 3 df read as a good fit."""
        budget = error_budget(self.NUMBERS["fit"], self.NUMBERS["folds"],
                              self.NUMBERS["weighting"])
        self.assertTrue(budget["the_fit_is_itself_rejected"])


class CandidateSizes(unittest.TestCase):
    def test_338_is_unusable_because_it_has_one_primitive_representative(self):
        """Stops us proposing a production the pipeline cannot run.

        ``338 = 2 * 13^2`` looked like the cheapest way to reach a m=2.5
        transition with a 5-adic valuation of one, but ``(17, 7)`` is its only
        primitive representative, so it has no second orientation and no spin-0
        combination at all.  The wrong number is a design cost quoted for a size
        that yields zero usable transitions.
        """
        self.assertEqual(primitive_representatives(338), [(17, 7)])
        entry = next(item for item in candidate_sizes() if item["size"] == 338)
        self.assertFalse(entry["usable_in_this_pipeline"])

    def test_725_breaks_the_valuation_interpolation_degeneracy(self):
        """Stops us believing the existing five transitions can ever separate labels.

        On the committed sizes, a 5-adic valuation of two and an *absent*
        interpolating spin-0 combination coincide exactly (325 and 425 have
        both).  ``725 = 5^2 * 29`` has valuation two and *does* admit an
        interpolating combination, because ``(26, 7)`` and ``(23, 14)`` straddle
        zero in ``cos 4 theta``.  The wrong number is a label screen p-value
        quoted as if the two labels were distinguishable.
        """
        entry = next(item for item in candidate_sizes() if item["size"] == 725)
        self.assertTrue(entry["usable_in_this_pipeline"])
        self.assertEqual(entry["five_adic_valuation"], 2)
        self.assertTrue(entry["an_interpolating_spin_zero_combination_exists"])
        self.assertTrue(entry["breaks_the_valuation_interpolation_degeneracy"])
        self.assertIn([26, 7], entry["primitive_representatives"])
        self.assertIn([23, 14], entry["primitive_representatives"])

    def test_the_committed_extrapolating_sizes_do_not_break_it(self):
        """Stops us believing 325/425/650 already supply the crossing."""
        for size in (325, 425, 650):
            self.assertEqual(five_adic_valuation(size), 2)
            entry = {"size": size}
            representatives = primitive_representatives(size)
            self.assertGreaterEqual(len(representatives), 2, entry)

    def test_725_continues_the_two_size_lineage_at_the_declared_multiplier(self):
        """Stops us believing a third curvature that is not on the same lineage.

        ``725 = 2.5 * 290`` and shares the parent prime 29 with 145 and 290, so
        it extends p50 rather than starting a fourth lineage; a fourth lineage
        would supply no second difference at all.
        """
        self.assertAlmostEqual(725 / 290, 2.5, places=12)
        self.assertEqual(725 % 29, 0)
        self.assertEqual(145 % 29, 0)

    def test_every_declared_candidate_is_a_sum_of_two_squares(self):
        """Stops us shipping a candidate list with an unrepresentable size."""
        for size in CANDIDATE_SIZES:
            self.assertTrue(primitive_representatives(size), size)


class Verdict(unittest.TestCase):
    GOOD_CONTROL = {"every_statistic_matches_exactly": True}

    def test_reproduction_failure_short_circuits_every_other_finding(self):
        """Stops us publishing an exponent fitted to the wrong reanalysis."""
        outcome = decide({"every_statistic_matches_exactly": False},
                         {"statistic": 1.0}, {"largest_relative_prediction_error": 0.0},
                         {"unit_exponent_is_excluded": False},
                         {"the_discrepancy_reproduces_across_lineages": True},
                         {"the_discrepancy_survives_the_weighting_change": True})
        self.assertEqual(outcome["verdict"], "REPRODUCTION_CONTROL_FAILED")

    def test_the_headline_verdict_needs_both_halves(self):
        """Stops us calling the law established on the tangent alone.

        The verdict names two things at once: the law predicts a held-out
        transition, and the same law is falsified by the curvature.  Dropping
        either half changes what the reader is entitled to conclude.
        """
        both = decide(self.GOOD_CONTROL, {"statistic": 277.0},
                      {"largest_relative_prediction_error": 0.0465},
                      {"unit_exponent_is_excluded": False},
                      {"the_discrepancy_reproduces_across_lineages": True},
                      {"the_discrepancy_survives_the_weighting_change": True})
        self.assertEqual(both["verdict"],
                         "ONE_EXPONENT_DESCRIBES_THE_TANGENT_AND_FAILS_THE_CURVATURE")
        without_curvature = decide(self.GOOD_CONTROL, {"statistic": 277.0},
                                   {"largest_relative_prediction_error": 0.0465},
                                   {"unit_exponent_is_excluded": False},
                                   {"the_discrepancy_reproduces_across_lineages": False},
                                   {"the_discrepancy_survives_the_weighting_change": True})
        self.assertEqual(without_curvature["verdict"], "AMPLITUDE_LAW_NOT_ESTABLISHED")
        without_prediction = decide(self.GOOD_CONTROL, {"statistic": 277.0},
                                    {"largest_relative_prediction_error": 0.9},
                                    {"unit_exponent_is_excluded": False},
                                    {"the_discrepancy_reproduces_across_lineages": True},
                                    {"the_discrepancy_survives_the_weighting_change": True})
        self.assertEqual(without_prediction["verdict"], "AMPLITUDE_LAW_NOT_ESTABLISHED")
        without_weighting = decide(
            self.GOOD_CONTROL, {"statistic": 277.0},
            {"largest_relative_prediction_error": 0.0465},
            {"unit_exponent_is_excluded": False},
            {"the_discrepancy_reproduces_across_lineages": True},
            {"the_discrepancy_survives_the_weighting_change": False})
        self.assertEqual(without_weighting["verdict"], "AMPLITUDE_LAW_NOT_ESTABLISHED")


class CommittedArtifact(unittest.TestCase):
    """Pins the published numbers, so a silent drift in the pipeline is visible."""

    @classmethod
    def setUpClass(cls):
        if not DEFAULT_OUTPUT.exists():
            raise unittest.SkipTest("run scripts/p582_amplitude_law.py first")
        cls.report = json.loads(Path(DEFAULT_OUTPUT).read_text())

    def test_the_affine_statistics_still_reproduce_582_exactly(self):
        """Stops us believing a law fitted to a drifted reanalysis of #582."""
        self.assertTrue(
            self.report["reproduction_control"]["every_statistic_matches_exactly"])

    def test_the_second_difference_discrepancy_reproduces_in_both_lineages(self):
        """Stops us believing the curvature failure is noise in one lineage.

        The wrong number is a single lineage's ratio quoted alone.  Two
        independent lineages -- different parent primes, different productions,
        different seeds -- give 1.54 and 1.56, agreeing with each other far more
        closely than either agrees with one.
        """
        check = self.report["second_difference_check"]
        self.assertEqual(check["lineages_checked"], 2)
        low, high = check["ratio_range"]
        self.assertGreater(low, 1.4)
        self.assertLess(high, 1.7)
        self.assertLess(check["relative_spread_between_lineages"], 0.05)
        self.assertTrue(check["the_discrepancy_reproduces_across_lineages"])
        for item in check["per_lineage"]:
            if item["usable"]:
                self.assertTrue(item["weights_annihilate_constants"])

    def test_the_two_lineages_share_their_third_rung_structure(self):
        """Stops us believing the 1.3% agreement is independent evidence.

        Under the primary weighting, 325 and 425 are the only sizes in the tree
        whose spin-0 combination extrapolates (weights 1.278/-0.278 and
        -0.026/1.026), and they are also the only sizes run at 5M per batch
        rather than 1M.  They are rung 3 of gaussian_13 and rung 3 of
        gaussian_17 respectively, so the two "independent" lineages share that
        structure exactly.  The wrong number is 1.3% quoted as cross-lineage
        replication when the shared part has not been varied.
        """
        sizes = {item["size"]: item for item in self.report["candidate_sizes"]}
        del sizes  # the committed sizes are pinned by the flow artifact, below
        flow = json.loads(
            (Path(DEFAULT_OUTPUT).parent.parent
             / "wasserstein-shape-flow" / "latest.json").read_text())["sizes"]
        extrapolating = {int(size) for size, block in flow.items()
                         if not block["spin4_correction_is_an_interpolation"]}
        self.assertEqual(extrapolating, {325, 425})
        deep = {int(size) for size, block in flow.items()
                if "500m" in block["source"] or "_500m" in block["source"]}
        self.assertEqual(deep, {325, 425})

    def test_the_discrepancy_survives_the_weighting_that_never_extrapolates(self):
        """Stops us believing a curvature failure manufactured by extrapolation.

        The equal weighting is 0.5/0.5 at every size, so it never extrapolates
        and carries none of the shared rung-3 structure above.  The ratio there
        is 1.58 and 1.44 -- still nowhere near one.  The wrong number this
        stops us believing is 1.55 read as a systematic of the spin-0
        combination, and equally 1.3% read as the honest spread: over all four
        weighting-by-lineage cells it is 8.8%.
        """
        weighting = self.report["weighting_systematic"]
        cells = weighting["second_difference_ratio_cells"]
        self.assertEqual(len(cells), 4)
        low, high = weighting["second_difference_ratio_range"]
        self.assertGreater(low, 1.4)
        self.assertLess(high, 1.7)
        self.assertTrue(weighting["the_discrepancy_survives_the_weighting_change"])
        self.assertGreater(
            weighting["second_difference_ratio_spread_over_all_cells"], 0.05)

    def test_the_held_out_amplitude_prediction_stays_under_five_percent(self):
        """Stops us believing the exponent is only an in-sample description."""
        folds = self.report["leave_one_transition_out"]
        self.assertLess(folds["largest_relative_prediction_error"], 0.05)

    def test_the_acquisition_is_a_single_size(self):
        """Stops us shipping a design that asks for more than the one crossing."""
        self.assertEqual(self.report["acquisition"]["preferred"], [725])


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
