#!/usr/bin/env python3
"""Lock the N=580 ladder runner against its own frozen design.

No engine runs here. What these tests exist to stop is a five-hour scoring run
that executes a geometry the frozen design does not name, or scores it against
competitor numbers nobody generated -- both of which produce a plausible
artifact rather than an error.
"""

from __future__ import annotations

from fractions import Fraction
import json
import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_aspect_ladder_n580 as ladder  # noqa: E402
from score_aspect_ladder_n580 import (  # noqa: E402
    COMPETING,
    SITE_COUNT,
    jackknife_se,
    load_design,
)


class ScoreAspectLadderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.rungs = load_design()

    def test_the_runner_reads_the_geometry_it_claims_to_execute(self) -> None:
        """The runner carries no copy of the period matrices.

        If it did, an edit to the design would leave the runner executing the
        old tori while the artifact cited the new design. It reads them, and
        `load_design` refuses anything whose determinant is not the site count.
        """
        self.assertEqual(set(self.rungs), {1, 2, 4})
        for aspect, rung in self.rungs.items():
            with self.subTest(aspect=aspect):
                self.assertEqual(rung["gaussian_norm"] * aspect, SITE_COUNT)
                for key in ("first_matrix_row_major", "second_matrix_row_major"):
                    a, b, c, d = rung[key]
                    self.assertEqual(a * d - b * c, SITE_COUNT)

    def test_all_three_rungs_share_the_leverage_the_estimator_divides_by(self) -> None:
        leverages = {rung["delta_cos4"] for rung in self.rungs.values()}
        self.assertEqual(leverages, {"8064/4205"})

    def test_the_discriminating_ratio_is_the_one_where_spin8_cancels(self) -> None:
        """r4_over_r1, not r2_over_r1, decides exclusion.

        The r=2 rung carries the opposite leakage sign, so its ratio keeps a
        systematic that does not shrink with samples. Excluding a competitor on
        that entry would be excluding it on a bias.
        """
        self.assertEqual(
            Fraction(self.rungs[1]["spin8_leakage"]),
            Fraction(self.rungs[4]["spin8_leakage"]),
        )
        self.assertEqual(
            Fraction(self.rungs[2]["spin8_leakage"]),
            -Fraction(self.rungs[1]["spin8_leakage"]),
        )

    def test_the_competitor_table_matches_the_frozen_prediction_file(self) -> None:
        """The run is scored against these numbers, so they may not be retyped.

        A competitor value that drifts from the frozen file turns the run into a
        comparison against something nobody predicted.
        """
        import yaml

        frozen = yaml.safe_load(
            (ROOT / "predictions" / "aspect_ladder_n580_20260905.yaml").read_text(encoding="utf-8")
        )
        predictions = frozen["competing_predictions"]
        self.assertEqual(set(predictions), set(COMPETING))
        for name, (at_two, at_four) in COMPETING.items():
            with self.subTest(prediction=name):
                self.assertAlmostEqual(float(predictions[name][0]), at_two, places=10)
                self.assertAlmostEqual(float(predictions[name][1]), at_four, places=10)

    def test_the_two_live_hypotheses_are_far_apart_at_r4(self) -> None:
        """If they were not, the run would be a null by construction."""
        weight4 = COMPETING["q4_jordan_weight4"]
        linear = COMPETING["bare_aspect_ratio"]
        self.assertLess(abs(weight4[0] - linear[0]), 1.0)
        self.assertGreater(abs(weight4[1] - linear[1]), 6.0)

    def test_jackknife_matches_the_closed_form_on_a_known_case(self) -> None:
        """Delete-one pseudovalues of a mean reduce to the ordinary SE."""
        values = [1.0, 2.0, 3.0, 4.0]
        full = sum(values) / len(values)
        deleted = [
            (sum(values) - value) / (len(values) - 1) for value in values
        ]
        mean = sum(values) / len(values)
        variance = sum((v - mean) ** 2 for v in values) / (len(values) * (len(values) - 1))
        self.assertAlmostEqual(jackknife_se(full, deleted), variance ** 0.5, places=12)
        with self.assertRaises(ValueError):
            jackknife_se(1.0, [1.0])



class FiellerContrastTests(unittest.TestCase):
    """The statistic that decides the ladder, tested on the numbers it decided.

    The run measured A4 at r=1 as 9.016e-04 +- 2.491e-04 -- 3.6 sigma from zero --
    against 4.132e-03 +- 1.951e-04 at r=4.  Dividing first turns that into
    4.58 +- 1.32 and then reads competitors off it, which is where the error was.
    """

    # The committed N=580 run, results/aspect-ladder-n580/latest.json.
    Y, SY = 4.131807633330795e-03, 1.9512853917124002e-04
    X, SX = 9.016433036753112e-04, 2.4913893253110386e-04
    RHO = -0.1526
    COV = RHO * SX * SY

    def test_the_two_statistics_agree_when_the_denominator_is_sharp(self) -> None:
        """The anchor: with a well-measured denominator there is nothing to fix.

        The wrong number this stops us believing is a Fieller implementation
        that is simply a different formula rather than a better one.  Shrink the
        denominator's error and it must converge on the naive ratio z.
        """
        var_x = (self.X * 1e-6) ** 2
        ratio = self.Y / self.X
        se_ratio = math.sqrt(self.SY ** 2) / self.X  # denominator noise negligible
        for predicted in (2.0, 4.0, 10.9908008589, 16.0):
            naive = (ratio - predicted) / se_ratio
            exact = ladder.fieller_z(
                self.Y, self.SY ** 2, self.X, var_x, 0.0, predicted
            )
            self.assertAlmostEqual(naive, exact, delta=abs(naive) * 1e-3 + 1e-9)

    def test_the_weak_denominator_is_where_the_two_disagree(self) -> None:
        """Pins the finding that changed the run's verdicts.

        The wrong numbers here are the pre-registered ones.  On this data the
        ratio z-test excluded the weight-4 modular shape at 4.9 sigma and plain
        area scaling at 8.7, and let "no modulus dependence" survive at 2.7.
        All three are artifacts of dividing by a 3.6-sigma denominator.
        """
        cases = {          # predicted r4/r1 -> (frozen ratio z, correct z)
            1.0:            (2.72, 9.53),
            4.0:            (0.44, 0.50),
            10.9908008589:  (-4.87, -2.08),
            16.0:           (-8.67, -2.56),
            120.79770352:   (-88.26, -3.48),
        }
        ratio = self.Y / self.X
        se_ratio = 1.3167396120465709
        for predicted, (frozen, correct) in cases.items():
            self.assertAlmostEqual((ratio - predicted) / se_ratio, frozen, delta=0.02)
            self.assertAlmostEqual(
                ladder.fieller_z(self.Y, self.SY ** 2, self.X, self.SX ** 2,
                                 self.COV, predicted),
                correct, delta=0.02,
            )

    def test_the_verdict_flips_in_both_directions(self) -> None:
        """Guards against a statistic chosen because it helped.

        A correction that only ever loosened, or only ever tightened, would be
        suspect.  This one excludes a competitor the frozen test let live and
        revives two it had killed, which is what an unmotivated fix looks like.
        """
        def decide(predicted: float) -> bool:
            return abs(ladder.fieller_z(self.Y, self.SY ** 2, self.X, self.SX ** 2,
                                        self.COV, predicted)) >= 3.0
        self.assertTrue(decide(1.0))            # survived the frozen test
        self.assertFalse(decide(10.9908008589))  # the frozen test excluded it
        self.assertFalse(decide(16.0))           # so did this one

    def test_the_conclusions_do_not_rest_on_the_correlation(self) -> None:
        """The covariance was reconstructed, not measured, for this run.

        The committed artifact dropped the delete-one replicates, so rho was
        backed out of the reported ratio standard error.  If the verdicts moved
        with rho, that reconstruction would be load-bearing and unusable.
        """
        for rho in (-0.5, -0.1526, 0.0, 0.5):
            cov = rho * self.SX * self.SY
            for predicted, expected in ((1.0, True), (10.9908008589, False), (16.0, False)):
                z = ladder.fieller_z(self.Y, self.SY ** 2, self.X, self.SX ** 2,
                                     cov, predicted)
                self.assertEqual(abs(z) >= 3.0, expected, f"rho={rho} R0={predicted}")

    def test_the_interval_is_unbounded_when_the_denominator_is_not_resolved(self) -> None:
        """Stops a wide interval being reported where there is no interval.

        Below 3 sigma on the denominator every large ratio is compatible with
        the data.  Reporting some finite upper limit there would be the same
        error as the ratio z-test, one level up.
        """
        bounds = ladder.fieller_interval(self.Y, self.SY ** 2, self.X, self.SX ** 2, self.COV)
        self.assertIsNotNone(bounds)
        self.assertAlmostEqual(bounds[0], 2.399, delta=0.01)
        self.assertAlmostEqual(bounds[1], 27.42, delta=0.05)
        # Same data, denominator error widened until it no longer clears 3 sigma.
        weak = ladder.fieller_interval(self.Y, self.SY ** 2, self.X, (self.X / 2.9) ** 2, 0.0)
        self.assertIsNone(weak)

    def test_paired_covariance_sees_a_correlation_a_separate_one_would_miss(self) -> None:
        """The batches must be deleted together or this is not a covariance.

        Two channels built from the same random blocks are correlated; deleting
        batch b from one and batch b from the other is what measures it.  A
        perfectly proportional pair must come back with the covariance implied
        by its own variances, and an antiproportional one with the sign flipped.
        """
        deleted_x = [1.0 + 0.01 * b for b in range(20)]
        full_x = sum(deleted_x) / 20
        for scale in (3.0, -3.0):
            deleted_y = [scale * v for v in deleted_x]
            full_y = scale * full_x
            cov = ladder.jackknife_covariance(full_x, deleted_x, full_y, deleted_y)
            sx = ladder.jackknife_se(full_x, deleted_x)
            sy = ladder.jackknife_se(full_y, deleted_y)
            self.assertAlmostEqual(cov / (sx * sy), 1.0 if scale > 0 else -1.0, places=9)

if __name__ == "__main__":
    unittest.main()
