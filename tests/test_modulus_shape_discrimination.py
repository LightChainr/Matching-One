#!/usr/bin/env python3
"""Lock the two facts the rectangular-torus design rests on.

If either is wrong the design is void rather than merely weaker, so both are
worth one assertion each: the weight-4 amplitude must flip sign when the
lattice turns 45 degrees (otherwise the difference channel is identically zero
and there is nothing to measure), and weight 8 must not (otherwise the channel
does not remove the scalar corrections it is there to remove).
"""

from __future__ import annotations

from pathlib import Path
import sys
import unittest

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from modulus_shape_discrimination import (  # noqa: E402
    axis_and_diagonal,
    eisenstein,
    production_design,
    render,
)


class SpinSelectionTests(unittest.TestCase):
    def test_weight_four_flips_sign_under_a_45_degree_turn(self) -> None:
        with mp.workdps(40):
            axis, diagonal = axis_and_diagonal(1, 4, lambda t: eisenstein(t, 4))
            self.assertLess(abs(mp.re(diagonal / axis) + 1), mp.mpf(10) ** -30)

    def test_weight_eight_does_not_flip_and_so_cancels(self) -> None:
        with mp.workdps(40):
            axis, diagonal = axis_and_diagonal(1, 8, lambda t: eisenstein(t, 8))
            self.assertLess(abs(mp.re(diagonal / axis) - 1), mp.mpf(10) ** -30)


class PredictionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = render(dps=40)

    def test_the_predicted_ratio_is_eleven_quarters(self) -> None:
        q4 = next(
            row for row in self.result["candidate_shapes"] if row["shape"] == "E4"
        )
        self.assertEqual(self.result["q4_jordan_prediction_at_2i"], "11/4")
        self.assertLess(abs(float(q4["difference_ratio"]["2"]) - 2.75), 1e-12)

    def test_the_area_competitor_is_not_close(self) -> None:
        """11/4 against 4.  If these were within a few percent the run is pointless."""
        row = self.result["separation_from_area_scaling"][0]
        self.assertEqual(row["area_scaling"], "4")
        self.assertGreater(float(row["relative_distance"]), 0.30)

    def test_longer_rectangles_add_nothing_against_that_competitor(self) -> None:
        """The claim that one rectangle is the whole experiment."""
        for row in self.result["marginal_value_of_extra_aspect_ratios"]:
            with self.subTest(aspect_ratio=row["aspect_ratio"]):
                self.assertLess(float(row["distance_from_the_r2_value"]), 0.01)


class ProductionDesignTests(unittest.TestCase):
    def test_spin8_leaks_equal_and_opposite(self) -> None:
        """The size of the one systematic the score cannot outspend."""
        leakage = production_design()["spin8_leakage"]
        self.assertTrue(leakage["equal_and_opposite"])
        self.assertEqual(leakage["square"], "1148/21025")
        self.assertEqual(leakage["rectangular"], "-1148/21025")

    def test_both_families_have_the_same_site_count(self) -> None:
        design = production_design()
        for family in ("square_family", "rectangular_family"):
            for member in design[family]["members"]:
                with self.subTest(family=family, w=member["gaussian_integer"]):
                    self.assertEqual(member["sites"], design["site_count"])

    def test_the_two_families_have_equal_angular_leverage(self) -> None:
        """If they differed the ratio estimator would pay for it in variance."""
        design = production_design()
        self.assertTrue(design["leverages_are_equal"])
        self.assertEqual(
            design["square_family"]["angular_leverage"],
            design["rectangular_family"]["angular_leverage"],
        )

    def test_period_vectors_span_the_stated_lattice(self) -> None:
        """The determinant is the site count; a wrong vector would be a wrong torus."""
        design = production_design()
        for family in ("square_family", "rectangular_family"):
            for member in design[family]["members"]:
                (a, b), (c, d) = member["period_vectors"]
                with self.subTest(family=family, w=member["gaussian_integer"]):
                    self.assertEqual(abs(a * d - b * c), member["sites"])


if __name__ == "__main__":
    unittest.main()
