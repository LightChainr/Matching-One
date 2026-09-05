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
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

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


if __name__ == "__main__":
    unittest.main()
