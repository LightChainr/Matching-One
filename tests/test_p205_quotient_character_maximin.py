#!/usr/bin/env python3

from __future__ import annotations

import sys
import unittest
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from design_p205_quotient_character_maximin import (  # noqa: E402
    QuotientPair,
    calibrate,
    campaign_score,
    character_difference,
    determinant,
    enumerate_pairs,
    gaussian_matrix,
    rank_campaigns,
)


SCORE = (
    ROOT
    / "results/server-20260829/P205-norm5-conjugate-coalescence/analysis/score.json"
)
PRISM = (
    QuotientPair(25, (5, 0), (4, 3)),
    QuotientPair(50, (7, 1), (5, 5)),
    QuotientPair(125, (11, 2), (10, 5)),
)


class P205QuotientCharacterMaximinTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.calibration = calibrate(SCORE, ROOT)

    def test_prism_has_declared_determinants_and_cross_smith_types(self) -> None:
        self.assertEqual(
            [(pair.first_smith, pair.second_smith) for pair in PRISM],
            [((5, 5), (1, 25)), ((1, 50), (5, 10)), ((1, 125), (5, 25))],
        )
        for pair in PRISM:
            self.assertEqual(determinant(gaussian_matrix(pair.first)), pair.n)
            self.assertEqual(determinant(gaussian_matrix(pair.second)), pair.n)
            self.assertNotEqual(pair.first_smith, pair.second_smith)

    def test_exact_character_prism_is_hadamard_like(self) -> None:
        signs = {
            harmonic: [
                "+" if character_difference(pair, harmonic) > 0 else "-"
                for pair in PRISM
            ]
            for harmonic in (4, 8, 12)
        }
        self.assertEqual(
            signs,
            {4: ["+", "+", "+"], 8: ["+", "-", "+"], 12: ["+", "+", "-"]},
        )
        self.assertEqual(character_difference(PRISM[0], 4), Fraction(1152, 625))
        self.assertEqual(character_difference(PRISM[1], 8), Fraction(-225792, 390625))
        self.assertEqual(
            character_difference(PRISM[2], 12),
            Fraction(-4983630043392, 3814697265625),
        )

    def test_calibration_is_the_completed_p205_noise_scale(self) -> None:
        self.assertAlmostEqual(
            self.calibration["mean_pair_contrast_variance"],
            5.314603523220254e-08,
            places=20,
        )
        self.assertAlmostEqual(
            self.calibration["cpu_seconds_per_site_update"],
            2.0709907904351133e-07,
            places=20,
        )

    def test_prism_beats_any_two_coordinate_campaign(self) -> None:
        prism = campaign_score(PRISM, self.calibration)
        self.assertAlmostEqual(prism["maximin_noncentrality"], 5.748455127225554)
        candidates = enumerate_pairs(25, 200)
        top_two = rank_campaigns(candidates, 2, self.calibration, top_n=1)[0]
        top_three = rank_campaigns(candidates, 3, self.calibration, top_n=1)[0]
        self.assertEqual(
            [(row["N"], row["first"], row["second"]) for row in top_three["pairs"]],
            [(25, [5, 0], [4, 3]), (50, [7, 1], [5, 5]), (125, [11, 2], [10, 5])],
        )
        self.assertLess(top_two["maximin_noncentrality"], 1.0)
        self.assertGreater(top_three["maximin_noncentrality"], 5.7)


if __name__ == "__main__":
    unittest.main()
