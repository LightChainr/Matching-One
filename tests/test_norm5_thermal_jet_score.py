#!/usr/bin/env python3


from __future__ import annotations
import math
from pathlib import Path
import sys
import unittest

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_norm5_thermal_jet import (  # noqa: E402
    State,
    generalized_chi_square,
    jackknife_covariance,
    load_prediction,
    multiplier_residuals,
    width_residuals,
)


def make_state(n: int, width: float, normalized_jet: list[float]) -> State:
    jet = tuple(
        mp.mpf(value) / mp.power(width, order)
        for order, value in enumerate(normalized_jet)
    )
    return State(n, mp.mpf("0.59"), jet, mp.mpf(width), mp.mpf("5"))


class Norm5ThermalJetScoreTests(unittest.TestCase):
    def test_frozen_prediction_has_required_scoring_order(self) -> None:
        payload = load_prediction(ROOT / "predictions/hermite_krawtchouk_jet_20260829.yaml")
        self.assertEqual(payload["scoring_order"][0], "rank_gap_width_collapse")

    def test_width_residual_cancels_shape_with_free_amplitude_and_width(self) -> None:
        shape = [2.0, -1.0, 3.0, -4.0, 5.0, -6.0, 7.0]
        states = {
            10: make_state(10, 0.8, shape),
            20: make_state(20, 0.9, shape),
            50: make_state(50, 1.3, [3 * value for value in shape]),
        }
        residual = width_residuals(states, [(10, 20, 50)])
        self.assertTrue(all(abs(value) < mp.mpf("1e-40") for value in residual))

    def test_multiplier_residual_recovers_declared_cocycle(self) -> None:
        multiplier = mp.mpf(8) / 5
        parent = [1.0, 2.0, -1.0, 4.0, -2.0, 3.0, 5.0]
        norm2 = [2.0, -1.0, 3.0, 1.0, 4.0, -2.0, 6.0]
        norm5 = [
            float(multiplier * b - (multiplier - 1) * a)
            for a, b in zip(parent, norm2)
        ]
        states = {
            10: make_state(10, 0.7, parent),
            20: make_state(20, 0.9, norm2),
            50: make_state(50, 1.2, norm5),
        }
        residual = multiplier_residuals(states, [(10, 20, 50)], multiplier)
        self.assertTrue(all(abs(value) < mp.mpf("1e-14") for value in residual))

    def test_jackknife_covariance_uses_delete_one_normalization(self) -> None:
        covariance = jackknife_covariance([[0.0, 0.0], [1.0, 2.0], [2.0, 4.0]])
        self.assertAlmostEqual(float(covariance[0][0]), 4.0 / 3.0)
        self.assertAlmostEqual(float(covariance[0][1]), 8.0 / 3.0)
        self.assertAlmostEqual(float(covariance[1][1]), 16.0 / 3.0)

    def test_generalized_score_whitens_diagonal_covariance(self) -> None:
        score = generalized_chi_square(
            [mp.mpf(1), mp.mpf(2)],
            [[mp.mpf(1), mp.mpf(0)], [mp.mpf(0), mp.mpf(4)]],
        )
        self.assertAlmostEqual(float(score["chi_square"]), 2.0)
        self.assertEqual(score["degrees_of_freedom"], 2)
        self.assertAlmostEqual(float(score["chi_square_survival"]), math.exp(-1.0))


if __name__ == "__main__":
    unittest.main()
