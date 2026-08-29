import math
import sys
import unittest
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import yaml


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from hermite_krawtchouk_scaling_jet import (  # noqa: E402
    cocycle_residual,
    dilate_jet,
    pooled_gap_convention_shift,
    response_from_modes,
    translate_jet,
    width_cross_residual,
)
from threshold_score_modes import binomial_weights, krawtchouk_mode  # noqa: E402


class HermiteKrawtchoukScalingJetTests(unittest.TestCase):
    def test_exact_finite_n_generating_function(self) -> None:
        mp.mp.dps = 60
        n = 12
        p0 = mp.mpf("0.592746050790")
        p = mp.mpf("0.631")
        response = [mp.mpf((x * x + 3 * x + 1) % 11) / 10 for x in range(n + 1)]
        weights0 = binomial_weights(n, p0)
        coefficients = [
            mp.fsum(
                weights0[x] * response[x] * krawtchouk_mode(n, x, order, p0)
                for x in range(n + 1)
            )
            for order in range(n + 1)
        ]
        direct = mp.fsum(
            weight * value for weight, value in zip(binomial_weights(n, p), response)
        )
        reconstructed = response_from_modes(coefficients, n, p0, p)
        self.assertLess(abs(direct - reconstructed), mp.mpf("1e-48"))

    def test_translation_and_width_operators_on_polynomial_jet(self) -> None:
        # F(z)=1+2z+3z^2+4z^3 has derivative jet [1,2,6,24].
        jet = [mp.mpf(1), mp.mpf(2), mp.mpf(6), mp.mpf(24)]
        delta = mp.mpf("0.125")
        translated = translate_jet(jet, delta)
        expected0 = 1 + 2 * delta + 3 * delta**2 + 4 * delta**3
        expected1 = 2 + 6 * delta + 12 * delta**2
        self.assertLess(abs(translated[0] - expected0), mp.mpf("1e-50"))
        self.assertLess(abs(translated[1] - expected1), mp.mpf("1e-50"))
        self.assertEqual(dilate_jet(jet, mp.mpf(3)), [1, 6, 54, 648])

    def test_width_cross_residual_cancels_width_and_amplitude(self) -> None:
        shape = [mp.mpf(value) for value in (2, 3, -5, 7, -11)]
        parent_width, child_width = mp.mpf("0.41"), mp.mpf("0.43")
        parent_amplitude, child_amplitude = mp.mpf("1.7"), mp.mpf("2.2")
        parent = [
            parent_amplitude * value / parent_width**order
            for order, value in enumerate(shape)
        ]
        child = [
            child_amplitude * value / child_width**order
            for order, value in enumerate(shape)
        ]
        residual = width_cross_residual(
            parent, child, parent_width, child_width
        )
        self.assertTrue(all(abs(value) < mp.mpf("1e-48") for value in residual))

    def test_q2_and_jordan_cocycles(self) -> None:
        base = [mp.mpf(2), mp.mpf(-3), mp.mpf(5)]
        correction = [mp.mpf(7), mp.mpf(11), mp.mpf(-13)]
        n = mp.mpf(65)
        q2 = [base[i] + correction[i] / n for i in range(3)]
        q2_2 = [base[i] + correction[i] / (2 * n) for i in range(3)]
        q2_5 = [base[i] + correction[i] / (5 * n) for i in range(3)]
        self.assertTrue(
            all(abs(value) < mp.mpf("1e-48") for value in cocycle_residual(q2, q2_2, q2_5, mp.mpf(8) / 5))
        )
        jordan = [base[i] + correction[i] * mp.log(n) for i in range(3)]
        jordan_2 = [base[i] + correction[i] * mp.log(2 * n) for i in range(3)]
        jordan_5 = [base[i] + correction[i] * mp.log(5 * n) for i in range(3)]
        self.assertTrue(
            all(
                abs(value) < mp.mpf("1e-48")
                for value in cocycle_residual(
                    jordan, jordan_2, jordan_5, mp.log(5) / mp.log(2)
                )
            )
        )

    def test_quarter_is_not_forced_by_two_orientation_rank_conventions(self) -> None:
        self.assertEqual(pooled_gap_convention_shift([(1, 1), (1, 1)]), 0)
        self.assertEqual(pooled_gap_convention_shift([(0, 1), (0, 0)]), Fraction(1, 2))
        shifts = {
            pooled_gap_convention_shift([(dm1, dp1), (dm2, dp2)])
            for dm1 in (-1, 0, 1)
            for dp1 in (-1, 0, 1)
            for dm2 in (-1, 0, 1)
            for dp2 in (-1, 0, 1)
        }
        self.assertNotIn(Fraction(-1, 4), shifts)

    def test_prediction_artifact_freezes_width_first(self) -> None:
        path = ROOT / "predictions" / "hermite_krawtchouk_jet_20260829.yaml"
        artifact = yaml.safe_load(path.read_text(encoding="utf-8"))
        self.assertEqual(artifact["status"], "source_theory_frozen_before_norm5_reveal")
        self.assertEqual(artifact["scoring_order"][0], "rank_gap_width_collapse")
        self.assertEqual(artifact["rank_gap_bridge"]["correction_exponent_in_N"], "5/8")


if __name__ == "__main__":
    unittest.main()
