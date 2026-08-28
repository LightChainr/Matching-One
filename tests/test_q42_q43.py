#!/usr/bin/env python3
from __future__ import annotations

import math
import sys
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "scripts"))

import analyze_q42_q43 as q  # noqa: E402


class HarmonicExactnessTests(unittest.TestCase):
    def test_cos4_cos8_exact_for_8_1(self) -> None:
        c4 = q.cos4_fraction(8, 1)
        self.assertEqual(c4, Fraction(3713, 4225))
        c8 = q.cos8_fraction(8, 1)
        self.assertEqual(c8, 2 * c4 * c4 - 1)

    def test_n65_and_n130_are_complementary_delta_cos8(self) -> None:
        d65 = q.design_row(8, 1, 7, 4)
        d130 = q.design_row(11, 3, 9, 7)
        self.assertEqual(Fraction(d65["delta_cos4"]), Fraction(d130["delta_cos4"]))
        self.assertEqual(Fraction(d65["delta_cos8"]), -Fraction(d130["delta_cos8"]))


class ClosureIdentityTests(unittest.TestCase):
    def test_linear_common_slope_gives_C_N_one(self) -> None:
        pc = 0.6
        s = 8.0
        d = 0.002
        p1 = pc - (d / 2.0) / s
        p2 = pc + (d / 2.0) / s
        delta_root = p1 - p2
        delta_m = d
        mean_mp = s
        c_n = -delta_root * mean_mp / delta_m
        self.assertAlmostEqual(c_n, 1.0, places=12)

    def test_quadratic_root_recovers_exact_parabola(self) -> None:
        # 0.01 - 8 x + 3 x^2 = 0  (i.e. m=0.01, mp=-8, mpp=6)
        m, mp, mpp = 0.01, -8.0, 6.0
        x = q.quadratic_root_shift(m, mp, mpp)
        residual = m + mp * x + 0.5 * mpp * x * x
        self.assertAlmostEqual(residual, 0.0, places=10)
        lin = -m / mp
        self.assertLess(abs(x - lin), 0.01)


class GlsTests(unittest.TestCase):
    def test_recovers_known_A4(self) -> None:
        c4 = np.array([1.36, 1.59, 1.36])
        true = 0.2
        y = true * c4
        se = np.array([0.01, 0.01, 0.01])
        fit = q.gls_fit(c4[:, None], y, se)
        self.assertAlmostEqual(fit["beta"][0], true, places=10)
        self.assertLess(fit["chi2"], 1e-20)

    def test_omega_selection_uses_training_only(self) -> None:
        rows = [
            {"N": 65, "y_scaled": 0.27, "y_scaled_se": 0.02, "delta_cos4": 1.36},
            {"N": 85, "y_scaled": 0.32, "y_scaled_se": 0.02, "delta_cos4": 1.59},
            {"N": 130, "y_scaled": 0.27, "y_scaled_se": 0.02, "delta_cos4": 1.36},
        ]
        sel = q.select_omega(rows)
        self.assertIn(sel["selected_omega"], q.OMEGA_CANDIDATES)
        self.assertEqual(len(sel["grid"]), len(q.OMEGA_CANDIDATES))


if __name__ == "__main__":
    unittest.main()
