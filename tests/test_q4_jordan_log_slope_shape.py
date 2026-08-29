import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from derive_q4_jordan_log_slope_shape import render  # noqa: E402


class Q4JordanLogSlopeShapeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = render(90)
        cls.artifact = json.loads(
            (ROOT / "predictions" / "q4_jordan_log_slope_shape_20260829.json")
            .read_text(encoding="utf-8")
        )

    def test_exact_module_slope_coefficients(self) -> None:
        exact = self.payload["exact_inputs"]
        self.assertEqual(Fraction(exact["Q4_Ward_coefficient"]), Fraction(493, 96))
        self.assertEqual(
            Fraction(exact["logN_coefficient_after_Ward"]),
            Fraction(-493, 192),
        )
        self.assertEqual(
            Fraction(exact["area_normalized_E4hat_2i_over_i"]),
            Fraction(11, 4),
        )

    def test_high_precision_cm_targets(self) -> None:
        checks = self.payload["numerical_qseries"]["checks"]
        self.assertLess(mp.mpf(checks["hexagonal_zero_abs"]), mp.mpf("1e-80"))
        self.assertLess(mp.mpf(checks["ratio_error_abs"]), mp.mpf("1e-68"))

    def test_interpretation_boundary_is_explicit(self) -> None:
        boundaries = self.payload["interpretation_boundary"]
        self.assertIn(
            "the_norm4_scale_cocycle_alone_cannot_identify_the_Q4_Jordan_module",
            boundaries,
        )
        self.assertIn(
            "the_full_top_intercept_shape_A_tilde_tau_has_not_been_derived",
            boundaries,
        )

    def test_frozen_artifact_matches_derivation(self) -> None:
        self.assertEqual(self.artifact, self.payload)


if __name__ == "__main__":
    unittest.main()
