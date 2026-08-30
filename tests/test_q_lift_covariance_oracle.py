from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from q_lift_covariance_oracle import (  # noqa: E402
    add,
    analyze_length,
    multiply_q,
    normalized_derivative_at_one,
    path_derivative_at_one,
    scale,
)


class QLiftCovarianceOracleTests(unittest.TestCase):
    def test_abstract_polynomial_identity(self) -> None:
        w0 = {(0, 0): 3, (2, 1): 5}
        w2 = {(1, 0): 7, (3, 2): 11}
        h = add(w2, scale(w0, -1))
        c = add(w2, scale(multiply_q(w0), -1))
        self.assertEqual(add(c, scale(h, -1)), add(w0, scale(multiply_q(w0), -1)))

    def test_paths_have_same_counterterm_but_different_individual_tangents(self) -> None:
        numerator = {(1, 0): 2, (0, 1): 3}
        denominator = {(0, 0): 5, (1, 1): 7}
        self.assertNotEqual(
            path_derivative_at_one(numerator, "fixed_v_1"),
            path_derivative_at_one(numerator, "critical_square_bond_v_sqrt_Q"),
        )
        self.assertIsInstance(
            normalized_derivative_at_one(numerator, denominator, "fixed_v_1"),
            Fraction,
        )

    def test_L2_exact_oracle(self) -> None:
        result = analyze_length(2)
        self.assertEqual(result["configurations"], 256)
        self.assertTrue(result["endpoint_H_equals_C"])
        self.assertTrue(result["exact_polynomial_difference_passed"])
        self.assertTrue(result["passed"])
        for row in result["paths"].values():
            self.assertEqual(
                row["unnormalized_dC_minus_dH"],
                row["expected_unnormalized_counterterm"],
            )
            self.assertEqual(
                row["normalized_dc_minus_dh"],
                row["expected_normalized_counterterm"],
            )


if __name__ == "__main__":
    unittest.main()
