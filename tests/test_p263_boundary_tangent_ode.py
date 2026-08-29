from __future__ import annotations

import importlib.util
import json
from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "p263_boundary_tangent_ode.py"
SPEC = importlib.util.spec_from_file_location("p263_boundary_tangent", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class BoundaryTangentOdeTests(unittest.TestCase):
    def test_q_kappa_and_high_branch_exponent_derivatives(self) -> None:
        exponent, _ = MODULE.frobenius_high_branch(0)
        self.assertEqual(exponent.value, Fraction(2))
        self.assertEqual(exponent.derivative, Fraction(-2, 3))
        payload = MODULE.build_payload(0)
        self.assertEqual(
            payload["fk_cle_map"]["d_kappa_d_Q_at_Q1"],
            "-3 sqrt(3)/(2 pi)",
        )
        self.assertEqual(
            payload["high_frobenius_branch"]["d_Q_r"], "sqrt(3)/pi"
        )

    def test_differentiated_operator_coefficients(self) -> None:
        payload = MODULE.build_payload(0)
        polynomials = payload["inhomogeneous_tangent_equation"][
            "operator_polynomials_ascending_lambda_power"
        ]
        self.assertEqual(polynomials["D3"]["L6"], ["0", "0", "9", "-18", "9"])
        self.assertEqual(
            polynomials["D3"]["d_kappa_L_at_6"],
            ["0", "0", "9/2", "-9", "9/2"],
        )
        self.assertEqual(polynomials["D1"]["L6"], ["-6", "8", "-8"])
        self.assertEqual(
            polynomials["D1"]["d_kappa_L_at_6"],
            ["-1", "-2/3", "2/3"],
        )

    def test_high_branch_matches_known_percolation_solution(self) -> None:
        _, coefficients = MODULE.frobenius_high_branch(3)
        # These are the first coefficients of Cai's
        # lambda^2 (1-lambda)^2 3F2(...;4 lambda(1-lambda)).
        self.assertEqual(
            [coefficient.value for coefficient in coefficients],
            [Fraction(1), Fraction(1, 3), Fraction(37, 198), Fraction(112, 891)],
        )
        self.assertEqual(coefficients[1].derivative, Fraction(1, 9))
        self.assertEqual(coefficients[2].derivative, Fraction(185, 2178))

    def test_dual_frobenius_recurrence_has_zero_residual(self) -> None:
        exponent, coefficients = MODULE.frobenius_high_branch(6)
        operator = MODULE.cle_operator(MODULE.Dual(Fraction(6), Fraction(1)))
        for m in range(7):
            residual = MODULE.Dual(Fraction(0))
            for derivative_order, polynomial in operator.items():
                for power, polynomial_coefficient in enumerate(polynomial):
                    n = m - 1 + derivative_order - power
                    if 0 <= n < len(coefficients):
                        residual += (
                            polynomial_coefficient
                            * MODULE.falling(exponent + n, derivative_order)
                            * coefficients[n]
                        )
            self.assertEqual(residual.value, 0, msg=f"value residual m={m}")
            self.assertEqual(
                residual.derivative, 0, msg=f"derivative residual m={m}"
            )

    def test_artifact_is_reproducible(self) -> None:
        expected = json.loads(
            (
                ROOT
                / "predictions"
                / "p263_boundary_tangent_ode_20260829.json"
            ).read_text()
        )
        self.assertEqual(MODULE.build_payload(6), expected)


if __name__ == "__main__":
    unittest.main()
