from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_polynomial_root_certificate import (  # noqa: E402
    build_artifact,
    classify_stationary_points,
    isolate_roots,
    polynomial_divmod,
    square_free_part,
)


class ExactPolynomialRootCertificateTests(unittest.TestCase):
    def test_checked_artifact_reproduces(self) -> None:
        artifact = build_artifact()
        checked = json.loads(
            (ROOT / "analysis" / "exact_polynomial_root_certificate.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(artifact, checked)

    def test_frozen_density_has_one_exact_maximum(self) -> None:
        self.assertEqual(
            classify_stationary_points([Fraction(1, 2), 3, -3, 0]),
            [{"left": "1/2", "right": "1/2", "classification": "strict_maximum"}],
        )

    def test_repeated_and_endpoint_roots_are_distinctly_isolated(self) -> None:
        # x (x-1) (x-1/3)^2
        polynomial = [Fraction(0), Fraction(-1, 9), Fraction(7, 9), Fraction(-5, 3), 1]
        roots = isolate_roots(polynomial, bits=12)
        self.assertEqual(roots[0], (Fraction(0), Fraction(0)))
        self.assertEqual(roots[-1], (Fraction(1), Fraction(1)))
        middle = roots[1]
        self.assertLessEqual(middle[0], Fraction(1, 3))
        self.assertGreaterEqual(middle[1], Fraction(1, 3))
        self.assertLessEqual(middle[1] - middle[0], Fraction(1, 1 << 12))
        self.assertEqual(len(roots), 3)

    def test_multiple_stationary_points_are_classified(self) -> None:
        # derivative = -(x-1/4)(x-1/2)(x-3/4)
        density = [
            Fraction(1),
            Fraction(3, 32),
            Fraction(-11, 32),
            Fraction(1, 2),
            Fraction(-1, 4),
        ]
        classes = classify_stationary_points(density, bits=12)
        self.assertEqual(
            [item["classification"] for item in classes],
            ["strict_maximum", "strict_minimum", "strict_maximum"],
        )

    def test_no_root_and_invalid_inputs_fail_closed(self) -> None:
        self.assertEqual(isolate_roots([1, 1], bits=8), [])
        with self.assertRaisesRegex(ValueError, "bits"):
            isolate_roots([0, 1], bits=0)
        with self.assertRaisesRegex(ValueError, "ordered"):
            from exact_polynomial_root_certificate import open_root_count, sturm_sequence

            open_root_count(sturm_sequence([0, 1]), Fraction(1), Fraction(0))

    def test_polynomial_arithmetic_is_exact(self) -> None:
        quotient, remainder = polynomial_divmod([-1, 0, 1], [-1, 1])
        self.assertEqual(quotient, [1, 1])
        self.assertEqual(remainder, [])
        self.assertEqual(square_free_part([0, 0, 1]), [0, 1])


if __name__ == "__main__":
    unittest.main()

