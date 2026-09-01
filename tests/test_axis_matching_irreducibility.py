from __future__ import annotations

import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from certify_axis_matching_irreducibility import (  # noqa: E402
    CERTIFICATE_PRIMES as AXIS_PRIMES,
    load_axis_coefficients,
    rabin_irreducible,
)
from certify_diamond_matching_irreducibility import (  # noqa: E402
    CERTIFICATE_PRIMES as DIAMOND_PRIMES,
    load_diamond_coefficients,
)


class ExactMatchingIrreducibilityTests(unittest.TestCase):
    def test_axis_l2_through_l5_certificates(self) -> None:
        coefficients = load_axis_coefficients()
        self.assertEqual(sorted(coefficients), [2, 3, 4, 5])
        self.assertEqual(
            {L: len(values) - 1 for L, values in coefficients.items()},
            {2: 4, 3: 9, 4: 16, 5: 25},
        )
        for L, values in sorted(coefficients.items()):
            result = rabin_irreducible(values, AXIS_PRIMES[L])
            self.assertTrue(result["x_qn_equals_x"], msg=f"axis L={L}")
            self.assertTrue(
                all(row["gcd_is_one"] for row in result["gcd_checks"]),
                msg=f"axis L={L}",
            )
            self.assertTrue(result["irreducible"], msg=f"axis L={L}")

    def test_diamond_l1_through_l3_certificates(self) -> None:
        coefficients = load_diamond_coefficients()
        self.assertEqual(sorted(coefficients), [1, 2, 3])
        self.assertEqual(
            {L: len(values) - 1 for L, values in coefficients.items()},
            {1: 2, 2: 8, 3: 18},
        )
        for L, values in sorted(coefficients.items()):
            result = rabin_irreducible(values, DIAMOND_PRIMES[L])
            self.assertTrue(result["x_qn_equals_x"], msg=f"diamond L={L}")
            self.assertTrue(
                all(row["gcd_is_one"] for row in result["gcd_checks"]),
                msg=f"diamond L={L}",
            )
            self.assertTrue(result["irreducible"], msg=f"diamond L={L}")

    def test_reducible_polynomial_is_rejected(self) -> None:
        # (x+1)(x+2)=x^2+3x+2 over F5.
        result = rabin_irreducible([2, 3, 1], 5)
        self.assertFalse(result["irreducible"])

    def test_linear_polynomial_is_trivially_irreducible(self) -> None:
        result = rabin_irreducible([86, 1], 331)
        self.assertTrue(result["x_qn_equals_x"])
        self.assertTrue(result["irreducible"])


if __name__ == "__main__":
    unittest.main()
