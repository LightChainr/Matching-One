
from __future__ import annotations
from fractions import Fraction
import unittest

from scripts.two_sided_hexagonal_pell import (
    Quadratic,
    build_contract,
    exact_limit_error_identity,
    pell_family,
    pell_residual,
    scaled_shape_defect,
    site_count,
    validate_contract,
)


class TwoSidedHexagonalPellTests(unittest.TestCase):
    def test_checked_in_contract_closes_exactly(self) -> None:
        self.assertEqual(validate_contract(), build_contract())

    def test_both_families_and_recurrence_are_exact(self) -> None:
        self.assertEqual(pell_family(1, 4), [(2, 1), (7, 4), (26, 15), (97, 56)])
        self.assertEqual(pell_family(-2, 4), [(1, 1), (5, 3), (19, 11), (71, 41)])
        for eta in (1, -2):
            for p, q in pell_family(eta, 8):
                self.assertEqual(pell_residual(p, q), eta)
                self.assertEqual(site_count(p, q), 2 * p * q)

    def test_scaled_defect_rationalization_identity(self) -> None:
        for eta in (1, -2):
            for p, q in pell_family(eta, 8):
                value = scaled_shape_defect(p, q)
                denominator = Quadratic(Fraction(p), Fraction(q))
                self.assertEqual(value * denominator, Quadratic(Fraction(eta * p)))

    def test_both_limits_are_approached_from_above_exactly(self) -> None:
        for eta in (1, -2):
            for p, q in pell_family(eta, 8):
                self.assertTrue(exact_limit_error_identity(p, q))
        self.assertEqual(Fraction(-2, 2) / Fraction(1, 2), -2)


if __name__ == "__main__":
    unittest.main()
