
from __future__ import annotations
from fractions import Fraction
import unittest

from scripts.landing_registry_spin_separation import (
    build_contract,
    cosine_alias_at_rational_angle,
    matrix_rank,
    response_columns,
    validate_contract,
)


class LandingRegistrySpinSeparationTests(unittest.TestCase):
    def test_checked_in_contract_closes_exactly(self) -> None:
        self.assertEqual(validate_contract(), build_contract())

    def test_pi_over_eight_multiples_are_exactly_the_alias_set(self) -> None:
        aliases = [
            numerator
            for numerator in range(13)
            if cosine_alias_at_rational_angle(4, 12, numerator, 24)
        ]
        self.assertEqual(aliases, [0, 3, 6, 9, 12])

    def test_axis_diagonal_pair_has_rank_one(self) -> None:
        columns = response_columns(False)
        self.assertEqual(columns[0], columns[1])
        self.assertEqual(matrix_rank(columns), 1)

    def test_pi_over_twelve_extension_has_rank_two(self) -> None:
        columns = response_columns(True)
        self.assertEqual(columns[0][-1], Fraction(1, 2))
        self.assertEqual(columns[1][-1], Fraction(-1))
        self.assertEqual(matrix_rank(columns), 2)
        with self.assertRaisesRegex(ValueError, "denominator must be positive"):
            cosine_alias_at_rational_angle(4, 12, 1, 0)


if __name__ == "__main__":
    unittest.main()
