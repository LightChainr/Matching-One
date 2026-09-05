
from __future__ import annotations
from fractions import Fraction
from itertools import product
import unittest

from scripts.descendant_jordan_rank_survival import jordan_nilpotent, multiply
from scripts.jordan_commutant_valuation import (
    build_contract,
    coefficient_valuation,
    commutant_coefficients,
    image_chain_rank,
    polynomial_in_nilpotent,
    validate_contract,
)


class JordanCommutantValuationTests(unittest.TestCase):
    def test_checked_in_contract_closes_exactly(self) -> None:
        self.assertEqual(validate_contract(), build_contract())

    def test_valuation_gives_every_possible_image_chain_rank(self) -> None:
        for rank in range(1, 8):
            for valuation in range(rank):
                coefficients = [Fraction()] * valuation + [Fraction(2)]
                coefficients += [Fraction()] * (rank - len(coefficients))
                self.assertEqual(coefficient_valuation(coefficients), valuation)
                self.assertEqual(image_chain_rank(rank, coefficients), rank - valuation)
            self.assertEqual(image_chain_rank(rank, [Fraction()] * rank), 0)

    def test_small_integer_commutants_are_exactly_upper_toeplitz(self) -> None:
        for rank in (2, 3):
            nilpotent = jordan_nilpotent(rank)
            for entries in product((-1, 0, 1), repeat=rank * rank):
                matrix = [
                    [Fraction(entries[row * rank + column]) for column in range(rank)]
                    for row in range(rank)
                ]
                if multiply(matrix, nilpotent) == multiply(nilpotent, matrix):
                    coefficients = commutant_coefficients(matrix)
                    self.assertEqual(polynomial_in_nilpotent(rank, coefficients), matrix)

    def test_noncommuting_and_oversized_inputs_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "does not commute"):
            commutant_coefficients(
                [[Fraction(1), Fraction()], [Fraction(1), Fraction(1)]]
            )
        with self.assertRaisesRegex(ValueError, "must fit"):
            polynomial_in_nilpotent(2, [Fraction(1)] * 3)


if __name__ == "__main__":
    unittest.main()
