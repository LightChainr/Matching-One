
from __future__ import annotations
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from formal_jordan_multiplier_semigroup import (  # noqa: E402
    ComplexFraction,
    gaussian_multiply,
    jordan_transfer,
    matrix_multiply,
    quartic_character,
    validate_contract,
)


class FormalJordanMultiplierSemigroupTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = json.loads(
            (
                ROOT
                / "analysis"
                / "formal_jordan_multiplier_semigroup_contract.json"
            ).read_text(encoding="utf-8")
        )

    def test_checked_in_contract_closes_exactly(self) -> None:
        result = validate_contract(self.contract)
        self.assertTrue(result["gaussian_products_close_both_orders"])
        self.assertTrue(result["quartic_character_multiplicative"])
        self.assertTrue(result["common_nilpotent_squares_to_zero"])
        self.assertTrue(result["norm2_norm5_matrix_paths_close"])
        self.assertFalse(result["contains_physical_jordan_claim"])
        self.assertFalse(result["contains_target_data"])

    def test_exact_quartic_character_values(self) -> None:
        self.assertEqual(quartic_character((1, 1)), ComplexFraction(Fraction(-1, 4), Fraction(0)))
        self.assertEqual(
            quartic_character((2, -1)),
            ComplexFraction(Fraction(-7, 625), Fraction(24, 625)),
        )
        self.assertEqual(
            quartic_character((3, 1)),
            ComplexFraction(Fraction(7, 2500), Fraction(-6, 625)),
        )

    def test_direct_and_both_composed_matrix_paths_agree(self) -> None:
        norm2 = (1, 1)
        norm5 = (2, -1)
        norm10 = gaussian_multiply(norm2, norm5)
        direct = jordan_transfer(norm10)
        self.assertEqual(direct, matrix_multiply(jordan_transfer(norm2), jordan_transfer(norm5)))
        self.assertEqual(direct, matrix_multiply(jordan_transfer(norm5), jordan_transfer(norm2)))

    def test_norm10_carries_both_formal_log_generators(self) -> None:
        upper = jordan_transfer((3, 1))[0][1]
        self.assertNotEqual(upper.log2, ComplexFraction())
        self.assertNotEqual(upper.log5, ComplexFraction())
        self.assertEqual(upper.constant, ComplexFraction())

    def test_zero_multiplier_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "zero Gaussian multiplier"):
            quartic_character((0, 0))


if __name__ == "__main__":
    unittest.main()
