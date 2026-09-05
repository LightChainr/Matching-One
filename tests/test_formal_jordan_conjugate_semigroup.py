
from __future__ import annotations
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from formal_jordan_conjugate_semigroup import (  # noqa: E402
    certify_multiplier,
    complex_conjugate,
    gaussian_conjugate,
    matrix_conjugate,
    validate_contract,
)
from formal_jordan_multiplier_semigroup import (  # noqa: E402
    ComplexFraction,
    jordan_transfer,
    matrix_multiply,
)


class FormalJordanConjugateSemigroupTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = json.loads(
            (
                ROOT
                / "analysis"
                / "formal_jordan_conjugate_semigroup_contract.json"
            ).read_text(encoding="utf-8")
        )

    def test_checked_in_contract_closes_exactly(self) -> None:
        result = validate_contract(self.contract)
        self.assertTrue(result["conjugation_involutive"])
        self.assertTrue(result["quartic_character_star_compatible"])
        self.assertTrue(result["jordan_transfer_star_compatible"])
        self.assertTrue(result["norm_factorization_exact"])
        self.assertTrue(result["composite_path_star_compatible"])
        self.assertFalse(result["constructs_group_inverse"])
        self.assertFalse(result["contains_physical_jordan_claim"])

    def test_complex_and_gaussian_conjugations_are_involutions(self) -> None:
        gaussian = (2, 1)
        self.assertEqual(gaussian_conjugate(gaussian_conjugate(gaussian)), gaussian)
        value = ComplexFraction(Fraction(3, 7), Fraction(-5, 11))
        self.assertEqual(complex_conjugate(complex_conjugate(value)), value)

    def test_norm5_transfer_factorizes_through_conjugate(self) -> None:
        multiplier = (2, 1)
        conjugate = gaussian_conjugate(multiplier)
        norm_transfer = jordan_transfer((5, 0))
        self.assertEqual(
            matrix_multiply(jordan_transfer(multiplier), jordan_transfer(conjugate)),
            norm_transfer,
        )
        self.assertEqual(
            jordan_transfer(conjugate),
            matrix_conjugate(jordan_transfer(multiplier)),
        )

    def test_all_declared_certificates_close_both_orders(self) -> None:
        for multiplier in ((1, 1), (2, 1), (1, 3)):
            row = certify_multiplier(multiplier)
            self.assertTrue(row["conjugation_intertwines_character"])
            self.assertTrue(row["conjugation_intertwines_transfer"])
            self.assertTrue(row["norm_factorization_closes_both_orders"])


if __name__ == "__main__":
    unittest.main()
