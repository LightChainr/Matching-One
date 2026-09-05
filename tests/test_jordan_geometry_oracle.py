
from __future__ import annotations
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from jordan_geometry_oracle import (  # noqa: E402
    discriminant,
    jordan_certificate,
    similarity_transform,
    validate_contract,
    validate_nilpotent_cocycle,
)


class JordanGeometryOracleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = json.loads(
            (ROOT / "analysis" / "jordan_geometry_contract.json").read_text(
                encoding="utf-8"
            )
        )

    def test_checked_in_contract_closes(self) -> None:
        result = validate_contract(self.contract)
        self.assertTrue(result["jordan_certificate"]["is_size_two_jordan_block"])
        self.assertTrue(result["similarity_invariants_preserved"])
        self.assertTrue(result["nilpotent_cocycle_composes"])
        self.assertTrue(result["coalescence_monotone"])
        self.assertFalse(result["contains_empirical_transfer_matrix_claim"])

    def test_scalar_repeated_eigenvalue_is_not_a_jordan_block(self) -> None:
        scalar = ((Fraction(3, 2), Fraction(0)), (Fraction(0), Fraction(3, 2)))
        certificate = jordan_certificate(scalar)
        self.assertTrue(certificate["repeated_eigenvalue"])
        self.assertFalse(certificate["nonzero_nilpotent_part"])
        self.assertFalse(certificate["is_size_two_jordan_block"])

    def test_similarity_preserves_exact_discriminant(self) -> None:
        matrix = ((Fraction(3, 2), Fraction(1)), (Fraction(0), Fraction(3, 2)))
        similarity = ((Fraction(2), Fraction(1)), (Fraction(1), Fraction(1)))
        transformed = similarity_transform(matrix, similarity)
        self.assertEqual(discriminant(transformed), discriminant(matrix))
        self.assertTrue(jordan_certificate(transformed)["is_size_two_jordan_block"])

    def test_coalescing_eigenvector_angles_strictly_decrease(self) -> None:
        result = validate_contract(self.contract)
        records = result["coalescence_family"]
        discriminants = [Fraction(record["discriminant"]) for record in records]
        angles = [Fraction(record["squared_sine_angle"]) for record in records]
        self.assertTrue(all(left > right for left, right in zip(discriminants, discriminants[1:])))
        self.assertTrue(all(left > right for left, right in zip(angles, angles[1:])))

    def test_non_nilpotent_generator_is_rejected(self) -> None:
        generator = ((Fraction(1), Fraction(0)), (Fraction(0), Fraction(0)))
        with self.assertRaisesRegex(ValueError, "square to zero"):
            validate_nilpotent_cocycle(generator, Fraction(2), Fraction(3))


if __name__ == "__main__":
    unittest.main()
