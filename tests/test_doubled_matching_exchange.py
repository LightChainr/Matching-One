from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from doubled_matching_exchange import (  # noqa: E402
    exchange_matrix,
    identity,
    inverse_2x2,
    lift_exchange_vector,
    matmul,
    matvec,
    validate_contract,
)


class DoubledMatchingExchangeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = json.loads(
            (ROOT / "analysis" / "doubled_matching_exchange_contract.json").read_text(
                encoding="utf-8"
            )
        )

    def test_checked_in_doubled_exchange_closes(self) -> None:
        result = validate_contract(self.contract)
        self.assertTrue(result["matching_map_invertible"])
        self.assertTrue(result["matching_map_intertwines_rg"])
        self.assertTrue(result["doubled_exchange_involutive"])
        self.assertTrue(result["doubled_exchange_commutes_with_rg"])
        self.assertTrue(result["degenerate_rg_does_not_fix_identification"])
        self.assertFalse(result["contains_physical_parity_assignment"])

    def test_exchange_eigenvectors_require_an_explicit_identification(self) -> None:
        a = ((Fraction(1), Fraction(1)), (Fraction(0), Fraction(1)))
        exchange = exchange_matrix(a, inverse_2x2(a))
        vector = (Fraction(2), Fraction(-1))
        for parity in (-1, 1):
            lifted = lift_exchange_vector(vector, a, parity)
            self.assertEqual(matvec(exchange, lifted), tuple(parity * value for value in lifted))

    def test_inverse_and_involution_are_exact(self) -> None:
        a = ((Fraction(2), Fraction(1)), (Fraction(1), Fraction(1)))
        a_inverse = inverse_2x2(a)
        self.assertEqual(matmul(a, a_inverse), identity(2))
        exchange = exchange_matrix(a, a_inverse)
        self.assertEqual(matmul(exchange, exchange), identity(4))

    def test_singular_identification_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "singular"):
            inverse_2x2(((Fraction(1), Fraction(2)), (Fraction(2), Fraction(4))))

    def test_intertwiner_drift_fails_closed(self) -> None:
        changed = deepcopy(self.contract)
        changed["RG_G_star"][0][1] = "3"
        with self.assertRaisesRegex(ValueError, "does not intertwine"):
            validate_contract(changed)

    def test_declared_inverse_drift_fails_closed(self) -> None:
        changed = deepcopy(self.contract)
        changed["matching_map_A_inverse"][0][1] = "0"
        with self.assertRaisesRegex(ValueError, "declared inverse drifted"):
            validate_contract(changed)


if __name__ == "__main__":
    unittest.main()
