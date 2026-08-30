from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from square_bond_transport_parity_theorem import (  # noqa: E402
    certify_centered_parity,
    certify_length,
    dual_bond_permutation,
    reversal_fixture,
    validate_contract,
)


class SquareBondTransportParityTheoremTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = json.loads(
            (ROOT / "analysis" / "square_bond_transport_parity_contract.json").read_text(
                encoding="utf-8"
            )
        )

    def test_checked_in_contract_closes_exactly(self) -> None:
        result = validate_contract(self.contract)
        self.assertEqual(result["lengths_checked"], list(range(2, 9)))
        self.assertTrue(result["transport_is_complement_plus_permutation"])
        self.assertTrue(result["transport_even_implies_centered_even"])
        self.assertTrue(result["transport_odd_implies_centered_odd"])
        self.assertFalse(result["enumerates_configurations"])
        self.assertFalse(result["contains_continuum_amplitude_claim"])

    def test_dual_permutation_square_is_translation(self) -> None:
        for length in range(2, 9):
            row = certify_length(length)
            self.assertEqual(len(dual_bond_permutation(length)), 2 * length * length)
            self.assertTrue(row["permutation_is_bijective"])
            self.assertEqual(row["permutation_square_translation"], [-1, -1])

    def test_coefficient_reversal_gives_both_centered_parities(self) -> None:
        for degree in (3, 4, 8, 18, 32):
            for sign, parity in ((1, "even"), (-1, "odd")):
                coefficients = reversal_fixture(degree, sign)
                self.assertTrue(
                    all(
                        coefficients[degree - k] == sign * coefficients[k]
                        for k in range(degree + 1)
                    )
                )
                self.assertEqual(
                    certify_centered_parity(degree, sign)["centered_polynomial_parity"],
                    parity,
                )
        self.assertEqual(reversal_fixture(4, -1)[2], Fraction(0))

    def test_contract_drift_fails_closed(self) -> None:
        changed = deepcopy(self.contract)
        changed["lengths"] = [3, 2]
        with self.assertRaisesRegex(ValueError, "lengths must be sorted"):
            validate_contract(changed)

        changed = deepcopy(self.contract)
        changed["permutation_square_translation"] = [0, 0]
        with self.assertRaisesRegex(ValueError, "translation drifted"):
            validate_contract(changed)


if __name__ == "__main__":
    unittest.main()
