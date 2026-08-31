from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from scripts.cyclic_deck_charge_selection import (
    build_contract,
    invariant_allowed,
    minimal_neutral_tensor_power,
    orbit_cancels,
    residue_counts,
    selection_matrix,
    total_charge,
    validate_contract,
)


class CyclicDeckChargeSelectionTests(unittest.TestCase):
    def test_checked_in_contract_closes_exactly(self) -> None:
        self.assertEqual(validate_contract(), build_contract())

    def test_linear_nontrivial_characters_cancel_exactly(self) -> None:
        for order in (2, 5):
            for charge in range(1, order):
                self.assertFalse(invariant_allowed(order, (charge,)))
                self.assertTrue(orbit_cancels(order, (charge,)))
                self.assertEqual(residue_counts(order, (charge,)), [1] * order)

    def test_only_conjugate_charge_pairs_enter_a_scalar(self) -> None:
        for order in (2, 5):
            matrix = selection_matrix(order)
            for left in range(order):
                for right in range(order):
                    self.assertEqual(matrix[left][right], (left + right) % order == 0)
        self.assertTrue(invariant_allowed(5, (4, 1)))
        self.assertFalse(invariant_allowed(5, (1, 1)))

    def test_minimal_neutral_tensor_powers_are_exact(self) -> None:
        self.assertEqual(minimal_neutral_tensor_power(2, 1), 2)
        for charge in range(1, 5):
            power = minimal_neutral_tensor_power(5, charge)
            self.assertEqual(power, 5)
            self.assertEqual(total_charge(5, (charge,) * power), 0)

    def test_contract_drift_fails_closed(self) -> None:
        frozen = build_contract()
        frozen["C5"]["minimal_tensor_power"]["1"] = 4
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "drift.json"
            path.write_text(json.dumps(frozen), encoding="utf-8")
            with self.assertRaisesRegex(AssertionError, "contract drifted"):
                validate_contract(path)


if __name__ == "__main__":
    unittest.main()
