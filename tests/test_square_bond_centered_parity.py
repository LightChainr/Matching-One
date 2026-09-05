
from __future__ import annotations
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from square_bond_centered_parity import (  # noqa: E402
    centered_from_bernstein_sums,
    enumerate_centered_parity,
    validate_contract,
)


class SquareBondCenteredParityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = json.loads(
            (ROOT / "analysis" / "square_bond_centered_parity_contract.json").read_text(
                encoding="utf-8"
            )
        )

    def test_checked_in_contract_closes_exactly(self) -> None:
        result = validate_contract(self.contract)
        self.assertTrue(result["S_is_centered_even_in_every_channel"])
        self.assertTrue(result["D_is_centered_odd_in_every_channel"])
        self.assertEqual(result["common_D_first_derivative_at_half"], "27/8")
        self.assertTrue(result["D_center_zero_does_not_force_D_derivative_zero"])
        self.assertFalse(result["contains_orientation_production_result"])

    def test_all_five_channels_obey_complement_and_taylor_parity(self) -> None:
        result = enumerate_centered_parity(2)
        self.assertTrue(result["passed"])
        self.assertEqual(result["configurations"], 256)
        self.assertEqual(result["distinct_D_centered_polynomials"], 1)
        for row in result["channels"].values():
            self.assertTrue(row["S_bernstein_complement_even"])
            self.assertTrue(row["D_bernstein_complement_odd"])
            self.assertTrue(row["S_all_odd_derivatives_vanish_at_half"])
            self.assertTrue(row["D_all_even_derivatives_vanish_at_half"])

    def test_centered_bernstein_conversion_is_exact(self) -> None:
        # p^2-(1-p)^2 = 2t after p=1/2+t.
        converted = centered_from_bernstein_sums(
            (Fraction(-1), Fraction(0), Fraction(1))
        )
        self.assertEqual(converted, (Fraction(0), Fraction(2), Fraction(0)))

    def test_center_zero_and_nonzero_slope_are_distinct_claims(self) -> None:
        result = enumerate_centered_parity(2)
        for row in result["channels"].values():
            self.assertEqual(row["D_at_half"], "0")
            self.assertEqual(row["D_first_derivative_at_half"], "27/8")


if __name__ == "__main__":
    unittest.main()
