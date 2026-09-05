
from __future__ import annotations
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from logodds_derivative_chain import (  # noqa: E402
    complement_values,
    eta_to_p_derivatives,
    p_to_eta_derivatives,
    response_chain,
    validate_contract,
)


class LogoddsDerivativeChainTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = json.loads(
            (ROOT / "analysis" / "logodds_derivative_contract.json").read_text(
                encoding="utf-8"
            )
        )

    def test_checked_in_chain_round_trips_and_obeys_complement_parity(self) -> None:
        result = validate_contract(self.contract)
        self.assertEqual(result["center_probability"], "2/5")
        self.assertEqual(result["complement_center_probability"], "3/5")
        self.assertTrue(result["round_trip_exact"])
        self.assertTrue(result["general_complement_parity_exact"])
        self.assertFalse(result["contains_empirical_result"])

    def test_arbitrary_derivative_jets_round_trip_at_rational_centers(self) -> None:
        p_jet = [Fraction(value) for value in (7, -3, 5, 11, -13, 17, 19)]
        for center in (Fraction(1, 3), Fraction(2, 5), Fraction(5, 7)):
            with self.subTest(center=center):
                eta_jet = p_to_eta_derivatives(p_jet, center, 6)
                self.assertEqual(eta_to_p_derivatives(eta_jet, center, 6), p_jet)

    def test_general_even_and_odd_complement_parity(self) -> None:
        conditional = [Fraction(value) for value in (0, 1, -2, 3, 5, -4, 2)]
        center = Fraction(2, 5)
        _p, eta = response_chain(conditional, center, 6)
        for parity in (-1, 1):
            with self.subTest(parity=parity):
                _cp, complement_eta = response_chain(
                    complement_values(conditional, parity), 1 - center, 6
                )
                self.assertEqual(
                    complement_eta,
                    [parity * (-1) ** order * eta[order] for order in range(7)],
                )

    def test_invalid_centers_parity_and_short_jets_fail_closed(self) -> None:
        for center in (Fraction(0), Fraction(1)):
            with self.subTest(center=center):
                with self.assertRaisesRegex(ValueError, "strictly between"):
                    p_to_eta_derivatives([Fraction(0), Fraction(1)], center, 1)
        with self.assertRaisesRegex(ValueError, "parity must be"):
            complement_values([Fraction(0)], 0)
        with self.assertRaisesRegex(ValueError, "jet is too short"):
            p_to_eta_derivatives([Fraction(0)], Fraction(1, 2), 1)


if __name__ == "__main__":
    unittest.main()
