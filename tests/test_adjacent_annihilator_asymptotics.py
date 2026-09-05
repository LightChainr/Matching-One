
from __future__ import annotations
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from adjacent_annihilator_asymptotics import (  # noqa: E402
    annihilator_profile,
    generalized_binomial,
    one_minus_power,
    validate_contract,
)


class AdjacentAnnihilatorAsymptoticsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = json.loads(
            (ROOT / "analysis" / "adjacent_annihilator_asymptotics_contract.json").read_text(
                encoding="utf-8"
            )
        )

    def test_checked_in_contract_closes_exactly(self) -> None:
        result = validate_contract(self.contract)
        self.assertTrue(result["maps_q_to_w_as_4_plus_q"])
        self.assertFalse(result["uses_root_or_production_data"])
        self.assertFalse(result["selects_winning_exponent"])
        self.assertEqual(
            [row["root_power_w"] for row in result["profiles"]],
            ["11/2", "6", "7", "8", "10"],
        )

    def test_generalized_binomial_supports_half_integer_q(self) -> None:
        self.assertEqual(generalized_binomial(Fraction(-3, 2), 2), Fraction(15, 8))
        series = one_minus_power(Fraction(-3, 2), 3)
        self.assertEqual(series[:3], (Fraction(0), Fraction(-3, 2), Fraction(-15, 8)))

    def test_leading_coefficient_is_q_over_four(self) -> None:
        for q in (Fraction(1, 7), Fraction(3, 2), Fraction(9, 4), Fraction(6)):
            profile = annihilator_profile(q)
            self.assertEqual(Fraction(profile["leading_coefficient"]), q / 4)
            self.assertEqual(Fraction(profile["root_power_w"]), q + 4)


if __name__ == "__main__":
    unittest.main()
