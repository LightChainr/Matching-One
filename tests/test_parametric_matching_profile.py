
from __future__ import annotations
from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_jet_algebra import (  # noqa: E402
    compose_derivatives,
    parametric_derivatives,
)
from parametric_matching_profile import validate_contract  # noqa: E402


class ParametricMatchingProfileTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = json.loads(
            (ROOT / "analysis" / "parametric_matching_profile_contract.json").read_text(
                encoding="utf-8"
            )
        )

    def test_checked_in_contract_eliminates_both_bare_coordinates(self) -> None:
        result = validate_contract(self.contract)
        self.assertEqual(result["common_parametric_derivatives"], ["0", "2", "0", "-16", "0", "384"])
        self.assertEqual(result["normalized_odd_invariants"], {"3": "-2", "5": "12"})
        self.assertEqual(len(result["representations_verified"]), 2)
        self.assertFalse(result["contains_empirical_result"])

    def test_exact_composition_and_elimination_round_trip(self) -> None:
        profile = [Fraction(value) for value in (3, -2, 5, 7, -11, 13)]
        coordinate = [Fraction(value) for value in (0, -3, 4, 2, -5, 6)]
        matching = compose_derivatives(profile, coordinate, 5)
        self.assertEqual(parametric_derivatives(matching, coordinate, 5), profile)

    def test_singular_reference_and_short_jets_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "first derivative must be nonzero"):
            parametric_derivatives(
                [Fraction(0), Fraction(1), Fraction(0)],
                [Fraction(0), Fraction(0), Fraction(1)],
                2,
            )
        with self.assertRaisesRegex(ValueError, "cover the requested order"):
            parametric_derivatives([Fraction(0), Fraction(1)], [Fraction(0), Fraction(1)], 2)

    def test_noncanonical_fraction_is_rejected(self) -> None:
        changed = deepcopy(self.contract)
        changed["expected_parametric_derivatives"][3] = "-32/2"
        with self.assertRaisesRegex(ValueError, "must be reduced"):
            validate_contract(changed)


if __name__ == "__main__":
    unittest.main()
