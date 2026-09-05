
from __future__ import annotations
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from angular_radial_projector_oracle import (  # noqa: E402
    angular_filter,
    compose_angular_then_radial,
    compose_radial_then_angular,
    outer,
    radial_filter,
    validate_contract,
)


class AngularRadialProjectorOracleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = json.loads(
            (ROOT / "analysis" / "angular_radial_projector_contract.json").read_text(
                encoding="utf-8"
            )
        )

    def test_checked_in_contract_closes_exactly(self) -> None:
        result = validate_contract(self.contract)
        self.assertTrue(result["angular_and_radial_filters_commute"])
        self.assertTrue(result["pure_h4_killed_by_each_filter"])
        self.assertTrue(result["pure_h4_double_filter_is_redundant"])
        self.assertTrue(result["mixed_sector_retains_scalar_residue"])
        self.assertFalse(result["contains_production_data"])
        self.assertFalse(result["identifies_h4_exponent"])

    def test_tensor_factor_maps_commute_for_a_generic_grid(self) -> None:
        grid = ((Fraction(1), Fraction(2)), (Fraction(3), Fraction(5)))
        angular = (Fraction(2, 3), Fraction(1, 3))
        radial = (Fraction(-2), Fraction(7))
        self.assertEqual(
            compose_angular_then_radial(grid, angular, radial),
            compose_radial_then_angular(grid, angular, radial),
        )

    def test_each_declared_filter_kills_pure_h4(self) -> None:
        pure_h4 = outer((Fraction(1), Fraction(-1)), (Fraction(8), Fraction(1)))
        self.assertEqual(
            angular_filter(pure_h4, (Fraction(1, 2), Fraction(1, 2))),
            (Fraction(0), Fraction(0)),
        )
        self.assertEqual(
            radial_filter(pure_h4, (Fraction(-1, 8), Fraction(1))),
            (Fraction(0), Fraction(0)),
        )


if __name__ == "__main__":
    unittest.main()
