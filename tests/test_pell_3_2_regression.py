
from __future__ import annotations
from copy import deepcopy
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from validate_pell_3_2_regression import validate_contract  # noqa: E402


class PellThreeTwoRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = json.loads(
            (ROOT / "analysis" / "pell_3_2_exact_regression.json").read_text(encoding="utf-8")
        )

    def test_checked_in_pair_regenerates_exactly(self) -> None:
        result = validate_contract(self.contract)
        self.assertEqual(result["pell_residual"], 1)
        self.assertEqual(result["site_difference"], 1)
        self.assertEqual(result["squared_length_ratio"], "9/8")
        self.assertEqual([row["degree"] for row in result["geometries"]], [9, 8])
        self.assertTrue(all(row["primitive"] for row in result["geometries"]))
        self.assertFalse(result["numerical_roots_evaluated"])
        self.assertFalse(result["contains_monte_carlo_result"])

    def test_axis_and_diamond_coefficient_drift_fail_closed(self) -> None:
        for index, label in ((0, "axis_L3"), (1, "diamond_L2")):
            changed = deepcopy(self.contract)
            changed["geometries"][index]["power_coefficients_ascending"][-1] += 1
            with self.subTest(geometry=label):
                with self.assertRaisesRegex(ValueError, f"{label} power coefficients drifted"):
                    validate_contract(changed)

    def test_numerical_or_monte_carlo_result_fields_are_forbidden(self) -> None:
        for key in ("physical_root", "root_gap", "samples"):
            changed = deepcopy(self.contract)
            changed[key] = "forbidden"
            with self.subTest(key=key):
                with self.assertRaisesRegex(ValueError, "forbidden result keys"):
                    validate_contract(changed)


if __name__ == "__main__":
    unittest.main()
