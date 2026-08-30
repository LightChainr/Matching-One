import sys
import unittest
from fractions import Fraction
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from derive_rectangular_thermal_q4_hecke import exact_ratios  # noqa: E402


class RectangularThermalQ4HeckeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.values = exact_ratios()
        artifact_path = (
            ROOT / "predictions" / "rectangular_thermal_q4_hecke_20260829.yaml"
        )
        cls.artifact = yaml.safe_load(artifact_path.read_text(encoding="utf-8"))

    def test_exact_holomorphic_ratios(self) -> None:
        self.assertEqual(self.values["T2_eigenvalue"], Fraction(9, 1))
        self.assertEqual(self.values["E4_2i_over_E4_i"], Fraction(11, 16))
        self.assertEqual(self.values["E4_i_over_2_over_E4_i"], Fraction(11, 1))
        self.assertEqual(self.values["E4_diagonal_over_E4_i"], Fraction(-4, 1))

    def test_exact_area_normalized_ratios(self) -> None:
        self.assertEqual(self.values["E4hat_2i_over_E4hat_i"], Fraction(11, 4))
        self.assertEqual(
            self.values["E4hat_i_over_2_over_E4hat_i"], Fraction(11, 4)
        )
        self.assertEqual(
            self.values["E4hat_diagonal_over_E4hat_i"], Fraction(-1, 1)
        )

    def test_exact_hecke_equation(self) -> None:
        x = self.values["E4_2i_over_E4_i"]
        diagonal = self.values["E4_diagonal_over_E4_i"]
        lhs = 8 * x + Fraction(1, 2) * (16 * x + diagonal)
        self.assertEqual(lhs, self.values["T2_eigenvalue"])

    def test_frozen_artifact_matches_derivation(self) -> None:
        ratios = self.artifact["exact_mathematical_layer"]["ratios"]
        self.assertEqual(ratios["E4_2i_over_E4_i"], "11/16")
        self.assertEqual(ratios["E4hat_2i_over_E4hat_i"], "11/4")
        self.assertEqual(ratios["E4hat_diagonal_over_E4hat_i"], "-1")
        self.assertEqual(
            self.artifact["lattice_observable_bridge"]["status"],
            "conditional_hypothesis",
        )


if __name__ == "__main__":
    unittest.main()
