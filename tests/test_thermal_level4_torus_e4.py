import sys
import unittest
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from derive_thermal_level4_torus_e4 import coefficients  # noqa: E402


class ThermalLevel4TorusE4Tests(unittest.TestCase):
    def test_exact_coefficients(self) -> None:
        values = coefficients()
        self.assertEqual(values["L31_over_L4"], Fraction(-2, 1))
        self.assertEqual(values["L22_over_L4"], Fraction(4, 3))
        self.assertEqual(values["Q4_over_L4"], Fraction(493, 3))
        self.assertEqual(values["Q4_over_primary_g2"], Fraction(493, 96))
        self.assertEqual(values["Q4_over_primary_pi4_E4"], Fraction(493, 72))


if __name__ == "__main__":
    unittest.main()
