import importlib.util
import pathlib
import unittest
from fractions import Fraction


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "verify_thermal_q4_torus_ward.py"
SPEC = importlib.util.spec_from_file_location("thermal_q4_torus_ward", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class ThermalQ4TorusWardTest(unittest.TestCase):
    def test_exact_rational_reduction(self):
        c = MODULE.self_test()
        self.assertEqual(c["L-3L-1_over_C"], Fraction(-2, 1))
        self.assertEqual(c["L-2^2_over_C"], Fraction(4, 3))
        self.assertEqual(c["Q4_over_C"], Fraction(493, 3))
        self.assertEqual(c["C_over_g2_primary"], Fraction(1, 32))
        self.assertEqual(c["Q4_over_g2_primary"], Fraction(493, 96))


if __name__ == "__main__":
    unittest.main()
