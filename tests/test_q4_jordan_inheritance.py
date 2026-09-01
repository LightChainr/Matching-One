import unittest
from fractions import Fraction

from scripts.verify_q4_jordan_inheritance import render


class Q4JordanInheritanceTest(unittest.TestCase):
    def test_exact_thermal_module_data(self) -> None:
        payload = render()
        self.assertEqual(payload["level2_null"]["norm"], "0")
        self.assertEqual(
            payload["level2_null"]["gram"],
            [["5/2", "15/4"], ["15/4", "45/8"]],
        )
        self.assertEqual(
            payload["q4"]["gram"],
            [["65/2", "75/4", "15"], ["75/4", "195/16", "35/4"], ["15", "35/4", "5"]],
        )
        self.assertEqual(payload["q4"]["norm"], "4930")

    def test_jordan_quantum_numbers(self) -> None:
        payload = render()
        jordan = payload["jordan_inheritance"]
        self.assertEqual(jordan["parent_weights"], ["5/8", "5/8"])
        self.assertEqual(jordan["parent_dimension"], "5/4")
        self.assertEqual(jordan["descendant_weights"], ["37/8", "5/8"])
        self.assertEqual(jordan["descendant_dimension"], "21/4")
        self.assertEqual(jordan["descendant_spin"], "4")
        self.assertEqual(jordan["rank"], 2)
        self.assertEqual(
            payload["he_normalization"]["unnormalized_descendant_cross_overlap_factor"],
            "-12325/2",
        )


if __name__ == "__main__":
    unittest.main()
