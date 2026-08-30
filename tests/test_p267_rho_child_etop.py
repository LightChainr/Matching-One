import pathlib
import sys
import unittest


sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))

from p267_rho_child_etop_mc import CHILD_ORDER, exact_gate
from score_p267_rho_child_etop import complex_zero_score, dft, real_zero_score


class RhoChildEtopTest(unittest.TestCase):
    def test_geometry_gate(self):
        gate = exact_gate()
        self.assertTrue(gate["passed"])
        self.assertEqual(gate["child_order"], list(CHILD_ORDER))
        self.assertEqual([row["N"] for row in gate["children"]], [112, 112, 112])

    def test_c3_constant_is_scalar(self):
        self.assertAlmostEqual(abs(complex(dft([1, 1, 1], 1))), 0.0, places=12)

    def test_complex_zero_score(self):
        values = [complex(-1, -1), complex(-1, 1), complex(1, -1), complex(1, 1)]
        score = complex_zero_score(values)
        self.assertAlmostEqual(score["chi_square"], 0.0)
        self.assertAlmostEqual(score["p"], 1.0)

    def test_scalar_gate_is_one_dimensional(self):
        score = real_zero_score([-1.0, -1.0, 1.0, 1.0])
        self.assertEqual(score["dof"], 1)
        self.assertAlmostEqual(score["chi_square"], 0.0)


if __name__ == "__main__":
    unittest.main()
