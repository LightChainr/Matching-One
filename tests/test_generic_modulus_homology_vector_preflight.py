import sys
import unittest
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from generic_modulus_homology_vector_preflight import (  # noqa: E402
    lattice_stabilizer,
    period_tau,
    primitive_lines,
    rank_readouts,
)


class GenericModulusHomologyVectorPreflightTest(unittest.TestCase):
    def test_base_geometry_is_generic(self):
        matrix = ((10, 3), (0, 10))
        self.assertEqual(period_tau(matrix), (Fraction(3, 10), Fraction(1, 1)))
        self.assertEqual(
            lattice_stabilizer(matrix),
            [((-1, 0), (0, -1)), ((1, 0), (0, 1))],
        )

    def test_primitive_lines_are_canonical_and_unique(self):
        lines = primitive_lines(3)
        self.assertEqual(len(lines), len(set(lines)))
        self.assertIn((1, -3), lines)
        self.assertIn((0, 1), lines)
        self.assertNotIn((-1, 0), lines)

    def test_rank_readout_contract(self):
        result = rank_readouts((0.2, 0.5, 0.3), (-0.1, 0.04, 0.06))
        self.assertAlmostEqual(result["q"], 0.1)
        self.assertAlmostEqual(result["E"], 0.5)
        self.assertAlmostEqual(result["dpq"], 0.16)
        self.assertAlmostEqual(result["dpE"], -0.04)


if __name__ == "__main__":
    unittest.main()
