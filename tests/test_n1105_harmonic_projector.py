import importlib.util
import unittest
from fractions import Fraction
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "n1105_harmonic_projector.py"
SPEC = importlib.util.spec_from_file_location("n1105_harmonic_projector", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class N1105HarmonicProjectorTest(unittest.TestCase):
    def test_exact_kronecker_projectors(self) -> None:
        for target in range(4):
            weights = MODULE.projector(target)
            for harmonic in range(4):
                self.assertEqual(
                    MODULE.response(weights, harmonic),
                    Fraction(int(harmonic == target)),
                )

    def test_h0_conditioning_and_declared_higher_leakage(self) -> None:
        weights = MODULE.projector(0)
        l1, l2 = MODULE.norms(weights)
        self.assertLess(l1, 2.0)
        self.assertLess(l2, 1.0)
        self.assertEqual(MODULE.response(weights, 1), 0)
        self.assertEqual(MODULE.response(weights, 2), 0)
        self.assertEqual(MODULE.response(weights, 3), 0)
        self.assertNotEqual(MODULE.response(weights, 4), 0)


if __name__ == "__main__":
    unittest.main()
