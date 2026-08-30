import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from diagnose_p28_birth_clock_mixture import log_mixture_decomposition  # noqa: E402


class P28BirthClockMixtureTests(unittest.TestCase):
    def test_log_mixture_identity(self):
        row = log_mixture_decomposition(0.3, 1.7)
        self.assertAlmostEqual(row["log_composite"], 0.0)
        self.assertAlmostEqual(row["reconstruction_error"], 0.0, places=15)
        self.assertAlmostEqual(row["responsibility_K1"] + row["responsibility_K2"], 1.0)

    def test_equal_components_have_zero_separation_entropy(self):
        row = log_mixture_decomposition(2.4, 2.4)
        self.assertAlmostEqual(row["separation_entropy"], 0.0, places=15)
        self.assertAlmostEqual(row["component_shape"], math.log(2.4), places=15)

    def test_unequal_components_have_negative_entropy_term(self):
        row = log_mixture_decomposition(0.1, 2.0)
        self.assertLess(row["separation_entropy"], 0.0)


if __name__ == "__main__":
    unittest.main()
