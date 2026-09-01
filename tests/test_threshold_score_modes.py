import math
import sys
import unittest
from pathlib import Path

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from threshold_score_modes import binomial_weights, krawtchouk_mode  # noqa: E402


class ThresholdScoreModeTests(unittest.TestCase):
    def test_orthonormality_and_positive_first_score(self) -> None:
        mp.mp.dps = 50
        n, p, maximum = 17, mp.mpf("0.592746050790"), 6
        weights = binomial_weights(n, p)
        modes = [[krawtchouk_mode(n, x, r, p) for x in range(n + 1)] for r in range(maximum + 1)]
        for r in range(maximum + 1):
            for s in range(maximum + 1):
                inner = mp.fsum(weights[x] * modes[r][x] * modes[s][x] for x in range(n + 1))
                self.assertLess(abs(inner - int(r == s)), mp.mpf("1e-42"))
        for x in range(n + 1):
            expected = (x - n * p) / mp.sqrt(n * p * (1 - p))
            self.assertLess(abs(modes[1][x] - expected), mp.mpf("1e-45"))

    def test_half_occupation_complement_parity(self) -> None:
        p = mp.mpf("0.5")
        for order in range(7):
            for x in range(15):
                left = krawtchouk_mode(14, 14 - x, order, p)
                right = (-1) ** order * krawtchouk_mode(14, x, order, p)
                self.assertLess(abs(left - right), mp.mpf("1e-40"))


if __name__ == "__main__":
    unittest.main()
