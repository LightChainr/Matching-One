from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"scripts"))
from p267_rank_clock_width import width_values


class RankClockWidthTests(unittest.TestCase):
    def test_exact_signed_moment_identity(self):
        # Includes a negative profile weight: the identity is algebraic.
        f = [F(0), F(3), F(-1), F(2), F(0)]
        n = len(f)-1
        canonical = [sum(v*F(1, n+1) for v in f),
                     sum(v*F(j+1, (n+1)*(n+2)) for j, v in enumerate(f)),
                     sum(v*F((j+1)*(j+2), (n+1)*(n+2)*(n+3)) for j, v in enumerate(f))]
        step = [sum(v*F(1, n) for v in f),
                sum(v*F(2*j+1, 2*n*n) for j, v in enumerate(f)),
                sum(v*F(3*j*j+3*j+1, 3*n**3) for j, v in enumerate(f))]
        cv = canonical[2]/canonical[0]-(canonical[1]/canonical[0])**2
        sv = step[2]/step[0]-(step[1]/step[0])**2
        out = width_values(np.array(f, dtype=float))
        np.testing.assert_allclose(out[:3], np.array([float(cv), float(sv), float(cv-sv)])*n**.75, rtol=1e-14)
        self.assertAlmostEqual(out[0], out[4]+out[5], places=14)


if __name__ == "__main__":
    unittest.main()
