"""Small exact-observer test, with no source replay or sampling."""
import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"scripts"))
from p334_full_global_conditional_clock import observation_fields


class CompleteObserverIdentity(unittest.TestCase):
    def test_strata_and_fixed_origin_add_back_pathwise(self):
        ranks = np.array([[0, 1], [1, 2], [2, 0], [1, 1]])
        f1 = np.array([[.2, .4, .3, .5], [.7, .8, .6, .7],
                       [.9, .3, .8, .4], [.8, .7, .7, .6]])
        f2 = f1*.4
        hybrid = f2+.01
        fields = observation_fields(ranks, f1, f2, hybrid, f2, -.7, -.185)
        for variant in ("baseline", "safe"):
            for ep in ("p_ref", "p_integral"):
                total = fields[f"{variant}.{ep}.A_H4"]
                np.testing.assert_allclose(sum(fields[f"{variant}.{ep}.A_R{r}_H4"] for r in range(3)), total, atol=1e-15)
            np.testing.assert_allclose(sum(fields[f"{variant}.p_integral.centered_A_R{r}_H4"] for r in range(3)), fields[f"{variant}.p_integral.A_H4"], atol=1e-15)


if __name__ == "__main__":
    unittest.main()
