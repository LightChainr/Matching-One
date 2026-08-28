from __future__ import annotations

from pathlib import Path
import sys
import unittest

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from finite_size_audit import Observation  # noqa: E402
from rational_finite_size_audit import (  # noqa: E402
    Family,
    denominator_poles,
    fit,
    pole_guard,
    predict,
    rolling_rows,
)


class RationalFiniteSizeAuditTests(unittest.TestCase):
    def setUp(self) -> None:
        mp.mp.dps = 80

    def test_recovers_synthetic_pade(self) -> None:
        family = Family("pade_1_1", 1, 1)
        truth = mp.matrix([mp.mpf("0.5927"), mp.mpf("-0.4"), mp.mpf("0.2"), mp.mpf("3")])
        observations = [Observation(n=n, value=predict(n, truth, family)) for n in range(4, 19)]
        estimated = fit(observations, family)
        for n in (7, 19, 30):
            self.assertLess(abs(predict(n, estimated, family) - predict(n, truth, family)), mp.mpf("1e-60"))

    def test_pole_guard_is_predeclared(self) -> None:
        family = Family("pade_1_1", 1, 1)
        # Q(x)=1-25x has a pole at x=0.04, inside [0,1/4^2].
        parameters = mp.matrix([mp.mpf("0.5"), 1, 1, -25])
        poles = denominator_poles(parameters, family)
        rejected, distance = pole_guard(poles, 4)
        self.assertTrue(rejected)
        self.assertEqual(distance, 0)

    def test_rolling_tests_never_exceed_cutoff(self) -> None:
        family = Family("poly_1", 1, 0)
        truth = mp.matrix([mp.mpf("0.59"), mp.mpf("-1"), mp.mpf("0.3")])
        observations = [Observation(n=n, value=predict(n, truth, family)) for n in range(4, 22)]
        rows = rolling_rows(observations, family, 4, 18, 2)
        self.assertTrue(rows)
        self.assertLessEqual(max(int(row["test_max"]) for row in rows), 18)


if __name__ == "__main__":
    unittest.main()
