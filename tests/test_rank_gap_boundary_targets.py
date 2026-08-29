from __future__ import annotations

from pathlib import Path
import sys
import unittest

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_rank_gap_boundary_targets import boundary_gls, score_targets  # noqa: E402


class RankGapBoundaryTargetTests(unittest.TestCase):
    def test_gls_recovers_fixed_exponent_plus_integer_offset(self) -> None:
        sizes = [16, 81, 256]
        points = [mp.mpf(2) * mp.power(n, mp.mpf(5) / 8) + 3 for n in sizes]
        covariance = mp.diag([mp.mpf("0.2"), mp.mpf("0.3"), mp.mpf("0.4")]).tolist()
        fit = boundary_gls(sizes, points, covariance, [625, 1296])
        self.assertAlmostEqual(float(fit["parameters"][0]), 2.0, places=12)
        self.assertAlmostEqual(float(fit["parameters"][1]), 3.0, places=12)
        self.assertAlmostEqual(float(fit["source_chi_square"]), 0.0, places=20)
        for n, prediction in zip((625, 1296), fit["target_prediction"]):
            self.assertAlmostEqual(
                float(prediction), float(2 * mp.power(n, mp.mpf(5) / 8) + 3), places=12
            )

    def test_target_score_adds_observation_and_fit_covariance(self) -> None:
        fit = {
            "target_prediction": mp.matrix([10, 20]),
            "target_prediction_covariance": mp.matrix([[1, mp.mpf("0.5")], [mp.mpf("0.5"), 4]]),
        }
        scored = score_targets([12, 17], mp.matrix([[3, 0], [0, 5]]), fit)
        self.assertEqual(scored["total_covariance"], mp.matrix([[4, mp.mpf("0.5")], [mp.mpf("0.5"), 9]]))
        expected = (mp.matrix([2, -3]).T * mp.matrix([[4, mp.mpf("0.5")], [mp.mpf("0.5"), 9]])**-1 * mp.matrix([2, -3]))[0]
        self.assertAlmostEqual(float(scored["joint_chi_square"]), float(expected), places=14)


if __name__ == "__main__":
    unittest.main()
