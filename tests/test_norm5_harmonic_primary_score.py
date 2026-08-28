import sys
import unittest
from pathlib import Path

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_norm5_harmonic_primary import (  # noqa: E402
    MODEL_ORDER,
    load_artifact,
    residual_covariance,
    score_from_summary,
)


class Norm5HarmonicPrimaryScoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        mp.mp.dps = 50
        cls.artifact = load_artifact(
            ROOT / "predictions/gaussian_norm5_harmonic_discrimination_20260828.yaml"
        )

    def test_exact_frozen_model_order_and_ratios(self) -> None:
        self.assertEqual(MODEL_ORDER, ("H4", "H12", "H8", "zero_effect"))
        ratios = self.artifact["ratios"]
        radial = mp.power(5, -mp.mpf(13) / 8)
        self.assertEqual(ratios["H4"], -mp.mpf(14) / 25 * radial)
        self.assertEqual(ratios["H12"], mp.mpf(23506) / 15625 * radial)
        self.assertEqual(ratios["H8"], -mp.mpf(1054) / 625 * radial)

    def test_exact_h4_synthetic_observation_selects_h4(self) -> None:
        ratios = self.artifact["ratios"]
        point = {65: mp.mpf("0.0012"), 85: mp.mpf("0.0008")}
        point[325] = ratios["H4"] * point[65]
        point[425] = ratios["H4"] * point[85]
        covariance = [[mp.mpf(0) for _ in range(4)] for _ in range(4)]
        for i, variance in enumerate(("1e-12", "2e-12", "3e-12", "4e-12")):
            covariance[i][i] = mp.mpf(variance)
        covariance[0][1] = covariance[1][0] = mp.mpf("2e-13")
        rows = score_from_summary(point, covariance, ratios)
        self.assertEqual([row["name"] for row in rows], list(MODEL_ORDER))
        self.assertEqual(mp.mpf(rows[0]["chi_square"]), 0)
        self.assertGreater(mp.mpf(rows[1]["chi_square"]), 0)

    def test_parent_covariance_propagates_with_signed_ratio(self) -> None:
        covariance = [[mp.mpf(0) for _ in range(4)] for _ in range(4)]
        covariance[0][0] = covariance[1][1] = mp.mpf(2)
        covariance[0][1] = covariance[1][0] = mp.mpf("0.5")
        covariance[2][2] = covariance[3][3] = mp.mpf(3)
        ratio = mp.mpf("-0.25")
        residual = residual_covariance(covariance, ratio)
        self.assertEqual(residual[0][0], mp.mpf("3.125"))
        self.assertEqual(residual[1][1], mp.mpf("3.125"))
        self.assertEqual(residual[0][1], mp.mpf("0.03125"))


if __name__ == "__main__":
    unittest.main()
