from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_n170_angle_flip import score  # noqa: E402


class N170AngleFlipScoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        raw = ROOT / "results/server-20260830/P337-N170-angle-flip/raw"
        cls.payload = score(
            ROOT / "analysis/p337_n170_angle_flip_preregistration.json",
            raw / "n170_8m.births.csv",
            raw / "n170_8m.metadata.json",
        )

    def test_production_matches_freeze(self) -> None:
        self.assertTrue(self.payload["freeze_gates"]["passed"])
        self.assertTrue(self.payload["exact_gates"]["passed"])

    def test_natural_coordinate(self) -> None:
        natural = self.payload["natural_coordinate"]
        self.assertAlmostEqual(natural["value"][0], 0.00643470859085471)
        self.assertAlmostEqual(natural["value"][1], -0.011282165118904897)
        self.assertAlmostEqual(natural["value"][2], -0.017716873709759606)
        self.assertAlmostEqual(natural["standard_error"][2], 0.0015556006552306509)

    def test_exact_geometry_sign_flip_is_resolved(self) -> None:
        pair = self.payload["primary_pair_contrast"]
        self.assertLess(pair["observed"], 0.0)
        self.assertLess(pair["z_vs_scalar_zero"], -11.0)
        self.assertLess(pair["z_to_H4"], -3.0)
        self.assertTrue(self.payload["reading"]["geometry_sign_flip_resolved"])

    def test_residual_is_curvature_not_projective_scalar(self) -> None:
        split = self.payload["curvature_projective_decomposition"]
        self.assertAlmostEqual(split["observed"][0], -0.011111494145226835)
        self.assertAlmostEqual(split["observed"][1], -0.0013809779393524303)
        self.assertLess(split["H4_curvature"]["z"], -3.0)
        self.assertLess(abs(split["A_projective_scalar"]["z"]), 0.8)
        self.assertTrue(self.payload["reading"]["H4_amplitude_curvature_resolved"])
        self.assertFalse(self.payload["reading"]["projective_scalar_resolved"])
        self.assertLess(split["basis_invariance_residual"], 1e-12)


if __name__ == "__main__":
    unittest.main()
