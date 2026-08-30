from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_n340_second_child import score  # noqa: E402


class N340SecondChildScoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        raw = ROOT / "results/server-20260830/P337-N340-second-child/raw"
        cls.payload = score(
            ROOT / "analysis/p337_n340_second_child_preregistration.json",
            raw / "n340_12m.births.csv",
            raw / "n340_12m.metadata.json",
        )

    def test_production_matches_freeze(self) -> None:
        self.assertTrue(self.payload["freeze_gates"]["passed"])
        self.assertTrue(self.payload["exact_gates"]["passed"])

    def test_heldout_coordinate(self) -> None:
        natural = self.payload["natural_coordinate"]
        self.assertAlmostEqual(natural["value"][0], -0.0025808238159022294)
        self.assertAlmostEqual(natural["value"][1], 0.005163900638930807)
        self.assertAlmostEqual(natural["value"][2], 0.007744724454833037)
        split = self.payload["decomposition"]
        self.assertAlmostEqual(split["H4_amplitude"]["value"], -0.0048572599119938095)
        self.assertAlmostEqual(split["A_projective_scalar"]["value"], 0.0010677461875082945)

    def test_fixed_target_discriminator(self) -> None:
        scores = {row["name"]: row for row in self.payload["fixed_model_scores"]}
        self.assertLess(abs(scores["nominal_area_H4"]["measurement_only_z"]), 1.01)
        self.assertGreater(scores["scale_neutral"]["measurement_only_z"], 5.0)
        self.assertGreater(scores["scale_neutral"]["predictive_z"], 3.9)
        self.assertEqual(
            self.payload["reading"]["closest_fixed_target_by_measurement_residual"],
            "nominal_area_H4",
        )

    def test_scalar_remains_null_and_sign_is_directional(self) -> None:
        scalar = self.payload["decomposition"]["A_projective_scalar"]
        self.assertLess(abs(scalar["z"]), 0.51)
        pair = self.payload["primary_pair_sign_flip"]
        self.assertGreater(pair["observed"], 0.0)
        self.assertGreater(pair["z_vs_scalar_zero"], 3.8)
        self.assertFalse(pair["resolved"])


if __name__ == "__main__":
    unittest.main()
