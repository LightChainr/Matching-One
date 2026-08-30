from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_n680_heldout import score  # noqa: E402


class N680HeldoutScoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        raw = ROOT / "results/server-20260830/P337-N680-heldout/raw"
        cls.payload = score(
            ROOT / "analysis/p337_n680_heldout_preregistration.json",
            raw / "n680_120m.births.csv.gz",
            raw / "n680_120m.metadata.json",
        )

    def test_production_matches_freeze(self) -> None:
        self.assertTrue(self.payload["freeze_gates"]["passed"])
        self.assertTrue(self.payload["exact_gates"]["passed"])
        self.assertEqual(self.payload["source"]["births_compression"], "gzip")
        self.assertEqual(
            self.payload["source"]["births_uncompressed_sha256"],
            "037575489db0b7139c9cc4af95aafa200aa976a7d67e1cad16ddbf198e54f4c6",
        )

    def test_heldout_coordinate(self) -> None:
        split = self.payload["decomposition"]
        self.assertAlmostEqual(split["H4_amplitude"]["value"], -0.0021675562817575636)
        self.assertAlmostEqual(split["H4_amplitude"]["standard_error"], 0.00055692853492321)
        self.assertAlmostEqual(split["A_projective_scalar"]["value"], -0.0003116496715883598)

    def test_frozen_forecast_ranking(self) -> None:
        scores = {row["name"]: row for row in self.payload["fixed_model_scores"]}
        self.assertLess(abs(scores["two_mode_recurrence"]["measurement_only_z"]), 0.59)
        self.assertGreater(scores["single_free_lambda"]["measurement_only_z"], 1.56)
        self.assertLess(scores["single_frozen_lambda0"]["measurement_only_z"], -2.04)
        self.assertGreater(scores["scale_neutral"]["measurement_only_z"], 12.1)
        self.assertEqual(
            self.payload["reading"]["closest_fixed_target_by_measurement_residual"],
            "two_mode_recurrence",
        )

    def test_controls_stay_directional_not_overclaimed(self) -> None:
        self.assertLess(abs(self.payload["decomposition"]["A_projective_scalar"]["z"]), 0.4)
        pair = self.payload["primary_pair_sign_flip"]
        self.assertLess(pair["observed"], 0.0)
        self.assertLess(pair["z_vs_scalar_zero"], -3.89)
        self.assertFalse(pair["resolved"])


if __name__ == "__main__":
    unittest.main()
