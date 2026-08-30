from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from fit_four_generation_recurrence import analyze  # noqa: E402


class FourGenerationRecurrenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = analyze(
            ROOT / "results/server-20260830/P337-natural-current-scale-N85/score.json",
            ROOT / "results/server-20260830/P337-N170-angle-flip/score.json",
            ROOT / "results/server-20260830/P337-N340-second-child/score.json",
            ROOT / "results/server-20260830/P337-N680-heldout/score.json",
        )

    def test_four_independent_generations(self) -> None:
        data = self.payload["data"]
        self.assertEqual(data["generation_order"], ["N85", "N170", "N340", "N680"])
        self.assertTrue(all(data["covariance"][i][j] == 0.0
                            for i in range(4) for j in range(4) if i != j))

    def test_recurrence_has_one_df_and_passes(self) -> None:
        model = self.payload["models"]["fixed_lambda0_plus_correction"]
        self.assertEqual(model["df"], 1)
        self.assertAlmostEqual(model["lambda1"], 0.27068136054908193, places=6)
        self.assertLess(model["quadratic"], 0.078)
        self.assertGreater(model["gof_p"], 0.78)
        self.assertGreater(model["lambda1_95pct_profile_interval"][0], 0.0)
        self.assertLess(model["lambda1_95pct_profile_interval"][0], 0.0001)

    def test_comparator_gof_and_aic_boundary(self) -> None:
        models = self.payload["models"]
        self.assertLess(models["free_single_lambda"]["quadratic"], 2.0)
        self.assertGreater(models["fixed_lambda0_single"]["quadratic"], 15.8)
        self.assertGreater(models["scale_neutral"]["quadratic"], 68.9)
        self.assertLess(models["free_single_lambda"]["delta_AIC_descriptive"], 0.001)
        self.assertLess(models["fixed_lambda0_plus_correction"]["delta_AIC_descriptive"], 0.1)

    def test_N1360_freeze(self) -> None:
        freeze = self.payload["N1360_forecast_freeze"]
        self.assertEqual(freeze["Smith_classes"], [[4, 340], [4, 340]])
        self.assertEqual(freeze["H4_covectors_exact"], ["4633/7225", "-6887/7225"])
        self.assertAlmostEqual(
            freeze["predictions"]["fixed_lambda0_plus_correction"]["H4_amplitude"],
            -0.0007955691743768322,
            places=6,
        )
        self.assertIn("not authorized", freeze["status"])


if __name__ == "__main__":
    unittest.main()
