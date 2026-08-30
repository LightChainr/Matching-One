from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analyze_theory_fixed_dressing_adversary import analyze  # noqa: E402


class TheoryFixedDressingAdversaryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = analyze(
            ROOT / "results/server-20260830/P337-natural-current-scale-N85/score.json",
            ROOT / "results/server-20260830/P337-N170-angle-flip/score.json",
            ROOT / "results/server-20260830/P337-N340-second-child/score.json",
            ROOT / "results/server-20260830/P337-N680-heldout/score.json",
        )

    def test_fixed_identity_dressing_passes(self) -> None:
        model = self.payload["models"]["fixed_identity_dressing"]
        self.assertEqual(model["df"], 2)
        self.assertAlmostEqual(model["lambda_id4"], 2.0 ** (-21.0 / 8.0))
        self.assertLess(model["quadratic"], 0.81)
        self.assertGreater(model["gof_p"], 0.66)
        self.assertLess(model["delta_AIC_descriptive"], 1e-12)

    def test_rank3_is_minimal_and_passes(self) -> None:
        model = self.payload["models"]["rank3_same_base_jordan"]
        self.assertEqual(model["rank"], 3)
        self.assertEqual(model["df"], 1)
        self.assertLess(model["quadratic"], 0.084)
        self.assertGreater(model["gof_p"], 0.77)

    def test_full_forecast_covariance_and_identifiability(self) -> None:
        joint = self.payload["N1360_forecast_joint"]
        self.assertEqual(len(joint["covariance"]), 6)
        rows = {row["right"]: row for row in joint["pairwise_identity_dressing_separation"]}
        self.assertLess(rows["free_lambda_recurrence"]["maximum_source_limited_z"], 0.7)
        self.assertLess(rows["rank3_same_base_jordan"]["maximum_source_limited_z"], 0.6)
        self.assertFalse(rows["free_single_lambda"]["three_sigma_possible_without_refitting_sources"])
        self.assertTrue(rows["fixed_single_lambda0"]["three_sigma_possible_without_refitting_sources"])

    def test_no_production(self) -> None:
        self.assertEqual(
            self.payload["N1360_forecast_joint"]["production_status"],
            "not authorized or started",
        )


if __name__ == "__main__":
    unittest.main()
