from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_natural_current_scale import score  # noqa: E402


class NaturalCurrentScaleScoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        raw = ROOT / "results/server-20260830/P337-natural-current-scale-N85/raw"
        cls.payload = score(
            ROOT / "analysis/p337_natural_current_scale_preregistration.json",
            raw / "n85_200k.births.csv",
            raw / "n85_200k.metadata.json",
        )

    def test_production_matches_freeze(self) -> None:
        self.assertTrue(self.payload["freeze_gates"]["passed"])
        self.assertTrue(self.payload["exact_gates"]["passed"])

    def test_independent_natural_coordinate(self) -> None:
        natural = self.payload["natural_coordinate"]
        self.assertAlmostEqual(natural["value"][0], -0.010587597656159131)
        self.assertAlmostEqual(natural["value"][1], 0.01500501522754386)
        self.assertAlmostEqual(natural["value"][2], 0.025592612883702993)
        self.assertAlmostEqual(natural["standard_error"][2], 0.007926872163140496)

    def test_frozen_H4_is_closest_predictive_target(self) -> None:
        self.assertEqual(
            self.payload["reading"]["closest_frozen_predictive_target"],
            "source_fitted_project_H4",
        )
        scores = {row["name"]: row for row in self.payload["target_comparison"]}
        self.assertGreater(scores["zero"]["predictive_quadratic"], 10.0)
        self.assertLess(scores["source_fitted_project_H4"]["predictive_quadratic"], 1.6)
        self.assertGreater(scores["source_fitted_scale_neutral"]["measurement_only_quadratic"], 30.0)


if __name__ == "__main__":
    unittest.main()
