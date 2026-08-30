from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_natural_current_third_scale import score  # noqa: E402


class NaturalCurrentThirdScaleScoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        raw = ROOT / "results/server-20260830/P337-natural-current-third-scale-N145/raw"
        cls.payload = score(
            ROOT / "analysis/p337_natural_current_third_scale_preregistration.json",
            raw / "n145_2p4m.births.csv",
            raw / "n145_2p4m.metadata.json",
        )

    def test_production_matches_freeze(self) -> None:
        self.assertTrue(self.payload["freeze_gates"]["passed"])
        self.assertTrue(self.payload["exact_gates"]["passed"])

    def test_heldout_natural_coordinate(self) -> None:
        natural = self.payload["natural_coordinate"]
        self.assertAlmostEqual(natural["value"][0], -0.010723781470020469)
        self.assertAlmostEqual(natural["value"][1], 0.009587692386372628)
        self.assertAlmostEqual(natural["value"][2], 0.0203114738563931)
        self.assertAlmostEqual(natural["standard_error"][2], 0.003136452035000094)

    def test_scale_neutral_is_closest_primary_target(self) -> None:
        self.assertEqual(
            self.payload["reading"]["closest_primary_predictive_target"],
            "source_fitted_scale_neutral",
        )
        scores = {
            row["name"]: row for row in self.payload["primary_target_comparison"]
        }
        self.assertGreater(scores["zero"]["predictive_quadratic"], 40.0)
        self.assertLess(
            scores["source_fitted_scale_neutral"]["predictive_quadratic"], 0.4
        )
        self.assertGreater(
            scores["source_fitted_project_H4"]["predictive_quadratic"], 4.0
        )

    def test_fast_secondary_transfer_is_disfavored_without_refit(self) -> None:
        secondary = self.payload["secondary_target_comparison"]
        self.assertEqual(len(secondary), 1)
        self.assertEqual(
            secondary[0]["name"], "secondary_post_reveal_effective_transfer"
        )
        self.assertGreater(secondary[0]["predictive_quadratic"], 11.0)
        self.assertAlmostEqual(self.payload["reading"]["ratio_to_N85"], 0.793645961383144)


if __name__ == "__main__":
    unittest.main()
