from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analyze_three_generation_recurrence import analyze  # noqa: E402


class ThreeGenerationRecurrenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = analyze(
            ROOT / "results/server-20260830/P337-natural-current-scale-N85/score.json",
            ROOT / "results/server-20260830/P337-N170-angle-flip/score.json",
            ROOT / "results/server-20260830/P337-N340-second-child/score.json",
        )

    def test_input_covariance_is_independent_blocks(self) -> None:
        covariance = self.payload["data"]["covariance"]
        self.assertTrue(all(covariance[i][j] == 0.0 for i in range(3) for j in range(3) if i != j))
        self.assertAlmostEqual(self.payload["data"]["H4_amplitude"][0], -0.01605092257680157)

    def test_exact_recurrence_solution(self) -> None:
        model = self.payload["two_mode_recurrence"]
        self.assertAlmostEqual(model["lambda1"], 0.21240402888481813)
        self.assertAlmostEqual(model["c0"], -0.06888917573635066)
        self.assertAlmostEqual(model["c1"], 0.052838253159549084)
        self.assertTrue(model["point_in_open_unit_interval"])
        self.assertLess(model["z_above_zero"], 1.0)

    def test_correction_decays_and_predicts_N680(self) -> None:
        model = self.payload["two_mode_recurrence"]
        ratios = model["absolute_correction_to_leading_ratio"]
        self.assertTrue(all(left > right for left, right in zip(ratios, ratios[1:])))
        self.assertAlmostEqual(model["N680_H4_amplitude_prediction"], -0.001841297042231899)
        self.assertGreater(model["N680_prediction_standard_error_delta"], 0.0009)

    def test_comparator_order(self) -> None:
        comparators = self.payload["comparators"]
        self.assertLess(comparators["single_free_lambda"]["quadratic"], comparators["single_frozen_lambda0"]["quadratic"])
        self.assertLess(comparators["single_frozen_lambda0"]["quadratic"], comparators["scale_neutral"]["quadratic"])
        self.assertAlmostEqual(comparators["single_free_lambda"]["lambda"], 0.5358337477920293, places=6)


if __name__ == "__main__":
    unittest.main()
