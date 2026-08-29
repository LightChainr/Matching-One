from __future__ import annotations

import sys
import unittest
from pathlib import Path

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_p205_quotient_character_prism import (  # noqa: E402
    MODEL_ORDER,
    SIZES,
    load_contract,
    model_vector,
    score_model,
)


DESIGN = ROOT / "predictions/p205_quotient_character_prism_20260829.yaml"


class P205QuotientCharacterPrismScoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        mp.mp.dps = 60
        cls.contract = load_contract(DESIGN)

    def test_contract_keeps_frozen_order_and_sign_code(self) -> None:
        self.assertEqual(MODEL_ORDER, ("H4", "H8", "H12"))
        self.assertEqual(SIZES, (25, 50, 125))
        signs = {
            model: [
                "+" if self.contract["characters"][model][n] > 0 else "-"
                for n in SIZES
            ]
            for model in MODEL_ORDER
        }
        self.assertEqual(
            signs,
            {"H4": ["+", "+", "+"], "H8": ["+", "-", "+"], "H12": ["+", "+", "-"]},
        )

    def test_exact_h8_line_scores_zero_for_synthetic_h8_data(self) -> None:
        h8 = model_vector("H8", self.contract["characters"])
        points = [mp.mpf("0.75") * value for value in h8]
        variances = [mp.mpf("1e-10"), mp.mpf("2e-10"), mp.mpf("3e-10")]
        score = score_model(points, variances, h8)
        self.assertLess(abs(mp.mpf(score["chi_square"])), mp.mpf("1e-50"))
        self.assertLess(
            mp.mpf(score["chi_square"]),
            mp.mpf(score_model(points, variances, model_vector("H4", self.contract["characters"]))["chi_square"]),
        )

    def test_one_amplitude_fit_has_two_degrees_of_freedom(self) -> None:
        vector = model_vector("H12", self.contract["characters"])
        score = score_model(
            [mp.mpf("1e-4"), mp.mpf("2e-4"), mp.mpf("-1e-4")],
            [mp.mpf("1e-8")] * 3,
            vector,
        )
        self.assertEqual(score["degrees_of_freedom"], 2)
        self.assertEqual(len(score["residuals"]), 3)


if __name__ == "__main__":
    unittest.main()
