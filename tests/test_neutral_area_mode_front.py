from __future__ import annotations

import math
import sys
import unittest
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import yaml


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analyze_neutral_area_mode_front import (  # noqa: E402
    area_covector,
    exact_n10_oracle,
    first_mode_within_tolerance,
    mode_front,
)
from threshold_score_modes import binomial_weights, krawtchouk_mode  # noqa: E402


class NeutralAreaModeFrontTests(unittest.TestCase):
    def test_exact_n10_half_filling_oracle(self) -> None:
        oracle = exact_n10_oracle()
        self.assertEqual(oracle["N"], 10)
        self.assertEqual(oracle["canonical_neutral_area"], "1/7")
        self.assertEqual(oracle["mean_rank_gap"], "11/7")
        self.assertTrue(oracle["odd_modes_zero"])
        self.assertEqual(oracle["R_epsilon"], {"0.05": 6, "0.10": 6})
        self.assertEqual(
            oracle["area_contributions"],
            [
                "5/16", "0", "-5/16", "0", "3/16", "0",
                "-5/112", "0", "0", "0", "0",
            ],
        )

    def test_general_covector_reconstructs_arbitrary_finite_curve(self) -> None:
        mp.mp.dps = 70
        n = 12
        p0 = mp.mpf("0.592746050790")
        curve = [mp.mpf((3 * k * k + 2 * k + 7) % 17) / 19 for k in range(n + 1)]
        weights = binomial_weights(n, p0)
        coefficients = []
        for order in range(n + 1):
            coefficients.append(
                mp.fsum(
                    weights[k]
                    * curve[k]
                    * krawtchouk_mode(n, k, order, p0)
                    for k in range(n + 1)
                )
            )
        reconstructed = mp.fsum(
            area_covector(n, p0, order) * coefficients[order]
            for order in range(n + 1)
        )
        direct = mp.fsum(curve) / (n + 1)
        self.assertLess(abs(reconstructed - direct), mp.mpf("1e-55"))

    def test_half_filling_covector_kills_odd_orders(self) -> None:
        mp.mp.dps = 50
        n = 16
        p0 = mp.mpf(1) / 2
        for order in range(n + 1):
            expected = (
                mp.mpf(0)
                if order % 2
                else mp.sqrt(math.comb(n, order)) / (order + 1)
            )
            self.assertLess(
                abs(area_covector(n, p0, order) - expected), mp.mpf("1e-45")
            )

    def test_epsilon_rule_and_exact_tail(self) -> None:
        mp.mp.dps = 50
        curve = [mp.mpf(value) for value in (0, 1, 1, 0)]
        area = mp.fsum(curve) / 4
        result = mode_front(
            curve, 3, mp.mpf(1) / 2, area, 3, (mp.mpf("0.05"),)
        )
        final_tail = mp.mpf(result["tail_after_order"][-1])
        self.assertLess(abs(final_tail), mp.mpf("1e-45"))
        self.assertIsNotNone(
            result["epsilon_diagnostics"]["0.05"]["R_epsilon"]
        )
        self.assertEqual(
            first_mode_within_tolerance(
                [mp.mpf("0.2"), mp.mpf("0.09"), mp.mpf("0.04")],
                mp.mpf("0.05"),
            ),
            2,
        )

    def test_manifest_freezes_bounded_diagnostic(self) -> None:
        manifest = yaml.safe_load(
            (ROOT / "analysis" / "neutral_area_mode_front_manifest.yaml").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(manifest["mode_front_diagnostic"]["epsilons"], ["0.05", "0.10"])
        self.assertEqual(manifest["score_basis"]["maximum_order"], 16)
        self.assertEqual(manifest["issues"], [182, 180])


if __name__ == "__main__":
    unittest.main()
