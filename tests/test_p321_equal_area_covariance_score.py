#!/usr/bin/env python3
"""Tests for the P321 equal-area covariance scorer and fixed-power fit."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import run_p321_equal_area_rectangles as runner  # noqa: E402
import score_p321_equal_area_rectangles as scorer  # noqa: E402


DESIGN = json.loads(scorer.DEFAULT_DESIGN.read_text(encoding="utf-8"))
ORACLE = json.loads(scorer.DEFAULT_ORACLE.read_text(encoding="utf-8"))
SMOKE = ROOT / "results/local-20260830/P321-equal-area-20k"


class P321EqualAreaCovarianceTests(unittest.TestCase):
    def test_geometry_router_has_one_square_and_four_rectangles(self) -> None:
        for n in (144, 576, 1296):
            square, rectangles = runner.rows_for_n(DESIGN, n)
            self.assertEqual(square["aspect_ratio"], "1/1")
            self.assertEqual(len(rectangles), 4)
            self.assertTrue(all(row["width"] * row["height"] == n for row in rectangles))

    def test_checked_in_smoke_has_full_covariance_and_exact_square_gate(self) -> None:
        result = json.loads((SMOKE / "score.json").read_text(encoding="utf-8"))
        campaign = result["campaigns"][0]
        self.assertEqual(campaign["root_order"], list(scorer.RHO_ORDER))
        self.assertTrue(campaign["square_histograms_byte_identical"])
        self.assertTrue(campaign["square_moments_byte_identical"])
        self.assertEqual(len(campaign["root_covariance"]), 5)
        self.assertTrue(all(len(row) == 5 for row in campaign["root_covariance"]))
        for i in range(5):
            for j in range(5):
                self.assertAlmostEqual(
                    campaign["root_covariance"][i][j],
                    campaign["root_covariance"][j][i],
                )
        self.assertEqual(
            result["scale_fit"]["status"],
            "insufficient_scales_for_fixed_N^-2_N^-3_fit",
        )

    def test_repeated_square_bytes_are_exactly_equal(self) -> None:
        manifest = json.loads((SMOKE / "N144/campaign.json").read_text(encoding="utf-8"))
        signatures = [
            scorer._orientation_bytes(SMOKE / "N144" / run["histogram"], "first")
            for run in manifest["runs"]
        ]
        self.assertTrue(all(value == signatures[0] for value in signatures[1:]))
        second_signatures = [
            scorer._orientation_bytes(SMOKE / "N144" / run["histogram"], "second")
            for run in manifest["runs"]
        ]
        self.assertGreater(len(set(second_signatures)), 1)

    def test_synthetic_fixed_power_fit_recovers_width_conversion_and_e4_curve(self) -> None:
        mp.mp.dps = 80
        oracle_by_rho = {row["rho"]: row for row in ORACLE["records"]}
        rho_values = [1.0, 16 / 9, 9 / 4, 4.0, 9.0]
        amplitude = 0.125
        c_width = [
            amplitude * float(oracle_by_rho[rho]["width_C_over_square_C"])
            for rho in scorer.RHO_ORDER
        ]
        c_n = [value * rho**2 for value, rho in zip(c_width, rho_values)]
        d_n = [0.01 * (index + 1) for index in range(5)]
        pc = 0.59274605079
        campaigns = []
        for n in (144, 576, 1296):
            roots = [pc + c_n[i] * n**-2 + d_n[i] * n**-3 for i in range(5)]
            covariance = [[1e-20 if i == j else 0.0 for j in range(5)] for i in range(5)]
            campaigns.append({"N": n, "roots": roots, "root_covariance": covariance})
        fit = scorer.fit_fixed_model(campaigns, ORACLE)
        self.assertEqual(fit["status"], "fixed_N^-2_N^-3_gls_fit")
        self.assertFalse(fit["free_exponent_fit"])
        self.assertAlmostEqual(fit["pc"], pc, places=13)
        for observed, expected in zip(fit["C_width_equals_C_N_over_rho_squared"], c_width):
            self.assertAlmostEqual(observed, expected, places=6)
        self.assertLess(fit["conditional_thermal_Q4_E4_score"]["chi_square"], 1e-3)
        self.assertEqual(fit["fit_degrees_of_freedom"], 4)


if __name__ == "__main__":
    unittest.main()
