#!/usr/bin/env python3
"""Fast algebraic gates for #633. Does not invert the N725 histogram."""
from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import threshold_quantile_lineage as L  # noqa: E402


class TestN725Spin0Algebra(unittest.TestCase):
    def test_lineage_cos4_and_weights_match_issue_633(self):
        c1 = L.cos_four_theta(26, 7)
        c2 = L.cos_four_theta(23, 14)
        self.assertAlmostEqual(c1, 0.4958535077288942, places=15)
        self.assertAlmostEqual(c2, -0.5780680142687277, places=15)
        w = L.spin_zero_weights({"first": c1, "second": c2})
        self.assertAlmostEqual(w["first"], 0.5382777069160998, places=15)
        self.assertAlmostEqual(w["second"], 0.4617222930839002, places=15)
        self.assertAlmostEqual(w["first"] + w["second"], 1.0, places=15)
        self.assertAlmostEqual(w["first"] * c1 + w["second"] * c2, 0.0, places=15)
        self.assertTrue(L.is_interpolation(w))

    def test_buggy_n725_zflow_formula_is_cos2(self):
        # The withdrawn formula, not imported from n725_zflow.py.
        c1 = (26 * 26 - 7 * 7) / 725
        c2 = (23 * 23 - 14 * 14) / 725
        w1 = -c2 / (c1 - c2)
        w2 = c1 / (c1 - c2)
        self.assertAlmostEqual(w1, -1.1326530612244896, places=12)
        self.assertAlmostEqual(w2, 2.13265306122449, places=12)
        true_c1 = L.cos_four_theta(26, 7)
        true_c2 = L.cos_four_theta(23, 14)
        leak = w1 * true_c1 + w2 * true_c2
        self.assertAlmostEqual(leak, -1.7944485136741972, places=12)

    def test_corrected_artifact_if_present(self):
        path = ROOT / "results/probe-invariant-shape/n725-zflow-corrected.json"
        if not path.exists():
            self.skipTest("corrected artifact not generated")
        data = json.loads(path.read_text())
        self.assertEqual(data["issue"], 633)
        self.assertEqual(len(data["levels_11"]), 11)
        self.assertIn(0.4, data["levels_11"])
        self.assertIn(0.6, data["levels_11"])
        self.assertIn(0.25, data["levels_11"])
        w = data["spin0_weights"]
        self.assertAlmostEqual(w["first"] + w["second"], 1.0, places=12)
        self.assertEqual(data["weight_checks"]["sum_w_cos4"], 0.0)
        self.assertAlmostEqual(
            data["tiny_torus_g_vs_deltaZ"]["unoriented_deg"], 20.375, places=2
        )


if __name__ == "__main__":
    unittest.main()
