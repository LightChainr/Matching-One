from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_intrinsic_quantile_center_n145_n290 import (  # noqa: E402
    FEATURE_ORDER,
    RESIDUAL_ORDER,
    pseudovalue_vectors,
)


class IntrinsicQuantileCenterN145N290Tests(unittest.TestCase):
    def test_feature_and_residual_contracts_are_frozen(self) -> None:
        self.assertEqual(
            FEATURE_ORDER,
            ("Q", "w_0.025_scaled", "w_0.05_scaled", "c_0.025", "c_0.05"),
        )
        self.assertEqual(
            RESIDUAL_ORDER,
            (
                "Q290_minus_frozen_ratio_Q145",
                "scaled_width_drift_0.025",
                "scaled_width_drift_0.05",
            ),
        )

    def test_delete_one_pseudovalues_are_coordinatewise(self) -> None:
        self.assertEqual(
            pseudovalue_vectors([10.0, 20.0], [[9.0, 21.0], [11.0, 19.0]]),
            [[11.0, 19.0], [9.0, 21.0]],
        )

    def test_frozen_primary_score_order_is_unchanged(self) -> None:
        plan = yaml.safe_load(
            (ROOT / "experiments/p50_n145_n290_fullcurve_20260829.yaml").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            plan["frozen_scoring"]["order"],
            [
                "intrinsic_and_thermal_even_DeltaM_transfer",
                "raw_asymptotic_slope_baseline",
                "frozen_scalar_plus_H4_slope_correction",
                "raw_root_ratio_baseline",
                "frozen_induced_root_ratio",
            ],
        )

    def test_committed_score_reaches_all_frozen_levels(self) -> None:
        path = (
            ROOT
            / "results/server-20260829/P50-n145-n290-fullcurve/analysis"
            / "intrinsic_quantile_center_score.json"
        )
        if not path.exists():
            self.skipTest("score artifact is generated after the focused unit test")
        payload = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(payload["frozen"]["u"], [0.025, 0.05])
        self.assertIn("delete-one", payload["covariance_rule"])
        self.assertIn("untouched", payload["primary_score_unchanged"])
        self.assertEqual(
            payload["provenance"]["chronology"],
            "freeze precedes first target-result commit",
        )
        for size in ("N145", "N290"):
            for u in ("0.025", "0.05"):
                level = payload["observations"][size]["levels"][u]
                self.assertAlmostEqual(level["M_minus"], -float(u), places=10)
                self.assertAlmostEqual(level["M_plus"], float(u), places=10)


if __name__ == "__main__":
    unittest.main()
