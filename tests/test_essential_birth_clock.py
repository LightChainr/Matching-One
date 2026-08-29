from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analyze_essential_birth_clock import (  # noqa: E402
    DEFAULT_INPUTS,
    build_report,
    cos4,
    row_coordinates,
)


class EssentialBirthClockTests(unittest.TestCase):
    def test_axis_and_diagonal_have_opposite_h4_phase(self) -> None:
        self.assertEqual(cos4(1, 0), 1.0)
        self.assertEqual(cos4(1, 1), -1.0)

    def test_coordinate_moment_bridge(self) -> None:
        # Two samples: (K1,K2)=(2,4),(3,3), with N=4.
        row = {
            "samples": "2",
            "n": "4",
            "sum_kminus": "5",
            "sum_kplus": "7",
            "sum_kminus2": "13",
            "sum_kplus2": "25",
            "sum_product": "17",
            "sum_gap": "2",
            "sum_gap2": "4",
        }
        result = row_coordinates(row)
        self.assertAlmostEqual(result["C"], 0.6)
        self.assertAlmostEqual(result["W"], 0.2)
        self.assertAlmostEqual(result["var_C"], 0.0)
        self.assertAlmostEqual(result["var_W"], 0.04)
        self.assertAlmostEqual(result["cov_CW"], 0.0)

    def test_canonical_archives_have_aligned_pair_batches(self) -> None:
        report = build_report([ROOT / path for path in DEFAULT_INPUTS])
        self.assertEqual(len(report["datasets"]), len(DEFAULT_INPUTS))
        self.assertTrue(all(row["batches"] >= 10 for row in report["datasets"]))
        self.assertTrue(
            all(math.isfinite(row["H4_contrast"]["delta_C"]) for row in report["datasets"])
        )

    def test_report_keeps_scaling_claim_exploratory(self) -> None:
        report = build_report([ROOT / path for path in DEFAULT_INPUTS])
        for coordinate in ("C", "W"):
            self.assertIn(
                "retrospective",
                report["retrospective_scaling_diagnostics"][coordinate]["boundary"],
            )

    def test_topological_clock_split_has_finite_fixed_power_scores(self) -> None:
        report = build_report([ROOT / path for path in DEFAULT_INPUTS])
        scores = report["retrospective_scaling_diagnostics"][
            "high_statistics_fixed_power_scores"
        ]
        self.assertEqual(scores["C_N_minus_13_over_8"]["dof"], 7)
        self.assertEqual(scores["W_N_minus_11_over_8"]["dof"], 7)
        self.assertLess(scores["C_N_minus_13_over_8"]["common_scaled_amplitude"], 0)
        self.assertGreater(scores["W_N_minus_11_over_8"]["common_scaled_amplitude"], 0)


if __name__ == "__main__":
    unittest.main()
