from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import rank_plane_crosswalk as crosswalk  # noqa: E402


class RankPlaneCrosswalkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = crosswalk.build_report()

    def test_all_three_rank_planes_close_exactly(self) -> None:
        for dataset in self.report["datasets"]:
            self.assertTrue(
                all(
                    abs(value) < 1e-24
                    for value in dataset["exact_basis_identity_residuals"].values()
                )
            )

    def test_full_common_batch_covariance_is_present_and_symmetric(self) -> None:
        dimension = len(crosswalk.METRICS)
        for dataset in self.report["datasets"]:
            for key in (
                "covariance_fixed_center_exact_batch_estimator",
                "covariance_intrinsic_center_first_order_influence",
            ):
                covariance = dataset[key]
                self.assertEqual(len(covariance), dimension)
                self.assertTrue(all(len(row) == dimension for row in covariance))
                for i in range(dimension):
                    self.assertGreaterEqual(covariance[i][i], 0)
                    for j in range(dimension):
                        self.assertAlmostEqual(covariance[i][j], covariance[j][i])

    def test_event_level_clock_covariance_is_reconstructed(self) -> None:
        for dataset in self.report["datasets"]:
            for row in dataset["event_clock_covariance_by_orientation"].values():
                covariance = row["covariance"]
                self.assertEqual(row["metric_order"], ["C", "W"])
                self.assertGreater(covariance[0][0], 0)
                self.assertGreater(covariance[1][1], 0)
                self.assertEqual(covariance[0][1], covariance[1][0])

    def test_high_statistics_h4_has_coherent_canonical_signs(self) -> None:
        counts = self.report["high_statistics_summary"]["sign_counts"]
        self.assertEqual(counts["P4_C"]["negative"], 8)
        self.assertEqual(counts["P4_W"]["positive"], 8)
        self.assertEqual(counts["P4_A_top"]["positive"], 8)
        self.assertEqual(counts["P4_E_top"]["negative"], 8)
        self.assertEqual(counts["P4_D_birth"]["positive"], 8)
        self.assertEqual(counts["P4_S_birth"]["negative"], 8)

    def test_k2_cancellation_is_state_localized_not_even_projection(self) -> None:
        summary = self.report["high_statistics_summary"]
        self.assertTrue(summary["K2_cancellation_answer"].startswith("No."))
        diagnostics = summary["median_direction_diagnostics"]
        self.assertLess(
            diagnostics["K2_state_cancellation_fraction"],
            diagnostics["K2_clock_cancellation_fraction"],
        )
        self.assertLess(
            diagnostics["K2_state_cancellation_fraction"],
            diagnostics["K2_density_cancellation_fraction"],
        )

    def test_missing_marked_cross_terms_are_not_claimed_recoverable(self) -> None:
        unavailable = self.report["archive_recoverability"]["not_recoverable"]
        self.assertTrue(any("ell/iota" in row for row in unavailable))
        self.assertTrue(any("A_top times J_D4" in row for row in unavailable))

    def test_maximin_mark_acquisition_is_frozen_without_marked_outcomes(self) -> None:
        design = self.report["maximin_next_acquisition"]
        self.assertEqual(
            design["single_size_mark_acquisition"]["selected"]["id"],
            "P50-N145",
        )
        self.assertEqual(
            design["minimal_radial_campaign"]["selected"], "N65_to_N130_q2"
        )

    def test_checked_in_artifacts_reproduce(self) -> None:
        checked = json.loads(
            (ROOT / "results/rank-plane-crosswalk/latest.json").read_text(
                encoding="utf-8"
            )
        )
        note = (ROOT / "notes/rank-plane-crosswalk.md").read_text(encoding="utf-8")
        self.assertEqual(checked, self.report)
        self.assertEqual(note, crosswalk.render_markdown(self.report))


if __name__ == "__main__":
    unittest.main()
