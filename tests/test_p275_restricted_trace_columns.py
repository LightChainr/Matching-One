from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import analyze_p275_restricted_trace_columns as audit  # noqa: E402


class P275RestrictedTraceColumnTests(unittest.TestCase):
    def test_common_normalizer_is_exact_kernel(self) -> None:
        result = audit.restricted_jet_audit()
        self.assertEqual(result["transport_rank"], 4)
        self.assertEqual(result["common_normalizer_kernel_dimension"], 2)
        self.assertTrue(result["common_normalizer_is_annihilated"])

    def test_current_candidate_family_ranks(self) -> None:
        result = audit.restricted_jet_audit()
        families = result["candidate_families"]
        self.assertEqual(families["vacuum_Ward_critical_surface"]["normalized_image_rank"], 3)
        self.assertEqual(
            families["vacuum_Ward_neighbourhood_stronger_case"]["normalized_image_rank"], 2
        )
        self.assertEqual(
            families["thermal_Q4_Jordan_current_envelope"]["normalized_image_rank"], 4
        )
        comparison = result["column_space_comparison"]
        self.assertTrue(comparison["critical_Ward_image_is_contained_in_thermal_envelope"])
        self.assertEqual(comparison["thermal_extra_directions_beyond_critical_Ward"], 1)

    def test_proxy_design_has_fixed_intersection(self) -> None:
        semisimple = audit.block_design(0.5)
        jordan = audit.jordan_design()
        self.assertEqual(audit.matrix_rank(semisimple), 8)
        self.assertEqual(audit.matrix_rank(jordan), 8)
        self.assertEqual(audit.image_intersection_dimension(semisimple, jordan), 4)
        self.assertEqual(audit.matrix_rank(np.column_stack([semisimple, jordan])), 12)

    def test_checked_assets_build_and_replay_archive_score(self) -> None:
        result = audit.build(audit.DEFAULT_MANIFEST)
        self.assertEqual(result["runtime"]["new_random_samples"], 0)
        self.assertEqual(
            result["status"],
            "PARTIAL_COLUMNS_COMPLETE_EXISTING_COVARIANCES_NOT_DIRECTLY_SCOREABLE",
        )
        proxy = result["existing_K1_K2_proxy_score"]
        self.assertAlmostEqual(
            proxy["semisimple_kappa_0p5"]["mahalanobis_chi_square"],
            proxy["archive_crosscheck"]["semisimple_chi_square"],
            places=8,
        )
        self.assertAlmostEqual(
            proxy["Jordan_kappa_1"]["mahalanobis_chi_square"],
            proxy["archive_crosscheck"]["Jordan_chi_square"],
            places=8,
        )


if __name__ == "__main__":
    unittest.main()
