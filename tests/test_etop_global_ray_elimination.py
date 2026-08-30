from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from etop_global_ray_elimination import fit_ray  # noqa: E402


def row(name, point, variance="0.01"):
    return {
        "id": name,
        "N": 1,
        "dependency_group": name,
        "estimate": [mp.mpf(str(value)) for value in point],
        "covariance": [
            [mp.mpf(variance), mp.mpf(0)],
            [mp.mpf(0), mp.mpf(variance)],
        ],
    }


class ETopGlobalRayEliminationTests(unittest.TestCase):
    def test_collinear_points_have_zero_profile_discrepancy(self) -> None:
        fit = fit_ray([row("a", (1, -2)), row("b", (3, -6))], grid_size=256)
        self.assertLess(fit["min_chi2"], 1e-20)
        self.assertAlmostEqual(fit["E_over_A_slope"], -2.0, places=10)

    def test_orthogonal_precise_points_reject_one_ray(self) -> None:
        fit = fit_ray(
            [row("a", (1, 0), "0.0001"), row("b", (0, 1), "0.0001")],
            grid_size=256,
        )
        self.assertGreater(fit["min_chi2"], 9000)

    def test_committed_result_has_nested_and_heldout_rows(self) -> None:
        report = json.loads(
            (ROOT / "results/etop-global-ray-elimination/latest.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            report["primary"]["global_vs_four_ray_likelihood_ratio"][
                "degrees_of_freedom"
            ],
            3,
        )
        self.assertEqual(
            len(report["primary"]["leave_one_lineage_out_predictive_profiles"]),
            4,
        )
        self.assertEqual(
            len(report["primary"]["pairwise_independent_geometry_determinants"]),
            6,
        )


if __name__ == "__main__":
    unittest.main()
