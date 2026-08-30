from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import etop_rank1_elimination as elimination  # noqa: E402


def row(estimate: tuple[str, str], covariance: tuple[tuple[str, str], tuple[str, str]]):
    return {
        "estimate": [mp.mpf(value) for value in estimate],
        "covariance": [[mp.mpf(value) for value in values] for values in covariance],
    }


class ETopRankOneEliminationTests(unittest.TestCase):
    def setUp(self) -> None:
        mp.mp.dps = 80

    def test_exact_common_ray_profiles_to_zero(self) -> None:
        parent = row(("2", "-1"), (("0.04", "0.006"), ("0.006", "0.02")))
        child = row(("1", "-0.5"), (("0.03", "-0.004"), ("-0.004", "0.01")))
        result = elimination.profile_rank_one(parent, child)
        self.assertAlmostEqual(result["lambda"], 0.5, places=12)
        self.assertLess(result["min_chi2"], 1e-25)

    def test_off_ray_has_positive_global_discrepancy(self) -> None:
        parent = row(("2", "0"), (("0.01", "0.002"), ("0.002", "0.02")))
        child = row(("0", "2"), (("0.02", "-0.003"), ("-0.003", "0.01")))
        result = elimination.profile_rank_one(parent, child)
        self.assertGreater(result["min_chi2"], 100)
        self.assertEqual(result["degrees_of_freedom"], 1)
        self.assertLess(
            result["optimizer_certificate"]["max_normalized_root_residual"],
            1e-40,
        )

    def test_production_json_has_four_named_dependency_groups(self) -> None:
        report = elimination.build_report(source="json")
        primary = report["primary"]
        self.assertEqual(len(primary["lineages"]), 4)
        self.assertEqual(
            primary["joint_rank_one"]["dependency_groups"],
            ["P49", "P43", "P50", "P57"],
        )
        self.assertEqual(primary["joint_rank_one"]["decision_at_alpha"], "survives")
        self.assertEqual(
            primary["joint_zero_even_baseline"]["decision_at_alpha"],
            "eliminated",
        )

    def test_checked_artifacts_reproduce_from_checked_crosswalk(self) -> None:
        report = elimination.build_report(source="json")
        checked = json.loads(
            (ROOT / "results/etop-rank1-elimination/latest.json").read_text(
                encoding="utf-8"
            )
        )
        # The committed artifact is generated from raw and therefore differs only
        # in its source-mode certificate, not in any scientific result.
        self.assertEqual(report["primary"], checked["primary"])
        self.assertEqual(
            report["fixed_center_sensitivity"],
            checked["fixed_center_sensitivity"],
        )


if __name__ == "__main__":
    unittest.main()
