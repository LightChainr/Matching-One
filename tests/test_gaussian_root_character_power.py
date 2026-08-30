#!/usr/bin/env python3
"""Tests for the source-only Issue 158 power gate."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import plan_gaussian_root_character_power as planner  # noqa: E402


class GaussianRootCharacterPowerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        p45 = json.loads(planner.DEFAULT_P45.read_text(encoding="utf-8"))
        p57 = json.loads(planner.DEFAULT_P57.read_text(encoding="utf-8"))
        metadata = {
            key: json.loads(path.read_text(encoding="utf-8"))
            for key, path in planner.DEFAULT_METADATA.items()
        }
        cls.result = planner.build_plan(p45, p57, metadata)

    def test_fixed_p_to_root_conversion_is_explicit(self) -> None:
        conversion = self.result["character_conversion"]
        self.assertEqual(conversion["delta_cos4_child_over_parent"], "14/25")
        self.assertEqual(conversion["root_child_over_parent"], "7/1250")
        self.assertAlmostEqual(
            conversion["fixed_p_delta_M_child_over_parent"]
            / conversion["slope_child_over_parent"],
            7 / 1250,
            places=14,
        )

    def test_leading_targets_and_root_linearization_close(self) -> None:
        expected = {
            650: (1.5091265192646443e-5, -7.595490736603176e-7),
            850: (1.0715235119122925e-5, -4.880960975681243e-7),
        }
        for row in self.result["lineages"]:
            delta_m, root = expected[row["target_N"]]
            self.assertAlmostEqual(row["leading_targets"]["delta_M_at_fixed_p"], delta_m)
            self.assertAlmostEqual(row["leading_targets"]["root_gap"], root)
            self.assertLess(row["leading_targets"]["root_linearization_relative_mismatch"], 2e-4)

    def test_no_crn_requirement_is_tens_of_billions(self) -> None:
        by_n = {row["target_N"]: row for row in self.result["lineages"]}
        self.assertAlmostEqual(
            by_n[650]["no_parent_child_crn"]["target_samples_from_delta_M"] / 1e9,
            27.2126228728,
            places=5,
        )
        self.assertAlmostEqual(
            by_n[850]["no_parent_child_crn"]["target_samples_from_delta_M"] / 1e9,
            55.1790543083,
            places=5,
        )
        for row in by_n.values():
            delta = row["no_parent_child_crn"]["target_samples_from_delta_M"]
            root = row["no_parent_child_crn"]["target_samples_from_root_gap"]
            self.assertLess(abs(delta / root - 1), 2e-4)
        self.assertGreater(
            by_n[850]["no_parent_child_crn"]["target_samples_for_relative_standard_error"]["0.1"],
            1e12,
        )
        self.assertEqual(
            by_n[650]["no_parent_child_crn"]["target_samples_for_relative_standard_error"]["0.05"],
            None,
        )

    def test_even_perfect_cover_crn_cannot_materially_help(self) -> None:
        for row in self.result["lineages"]:
            bound = row["unattainable_perfect_parent_child_crn_bound"]
            self.assertLess(bound["maximum_variance_reduction_fraction"], 0.024)
        self.assertEqual(
            self.result["crn_conclusion"]["classification"],
            "power_no_go_for_new_norm10_covering_crn",
        )

    def test_joint_gate_is_not_per_lineage_closure(self) -> None:
        joint = self.result["joint_equal_target_depth"]
        self.assertAlmostEqual(joint["samples_per_lineage"] / 1e9, 17.8779941135, places=5)
        self.assertIn("does not establish", joint["warning"])
        self.assertFalse(self.result["contains_target_data"])


if __name__ == "__main__":
    unittest.main()
