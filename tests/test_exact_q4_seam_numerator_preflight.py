#!/usr/bin/env python3

from __future__ import annotations

import importlib.util
import json
from fractions import Fraction
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "exact_q4_seam_numerator_preflight.py"
SPEC = importlib.util.spec_from_file_location("q4_seam_preflight", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class ExactQ4SeamNumeratorPreflightTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = MODULE.render()

    def test_transfer_and_independent_spin_sum_agree(self) -> None:
        self.assertEqual(self.payload["partition_sums"],
                         self.payload["direct_spin_enumeration_partition_sums"])
        self.assertEqual(self.payload["partition_sums"]["identity"]["text"], "33024")
        self.assertEqual(self.payload["partition_sums"]["transposition"]["text"], "21568")

    def test_two_nonzero_witnesses_give_one_vs_zero(self) -> None:
        sectors = self.payload["sector_numerator_oracles"]
        for witness in ("transfer_power_witness", "logv_inserted_witness"):
            singlet = sectors["singlet"][witness]
            charged = sectors["two_row_2"][witness]
            self.assertNotEqual(Fraction(singlet["unnormalized_numerator"]["identity"]["text"]), 0)
            self.assertNotEqual(Fraction(charged["unnormalized_numerator"]["identity"]["text"]), 0)
            self.assertEqual(singlet["unnormalized_twist_to_identity_ratio"]["text"], "1")
            self.assertEqual(charged["unnormalized_twist_to_identity_ratio"]["text"], "0")

    def test_normalized_singlet_ratio_exhibits_and_repairs_denominator_trap(self) -> None:
        row = self.payload["sector_numerator_oracles"]["singlet"]["transfer_power_witness"]
        self.assertEqual(row["normalized_expectation_twist_to_identity_ratio"]["text"], "516/337")
        self.assertEqual(row["partition_restored_ratio"]["text"], "1")

    def test_projectors_commute_and_character_targets_cross_validate_p257(self) -> None:
        self.assertTrue(self.payload["exact_checks"]["all_projectors_commute_with_transfer"])
        self.assertTrue(self.payload["exact_checks"]["both_witnesses_match_frozen_character_targets"])

    def test_geometry_is_frozen_not_silently_generalized(self) -> None:
        with self.assertRaises(ValueError):
            MODULE.render(width=3, height=2)

    def test_committed_artifact_is_reproducible(self) -> None:
        expected = json.loads((ROOT / "results/q4-seam-numerator-preflight/latest.json")
                              .read_text(encoding="utf-8"))
        self.assertEqual(self.payload, expected)


if __name__ == "__main__":
    unittest.main()

