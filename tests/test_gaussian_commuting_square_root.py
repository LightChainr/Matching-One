#!/usr/bin/env python3
"""Tests for the exact, data-free Issue 158 commuting-square contract."""

from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import validate_gaussian_commuting_square_root as validator  # noqa: E402


CONTRACT_PATH = ROOT / "analysis" / "gaussian_commuting_square_root_contract.json"
ARITHMETIC_PATH = ROOT / "scripts" / "gaussian_harmonic_arithmetic.py"


class GaussianCommutingSquareRootTests(unittest.TestCase):
    def setUp(self) -> None:
        self.contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
        self.arithmetic_bytes = ARITHMETIC_PATH.read_bytes()

    def validate(self, contract=None):
        return validator.validate_contract(
            self.contract if contract is None else contract,
            self.arithmetic_bytes,
        )

    def test_checked_in_contract_closes_both_squares_exactly(self) -> None:
        result = self.validate()
        self.assertEqual(result["status"], "valid_exact_design_only")
        self.assertEqual(result["common_root_factor"], "7/1250")
        self.assertTrue(all(row["paths_close"] for row in result["lineages"]))
        self.assertFalse(result["contains_target_data"])

    def test_multiplier_or_child_drift_fails(self) -> None:
        changed = copy.deepcopy(self.contract)
        changed["lineages"][0]["multipliers"]["norm5"] = [2, 1]
        with self.assertRaisesRegex(ValueError, "commuting product drift"):
            self.validate(changed)

        changed = copy.deepcopy(self.contract)
        changed["lineages"][1]["canonical_children"]["norm10"][0] = [29, 1]
        with self.assertRaisesRegex(ValueError, "canonical child drift"):
            self.validate(changed)

    def test_character_target_and_source_digest_drift_fail(self) -> None:
        changed = copy.deepcopy(self.contract)
        changed["character_contract"]["norm10_factor"] = "1/100"
        with self.assertRaisesRegex(ValueError, "norm10 factor drift"):
            self.validate(changed)

        changed = copy.deepcopy(self.contract)
        changed["arithmetic_source"]["sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
            self.validate(changed)

    def test_target_data_fields_fail_closed(self) -> None:
        for key in ("samples", "seed", "roots", "covariance", "score"):
            with self.subTest(key=key):
                changed = copy.deepcopy(self.contract)
                changed["lineages"][0][key] = "forbidden"
                with self.assertRaisesRegex(ValueError, "target-data fields"):
                    self.validate(changed)

    def test_scoring_protocol_cannot_move_after_target_read(self) -> None:
        changed = copy.deepcopy(self.contract)
        changed["scoring_protocol"]["freeze_corrections_from_source_data_before_target_read"] = False
        with self.assertRaisesRegex(ValueError, "scoring protocol drift"):
            self.validate(changed)


if __name__ == "__main__":
    unittest.main()