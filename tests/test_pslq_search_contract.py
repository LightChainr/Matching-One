#!/usr/bin/env python3
"""Tests for the protocol-only bounded integer-relation search contract."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import validate_pslq_search_contract as contract_validator  # noqa: E402


CONTRACT_PATH = ROOT / "analysis" / "pslq_search_contract.json"
PROVENANCE_PATH = ROOT / "data" / "literature_threshold_sources.json"


class PslqSearchContractTests(unittest.TestCase):
    def setUp(self) -> None:
        self.contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
        self.provenance_bytes = PROVENANCE_PATH.read_bytes()
        self.provenance = json.loads(self.provenance_bytes)

    def validate(self, contract=None):
        return contract_validator.validate_contract(
            self.contract if contract is None else contract,
            self.provenance,
            self.provenance_bytes,
        )

    def test_checked_in_contract_is_protocol_only_and_source_bound(self) -> None:
        result = self.validate()
        self.assertEqual(result["status"], "valid_protocol_only")
        self.assertEqual(result["method_specific_interval_count"], 4)
        self.assertFalse(result["contains_search_results"])
        self.assertEqual(result["provenance_sha256"], hashlib.sha256(self.provenance_bytes).hexdigest())

    def test_preferred_point_or_result_fields_fail_closed(self) -> None:
        for key in ("preferred_point", "combined_interval", "search_results", "near_relations"):
            with self.subTest(key=key):
                changed = copy.deepcopy(self.contract)
                changed[key] = "not allowed"
                with self.assertRaisesRegex(ValueError, "forbidden"):
                    self.validate(changed)

    def test_source_value_and_interval_endpoint_drift_fail(self) -> None:
        changed = copy.deepcopy(self.contract)
        changed["intervals"][0]["central_value"] = "0.59274605079211"
        with self.assertRaisesRegex(ValueError, "central value drift"):
            self.validate(changed)

        changed = copy.deepcopy(self.contract)
        changed["intervals"][0]["lower"] = "0.59274605079209"
        with self.assertRaisesRegex(ValueError, "lower endpoint drift"):
            self.validate(changed)

    def test_library_must_be_frozen_before_search(self) -> None:
        changed = copy.deepcopy(self.contract)
        changed["search_stages"]["standard_constant_pairwise"]["library_frozen_before_search"] = False
        with self.assertRaisesRegex(ValueError, "library is not frozen"):
            self.validate(changed)

    def test_binary_float_and_matching_partner_overclaims_fail(self) -> None:
        changed = copy.deepcopy(self.contract)
        changed["arithmetic"]["binary_float_exclusion_claims_allowed"] = True
        with self.assertRaisesRegex(ValueError, "binary-float"):
            self.validate(changed)

        changed = copy.deepcopy(self.contract)
        changed["false_positive_controls"]["matching_partner_is_independent_evidence"] = True
        with self.assertRaisesRegex(ValueError, "not independent evidence"):
            self.validate(changed)


if __name__ == "__main__":
    unittest.main()
