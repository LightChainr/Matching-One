#!/usr/bin/env python3
"""Tests for the exact, data-free Issue 55 design contract."""

from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import validate_h4_h12_orthogonal_design as validator  # noqa: E402


MANIFEST_PATH = ROOT / "analysis" / "h4_h12_orthogonal_design_manifest.json"
NOTE_PATH = ROOT / "notes" / "h4-h12-orthogonal-gaussian-design.md"


class H4H12OrthogonalDesignTests(unittest.TestCase):
    def setUp(self) -> None:
        self.manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
        self.note_bytes = NOTE_PATH.read_bytes()

    def validate(self, manifest=None):
        return validator.validate_manifest(
            self.manifest if manifest is None else manifest,
            self.note_bytes,
        )

    def test_checked_in_design_recomputes_exactly_without_target_data(self) -> None:
        result = self.validate()
        self.assertEqual(result["status"], "valid_design_only")
        self.assertEqual(result["design_count"], 2)
        self.assertTrue(result["opposite_alias_signs"])
        self.assertFalse(result["contains_target_data"])
        self.assertEqual(result["parent_issue"], "remain open")

    def test_orientation_swap_fails_signed_contrast(self) -> None:
        changed = copy.deepcopy(self.manifest)
        design = changed["designs"][0]
        design["first"], design["second"] = design["second"], design["first"]
        with self.assertRaisesRegex(ValueError, "signed delta_cos4 drift"):
            self.validate(changed)

    def test_exact_fraction_and_decimal_drift_fail(self) -> None:
        changed = copy.deepcopy(self.manifest)
        changed["designs"][1]["delta_cos12"] = "1"
        with self.assertRaisesRegex(ValueError, "signed delta_cos12 drift"):
            self.validate(changed)

        changed = copy.deepcopy(self.manifest)
        changed["designs"][0]["h4_only_target_mean"] = "4.9320781409e-5"
        with self.assertRaisesRegex(ValueError, "target mean decimal drift"):
            self.validate(changed)

        changed = copy.deepcopy(self.manifest)
        changed["designs"][1]["source_coefficient_only_se"] = "2.22590e-6"
        with self.assertRaisesRegex(ValueError, "source-only SE decimal drift"):
            self.validate(changed)

    def test_scoring_order_and_protocol_drift_fail(self) -> None:
        changed = copy.deepcopy(self.manifest)
        changed["scoring_order"][0], changed["scoring_order"][1] = (
            changed["scoring_order"][1],
            changed["scoring_order"][0],
        )
        with self.assertRaisesRegex(ValueError, "scoring order drift"):
            self.validate(changed)

        changed = copy.deepcopy(self.manifest)
        changed["protocol"]["inspect_target_before_sample_count_freeze"] = True
        with self.assertRaisesRegex(ValueError, "must remain false"):
            self.validate(changed)

    def test_target_data_fields_fail_closed(self) -> None:
        for key in ("samples", "seed", "observed_mean", "covariance", "chi2"):
            with self.subTest(key=key):
                changed = copy.deepcopy(self.manifest)
                changed["designs"][0][key] = "forbidden"
                with self.assertRaisesRegex(ValueError, "target-data fields"):
                    self.validate(changed)

    def test_source_note_digest_drift_fails(self) -> None:
        changed = copy.deepcopy(self.manifest)
        changed["source_note"]["sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
            self.validate(changed)


if __name__ == "__main__":
    unittest.main()
