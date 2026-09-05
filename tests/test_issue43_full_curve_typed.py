
from __future__ import annotations
import copy
from pathlib import Path
import sys
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_issue43_full_curve_locked_typed as locked_typed  # noqa: E402
import score_issue43_full_curve_typed as typed  # noqa: E402


def frozen_result() -> dict:
    return {
        "protocol": "Issue #43 prospective N=185/265 two-spin4 full-curve score",
        "status": "frozen primary score; no refit",
        "p_ref": typed.frozen.P_REF,
        "scores": {"DeltaM": {}, "DeltaS": {}},
    }


class TypedIssue43FullCurveTests(unittest.TestCase):
    def test_both_sector_maps_are_exact_identity(self) -> None:
        gate, validated = typed.load_semantic_gate(ROOT)
        self.assertEqual(gate["sector_order"], ["DeltaM", "DeltaS"])
        self.assertEqual(
            validated["DeltaM"]["source_descriptor"].combination.value, "odd"
        )
        self.assertEqual(
            validated["DeltaS"]["source_descriptor"].combination.value, "even"
        )
        for values in validated.values():
            self.assertEqual(
                (values["transform"].scale, values["transform"].offset),
                (1.0, 0.0),
            )

    def test_typed_analysis_preserves_frozen_payload(self) -> None:
        expected = frozen_result()
        with mock.patch.object(
            typed.frozen, "analyze", return_value=copy.deepcopy(expected)
        ) as kernel:
            result = typed.analyze_typed(ROOT, {}, Path("predictions.yaml"))
        kernel.assert_called_once_with({}, Path("predictions.yaml"))
        semantics = result.pop("observable_semantics")
        self.assertEqual(result, expected)
        self.assertEqual(semantics["sector_order"], ["DeltaM", "DeltaS"])
        self.assertEqual(
            semantics["validation_order"],
            "semantic_maps_before_frozen_full_curve_score",
        )

    def test_locked_entrypoint_freezes_production_and_activates_validators(self) -> None:
        gate, _ = locked_typed.load_semantic_gate(ROOT)
        self.assertEqual(
            gate["production_lock"]["source_commit"],
            locked_typed.locked.FROZEN_SOURCE_COMMIT,
        )
        original_metadata = locked_typed.typed.frozen.validate_metadata
        original_moments = locked_typed.typed.frozen.validate_moments
        self.addCleanup(
            setattr,
            locked_typed.typed.frozen,
            "validate_metadata",
            original_metadata,
        )
        self.addCleanup(
            setattr,
            locked_typed.typed.frozen,
            "validate_moments",
            original_moments,
        )
        with mock.patch.object(locked_typed.typed, "main", return_value=0) as main:
            self.assertEqual(locked_typed.main(), 0)
        main.assert_called_once_with(operational_entrypoint="production_lock")
        self.assertIs(
            locked_typed.typed.frozen.validate_metadata,
            locked_typed.locked.validate_metadata,
        )
        self.assertIs(
            locked_typed.typed.frozen.validate_moments,
            locked_typed.locked.validate_moments,
        )


if __name__ == "__main__":
    unittest.main()
