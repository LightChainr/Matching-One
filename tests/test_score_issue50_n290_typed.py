from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_issue50_n290_typed import (  # noqa: E402
    SEMANTIC_MANIFEST,
    load_semantic_gate,
    score_typed,
)


def fixture_rows() -> dict:
    rows = {}
    values = (
        (60, 40, 55, 45),
        (62, 38, 53, 47),
    )
    for batch, (fp, fm, sp, sm) in enumerate(values):
        rows[("either", batch)] = {
            "samples": 100,
            "first_primal": fp,
            "first_matching": fm,
            "second_primal": sp,
            "second_matching": sm,
        }
    return rows


def fixture_run() -> dict:
    return {
        "samples": 200,
        "batches": 2,
        "seed": 7,
        "counter_first": 0,
        "counter_last": 200,
        "p_ref": 0.5,
        "source_commit": "0" * 40,
    }


class TypedIssue50N290Tests(unittest.TestCase):
    def test_registered_map_is_exact_identity(self) -> None:
        _, source, target, transform = load_semantic_gate(ROOT)
        self.assertEqual(source.to_dict(), target.to_dict())
        self.assertEqual((transform.scale, transform.offset), (1.0, 0.0))

    def test_typed_score_preserves_kernel_value_and_records_semantics(self) -> None:
        result = score_typed(ROOT, fixture_rows(), fixture_run())
        self.assertAlmostEqual(result["child_delta_M"], 0.14)
        self.assertEqual(result["channel"], "either")
        self.assertEqual(result["sector"], "matching_function")
        semantics = result["observable_semantics"]
        self.assertEqual(semantics["applied_transform"]["scale"], 1.0)
        self.assertEqual(semantics["validation_order"], "semantic_map_before_frozen_kernel_score")

    def test_reversed_orientation_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / SEMANTIC_MANIFEST
            target.parent.mkdir(parents=True)
            payload = json.loads((ROOT / SEMANTIC_MANIFEST).read_text(encoding="utf-8"))
            payload["target_descriptor"]["orientation_order"] = "second_minus_first"
            target.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "differs from the semantic gate"):
                load_semantic_gate(root)

    def test_lineage_drift_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / SEMANTIC_MANIFEST
            target.parent.mkdir(parents=True)
            payload = json.loads((ROOT / SEMANTIC_MANIFEST).read_text(encoding="utf-8"))
            payload["lineage_first"] = [17, 1]
            target.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "first lineage"):
                load_semantic_gate(root)


if __name__ == "__main__":
    unittest.main()
