from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_p48_sprime_prospective_typed import (  # noqa: E402
    SEMANTIC_MANIFEST,
    annotate_output,
    find_output_path,
    load_semantic_gate,
)


class P48SPrimeSemanticGateTests(unittest.TestCase):
    def test_declared_map_is_exact_identity(self) -> None:
        manifest, source, target, transform = load_semantic_gate(ROOT)
        self.assertEqual(source, target)
        self.assertEqual((transform.scale, transform.offset), (1.0, 0.0))
        self.assertEqual(
            manifest["canonical_entrypoint"],
            "scripts/score_p48_sprime_prospective_typed.py",
        )
        self.assertEqual(
            manifest["kernel_scorer"], "scripts/score_p48_sprime_prospective.py"
        )

    def test_audit_row_points_to_typed_entrypoint_and_both_artifacts(self) -> None:
        audit = json.loads(
            (ROOT / "predictions/wrapping_channel_audit_20260828.json").read_text(
                encoding="utf-8"
            )
        )
        row = next(record for record in audit["records"] if record["id"] == "p48_cross_even_normalized_identity")
        self.assertEqual(
            row["scorer"], "scripts/score_p48_sprime_prospective_typed.py"
        )
        self.assertEqual(row["semantic_manifest"], SEMANTIC_MANIFEST)
        self.assertEqual(
            row["frozen_model_artifact"],
            "predictions/p48_sprime_correction_20260828.yaml",
        )

    def test_descriptor_drift_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            destination = root / SEMANTIC_MANIFEST
            destination.parent.mkdir(parents=True)
            payload = yaml.safe_load(
                (ROOT / SEMANTIC_MANIFEST).read_text(encoding="utf-8")
            )
            payload["source_descriptor"]["channel"] = "either"
            destination.write_text(yaml.safe_dump(payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "no longer matches"):
                load_semantic_gate(root)

    def test_missing_output_path_fails_before_scoring(self) -> None:
        with self.assertRaisesRegex(ValueError, "requires --output"):
            find_output_path(["--target", "target.json"])
        self.assertEqual(find_output_path(["--output=result.json"]), Path("result.json"))

    def test_annotation_records_complete_contract(self) -> None:
        manifest, source, target, transform = load_semantic_gate(ROOT)
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "score.json"
            output.write_text('{"status":"frozen score"}\n', encoding="utf-8")
            annotate_output(output, manifest, source, target, transform)
            payload = json.loads(output.read_text(encoding="utf-8"))
            semantics = payload["observable_semantics"]
            self.assertEqual(semantics["semantic_manifest"], SEMANTIC_MANIFEST)
            self.assertEqual(semantics["source_descriptor"], semantics["target_descriptor"])
            self.assertEqual(semantics["applied_transform"]["scale"], 1.0)
            self.assertEqual(
                semantics["validation_order"],
                "semantic_map_before_frozen_kernel_score",
            )

    def test_existing_annotation_is_not_overwritten(self) -> None:
        manifest, source, target, transform = load_semantic_gate(ROOT)
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "score.json"
            output.write_text(
                json.dumps({"observable_semantics": {"old": True}}),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "already contains"):
                annotate_output(output, manifest, source, target, transform)


if __name__ == "__main__":
    unittest.main()
