from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_norm5_harmonic_primary_typed import (  # noqa: E402
    annotate_output,
    find_output_path,
    load_semantic_gate,
)


class Norm5PrimarySemanticGateTests(unittest.TestCase):
    def test_exact_either_odd_to_cross_odd_map_is_registered(self) -> None:
        manifest, source, target, transform = load_semantic_gate(ROOT)
        self.assertEqual(source.channel.value, "either")
        self.assertEqual(target.channel.value, "cross")
        self.assertEqual(source.combination.value, "odd")
        self.assertEqual(target.combination.value, "odd")
        self.assertEqual(source.normalization.value, "raw")
        self.assertEqual(target.normalization.value, "raw")
        self.assertEqual(transform.scale, 1.0)
        self.assertEqual(transform.offset, 0.0)
        self.assertEqual(manifest["exact_registered_map"]["identity"], "D_either = D_cross")

    def test_wrapper_requires_output_path(self) -> None:
        self.assertEqual(find_output_path(["--output", "score.json"]), Path("score.json"))
        self.assertEqual(find_output_path(["--output=other.json"]), Path("other.json"))
        with self.assertRaises(ValueError):
            find_output_path(["--run", "65:x:y"])

    def test_semantics_are_appended_to_score_output(self) -> None:
        manifest, source, target, transform = load_semantic_gate(ROOT)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "score.json"
            path.write_text('{"schema":"synthetic norm5 score"}\n', encoding="utf-8")
            annotate_output(path, manifest, source, target, transform)
            payload = json.loads(path.read_text(encoding="utf-8"))
        semantics = payload["observable_semantics"]
        self.assertEqual(semantics["source_descriptor"]["channel"], "either")
        self.assertEqual(semantics["target_descriptor"]["channel"], "cross")
        self.assertEqual(semantics["applied_transform"]["scale"], 1.0)
        self.assertEqual(
            semantics["validation_order"], "semantic_map_before_frozen_kernel_score"
        )

    def test_channel_audit_marks_primary_norm5_score_registered(self) -> None:
        payload = json.loads(
            (ROOT / "predictions/wrapping_channel_audit_20260828.json").read_text(
                encoding="utf-8"
            )
        )
        row = next(
            item
            for item in payload["records"]
            if item["id"] == "issue57_norm5_fixed_p_raw_deltaM"
        )
        self.assertEqual(row["status"], "registered")
        self.assertEqual(row["expected_transform"], {"scale": 1.0, "offset": 0.0})


if __name__ == "__main__":
    unittest.main()
