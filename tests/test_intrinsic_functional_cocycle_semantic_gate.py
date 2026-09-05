
from __future__ import annotations
import json
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_intrinsic_functional_cocycle_typed import (  # noqa: E402
    annotate_output,
    find_json_path,
    load_semantic_gate,
)


class IntrinsicFunctionalSemanticGateTests(unittest.TestCase):
    def test_all_primitive_maps_are_identity_cross_channel_maps(self) -> None:
        manifest, descriptors, transforms = load_semantic_gate(ROOT)
        self.assertEqual(
            manifest["status"], "semantic_gate_added_before_norm5_target_reveal"
        )
        self.assertEqual(
            set(descriptors),
            {"Mbar_center", "P4_S", "P4_D", "P4_S_prime", "P4_D_prime"},
        )
        for name, descriptor in descriptors.items():
            with self.subTest(name=name):
                self.assertEqual(descriptor.channel.value, "cross")
                self.assertEqual(transforms[name].scale, 1.0)
                self.assertEqual(transforms[name].offset, 0.0)
        self.assertEqual(descriptors["P4_S"].combination.value, "even")
        self.assertEqual(descriptors["P4_D"].combination.value, "odd")
        self.assertEqual(descriptors["P4_S"].normalization.value, "angular_normalized")
        self.assertEqual(descriptors["P4_D"].normalization.value, "angular_normalized")

    def test_wrapper_requires_json_output(self) -> None:
        self.assertEqual(find_json_path(["--json", "score.json"]), Path("score.json"))
        self.assertEqual(find_json_path(["--json=other.json"]), Path("other.json"))
        with self.assertRaises(ValueError):
            find_json_path(["--histograms", "x"])

    def test_semantics_are_appended_to_score_output(self) -> None:
        manifest, descriptors, transforms = load_semantic_gate(ROOT)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "score.json"
            path.write_text('{"format_version":1}\n', encoding="utf-8")
            annotate_output(path, manifest, descriptors, transforms)
            payload = json.loads(path.read_text(encoding="utf-8"))
        semantics = payload["observable_semantics"]
        self.assertEqual(semantics["primitive_descriptors"]["P4_S"]["channel"], "cross")
        self.assertEqual(
            semantics["primitive_descriptors"]["P4_D"]["normalization"],
            "angular_normalized",
        )
        self.assertEqual(
            semantics["validation_order"],
            "semantic_identity_maps_before_frozen_kernel_score",
        )

    def test_channel_audit_marks_fullcurve_norm5_score_registered(self) -> None:
        payload = json.loads(
            (ROOT / "predictions/wrapping_channel_audit_20260828.json").read_text(
                encoding="utf-8"
            )
        )
        row = next(
            item
            for item in payload["records"]
            if item["id"] == "issue57_norm5_full_curve_projectors"
        )
        self.assertEqual(row["status"], "registered")
        self.assertEqual(row["scorer"], "scripts/score_intrinsic_functional_cocycle_typed.py")


if __name__ == "__main__":
    unittest.main()
