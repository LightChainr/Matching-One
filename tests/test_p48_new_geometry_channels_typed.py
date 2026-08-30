from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_p48_new_geometry_channels_typed as typed  # noqa: E402


def inputs() -> tuple[dict, dict[int, dict]]:
    gate, _ = typed.load_semantic_gate(ROOT)
    canonical = gate["canonical_inputs"]
    return typed.load_canonical_inputs(
        gate,
        ROOT / canonical["source"]["path"],
        ROOT / canonical["185"]["path"],
        ROOT / canonical["265"]["path"],
    )


class TypedP48NewGeometryChannelsTests(unittest.TestCase):
    def test_four_projector_maps_and_powers_are_frozen(self) -> None:
        gate, validated = typed.load_semantic_gate(ROOT)
        self.assertEqual(gate["projector_order"], list(typed.frozen.CHANNELS))
        for projector, values in validated.items():
            self.assertEqual(
                (values["transform"].scale, values["transform"].offset),
                (1.0, 0.0),
                projector,
            )
        self.assertEqual(
            validated["P4_S"]["source_descriptor"].combination.value, "even"
        )
        self.assertEqual(
            validated["P4_D"]["source_descriptor"].combination.value, "odd"
        )
        self.assertEqual(
            validated["P4_S_prime"]["normalization_power_in_N"],
            {"numerator": 5, "denominator": 4},
        )

    def test_typed_score_preserves_every_frozen_number(self) -> None:
        source, targets = inputs()
        expected = typed.frozen.score(copy.deepcopy(source), copy.deepcopy(targets))
        result = typed.score_typed(ROOT, source, targets)
        semantics = result.pop("observable_semantics")
        self.assertEqual(result, expected)
        self.assertEqual(
            semantics["validation_order"],
            "semantic_maps_before_frozen_four_projector_score",
        )
        self.assertEqual(
            semantics["projectors"]["P4_D_prime"]["response_coordinate"],
            "first_p_derivative",
        )

    def test_descriptor_drift_fails_before_frozen_score(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            destination = root / typed.SEMANTIC_MANIFEST
            destination.parent.mkdir(parents=True)
            payload = json.loads(
                (ROOT / typed.SEMANTIC_MANIFEST).read_text(encoding="utf-8")
            )
            payload["projectors"]["P4_S"]["source_descriptor"]["combination"] = "odd"
            destination.write_text(json.dumps(payload), encoding="utf-8")
            with mock.patch.object(typed.frozen, "score") as kernel:
                with self.assertRaisesRegex(ValueError, "P4_S source descriptor"):
                    typed.score_typed(root, {}, {})
                kernel.assert_not_called()

    def test_normalization_power_drift_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            destination = root / typed.SEMANTIC_MANIFEST
            destination.parent.mkdir(parents=True)
            payload = json.loads(
                (ROOT / typed.SEMANTIC_MANIFEST).read_text(encoding="utf-8")
            )
            payload["projectors"]["P4_D"]["normalization_power_in_N"]["numerator"] = 12
            destination.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "P4_D normalization power"):
                typed.load_semantic_gate(root)

    def test_canonical_input_hash_drift_fails_closed(self) -> None:
        gate, _ = typed.load_semantic_gate(ROOT)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "source.json"
            path.write_text("{}\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "canonical source input"):
                typed.load_canonical_inputs(
                    gate,
                    path,
                    ROOT / gate["canonical_inputs"]["185"]["path"],
                    ROOT / gate["canonical_inputs"]["265"]["path"],
                )


if __name__ == "__main__":
    unittest.main()
