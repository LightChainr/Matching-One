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

import score_c4_tangent_orthogonal_holdout_typed as typed  # noqa: E402


def frozen_result() -> dict:
    return {
        "schema": "matching-one/c4-selfmatching-orthogonal-holdout/v1",
        "status": "prospective_N170_score",
        "primary_channel": "cross",
        "source_N": 130,
        "target_N": 170,
        "point": {
            "source_t": 2.0,
            "source_lambda": 1.0,
            "source_c": 0.5,
            "target_t": 4.0,
            "target_lambda": 2.0,
            "target_c": 0.5,
            "orthogonal_residual": 0.0,
            "thermal_scaling_residual": 1.0,
        },
    }


class TypedC4TangentOrthogonalHoldoutTests(unittest.TestCase):
    def test_registered_size_map_is_exact_identity(self) -> None:
        gate, source, target, transform = typed.load_semantic_gate(ROOT)
        self.assertEqual(source.to_dict(), target.to_dict())
        self.assertEqual((transform.scale, transform.offset), (1.0, 0.0))
        self.assertEqual(gate["response_coordinates_in_order"], ["t", "lambda"])

    def test_typed_entrypoint_preserves_frozen_payload(self) -> None:
        expected = frozen_result()
        with mock.patch.object(typed.frozen_kernel, "render", return_value=copy.deepcopy(expected)):
            result = typed.render_typed(ROOT, Path("source.csv"), Path("target.csv"))
        semantics = result.pop("observable_semantics")
        self.assertEqual(result, expected)
        self.assertEqual(semantics["applied_transform"]["scale"], 1.0)
        self.assertEqual(
            semantics["validation_order"],
            "semantic_map_before_frozen_kernel_score",
        )

    def test_descriptor_drift_fails_before_kernel(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            destination = root / typed.SEMANTIC_MANIFEST
            destination.parent.mkdir(parents=True)
            payload = json.loads(
                (ROOT / typed.SEMANTIC_MANIFEST).read_text(encoding="utf-8")
            )
            payload["target_descriptor"]["channel"] = "either"
            destination.write_text(json.dumps(payload), encoding="utf-8")
            with mock.patch.object(typed.frozen_kernel, "render") as kernel:
                with self.assertRaisesRegex(ValueError, "no exact topology map"):
                    typed.render_typed(root, Path("source.csv"), Path("target.csv"))
                kernel.assert_not_called()

    def test_response_coordinate_drift_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            destination = root / typed.SEMANTIC_MANIFEST
            destination.parent.mkdir(parents=True)
            payload = json.loads(
                (ROOT / typed.SEMANTIC_MANIFEST).read_text(encoding="utf-8")
            )
            payload["response_coordinates_in_order"].reverse()
            destination.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "coordinate order"):
                typed.load_semantic_gate(root)

    def test_result_channel_drift_fails_closed(self) -> None:
        payload = frozen_result()
        payload["primary_channel"] = "either"
        with mock.patch.object(typed.frozen_kernel, "render", return_value=payload):
            with self.assertRaisesRegex(ValueError, "result channel"):
                typed.render_typed(ROOT, Path("source.csv"), Path("target.csv"))


if __name__ == "__main__":
    unittest.main()
