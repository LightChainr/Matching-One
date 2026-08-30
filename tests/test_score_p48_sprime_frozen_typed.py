from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_p48_sprime_frozen as frozen_kernel  # noqa: E402
from score_p48_sprime_frozen_typed import (  # noqa: E402
    SEMANTIC_MANIFEST,
    load_semantic_gate,
    score_typed,
)


def fixture_target() -> dict:
    return {
        "sizes": [185, 265],
        "independent_of_retrospective_source": True,
        "P4_S_prime": [0.002, 0.0015],
        "covariance_P4_S_prime": [[1e-7, 0.0], [0.0, 1e-7]],
    }


def fixture_manifest() -> dict:
    return {
        "status": "frozen_prospective_scoring_manifest",
        "target_sizes_N": [185, 265],
        "leading_power_in_N": 1.25,
        "models_in_scoring_order": [
            {
                "name": "pure_N^-5/4",
                "basis": "constant",
                "parameters": [1.0],
                "parameter_covariance": [[0.01]],
            },
            {
                "name": "zero_effect",
                "basis": "zero",
                "parameters": [],
                "parameter_covariance": [],
            },
            {
                "name": "q2_even_scalar",
                "basis": "constant_plus_inverse_N",
                "parameters": [1.0, 0.5],
                "parameter_covariance": [[0.01, 0.0], [0.0, 0.01]],
            },
            {
                "name": "rank2_jordan_log",
                "basis": "constant_plus_log_N",
                "parameters": [1.0, 0.05],
                "parameter_covariance": [[0.01, 0.0], [0.0, 0.001]],
            },
        ],
    }


class TypedFrozenP48SPrimeTests(unittest.TestCase):
    def test_registered_descriptor_is_exact_identity(self) -> None:
        _, source, target, transform = load_semantic_gate(ROOT)
        self.assertEqual(source.to_dict(), target.to_dict())
        self.assertEqual((transform.scale, transform.offset), (1.0, 0.0))

    def test_typed_score_preserves_frozen_numerics(self) -> None:
        target = fixture_target()
        manifest = fixture_manifest()
        expected = frozen_kernel.score(target, manifest)
        result = score_typed(ROOT, target, manifest)
        semantics = result.pop("observable_semantics")
        self.assertEqual(result, expected)
        self.assertEqual(
            semantics["validation_order"],
            "semantic_map_before_frozen_kernel_score",
        )
        self.assertEqual(semantics["applied_transform"]["scale"], 1.0)

    def test_descriptor_sector_drift_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / SEMANTIC_MANIFEST
            target.parent.mkdir(parents=True)
            payload = json.loads((ROOT / SEMANTIC_MANIFEST).read_text(encoding="utf-8"))
            payload["target_descriptor"]["combination"] = "odd"
            target.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaises(ValueError):
                load_semantic_gate(root)

    def test_model_order_drift_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / SEMANTIC_MANIFEST
            target.parent.mkdir(parents=True)
            payload = json.loads((ROOT / SEMANTIC_MANIFEST).read_text(encoding="utf-8"))
            payload["models_in_scoring_order"].reverse()
            target.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "model order"):
                load_semantic_gate(root)

    def test_target_independence_gate_is_preserved(self) -> None:
        target = copy.deepcopy(fixture_target())
        target["independent_of_retrospective_source"] = False
        with self.assertRaisesRegex(ValueError, "independent"):
            score_typed(ROOT, target, fixture_manifest())


if __name__ == "__main__":
    unittest.main()
