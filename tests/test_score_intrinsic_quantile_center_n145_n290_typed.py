from __future__ import annotations

import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_intrinsic_quantile_center_n145_n290_typed as typed  # noqa: E402


def frozen_payload() -> dict:
    return {
        "frozen": {"u": [0.025, 0.05]},
        "observations": {"N145": {}, "N290": {}},
        "size_local_feature_order": typed.FEATURE_ORDER,
        "joint_residual_order": typed.RESIDUAL_ORDER,
        "joint_residual": [0.0, 0.0, 0.0],
    }


class IntrinsicQuantileCenterTypedTests(unittest.TestCase):
    def copied_root(self) -> tuple[tempfile.TemporaryDirectory, Path]:
        directory = tempfile.TemporaryDirectory()
        root = Path(directory.name)
        destination = root / typed.SEMANTIC_GATE
        destination.parent.mkdir(parents=True)
        shutil.copy(ROOT / typed.SEMANTIC_GATE, destination)
        return directory, root

    def test_gate_registers_exact_identity(self) -> None:
        gate, source, target, transform = typed.load_semantic_gate(ROOT)
        self.assertEqual((gate["source_size"], gate["target_size"]), (145, 290))
        self.assertEqual(source, target)
        self.assertEqual((transform.scale, transform.offset), (1.0, 0.0))

    def test_delegation_preserves_frozen_payload(self) -> None:
        frozen = frozen_payload()
        result = typed.score_typed(
            ROOT, *(Path(name) for name in ("ph", "pm", "ch", "cm", "freeze")),
            runner=lambda *_: frozen.copy(),
        )
        semantics = result.pop("observable_semantics")
        self.assertEqual(result, frozen)
        self.assertEqual(semantics["cross_size_covariance"], "zero_by_independent_rng_domains")

    def test_descriptor_drift_fails_closed(self) -> None:
        directory, root = self.copied_root()
        self.addCleanup(directory.cleanup)
        path = root / typed.SEMANTIC_GATE
        gate = json.loads(path.read_text(encoding="utf-8"))
        gate["target_descriptor"]["channel"] = "either"
        path.write_text(json.dumps(gate), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "no exact topology map"):
            typed.load_semantic_gate(root)

    def test_feature_order_drift_fails_closed(self) -> None:
        directory, root = self.copied_root()
        self.addCleanup(directory.cleanup)
        path = root / typed.SEMANTIC_GATE
        gate = json.loads(path.read_text(encoding="utf-8"))
        gate["size_local_feature_order"].reverse()
        path.write_text(json.dumps(gate), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "feature order"):
            typed.load_semantic_gate(root)

    def test_kernel_contract_drift_fails_closed(self) -> None:
        bad = frozen_payload()
        bad["observations"] = {"N290": {}, "N145": {}}
        with self.assertRaisesRegex(ValueError, "size order"):
            typed.score_typed(
                ROOT, *(Path(name) for name in ("ph", "pm", "ch", "cm", "freeze")),
                runner=lambda *_: bad,
            )


if __name__ == "__main__":
    unittest.main()
