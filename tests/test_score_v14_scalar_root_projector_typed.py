from __future__ import annotations

import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_v14_scalar_root_projector_typed as typed  # noqa: E402


def payload() -> dict:
    return {
        "format_version": 1,
        "hypothesis": {"beta_in_N": 3.5},
        "sizes": {65: {}, 85: {}, 130: {}, 170: {}},
        "lineages": {"65->130": {}, "85->170": {}},
        "two_lineage_consistency": {"pc_difference": 0.0},
    }


class V14ScalarRootTypedTests(unittest.TestCase):
    def copied_root(self) -> tuple[tempfile.TemporaryDirectory, Path]:
        directory = tempfile.TemporaryDirectory()
        root = Path(directory.name)
        destination = root / typed.SEMANTIC_GATE
        destination.parent.mkdir(parents=True)
        shutil.copy(ROOT / typed.SEMANTIC_GATE, destination)
        return directory, root

    def test_gate_registers_exact_identity(self) -> None:
        gate, source, target, transform = typed.load_semantic_gate(ROOT)
        self.assertEqual(gate["response_coordinate"], "implicit_matching_root")
        self.assertEqual(source, target)
        self.assertEqual((transform.scale, transform.offset), (1.0, 0.0))

    def test_payload_is_unchanged_before_annotation(self) -> None:
        frozen = payload()
        result = typed.score_typed(
            ROOT, [Path("histograms.csv")], runner=lambda *_: frozen
        )
        semantics = result.pop("observable_semantics")
        self.assertIs(result, frozen)
        self.assertEqual(result, payload())
        self.assertEqual(semantics["fixed_beta_in_N"], {"numerator": 7, "denominator": 2})

    def test_descriptor_drift_fails_closed(self) -> None:
        directory, root = self.copied_root()
        self.addCleanup(directory.cleanup)
        path = root / typed.SEMANTIC_GATE
        gate = json.loads(path.read_text(encoding="utf-8"))
        gate["target_descriptor"]["channel"] = "direction_0"
        path.write_text(json.dumps(gate), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "no exact topology map"):
            typed.load_semantic_gate(root)

    def test_runtime_beta_drift_fails_before_kernel(self) -> None:
        called = []
        with self.assertRaisesRegex(ValueError, "runtime beta"):
            typed.score_typed(
                ROOT, [], beta=3.0, runner=lambda *_: called.append(True)
            )
        self.assertEqual(called, [])

    def test_runtime_lineage_drift_fails_before_kernel(self) -> None:
        called = []
        with self.assertRaisesRegex(ValueError, "runtime lineages"):
            typed.score_typed(
                ROOT, [], lineages=((65, 130),), runner=lambda *_: called.append(True)
            )
        self.assertEqual(called, [])

    def test_frozen_payload_drift_fails_closed(self) -> None:
        bad = payload()
        bad["two_lineage_consistency"] = None
        with self.assertRaisesRegex(ValueError, "consistency contract"):
            typed.score_typed(ROOT, [], runner=lambda *_: bad)


if __name__ == "__main__":
    unittest.main()
