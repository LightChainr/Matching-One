from __future__ import annotations

import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import threshold_score_modes_typed as typed  # noqa: E402


def payload(max_order: int = 2) -> dict:
    modes = [
        *(f"P4_S_mode_{order}" for order in range(max_order + 1)),
        *(f"P4_D_mode_{order}" for order in range(max_order + 1)),
    ]
    tower = {
        f"entry_{index}": str(index) for index in range(2 * (max_order + 1))
    }
    return {
        "schema": "matching-one/threshold-krawtchouk-score-modes/v1",
        "coordinate": "eta=log(p/(1-p))",
        "basis": "orthonormal Bin(N,p0) Krawtchouk; H1=(K-Np)/sqrt(Np(1-p))",
        "max_order": max_order,
        "inputs": [],
        "by_N": {
            "65": {
                "N": 65,
                "mode_order": modes,
                "exact_view_identities": {key: "0" for key in typed.EXACT_VIEW_KEYS},
                "parity_tower_scaled": tower,
            }
        },
        "evidence_guard": typed.EVIDENCE_GUARD,
    }


class ThresholdScoreModesTypedTests(unittest.TestCase):
    def copied_root(self) -> tuple[tempfile.TemporaryDirectory, Path]:
        directory = tempfile.TemporaryDirectory()
        root = Path(directory.name)
        destination = root / typed.SEMANTIC_GATE
        destination.parent.mkdir(parents=True)
        shutil.copy(ROOT / typed.SEMANTIC_GATE, destination)
        return directory, root

    def test_gate_registers_exact_sector_identities(self) -> None:
        _, validated = typed.load_semantic_gate(ROOT)
        self.assertEqual(list(validated), ["S", "D"])
        self.assertTrue(all(
            (item["transform"].scale, item["transform"].offset) == (1.0, 0.0)
            for item in validated.values()
        ))

    def test_payload_is_unchanged_before_annotation(self) -> None:
        frozen = payload()
        result = typed.score_typed(
            ROOT, [], max_order=2, runner=lambda *_: frozen
        )
        semantics = result.pop("observable_semantics")
        self.assertIs(result, frozen)
        self.assertEqual(result, payload())
        self.assertEqual(semantics["sector_order"], ["S", "D"])

    def test_descriptor_drift_fails_closed(self) -> None:
        directory, root = self.copied_root()
        self.addCleanup(directory.cleanup)
        path = root / typed.SEMANTIC_GATE
        gate = json.loads(path.read_text(encoding="utf-8"))
        gate["sector_descriptors"]["S"]["channel"] = "direction_0"
        path.write_text(json.dumps(gate), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "descriptor changed"):
            typed.load_semantic_gate(root)

    def test_runtime_order_fails_before_kernel(self) -> None:
        called = []
        with self.assertRaisesRegex(ValueError, "runtime order"):
            typed.score_typed(ROOT, [], max_order=13, runner=lambda *_: called.append(True))
        self.assertEqual(called, [])

    def test_mode_order_drift_fails_closed(self) -> None:
        bad = payload()
        bad["by_N"]["65"]["mode_order"].reverse()
        with self.assertRaisesRegex(ValueError, "mode order"):
            typed.score_typed(ROOT, [], max_order=2, runner=lambda *_: bad)

    def test_evidence_guard_drift_fails_closed(self) -> None:
        bad = payload()
        bad["evidence_guard"] = "independent"
        with self.assertRaisesRegex(ValueError, "evidence guard"):
            typed.score_typed(ROOT, [], max_order=2, runner=lambda *_: bad)


if __name__ == "__main__":
    unittest.main()
