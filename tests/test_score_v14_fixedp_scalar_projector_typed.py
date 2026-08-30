from __future__ import annotations

import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_v14_fixedp_scalar_projector_typed as typed  # noqa: E402


def row(n: int = 65) -> dict:
    return {
        "N": n,
        "first": [1, 0],
        "second": [1, 1],
        "cos4_first": 1.0,
        "cos4_second": -1.0,
        "delta_cos4": 2.0,
        "M_first": 0.1,
        "M_second": 0.2,
        "M_scalar_H4_null": 0.15,
        "M_scalar_se": 0.01,
        "M_scalar_z": 15.0,
        "within_pair_correlation": 0.0,
        "N25_8_scaled_scalar": 1.0,
        "N25_8_scaled_se": 0.1,
    }


class V14FixedPScalarTypedTests(unittest.TestCase):
    def copied_root(self) -> tuple[tempfile.TemporaryDirectory, Path]:
        directory = tempfile.TemporaryDirectory()
        root = Path(directory.name)
        destination = root / typed.SEMANTIC_GATE
        destination.parent.mkdir(parents=True)
        shutil.copy(ROOT / typed.SEMANTIC_GATE, destination)
        return directory, root

    def test_gate_registers_exact_identity(self) -> None:
        gate, source, target, transform = typed.load_semantic_gate(ROOT)
        self.assertEqual(gate["channel"], "direction_1")
        self.assertEqual(source, target)
        self.assertEqual((transform.scale, transform.offset), (1.0, 0.0))

    def test_payload_matches_frozen_contract_before_annotation(self) -> None:
        result = typed.score_typed(
            ROOT, Path("analysis.csv"), runner=lambda *_: [row()]
        )
        semantics = result.pop("observable_semantics")
        self.assertEqual(result["classification"], "retrospective discovery/power diagnostic")
        self.assertEqual(result["rows"], [row()])
        self.assertEqual(semantics["normalization_power_in_N"], {"numerator": 25, "denominator": 8})

    def test_descriptor_drift_fails_closed(self) -> None:
        directory, root = self.copied_root()
        self.addCleanup(directory.cleanup)
        path = root / typed.SEMANTIC_GATE
        gate = json.loads(path.read_text(encoding="utf-8"))
        gate["target_descriptor"]["channel"] = "direction_0"
        path.write_text(json.dumps(gate), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "no exact topology map"):
            typed.load_semantic_gate(root)

    def test_runtime_channel_drift_fails_before_rows(self) -> None:
        called = []
        with self.assertRaisesRegex(ValueError, "runtime channel"):
            typed.score_typed(
                ROOT, Path("analysis.csv"), channel="direction_0",
                runner=lambda *_: called.append(True),
            )
        self.assertEqual(called, [])

    def test_runtime_p_ref_drift_fails_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "runtime p_ref"):
            typed.score_typed(
                ROOT, Path("analysis.csv"), p_ref=0.6, runner=lambda *_: [row()]
            )

    def test_row_contract_drift_fails_closed(self) -> None:
        bad = row()
        bad.pop("N25_8_scaled_se")
        with self.assertRaisesRegex(ValueError, "row schema"):
            typed.score_typed(
                ROOT, Path("analysis.csv"), runner=lambda *_: [bad]
            )


if __name__ == "__main__":
    unittest.main()
