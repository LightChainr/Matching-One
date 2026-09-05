
from __future__ import annotations
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_p159_pell_hex_filter_typed as typed  # noqa: E402


class P159PellHexFilterTypedTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        gate = json.loads((ROOT / typed.SEMANTIC_GATE).read_text(encoding="utf-8"))
        cls.gate = gate
        cls.batches = ROOT / gate["canonical_inputs"]["batches"]["path"]
        cls.source = ROOT / gate["canonical_inputs"]["source_result"]["path"]
        cls.committed = ROOT / gate["canonical_inputs"]["committed_score"]["path"]

    def copied_root(self) -> tuple[tempfile.TemporaryDirectory, Path]:
        directory = tempfile.TemporaryDirectory()
        root = Path(directory.name)
        destination = root / typed.SEMANTIC_GATE
        destination.parent.mkdir(parents=True)
        shutil.copy(ROOT / typed.SEMANTIC_GATE, destination)
        return directory, root

    def test_gate_registers_envelope_identity_and_freezes_response_basis(self) -> None:
        gate, source, target, transform = typed.load_semantic_gate(ROOT)
        self.assertEqual(source, target)
        self.assertEqual((transform.scale, transform.offset), (1.0, 0.0))
        self.assertEqual(
            gate["primitive_line_contract"]["target_lines_in_order"],
            typed.TARGET_LINES,
        )
        self.assertIn("not registered topology channels", gate["semantic_boundary"])

    def test_canonical_typed_replay_preserves_frozen_score(self) -> None:
        result = typed.score_typed(ROOT, self.batches, self.source)
        semantics = result.pop("observable_semantics")
        frozen = json.loads(self.committed.read_text(encoding="utf-8"))
        self.assertEqual(result, frozen)
        self.assertEqual(
            semantics["typed_scope"],
            "full_configuration_rank_positive_topology_envelope_only",
        )
        self.assertEqual(
            semantics["validation_order"],
            "semantic_gate_and_canonical_inputs_before_frozen_score",
        )

    def test_noncanonical_input_fails_before_runner(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            altered = Path(directory) / "batches.csv"
            altered.write_bytes(self.batches.read_bytes() + b"\n")
            called = False

            def runner(*_args: Path) -> dict:
                nonlocal called
                called = True
                return {}

            with self.assertRaisesRegex(ValueError, "canonical batches git blob"):
                typed.score_typed(ROOT, altered, self.source, runner=runner)
            self.assertFalse(called)


if __name__ == "__main__":
    unittest.main()
