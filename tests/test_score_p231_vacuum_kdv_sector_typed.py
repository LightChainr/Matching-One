from __future__ import annotations

import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_p231_vacuum_kdv_sector_typed as typed  # noqa: E402


def frozen_payload() -> dict:
    return {
        "joint_order": typed.JOINT_ORDER,
        "observed": [1.0] * 6,
        "theory_vector_per_unit_g4": [2.0] * 6,
        "covariance": [[1.0 if i == j else 0.0 for j in range(6)] for i in range(6)],
        "governance": {"new_independent_evidence": False},
    }


class P231VacuumKdVSectorTypedTests(unittest.TestCase):
    def copied_root(self) -> tuple[tempfile.TemporaryDirectory, Path]:
        directory = tempfile.TemporaryDirectory()
        root = Path(directory.name)
        destination = root / typed.SEMANTIC_GATE
        destination.parent.mkdir(parents=True)
        shutil.copy(ROOT / typed.SEMANTIC_GATE, destination)
        return directory, root

    def test_gate_registers_exact_identity(self) -> None:
        gate, source, target, transform = typed.load_semantic_gate(ROOT)
        self.assertEqual(gate["coordinates_in_order"], typed.COORDINATES)
        self.assertEqual(source, target)
        self.assertEqual((transform.scale, transform.offset), (1.0, 0.0))

    def test_delegation_preserves_frozen_payload(self) -> None:
        frozen = frozen_payload()
        result = typed.score_typed(
            ROOT, Path("pilot.json"), Path("oracle.json"),
            runner=lambda *_: frozen.copy(),
        )
        semantics = result.pop("observable_semantics")
        self.assertEqual(result, frozen)
        self.assertEqual(semantics["non_scalar_C_only_indices"], [0, 3])

    def test_descriptor_drift_fails_closed(self) -> None:
        directory, root = self.copied_root()
        self.addCleanup(directory.cleanup)
        path = root / typed.SEMANTIC_GATE
        gate = json.loads(path.read_text(encoding="utf-8"))
        gate["target_descriptor"]["channel"] = "cross"
        path.write_text(json.dumps(gate), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "no exact topology map"):
            typed.load_semantic_gate(root)

    def test_coordinate_order_drift_fails_closed(self) -> None:
        directory, root = self.copied_root()
        self.addCleanup(directory.cleanup)
        path = root / typed.SEMANTIC_GATE
        gate = json.loads(path.read_text(encoding="utf-8"))
        gate["coordinates_in_order"].reverse()
        path.write_text(json.dumps(gate), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "coordinate order"):
            typed.load_semantic_gate(root)

    def test_result_vector_width_drift_fails_closed(self) -> None:
        bad = frozen_payload()
        bad["observed"] = [1.0] * 5
        with self.assertRaisesRegex(ValueError, "vector width"):
            typed.score_typed(
                ROOT, Path("pilot.json"), Path("oracle.json"),
                runner=lambda *_: bad,
            )

    def test_evidence_boundary_drift_fails_closed(self) -> None:
        bad = frozen_payload()
        bad["governance"] = {"new_independent_evidence": True}
        with self.assertRaisesRegex(ValueError, "evidence boundary"):
            typed.score_typed(
                ROOT, Path("pilot.json"), Path("oracle.json"),
                runner=lambda *_: bad,
            )


if __name__ == "__main__":
    unittest.main()
