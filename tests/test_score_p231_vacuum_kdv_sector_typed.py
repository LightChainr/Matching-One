
from __future__ import annotations
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


if __name__ == "__main__":
    unittest.main()
