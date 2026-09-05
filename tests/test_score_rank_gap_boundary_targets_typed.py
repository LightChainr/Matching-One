
from __future__ import annotations
from pathlib import Path
import shutil
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_rank_gap_boundary_targets_typed as typed  # noqa: E402


def frozen_payload() -> dict:
    return {
        "status": "frozen_source_fit_targets_unseen",
        "model": "E[G]=A*N^(5/8)+B; exponent fixed, not fitted",
        "target_order": [325, 425],
        "target_gap_mean_prediction": ["15.6", "18.5"],
    }


class RankGapBoundaryTypedTests(unittest.TestCase):
    def copied_root(self) -> tuple[tempfile.TemporaryDirectory, Path]:
        directory = tempfile.TemporaryDirectory()
        root = Path(directory.name)
        destination = root / typed.SEMANTIC_GATE
        destination.parent.mkdir(parents=True)
        shutil.copy(ROOT / typed.SEMANTIC_GATE, destination)
        return directory, root

    def test_gate_registers_exact_identity(self) -> None:
        gate, source, target, transform = typed.load_semantic_gate(ROOT)
        self.assertEqual(gate["paired_quantity"], "K_plus_minus_K_minus")
        self.assertEqual(source, target)
        self.assertEqual((transform.scale, transform.offset), (1.0, 0.0))

    def test_delegation_preserves_frozen_payload(self) -> None:
        frozen = frozen_payload()
        result = typed.score_typed(
            ROOT, Path("manifest"), Path("source"), runner=lambda *_: frozen.copy()
        )
        semantics = result.pop("observable_semantics")
        self.assertEqual(result, frozen)
        self.assertEqual(semantics["units"], "rank")


if __name__ == "__main__":
    unittest.main()
