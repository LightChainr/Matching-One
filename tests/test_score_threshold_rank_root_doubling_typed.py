from __future__ import annotations

import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_threshold_rank_root_doubling_typed as typed  # noqa: E402


def frozen_payload() -> dict:
    lineages = [
        {"parent_N": 65, "child_N": 130},
        {"parent_N": 85, "child_N": 170},
    ]
    block = {"target_ratio": -0.25, "lineages": lineages}
    return {
        "full_cross_size_covariance": dict(block),
        "diagonal_cross_size_covariance": dict(block),
    }


class ThresholdRankRootDoublingTypedTests(unittest.TestCase):
    def copied_root(self) -> tuple[tempfile.TemporaryDirectory, Path]:
        directory = tempfile.TemporaryDirectory()
        root = Path(directory.name)
        destination = root / typed.SEMANTIC_GATE
        destination.parent.mkdir(parents=True)
        shutil.copy(ROOT / typed.SEMANTIC_GATE, destination)
        return directory, root

    def test_gate_registers_exact_child_sign_map(self) -> None:
        gate, source, child, transform = typed.load_semantic_gate(ROOT)
        self.assertNotEqual(source.orientation_order, child.orientation_order)
        self.assertEqual((transform.scale, transform.offset), (-1.0, 0.0))
        self.assertEqual(gate["lineages_in_order"], [[65, 130], [85, 170]])

    def test_delegation_preserves_frozen_payload(self) -> None:
        frozen = frozen_payload()
        result = typed.score_typed(ROOT, {}, runner=lambda _: frozen.copy())
        semantics = result.pop("observable_semantics")
        self.assertEqual(result, frozen)
        self.assertEqual(semantics["applied_stored_child_transform"]["scale"], -1.0)

    def test_orientation_drift_fails_closed(self) -> None:
        directory, root = self.copied_root()
        self.addCleanup(directory.cleanup)
        path = root / typed.SEMANTIC_GATE
        gate = json.loads(path.read_text(encoding="utf-8"))
        gate["stored_child_descriptor"]["orientation_order"] = "first_minus_second"
        path.write_text(json.dumps(gate), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "child sign map changed"):
            typed.load_semantic_gate(root)

    def test_genealogy_drift_fails_closed(self) -> None:
        directory, root = self.copied_root()
        self.addCleanup(directory.cleanup)
        path = root / typed.SEMANTIC_GATE
        gate = json.loads(path.read_text(encoding="utf-8"))
        gate["lineages_in_order"].reverse()
        path.write_text(json.dumps(gate), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "lineage order"):
            typed.load_semantic_gate(root)

    def test_result_ratio_drift_fails_closed(self) -> None:
        bad = frozen_payload()
        bad["full_cross_size_covariance"] = {
            "target_ratio": 0.25,
            "lineages": bad["full_cross_size_covariance"]["lineages"],
        }
        with self.assertRaisesRegex(ValueError, "target ratio differs"):
            typed.score_typed(ROOT, {}, runner=lambda _: bad)

    def test_result_lineage_drift_fails_closed(self) -> None:
        bad = frozen_payload()
        bad["diagonal_cross_size_covariance"] = {
            "target_ratio": -0.25,
            "lineages": list(reversed(bad["diagonal_cross_size_covariance"]["lineages"])),
        }
        with self.assertRaisesRegex(ValueError, "lineage order differs"):
            typed.score_typed(ROOT, {}, runner=lambda _: bad)


if __name__ == "__main__":
    unittest.main()
