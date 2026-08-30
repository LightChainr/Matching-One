from __future__ import annotations

from argparse import Namespace
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from z5_projective_leg_radius6_flat_mc import (  # noqa: E402
    ACQUIRED_ENDPOINTS,
    DEGREE3,
    FIELD_ORDER,
    SOURCE_ENDPOINTS,
    TARGET_ENDPOINTS,
    add,
    exact_gate,
    selected_r2_alexander,
    validate_manifest,
)


class Z5ProjectiveLegRadius6FlatTests(unittest.TestCase):
    def test_minimal_endpoint_sumset_and_fixed_gauge(self) -> None:
        self.assertEqual({add(left, right) for left in DEGREE3 for right in DEGREE3}, set(SOURCE_ENDPOINTS))
        self.assertEqual({selected_r2_alexander(point) for point in SOURCE_ENDPOINTS}, set(TARGET_ENDPOINTS))
        self.assertEqual(len(ACQUIRED_ENDPOINTS), 13)

    def test_exact_gate(self) -> None:
        gate = exact_gate()
        self.assertTrue(gate["passed"])
        self.assertEqual(gate["acquired_spatial_points"], 13)
        self.assertEqual(gate["full_radius6_shell_points"], 24)
        self.assertEqual(gate["complex_coordinates_per_batch"], 28)
        self.assertEqual(len(FIELD_ORDER), 56)

    def test_frozen_manifest_cannot_execute(self) -> None:
        args = Namespace(
            samples=1_200_000,
            batches=400,
            workers=16,
            p=0.59274605079,
            seed=25060610120261250,
            replica_offset=0,
        )
        with self.assertRaisesRegex(ValueError, "not authorized"):
            validate_manifest(args, {"status": "frozen_unexecuted_radius6_flat_extension", "execution_authorized": False})

    def test_repository_authorization_preserves_frozen_geometry(self) -> None:
        manifest = json.loads((ROOT / "analysis/p250_projective_leg_radius6_flat_freeze.json").read_text())
        gate = json.loads((ROOT / "analysis/p250_projective_leg_radius6_flat_exact_gate.json").read_text())
        self.assertTrue(manifest["execution_authorized"])
        self.assertFalse(manifest["do_not_run_from_this_commit"])
        self.assertEqual(manifest["authorization_parent_commit"], "a6d3b7790fc92ee86ac157712905855e8e1f50bf")
        self.assertEqual(manifest["geometry"]["points_by_hand"], gate["hands"] if "hands" in gate else {
            "plus": gate["source_degree6_endpoints"],
            "minus": gate["target_degree6_endpoints"],
        })
        self.assertEqual(manifest["run"]["seed"], 25060610120261250)
        self.assertEqual(manifest["run"]["replica_last_exclusive"], 1_200_000)


if __name__ == "__main__":
    unittest.main()
