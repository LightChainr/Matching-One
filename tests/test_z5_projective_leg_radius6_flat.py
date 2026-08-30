from __future__ import annotations

from argparse import Namespace
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


if __name__ == "__main__":
    unittest.main()
