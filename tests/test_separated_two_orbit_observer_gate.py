from __future__ import annotations

import sys
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from separated_two_orbit_observer_gate import build_artifact  # noqa: E402


class SeparatedTwoOrbitObserverGateTests(unittest.TestCase):
    def test_production_mapping_and_rank_gate(self) -> None:
        artifact = build_artifact(6)
        self.assertEqual(artifact["response_determinant"], -2)
        self.assertEqual(artifact["response_rank"], 2)
        for gate in artifact["geometries"].values():
            self.assertEqual(gate["axis_orbit_size"], 4)
            self.assertEqual(gate["diagonal_orbit_size"], 4)
            self.assertTrue(gate["cross_orbit_disjoint"])
            self.assertTrue(gate["local_ring_injective"])
            self.assertTrue(gate["source_root_excluded_from_all_rings"])
            self.assertGreaterEqual(gate["minimum_axis_quotient_distance"], 6)


if __name__ == "__main__":
    unittest.main()
