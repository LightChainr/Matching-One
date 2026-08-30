from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import p418_anchor_paired_pilot as pilot  # noqa: E402


class P418AnchorPairedPilotTests(unittest.TestCase):
    def test_independent_anchor_stream_is_frozen_and_distinct(self) -> None:
        observed = [pilot.independent_anchor(41850510120260830, replica) for replica in range(5)]
        self.assertEqual(observed, [6, 43, 95, 3, 78])
        current = [pilot.translated_points(25050510120261130, replica)[0] for replica in range(5)]
        self.assertNotEqual(observed, current)

    def test_one_replica_full_is_mean_of_all_anchors(self) -> None:
        values, _, _, _ = pilot.replica_observables(25050510120261130, 41850510120260830, 0)
        for hand in pilot.HANDS:
            for charge in pilot.CHARGES:
                name = f"full__ap0_bp0_r{charge}_{hand}_im"
                self.assertLess(abs(values[name]), 1e-14)


if __name__ == "__main__":
    unittest.main()
