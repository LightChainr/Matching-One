from __future__ import annotations

import sys
from pathlib import Path
import unittest

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from rho_child_primitive_h4_mc import (  # noqa: E402
    N60_CHILD_DESIGNS,
    child_gate,
    run,
)
from score_rho_child_opposite_pell_phase import ray_score  # noqa: E402


class OppositePellPhaseTests(unittest.TestCase):
    def test_n60_frozen_geometry_and_alias_gate(self):
        self.assertEqual(
            [matrix for _, matrix in N60_CHILD_DESIGNS],
            [((6, 6), (0, 10)), ((12, 3), (0, 5)), ((12, 9), (0, 5))],
        )
        gate = child_gate("N60_Dminus2")
        self.assertEqual(gate["parent_N"], 30)
        self.assertTrue(gate["direction_alias_gate"]["all_rank_two"])
        self.assertTrue(all(row["N"] == 60 and row["bonds"] == 120 for row in gate["children"]))

    def test_n60_stream_counter_contract(self):
        rows, summary = run(
            20, 2, 1, 2672266001,
            family="N60_Dminus2", child="2omega", replica_offset=18000000000,
        )
        self.assertEqual(rows[0]["replica_first"], 18000000000)
        self.assertEqual(summary["primary_order"], ["2omega_re", "2omega_im"])
        self.assertTrue(summary["all_invariant_failures_zero"])

    def test_ray_score_separates_preserve_and_flip(self):
        mp.mp.dps = 50
        reference = mp.matrix([mp.mpf("-5e-5"), mp.mpf("2e-6")])
        covariance = mp.eye(2) * mp.mpf("1e-14")
        preserve = ray_score(reference, covariance, 2 * reference, covariance, +1)
        wrong_flip = ray_score(reference, covariance, 2 * reference, covariance, -1)
        self.assertLess(mp.mpf(preserve["chi_square"]), mp.mpf("1e-20"))
        self.assertGreater(mp.mpf(wrong_flip["chi_square"]), 100)
        flip = ray_score(reference, covariance, -2 * reference, covariance, -1)
        self.assertLess(mp.mpf(flip["chi_square"]), mp.mpf("1e-20"))


if __name__ == "__main__":
    unittest.main()
