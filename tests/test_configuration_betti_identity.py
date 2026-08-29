from __future__ import annotations

import sys
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from configuration_betti_identity import (  # noqa: E402
    Q_EQUALS_RANK_DIFFERENCE,
    betti_record,
    exhaustive_betti,
    run_betti_suite,
)
from euler_motif_controls import named_tiny_geometries  # noqa: E402
from integer_period_torus import axis_integer_torus  # noqa: E402


class ConfigurationBettiIdentityTests(unittest.TestCase):
    def test_suite_passes_on_named_tiny_geometries(self) -> None:
        payload = run_betti_suite()
        self.assertTrue(payload["passed"])
        self.assertEqual(payload["variance_reduction_branch"], "closed")
        self.assertEqual(payload["claim_level"], "C5")
        names = {row["name"] for row in payload["exhaustive"]}
        self.assertEqual(names, {"axis", "gaussian-2-1", "diamond"})

    def test_empty_axis_l2_is_the_rank_difference_counterexample(self) -> None:
        geometry = axis_integer_torus(2)
        empty = betti_record(geometry, [False] * 4, 0)
        self.assertEqual(empty.chi, 0)
        self.assertEqual(empty.q, -1)
        self.assertEqual(empty.beta0_black, 0)
        self.assertEqual(empty.beta0_white, 1)
        self.assertEqual(empty.r_black, 0)
        self.assertEqual(empty.r_white, 2)
        self.assertEqual(empty.residual, 0)
        self.assertFalse(empty.q_equals_rank_difference)
        self.assertGreaterEqual(empty.kappa_black, empty.r_black)
        self.assertGreaterEqual(empty.kappa_white, empty.r_white)

    def test_full_axis_l2_matches_euler_poincare(self) -> None:
        geometry = axis_integer_torus(2)
        full = betti_record(geometry, [True] * 4, 15)
        self.assertEqual(full.chi, 0)
        self.assertEqual(full.q, 1)
        self.assertEqual(full.beta0_black, 1)
        self.assertEqual(full.beta0_white, 0)
        self.assertEqual(full.r_black, 2)
        self.assertEqual(full.r_white, 0)
        self.assertEqual(full.residual, 0)
        self.assertFalse(full.q_equals_rank_difference)
        self.assertEqual(full.chi, full.beta0_black - full.beta0_white - full.q)

    def test_q_equals_rank_difference_locked_counts(self) -> None:
        for geometry in named_tiny_geometries():
            result = exhaustive_betti(geometry)
            self.assertTrue(result["passed"], result)
            locked = Q_EQUALS_RANK_DIFFERENCE[(geometry.name, geometry.L)]
            self.assertEqual(
                result["q_equals_rank_difference"],
                locked[0],
                geometry.name,
            )
            self.assertEqual(result["configurations"], locked[1])
            self.assertLess(locked[0], locked[1])
            self.assertEqual(result["q_values"], [-1, 0, 1])
            self.assertEqual(result["identity_failures"], 0)
            self.assertEqual(result["cyclomatic_failures"], 0)
            self.assertEqual(result["wrapping_not_identical"], 0)


if __name__ == "__main__":
    unittest.main()
