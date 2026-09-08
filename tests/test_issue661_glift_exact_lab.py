"""#661: checks for the compatible g-lift vs exact-lab DeltaZ comparison."""
from __future__ import annotations

import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "probe_invariant_shape"))

from issue661_glift_exact_lab import (  # noqa: E402
    A_G,
    B_G,
    LEVELS,
    dz_g_lift,
    main,
    q_on_deciles_from_census,
    unoriented_angle_deg,
)


class TestIssue661GliftExactLab(unittest.TestCase):
    def test_g_vectors_agree_between_the_two_freeze_files(self) -> None:
        t582 = json.loads((ROOT / "results/type582-residual/latest.json").read_text())
        p612 = json.loads((ROOT / "results/p612-n725-score/latest.json").read_text())
        self.assertEqual(t582["levels"], list(LEVELS))
        self.assertEqual(
            t582["consensus_g"],
            p612["weightings"]["spin0"][
                "consensus_direction_frozen_from_five_committed_transitions"])

    def test_repaired_census_in_use(self) -> None:
        census = json.loads(
            (ROOT / "results/probe-invariant-shape/census-exact.json").read_text())
        self.assertEqual(census["bond"]["3"]["dual_fail"], 0)

    def test_lift_is_affine_gauge_invariant(self) -> None:
        census = json.loads(
            (ROOT / "results/probe-invariant-shape/census-exact.json").read_text())
        t582 = json.loads((ROOT / "results/type582-residual/latest.json").read_text())
        g = t582["consensus_g"]
        z3 = census["site"]["3"]["Z_levels_float"]
        Q = q_on_deciles_from_census(
            z3, census["site"]["3"]["Q_quarter"], census["site"]["3"]["Q_threequarters"])
        base = dz_g_lift(Q, g)
        scaled = dz_g_lift([3.7 * q - 2.1 for q in Q], g)
        n1 = math.sqrt(sum(x * x for x in base))
        n2 = math.sqrt(sum(x * x for x in scaled))
        cos = sum(a * b for a, b in zip(base, scaled)) / (n1 * n2)
        self.assertAlmostEqual(abs(cos), 1.0, places=12)

    def test_lift_is_zero_at_anchors(self) -> None:
        t582 = json.loads((ROOT / "results/type582-residual/latest.json").read_text())
        g = t582["consensus_g"]
        Q = [0.1 * (i + 1) for i in range(len(LEVELS))]
        lifted = dz_g_lift(Q, g)
        ia, ib = LEVELS.index(A_G), LEVELS.index(B_G)
        self.assertEqual(lifted[ia], 0.0)
        self.assertEqual(lifted[ib], 0.0)

    def test_unoriented_angle_properties(self) -> None:
        a = unoriented_angle_deg([1.0, 0.0], [0.0, 1.0])
        self.assertEqual(a["unoriented_angle_deg"], 90.0)
        b = unoriented_angle_deg([1.0, 1.0], [-1.0, -1.0])
        self.assertAlmostEqual(b["unoriented_angle_deg"], 0.0, places=5)  # same line, opposite orientation
        z = unoriented_angle_deg([1.0], [0.0])
        self.assertTrue(z["undetermined"])

    def test_main_writes_expected_angles(self) -> None:
        main()
        out = json.loads(
            (ROOT / "results/probe-invariant-shape/issue661-glift-exact-lab.json").read_text())
        ang = out["angles"]["spin0"]["DZ_Q3_vs_deltaZ_site"]["angle"]
        self.assertFalse(ang["undetermined"])
        # 0 < angle < 90: the g-lift and the lab tangent are NOT the same line
        self.assertGreater(ang["unoriented_angle_deg"], 0.0)
        self.assertLess(ang["unoriented_angle_deg"], 90.0)
        self.assertAlmostEqual(ang["abs_cosine"], 0.9072881745290168, places=12)
        # equal-weighting sensitivity moves the angle by well under a degree
        self.assertLess(
            abs(out["sensitivity_equal_weighting"]["delta_unoriented_deg_primary"]), 1.0)


if __name__ == "__main__":
    unittest.main()
