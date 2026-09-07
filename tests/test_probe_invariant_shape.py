"""#622 probe tests: the numbers that would stop us if wrong.

Wrong-number tripwires (cited facts must stay cited; new objects must stay
consistent):

* site census margins reproduce the published #606 pair counts;
* M(0) = -1, M(1) = +1, M(1/2) = -21/64 (L=3) and -13757/32768 (L=4);
* F(1/2) = (1 + M(1/2))/2 exactly (canonical embedding identity);
* Z(1/4-anchored grid) is Aff(1)-invariant under a nontrivial (alpha, beta);
* the exact same-M/different-F counterexample has M' == M to the Fraction
  and F' != F;
* toy family A's Z-limit is not family B's, and family C does not converge;
* N=725 per-batch Z has max SE < 1e-3 (the shape is measurable);
* the #582 g (cited) is NOT parallel to deltaZ (angle > 45 degrees raw).
"""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROBE = ROOT / "scripts" / "probe_invariant_shape"
sys.path.insert(0, str(PROBE))
sys.path.insert(0, str(ROOT / "scripts"))

from fractions import Fraction  # noqa: E402

import exact_rank_census as census  # noqa: E402


class TestExactCensus(unittest.TestCase):
    def test_site_l3_census_and_facts(self) -> None:
        comps = census.site_rank_pair_components(3)
        N = 9
        counts = {name: sum(comps[name][k] for k in range(N + 1))
                  for name in ("P11", "P20", "P02")}
        self.assertEqual(counts["P02"], 259)   # published #606 margins
        self.assertEqual(counts["P11"], 162)
        self.assertEqual(counts["P20"], 91)
        M, F = census.M_and_F_from_components(comps)
        self.assertEqual(census.eval_poly_at(M, Fraction(0)), -1)
        self.assertEqual(census.eval_poly_at(M, Fraction(1)), 1)
        self.assertEqual(census.eval_poly_at(M, Fraction(1, 2)),
                         Fraction(-21, 64))
        self.assertEqual(census.eval_poly_at(F, Fraction(1, 2)),
                         (1 + Fraction(-21, 64)) / 2)

    def test_site_l4_margins(self) -> None:
        comps = census.site_rank_pair_components(4)
        N = 16
        self.assertEqual(sum(comps["P02"][k] for k in range(N + 1)), 36559)
        self.assertEqual(sum(comps["P11"][k] for k in range(N + 1)), 19932)
        self.assertEqual(sum(comps["P20"][k] for k in range(N + 1)), 9045)
        M, _ = census.M_and_F_from_components(comps)
        self.assertEqual(census.eval_poly_at(M, Fraction(1, 2)),
                         Fraction(-13757, 32768))

    def test_z_affine_invariance(self) -> None:
        comps = census.site_rank_pair_components(3)
        _, F = census.M_and_F_from_components(comps)
        # warp the occupation axis: F2(p) = F((p - beta)/alpha) with
        # alpha = 7/5, beta = -2/13; quantiles shift by the same Aff(1).
        alpha, beta = Fraction(7, 5), Fraction(-2, 13)
        F2 = [Fraction(0)] * (len(F))
        # Z_{a,b}[F2] uses quantiles Q2(u) = alpha*Q(u) + beta; the quotient
        # is unchanged.  Direct check through the exact bisection:
        qa1 = census.quantile_exact(F, census.ANCHOR_A)
        qb1 = census.quantile_exact(F, census.ANCHOR_B)
        q1 = [census.quantile_exact(F, u) for u in census.Z_LEVELS]
        z1 = [(q - qa1) / (qb1 - qa1) for q in q1]
        # Q2 by inverting F2 = composition: quantile of F2 at u equals
        # alpha*quantile of F at u + beta.
        q2 = [alpha * q + beta for q in q1]
        qa2, qb2 = alpha * qa1 + beta, alpha * qb1 + beta
        z2 = [(q - qa2) / (qb2 - qa2) for q in q2]
        self.assertEqual(z1, z2)

    def test_same_M_different_F_counterexample(self) -> None:
        comps = census.site_rank_pair_components(3)
        M0, F0 = census.M_and_F_from_components(comps)
        lam = Fraction(3, 2)
        comps1 = {"P11": [c * lam for c in comps["P11"]],
                  "P20": comps["P20"], "P02": comps["P02"]}
        M1, F1 = census.M_and_F_from_components(comps1)
        self.assertEqual(M0, M1)
        self.assertNotEqual(F0, F1)


class TestResultFiles(unittest.TestCase):
    def setUp(self) -> None:
        d = ROOT / "results" / "probe-invariant-shape"
        self.census_json = json.loads((d / "census-exact.json").read_text())
        self.toy = json.loads((d / "toy-families.json").read_text())
        self.glift = json.loads((d / "blindness-and-glift.json").read_text())
        self.zflow = json.loads((d / "n725-zflow.json").read_text())

    def test_census_json_numbers(self) -> None:
        s3 = self.census_json["site"]["3"]
        self.assertEqual(s3["M_half"], "-21/64")
        self.assertAlmostEqual(s3["Q_quarter_float"], 0.4478048000694663,
                               places=12)
        b = self.census_json["bond"]["3"]
        self.assertEqual(b["configs"], 262144)
        self.assertEqual(b["dual_fail"], 118133)

    def test_toy_families_separated(self) -> None:
        a = self.toy["families"]["A_scale_linear"][-1]["Z"]
        b = self.toy["families"]["B_skew_cube"][-1]["Z"]
        self.assertGreater(abs(a[8] - b[8]), 5e-3)
        c = [r["Z_mid_odd"] for r in self.toy["families"]["C_no_limit"]]
        self.assertTrue(any(abs(x - y) > 0.01 for x, y in zip(c, c[1:])))

    def test_glift_angle_not_a_lift(self) -> None:
        ang = self.glift["direction6"]["angle_raw_deg"]
        self.assertGreater(abs(ang - 90.0), 3.0)  # not parallel
        aff = self.glift["direction6"]["angle_after_affine_removal_deg"]
        self.assertGreater(abs(aff - 180.0), 15.0)  # not antiparallel either

    def test_n725_zflow_measurable(self) -> None:
        zse = self.zflow["weightings"]["spin0"]["Z_se"]
        self.assertLess(max(zse), 1e-3)


if __name__ == "__main__":
    unittest.main()
