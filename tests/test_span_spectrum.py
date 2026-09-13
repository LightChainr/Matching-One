"""Regression lock: the span-spectrum solver against committed exact controls.

The heavy C++ builder is NOT run here. The two committed validation tables
(record layout of scripts/span_spectrum_build.cpp) are solved by
scripts/span_spectrum_solve.py and must reproduce:

1. d_1 = p^w (1-p)^(2w) exactly at NN w=2, p=1/2 (the full-row cluster), and
2. sum_{h<=3} d_h = 9087/1048576 exactly at NN w=4, p=1/2, D_MAX=3 -- the
   section-4.1 finite-height control of the sewing-with-memory note, plus the
   full closure sum_{h<=3} d_h + tail = the certified nu_w(4, 1/2) = 323849/5576960.

Round 2 (erratum 2026-09-13, PR #739 issuecomment-5653178247) adds three locks:

3. the `light` chain is normalised -- it must agree with the exact solve, not
   return a non-stationary vector (the pre-erratum code raised
   RuntimeError: Factor is exactly singular on w2_nn.bin);
4. the headline moments of the published tables regenerate from the committed
   results JSON (they were previously produced by an uncommitted driver);
5. censoring is explicit -- the 4 p=1/2 configurations whose tail bin is not
   negligible are flagged instead of being read as complete moments.
"""
from __future__ import annotations
import json
import sys, unittest
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import span_spectrum_solve as S  # noqa: E402

TABLES = ROOT / "results" / "geometric-consistency" / "span-spectrum-20260913" / "validation-tables"
RESULTS = ROOT / "results" / "geometric-consistency" / "span-spectrum-20260913.json"


def _runs():
    with open(RESULTS) as fh:
        return json.load(fh)["runs"]


def _find(label, width, p):
    for r in _runs():
        if r["label"] == label and r["width"] == width and r["p"] == p:
            return r
    raise AssertionError(f"run {label} w={width} p={p} missing from the results JSON")


class SpanSpectrumTests(unittest.TestCase):
    def test_full_row_weight_is_exact(self):
        r = S.analyse(str(TABLES / "w2_nn.bin"), 1, 2, nu_ref="7/48")
        self.assertTrue(r["closure_exact"])
        self.assertEqual(Fraction(1, 64), Fraction(r["d_h_float"][0]))
        self.assertEqual(Fraction(7, 48), Fraction(r["sum_dh"]) + Fraction(r["tail_bin"]))

    def test_section_4_1_finite_height_control(self):
        r = S.analyse(str(TABLES / "w4_nn_D3.bin"), 1, 2, nu_ref="323849/5576960")
        self.assertEqual(Fraction("9087/1048576"), Fraction(r["sum_dh_float"]))
        self.assertEqual(Fraction(1, 4096), Fraction(r["d_h_float"][0]))
        self.assertTrue(r["closure_exact"])

    def test_certificate_bound_present(self):
        r = S.analyse(str(TABLES / "w2_nn.bin"), 1, 2)
        self.assertIsNotNone(r["certificate_bound"])

    # --- round 2: the erratum locks -------------------------------------

    def test_light_chain_is_normalised(self):
        """`light` must divide by p_den**w, exactly like the normal branch.

        Pre-erratum this raised RuntimeError: Factor is exactly singular (the
        unscaled (K^T - I) is nonsingular, so the splu path had no stationary
        solution to return).
        """
        ref = S.analyse(str(TABLES / "w2_nn.bin"), 1, 2)
        lit = S.analyse(str(TABLES / "w2_nn.bin"), 1, 2, light=True)
        self.assertEqual("exact-rational", ref["mode"])
        self.assertIn("light", lit["mode"])
        self.assertLess(abs(lit["sum_dh_float"] - ref["sum_dh_float"]), 1e-15)
        self.assertLess(lit["float_residual_l1"], 1e-12)

    def test_published_moments_regenerate(self):
        """The three headline numbers of the tables come from committed data."""
        r = _find("NN p=1/4", 6, "1/4")
        m = S.spectrum_moments(r["d_h_float"], r["tail_bin_float"], r["d_max"])
        self.assertAlmostEqual(m["E_L"], 3.9380183818, places=8)
        self.assertAlmostEqual(m["var_over_E2"], 0.1570569824, places=8)

        r = _find("NN p=1/8", 6, "1/8")
        m = S.spectrum_moments(r["d_h_float"], r["tail_bin_float"], r["d_max"])
        self.assertAlmostEqual(m["E_L"], 2.4488348685, places=8)
        self.assertAlmostEqual(m["var_over_E2"], 0.1707422693, places=8)

    def test_censoring_is_flagged(self):
        """The 4 non-negligible tails are flagged, not silently averaged in."""
        censored = []
        for r in _runs():
            m = S.spectrum_moments(r["d_h_float"], r["tail_bin_float"], r["d_max"])
            if m["censored"]:
                censored.append((r["label"], r["width"], m["tail_fraction"]))
        self.assertEqual(len(censored), 4)
        worst = max(censored, key=lambda t: t[2])
        self.assertEqual(("matching p=1/2", 4), worst[:2])
        self.assertAlmostEqual(worst[2], 2.411e-02, places=5)
        # every subcritical p in the tables is complete to < 1e-6
        for r in _runs():
            if r["p"] not in ("1/2",):
                m = S.spectrum_moments(r["d_h_float"], r["tail_bin_float"], r["d_max"])
                self.assertFalse(m["censored"], r["label"] + " w=%d" % r["width"])

    def test_moment_limit_reading(self):
        """`limit 0` is excluded; the conjectured constant is the better reading.

        Locks the erratum 5.1 table: for all four families the slope of
        (Var/(E L)^2 - (pi/3 - 1))*w is far from the -0.0472 per unit width that
        a decay to zero would require, and model A beats model B on rms.
        """
        sys.path.insert(0, str(ROOT / "scripts"))
        import span_moment_limits as M  # noqa: E402

        self.assertAlmostEqual(M.TARGET, 0.04719755119659763, places=15)
        rows = {r["family"]: r for r in M.report(str(RESULTS), verbose=False)}
        expected_c = {"NN p=1/4": 0.619, "NN p=1/8": 0.631,
                      "matching p=1/8": 0.495, "matching p=1/16": 0.257}
        for fam, want_c in expected_c.items():
            r = rows[fam]
            self.assertTrue(r["excludes_limit_zero"], fam)
            self.assertGreater(r["slope_Rw_w_ge_4"], -0.5 * M.TARGET, fam)
            self.assertLess(r["model_A"]["rms"], r["model_B"]["rms"], fam)
            self.assertAlmostEqual(r["model_A"]["c"], want_c, places=3)


if __name__ == "__main__":
    unittest.main()
