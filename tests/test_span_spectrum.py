"""Regression lock: the span-spectrum solver against committed exact controls.

The heavy C++ builder is NOT run here. The two committed validation tables
(record layout of scripts/span_spectrum_build.cpp) are solved by
scripts/span_spectrum_solve.py and must reproduce:

1. d_1 = p^w (1-p)^(2w) exactly at NN w=2, p=1/2 (the full-row cluster), and
2. sum_{h<=3} d_h = 9087/1048576 exactly at NN w=4, p=1/2, D_MAX=3 -- the
   section-4.1 finite-height control of the sewing-with-memory note, plus the
   full closure sum_{h<=3} d_h + tail = the certified nu_w(4, 1/2) = 323849/5576960.
"""
from __future__ import annotations
import sys, unittest
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import span_spectrum_solve as S  # noqa: E402

TABLES = ROOT / "results" / "geometric-consistency" / "span-spectrum-20260913" / "validation-tables"


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


if __name__ == "__main__":
    unittest.main()
