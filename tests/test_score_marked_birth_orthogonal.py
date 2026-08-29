from __future__ import annotations

from pathlib import Path
import sys
import unittest

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_marked_birth_orthogonal import build_report  # noqa: E402


class MarkedBirthOrthogonalTests(unittest.TestCase):
    def test_exact_contact_and_heldout_protocol(self) -> None:
        mp.mp.dps = 50
        raw = ROOT / "results" / "server-20260829" / "marked-birth-pilot" / "raw"
        report = build_report(
            raw / "q2_N65_20k",
            raw / "q2_N130_20k",
            raw / "P50_N145_10k",
        )
        self.assertEqual(report["training"], "N65 only")
        self.assertEqual(report["held_out"], "N130 q2 child")
        self.assertLess(mp.mpf(report["gram_normal_equation_residual"]["abs"]), mp.mpf("1e-40"))
        for score in report["scores"].values():
            self.assertLess(
                max(mp.mpf(value) for value in score["contact_identity_max_residual"].values()),
                mp.mpf("1e-18"),
            )
        # Frozen discovery result: estimator-level Gram projection leaves the
        # held-out connected thermal growth essentially unchanged.
        raw_ratio = mp.mpf(report["transfers"]["raw_connected_N130_over_conj_N65"]["abs"])
        residual_ratio = mp.mpf(
            report["transfers"]["orthogonal_connected_N130_over_conj_N65"]["abs"]
        )
        self.assertLess(abs(raw_ratio - residual_ratio), mp.mpf("0.001"))


if __name__ == "__main__":
    unittest.main()
