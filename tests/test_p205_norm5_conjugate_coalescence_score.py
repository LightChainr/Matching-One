from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_p205_norm5_conjugate_coalescence import (  # noqa: E402
    MODEL_ORDER,
    jackknife_covariance,
    load_contract,
    orientation_rows_sha256,
    residual_score,
)


class P205Norm5ConjugateCoalescenceScoreTests(unittest.TestCase):
    def setUp(self) -> None:
        mp.mp.dps = 50

    def test_contract_preserves_frozen_order_and_h0_cancellation(self) -> None:
        contract = load_contract(
            ROOT / "predictions/norm5_conjugate_coalescence_20260829.yaml",
            ROOT / "experiments/p205_norm5_conjugate_coalescence_20260829.yaml",
        )
        self.assertEqual(MODEL_ORDER, ("H4", "H12", "H8"))
        self.assertEqual(contract["weights"]["H4"][325], (5, -11, 6))
        self.assertEqual(contract["weights"]["H4"][425], (20, 13, -33))
        for model in MODEL_ORDER:
            for n in (325, 425):
                self.assertEqual(sum(contract["weights"][model][n]), 0)

    def test_exact_common_orientation_digest_is_byte_sensitive(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            first = Path(directory) / "first.csv"
            second = Path(directory) / "second.csv"
            content = (
                "n,a,b,orientation,batch,samples,kind,k,count\n"
                "5,2,1,first,0,2,minus,1,1\n"
                "5,2,1,second,0,2,minus,1,1\n"
            )
            first.write_text(content, encoding="utf-8")
            second.write_text(content, encoding="utf-8")
            self.assertEqual(
                orientation_rows_sha256(first, "first"),
                orientation_rows_sha256(second, "first"),
            )
            second.write_text(content.replace("first,0,2", "first,0,3"), encoding="utf-8")
            self.assertNotEqual(
                orientation_rows_sha256(first, "first"),
                orientation_rows_sha256(second, "first"),
            )

    def test_jackknife_covariance_keeps_cross_geometry_terms(self) -> None:
        covariance = jackknife_covariance([
            [mp.mpf(1), mp.mpf(2), mp.mpf(4)],
            [mp.mpf(2), mp.mpf(4), mp.mpf(8)],
            [mp.mpf(4), mp.mpf(8), mp.mpf(16)],
        ])
        self.assertGreater(covariance[0][1], 0)
        self.assertEqual(covariance[1][2], 8 * covariance[0][0])

    def test_frozen_h4_affine_law_has_zero_residual(self) -> None:
        # 5*C - 11*A + 6*B = 0 with C=(11*A-6*B)/5.
        point = {"C": mp.mpf(1), "A": mp.mpf(2), "B": mp.mpf("17") / 6}
        covariance = [
            [mp.mpf("0.04"), mp.mpf("0.01"), mp.mpf("0.00")],
            [mp.mpf("0.01"), mp.mpf("0.09"), mp.mpf("0.02")],
            [mp.mpf("0.00"), mp.mpf("0.02"), mp.mpf("0.16")],
        ]
        score = residual_score(point, covariance, (5, -11, 6))
        self.assertLess(abs(mp.mpf(score["residual"])), mp.mpf("1e-45"))
        self.assertLess(abs(mp.mpf(score["signed_z"])), mp.mpf("1e-45"))
        self.assertGreater(mp.mpf(score["variance"]), 0)


if __name__ == "__main__":
    unittest.main()
