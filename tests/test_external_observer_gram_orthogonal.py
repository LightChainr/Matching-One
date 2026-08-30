from __future__ import annotations

from pathlib import Path
import sys
import unittest

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_external_observer_gram_orthogonal import (  # noqa: E402
    METRICS,
    delete_one_covariance,
    orthogonalize,
)


class ExternalObserverGramOrthogonalTests(unittest.TestCase):
    def test_gram_rule_and_coupling_are_exact(self) -> None:
        values = {
            "Gram_abs_J_S2": mp.mpf(4),
            "Gram_J_D_conj_J_S_re": mp.mpf(2),
            "Gram_J_D_conj_J_S_im": mp.mpf(0),
            "connected_O_far_J_D_re": mp.mpf(3),
            "connected_O_far_J_D_im": mp.mpf(4),
            "connected_O_far_J_S_re": mp.mpf(1),
            "connected_O_far_J_S_im": mp.mpf(2),
        }
        result = orthogonalize(values)
        self.assertEqual(result["beta"], mp.mpf("0.5"))
        self.assertEqual(result["gram_residual"], 0)
        self.assertEqual(result["perp"], mp.mpc("2.5", "3"))

    def test_imaginary_gram_is_audit_not_fit_direction(self) -> None:
        values = {
            "Gram_abs_J_S2": mp.mpf(5),
            "Gram_J_D_conj_J_S_re": mp.mpf(-1),
            "Gram_J_D_conj_J_S_im": mp.mpf("1e-20"),
            "connected_O_far_J_D_re": mp.mpf(1),
            "connected_O_far_J_D_im": mp.mpf(0),
            "connected_O_far_J_S_re": mp.mpf(2),
            "connected_O_far_J_S_im": mp.mpf(0),
        }
        result = orthogonalize(values)
        self.assertEqual(result["beta"], mp.mpf("-0.2"))
        self.assertEqual(mp.im(result["gram_residual"]), mp.mpf("1e-20"))

    def test_delete_one_covariance_uses_local_metric_order(self) -> None:
        rows = [
            {name: mp.mpf(index + offset) for index, name in enumerate(METRICS)}
            for offset in (0, 1, 2)
        ]
        covariance = delete_one_covariance(rows)
        self.assertEqual(len(covariance), len(METRICS))
        self.assertEqual(covariance[0][0], mp.mpf(4) / 3)


if __name__ == "__main__":
    unittest.main()
