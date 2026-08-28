from __future__ import annotations

import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from audit_p43_even_channel_mismatch import P31_DEFAULT, P43_DEFAULT, audit  # noqa: E402


class P43EvenChannelAuditTests(unittest.TestCase):
    def test_wrong_channel_is_detected_and_cross_even_closes_posthoc(self) -> None:
        payload = audit(P31_DEFAULT, P43_DEFAULT)
        cross = payload["p31_common_amplitudes"]["cross_even"]
        either = payload["p31_common_amplitudes"]["either_even"]
        self.assertAlmostEqual(cross["mean"], -0.010603216462677733, places=14)
        self.assertAlmostEqual(either["mean"], +0.010603216462677735, places=14)
        self.assertAlmostEqual(cross["se"], 0.0009366870182463298, places=15)
        self.assertAlmostEqual(cross["chi_square"], 4.658014420830224, places=10)
        self.assertLess(
            payload["channel_relation"]["max_abs_cross_plus_either_scaled_amplitude"],
            2e-15,
        )

        corrected = payload["p43_posthoc_cross_even_score"]
        self.assertAlmostEqual(corrected["chi_square"], 0.5700315435551194, places=10)
        self.assertEqual(corrected["df"], 2)
        self.assertAlmostEqual(corrected["marginal_z"][0], 0.66723003, places=6)
        self.assertAlmostEqual(corrected["marginal_z"][1], -0.11888929, places=6)
        self.assertTrue(
            payload["governance"]["original_issue43_even_score_remains_failed"]
        )


if __name__ == "__main__":
    unittest.main()
