from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from gaussian_crt_mixed_join_phaseb import (  # noqa: E402
    incidence_cycle_rank,
    isolated_fiber_join,
    local_color_record,
    render,
)


class GaussianCrtMixedJoinPhaseBTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = render()
        cls.frozen = json.loads(
            (ROOT / "predictions" / "p200_n650_mixed_join_phaseB_20260829.json").read_text(
                encoding="utf-8"
            )
        )

    def test_all_isolated_fibers_match_incidence_cycle_rank(self) -> None:
        for bits in range(1 << 10):
            selected = [bool((bits >> index) & 1) for index in range(10)]
            self.assertEqual(isolated_fiber_join(selected), incidence_cycle_rank(selected))

    def test_two_row_closed_form_for_both_colors(self) -> None:
        for bits in range(1 << 10):
            black, white, odd, even, _ = local_color_record(bits)
            self.assertEqual(odd, black - white)
            self.assertEqual(even, black + white)

    def test_exact_half_normalization(self) -> None:
        half = self.payload["toy_exact_normalization"]["p_half"]
        self.assertEqual(half["mean_exact"], ["499/1024", "499/1024", "0", "499/512"])
        self.assertEqual(half["covariance_exact"][2][2], "681/512")
        self.assertEqual(
            self.payload["toy_exact_normalization"]["p_half_C_odd_distribution_counts"],
            {"-4": 1, "-3": 15, "-2": 80, "-1": 210, "0": 412, "1": 210, "2": 80, "3": 15, "4": 1},
        )

    def test_additive_null_and_connected_witness(self) -> None:
        null = self.payload["null_and_alternative"]
        self.assertEqual(null["exact_algebra"]["mixed_defect"], 0)
        self.assertEqual(null["synthetic_interaction"]["recovered_mixed_defect"], 4)
        witness = self.payload["connected_toy_witness"]
        self.assertEqual((witness["J_isolated"], witness["J_full"], witness["R_connected"]), (0, 1, 1))

    def test_frozen_protocol_matches(self) -> None:
        self.assertEqual(self.payload, self.frozen)


if __name__ == "__main__":
    unittest.main()
