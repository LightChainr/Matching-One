from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from select_z5_charged_radius import CANDIDATE_ORDER, read_batches, select  # noqa: E402


MULTI = ROOT / "results/local-20260830/P250-z5-multiseparation-smoke"
RESULT = ROOT / "results/local-20260830/P250-z5-radius-selector"


def checked_candidates() -> dict:
    output = {}
    for radius in range(1, 5):
        base = MULTI if radius == 1 else RESULT / f"R{radius}"
        output[f"R{radius}"] = {
            "radius": radius,
            "response": json.loads((base / "response_4k.json").read_text()),
            "batches": read_batches(base / "response_4k.batches.csv"),
            "score": json.loads((base / "score_4k.json").read_text()),
        }
    return output


class Z5ChargedRadiusSelectorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.observed = select(checked_candidates())
        cls.checked = json.loads((RESULT / "selection.json").read_text())

    def test_checked_result_is_reproducible(self) -> None:
        self.assertEqual(self.observed, self.checked)

    def test_common_counter_is_exactly_aligned(self) -> None:
        self.assertTrue(self.observed["common_counter_certificate"]["passed"])
        self.assertEqual(self.observed["common_counter_certificate"]["failures"], [])

    def test_no_radius_has_two_usable_separations(self) -> None:
        self.assertEqual(tuple(self.observed["candidates"]), CANDIDATE_ORDER)
        self.assertIsNone(self.observed["selected_candidate"])
        for candidate in self.observed["candidates"].values():
            self.assertEqual(candidate["usable_separation_count"], 1)
            self.assertFalse(candidate["candidate_pass_at_least_two"])

    def test_phase_is_never_a_selection_coordinate(self) -> None:
        for candidate in self.observed["candidates"].values():
            for separation in candidate["separations"].values():
                self.assertFalse(separation["phase_score_used_for_selection"])

    def test_large_R4_d3_cubic_score_cannot_bypass_pair_gate(self) -> None:
        row = self.observed["candidates"]["R4"]["separations"]["3"]
        self.assertLess(row["descriptive_cubic_support_p"], 1e-100)
        self.assertLess(row["minimum_two_point_abs_z"], 2.0)
        self.assertFalse(row["usable_for_production_selector"])

    def test_all_scanned_annuli_are_injective(self) -> None:
        for candidate in self.observed["candidates"].values():
            self.assertTrue(candidate["annulus_injectivity"]["passed"])


if __name__ == "__main__":
    unittest.main()
