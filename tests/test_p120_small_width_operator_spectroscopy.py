from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p120_small_width_operator_spectroscopy import build_oracle  # noqa: E402


class SmallWidthOperatorSpectroscopyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = build_oracle()

    def test_transfer_and_insertion_symmetries(self) -> None:
        checks = cls_checks = self.result["control"]["symmetry_checks"]
        self.assertTrue(all(cls_checks.values()), checks)

    def test_spatial_sector_is_live_but_aliased(self) -> None:
        spatial = self.result["spatial_sector"]
        self.assertEqual(spatial["combined_sector_ranks"], {"singlet_[4]": "7", "two_row_[2,2]": "12"})
        self.assertIn("cannot distinguish", spatial["conclusion"])

    def test_global_to_v22_matrix_element_is_exact_zero(self) -> None:
        matrix = self.result["exact_matrix_element"]
        self.assertTrue(matrix["is_zero"])
        self.assertEqual(matrix["max_absolute_coefficient_in_full_W4_basis"], "0")

    def test_zero_is_not_an_absence_claim(self) -> None:
        counterexample = self.result["minimal_counterexample_to_absence"]
        self.assertEqual(counterexample["trace_P_singlet_O"], "4")
        self.assertEqual(counterexample["trace_P_[2,2]_O"], "12")

    def test_colour_seam_is_one_versus_zero(self) -> None:
        seam = self.result["colour_seam_one_shot"]
        self.assertEqual(seam["singlet_ratio"], "1")
        self.assertEqual(seam["[2,2]_ratio"], "0")

    def test_committed_artifact_is_reproducible(self) -> None:
        committed = json.loads(
            (ROOT / "results" / "exact-small-width-operator-spectroscopy" / "latest.json").read_text()
        )
        self.assertEqual(committed, self.result)


if __name__ == "__main__":
    unittest.main()
