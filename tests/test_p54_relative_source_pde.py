from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p54_relative_source_pde import build_oracle  # noqa: E402


class RelativeSourcePDETests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = build_oracle()
        cls.tiny = cls.result["tiny_exact_oracle"]

    def test_all_mixed_raw_source_checks_close(self) -> None:
        self.assertTrue(all(row["equal"] for row in self.tiny["mixed_raw_PDE_checks"]))

    def test_logarithmic_closure_is_exact(self) -> None:
        self.assertEqual(
            self.tiny["third_q_cumulant"],
            self.tiny["third_q_cumulant_from_closure"],
        )
        self.assertEqual(self.tiny["third_q_cumulant"]["text"], "555/2048")
        self.assertTrue(self.tiny["first_thermal_log_closure"]["equal"])

    def test_first_two_cumulants_reconstruct_all_sectors(self) -> None:
        self.assertEqual(
            self.tiny["sector_weights"],
            self.tiny["sector_weights_from_mu_variance"],
        )
        self.assertEqual(
            {q: row["text"] for q, row in self.tiny["sector_weights"].items()},
            {"-1": "1/2", "0": "5/16", "1": "3/16"},
        )

    def test_thermal_kappa3_is_explicitly_not_claimed(self) -> None:
        boundary = self.result["claim_boundary"]
        self.assertTrue(any("thermal third-derivative" in row for row in boundary["not_proved"]))
        self.assertIn("F_suuu", boundary["critical_distinction"])

    def test_committed_artifact_is_reproducible(self) -> None:
        committed = json.loads(
            (ROOT / "results" / "exact-relative-source-pde" / "latest.json").read_text()
        )
        self.assertEqual(committed, self.result)


if __name__ == "__main__":
    unittest.main()
