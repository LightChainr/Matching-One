from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p250_pair_transfer_power_freeze import freeze  # noqa: E402


class P250PairTransferPowerFreezeTests(unittest.TestCase):
    def test_existing_10k_selects_40k(self) -> None:
        path = ROOT / "results/huawei-20260830/P250-z5-projective-leg-fresh-10k/response_10k.batches.csv"
        result = freeze(path)
        self.assertEqual(result["selected_samples"], 40_000)
        self.assertLess(result["source_resolution"]["3"]["minimum_real_abs_z"], 3.0)
        self.assertTrue(next(row for row in result["grid"] if row["samples"] == 40_000)["qualifies"])


if __name__ == "__main__":
    unittest.main()
