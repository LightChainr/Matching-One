import unittest
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "predictions" / "hexagonal_pell_spin_filter_20260828.yaml"


class HexagonalPellPredictionArtifactTests(unittest.TestCase):
    def test_frozen_targets_and_exponents(self) -> None:
        payload = yaml.safe_load(ARTIFACT.read_text(encoding="utf-8"))
        rows = payload["pell_family"]["targets"]
        self.assertEqual([row["N"] for row in rows], [56, 780, 10864])
        h4 = payload["frozen_hypotheses"]["H4_thermal_level4"]
        h12 = payload["frozen_hypotheses"]["H12_same_radial_alias"]
        self.assertEqual(h4["pell_root_bias_length_exponent"], "6")
        self.assertEqual(h12["pell_root_bias_length_exponent"], "4")
        self.assertEqual(h4["pell_matching_residual_length_exponent"], "21/4")
        self.assertEqual(h12["pell_matching_residual_length_exponent"], "13/4")


if __name__ == "__main__":
    unittest.main()
