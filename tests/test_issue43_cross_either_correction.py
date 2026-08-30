import csv
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_issue43_cross_either_correction import corrected_score  # noqa: E402


class Issue43CrossEitherCorrectionTests(unittest.TestCase):
    def test_p31_cross_and_either_even_are_opposite(self) -> None:
        path = ROOT / "results/server-20260828/P31/p31_confirmation_seed2026093001.analysis.csv"
        values = {}
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                if row["sector"] == "even" and row["channel"] in ("cross", "either"):
                    values[(int(row["N"]), row["channel"])] = float(row["hypothesis_scaled_amplitude"])
        sizes = sorted({n for n, channel in values if channel == "cross"} & {n for n, channel in values if channel == "either"})
        self.assertGreaterEqual(len(sizes), 5)
        for n in sizes:
            self.assertAlmostEqual(values[(n, "cross")], -values[(n, "either")], places=12)

    def test_current_p43_result_maps_to_small_residual(self) -> None:
        primary_path = ROOT / "results/server-20260828/P43-heldout-fullcurve-500m/analysis/primary_score.json"
        prediction_path = ROOT / "predictions/two_spin4_heldout_20260828.yaml"
        primary = json.loads(primary_path.read_text(encoding="utf-8"))
        result = corrected_score(primary, prediction_path.read_text(encoding="utf-8"))
        self.assertEqual(result["corrected_frozen_cross_mean"][0], -primary["scores"]["DeltaS"]["frozen_mean"][0])
        self.assertEqual(result["corrected_frozen_cross_mean"][1], -primary["scores"]["DeltaS"]["frozen_mean"][1])
        self.assertAlmostEqual(result["chi_square"], 0.5700315435551194, places=10)
        self.assertAlmostEqual(result["marginal_signed_z"][0], 0.66723003, places=6)
        self.assertAlmostEqual(result["marginal_signed_z"][1], -0.11888929, places=6)
        self.assertEqual(result["target_refit_parameters"], 0)
        self.assertTrue(math.isfinite(result["chi_square_survival_df2"]))


if __name__ == "__main__":
    unittest.main()
