from __future__ import annotations

import csv
import json
from pathlib import Path
import sys
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from audit_issue43_even_channel_map import _p31_even_rows, audit  # noqa: E402


class Issue43EvenChannelMapTests(unittest.TestCase):
    def test_p31_cross_and_either_even_are_equal_and_opposite(self) -> None:
        rows = _p31_even_rows(
            ROOT / "results/server-20260828/P31/p31_confirmation_seed2026093001.analysis.csv"
        )
        for n in (65, 85, 130, 145, 170):
            cross = rows[n]["cross"]
            either = rows[n]["either"]
            self.assertAlmostEqual(
                float(cross["hypothesis_scaled_amplitude"]),
                -float(either["hypothesis_scaled_amplitude"]),
                places=14,
            )
            self.assertAlmostEqual(
                float(cross["hypothesis_scaled_batch_se"]),
                float(either["hypothesis_scaled_batch_se"]),
                places=14,
            )

    def test_protocol_repair_recomputes_committed_diagnostic(self) -> None:
        with (
            ROOT / "predictions/two_spin4_heldout_20260828.yaml"
        ).open(encoding="utf-8") as handle:
            prediction = yaml.safe_load(handle)
        source_rows = _p31_even_rows(
            ROOT / "results/server-20260828/P31/p31_confirmation_seed2026093001.analysis.csv"
        )
        primary = json.loads(
            (
                ROOT
                / "results/server-20260828/P43-heldout-fullcurve-500m/analysis/primary_score.json"
            ).read_text(encoding="utf-8")
        )
        actual = audit(prediction, source_rows, primary)
        committed = json.loads(
            (
                ROOT
                / "results/server-20260828/P43-heldout-fullcurve-500m/analysis/channel_map_corrected_DeltaS.json"
            ).read_text(encoding="utf-8")
        )

        self.assertEqual(actual["target_refit_parameters"], 0)
        self.assertEqual(actual["source_channel"], "either/even")
        self.assertEqual(actual["target_channel"], "cross/even")
        self.assertAlmostEqual(actual["chi_square"], 0.5700315435551194, places=12)
        self.assertAlmostEqual(actual["chi_square"], committed["chi_square"], places=12)
        for left, right in zip(actual["corrected_frozen_cross_mean"], committed["corrected_frozen_cross_mean"]):
            self.assertAlmostEqual(left, right, places=18)


if __name__ == "__main__":
    unittest.main()
