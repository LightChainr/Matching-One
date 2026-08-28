#!/usr/bin/env python3
"""Tests for the frozen P31/P32 radial model challenge."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "challenge_orientation_radial_models.py"
SIZES = (65, 85, 130, 145, 170)
DELTA_COS4 = {65: 1.363, 85: 1.594, 130: 1.363, 145: 1.918, 170: 1.594}
DELTA_COS8_OVER_COS4 = {
    65: 0.788639, 85: -0.623945, 130: -0.788639,
    145: -0.054602, 170: 0.623945,
}


class OrientationRadialChallengeTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self) -> None:
        self.temp.cleanup()

    def write_data(self, path: Path, sizes, heldout_shift: float = 0.0) -> None:
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "row_id", "N", "seed", "delta_M", "delta_M_se",
                    "delta_cos4", "delta_cos8",
                ],
                lineterminator="\n",
            )
            writer.writeheader()
            for n in sizes:
                for seed_index, seed in enumerate(("seed-a", "seed-b")):
                    amplitude = 0.65 + 0.4 * n ** -0.5
                    signal = DELTA_COS4[n] * n ** (-13.0 / 8.0) * amplitude
                    noise = (-0.35 if seed_index == 0 else 0.35) * 2e-5
                    if n in (145, 170):
                        signal += heldout_shift
                    writer.writerow({
                        "row_id": "{}:{}".format(n, seed),
                        "N": n,
                        "seed": seed,
                        "delta_M": signal + noise,
                        "delta_M_se": 2e-5,
                        "delta_cos4": DELTA_COS4[n],
                        "delta_cos8": (
                            DELTA_COS4[n] * DELTA_COS8_OVER_COS4[n]
                        ),
                    })

    def write_covariance(self, path: Path) -> None:
        variance = (2e-5) ** 2
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=["row_id_i", "row_id_j", "covariance"],
                lineterminator="\n",
            )
            writer.writeheader()
            for seed in ("seed-a", "seed-b"):
                ids = ["{}:{}".format(n, seed) for n in SIZES]
                for i, first in enumerate(ids):
                    for second in ids[i + 1:]:
                        writer.writerow({
                            "row_id_i": first,
                            "row_id_j": second,
                            "covariance": 0.05 * variance,
                        })

    def run_challenge(self, data: Path, output: Path, covariance: Path = None):
        command = [sys.executable, str(SCRIPT), str(data), "--output-dir", str(output)]
        if covariance is not None:
            command += ["--covariance", str(covariance)]
        return subprocess.run(command, text=True, capture_output=True, check=True)

    def test_missing_frozen_sizes_is_explicitly_not_ready(self) -> None:
        data = self.root / "partial.csv"
        output = self.root / "partial-output"
        self.write_data(data, (65, 85, 145))
        completed = self.run_challenge(data, output)
        self.assertIn("NOT_READY", completed.stdout)
        payload = json.loads((output / "challenge.json").read_text(encoding="utf-8"))
        self.assertEqual(payload["status"], "NOT_READY")
        self.assertEqual(payload["missing_sizes"], [130, 170])
        self.assertFalse((output / "heldout_signed_errors.csv").exists())

    def test_ready_challenge_is_covariance_aware_and_heldout_blind(self) -> None:
        data = self.root / "complete.csv"
        shifted = self.root / "complete-shifted.csv"
        covariance = self.root / "covariance.csv"
        output = self.root / "output"
        shifted_output = self.root / "shifted-output"
        self.write_data(data, SIZES)
        self.write_data(shifted, SIZES, heldout_shift=0.01)
        self.write_covariance(covariance)

        self.run_challenge(data, output, covariance)
        self.run_challenge(shifted, shifted_output, covariance)
        payload = json.loads((output / "challenge.json").read_text(encoding="utf-8"))
        shifted_payload = json.loads(
            (shifted_output / "challenge.json").read_text(encoding="utf-8")
        )
        self.assertEqual(payload["status"], "READY")
        self.assertTrue(payload["covariance_source"].startswith("full_edge_list:"))
        self.assertEqual(payload["protocol"]["training_sizes"], [65, 85, 130])
        self.assertEqual(payload["protocol"]["heldout_sizes"], [145, 170])
        self.assertFalse(payload["protocol"]["selection_uses_heldout"])
        self.assertIn(payload["power_correction_selection"]["selected_omega"], (1, 2, 3))
        self.assertEqual(len(payload["angular_design"]), 5)
        self.assertIn("fixed_13_8_h4_h8", {
            model["model"] for model in payload["models"]
        })
        self.assertGreater(
            payload["zero_effect_benchmark"]["heldout_chi_square"], 0
        )

        # Changing only held-out values cannot alter any selected/fitted training object.
        self.assertEqual(
            payload["power_correction_selection"],
            shifted_payload["power_correction_selection"],
        )
        for original, changed in zip(payload["models"], shifted_payload["models"]):
            self.assertEqual(original["model"], changed["model"])
            self.assertEqual(original["parameters"], changed["parameters"])
            self.assertEqual(original["training_chi_square"], changed["training_chi_square"])
            self.assertTrue(math.isfinite(original["condition_number"]))
            self.assertTrue(math.isfinite(original["heldout"]["chi_square"]))
            self.assertIn(
                "chi_square_improvement_over_zero", original["heldout"]
            )
        self.assertNotEqual(
            payload["models"][0]["heldout"]["chi_square"],
            shifted_payload["models"][0]["heldout"]["chi_square"],
        )

        with (output / "heldout_signed_errors.csv").open(newline="", encoding="utf-8") as handle:
            heldout = list(csv.DictReader(handle))
        with (output / "amplitude_drift.csv").open(newline="", encoding="utf-8") as handle:
            amplitudes = list(csv.DictReader(handle))
        with (output / "loso_predictions.csv").open(newline="", encoding="utf-8") as handle:
            loso = list(csv.DictReader(handle))
        self.assertEqual(len(heldout), 5 * 2 * 2)
        self.assertEqual({int(row["N"]) for row in heldout}, {145, 170})
        self.assertEqual(len(amplitudes), 5)
        self.assertEqual(len(loso), 5 * 3 * 2)
        self.assertTrue(all(row["signed_error_observed_minus_predicted"] for row in heldout))


if __name__ == "__main__":
    unittest.main()
