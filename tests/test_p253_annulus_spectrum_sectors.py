#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from analyze_p253_annulus_spectrum_sectors import (  # noqa: E402
    CHANNELS,
    RADII,
    add_matrices,
    basis,
    combine_blocks,
    rank2_residual,
)
from score_p253_n365_heldout import read_n365  # noqa: E402


class P253AnnulusSpectrumSectorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result_path = ROOT / "results/p253-annulus-spectrum/latest.json"
        if cls.result_path.exists():
            cls.result = json.loads(cls.result_path.read_text(encoding="utf-8"))

    def test_rank_two_amplitude_calibration_predicts_exact_synthetic_curve(self) -> None:
        for model, parameters in (("J2", [-1.25]), ("R2_gap1", [0.4])):
            spectra = {channel: parameters for channel in CHANNELS}
            vector = []
            for channel_index, _ in enumerate(CHANNELS):
                amplitudes = [0.7 + channel_index, -0.2 + 0.3 * channel_index]
                for radius in RADII:
                    vector.append(sum(a * b for a, b in
                                      zip(basis(model, radius, parameters), amplitudes)))
            residual = rank2_residual(model, spectra, vector)
            self.assertLess(max(abs(value) for value in residual), 1e-12)

    def test_dependency_groups_are_block_diagonal_only_between_streams(self) -> None:
        old = json.loads((ROOT / "results/server-20260829/P225-norm5-multiradius/analysis.json")
                         .read_text(encoding="utf-8"))["contrast_vector"]
        new = read_n365(ROOT / "results/server-20260829/P253-n365-annulus/raw/n365_200k.metadata.json")
        labels, _, covariance = combine_blocks(old, new)
        for i, left in enumerate(labels):
            for j, right in enumerate(labels):
                left_old = not left.startswith("N365_")
                right_old = not right.startswith("N365_")
                if left_old != right_old:
                    self.assertEqual(covariance[i][j], 0.0)
        # The within-N365 plus/minus covariance is deliberately retained.
        plus = labels.index("N365_R2_Delta_A_plus")
        minus = labels.index("N365_R2_Delta_A_minus")
        self.assertNotEqual(covariance[plus][minus], 0.0)

    def test_committed_profile_counts_and_nested_order(self) -> None:
        if not self.result_path.exists():
            self.skipTest("generated spectrum result is not present")
        payload = self.result
        self.assertEqual(payload["schema"], "matching-one/p253-annulus-continuous-spectrum/v1")
        profiles = payload["continuous_spectrum_profiles"]
        self.assertEqual(profiles["per_sector"]["A_plus"]["J2"]["degrees_of_freedom"], 5)
        self.assertEqual(profiles["joint_shared_spectrum"]["J2"]["degrees_of_freedom"], 11)
        self.assertEqual(profiles["joint_sector_separated_spectra"]["J2"]["degrees_of_freedom"], 10)
        for model, record in profiles["sector_sharing_likelihood_ratio"].items():
            self.assertGreaterEqual(record["delta_chi_square"], -1e-8, model)

    def test_heldout_covariance_is_sum_of_independent_group_components(self) -> None:
        if not self.result_path.exists():
            self.skipTest("generated spectrum result is not present")
        for model in ("J2", "R2_gap1"):
            for scope in ("shared", "sector_separated"):
                record = self.result["heldout_N365_spectrum_transfer"][model][scope]
                components = record["covariance_components"]
                expected = add_matrices(
                    components["G_old_spectrum_training_jackknife"],
                    components["G_n365_calibration_and_target_jackknife"],
                )
                for row, target in zip(record["covariance"], expected):
                    for value, wanted in zip(row, target):
                        self.assertAlmostEqual(value, wanted, places=15)
                self.assertEqual(record["degrees_of_freedom"], 4)

    def test_ineligible_archives_are_not_silent_evidence_rows(self) -> None:
        if not self.result_path.exists():
            self.skipTest("generated spectrum result is not present")
        inventory = self.result["archive_eligibility"]
        self.assertEqual(len(inventory), 5)
        self.assertTrue(all(not row["eligible_for_numeric_merge"] for row in inventory))


if __name__ == "__main__":
    unittest.main()
