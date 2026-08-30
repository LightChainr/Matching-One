#!/usr/bin/env python3
"""Regression tests for the frozen Issue #55 acquisition and score."""

from __future__ import annotations

import copy
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import run_issue55_h4_h12_acquisition as runner  # noqa: E402
import score_issue55_h4_h12_orthogonal as scorer  # noqa: E402


MANIFEST = ROOT / "experiments" / "issue55_h4_h12_orthogonal_acquisition_20260830.json"
COMMIT = "c1a353a0718d86894ebf49f7b7200152e402ad09"


class Issue55AcquisitionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.manifest = runner.load_and_validate(MANIFEST)

    def test_frozen_signed_rows_and_rng_domains_emit_exact_smoke_commands(self) -> None:
        commands = runner.commands(
            self.manifest, "smoke", 0, Path("runner"), Path("out"), 1, COMMIT
        )
        self.assertEqual([item["N"] for item in commands], [305, 325])
        self.assertEqual(commands[0]["counter_interval"], [15550000000, 15550020000])
        self.assertEqual(commands[1]["seed"], 859043225199512)
        first = commands[0]["argv"]
        matrix_index = first.index("--first-matrix")
        self.assertEqual(first[matrix_index + 1:matrix_index + 5], ["17", "-4", "4", "17"])
        second_index = first.index("--second-matrix")
        self.assertEqual(first[second_index + 1:second_index + 5], ["16", "-7", "7", "16"])

    def test_three_production_shards_are_disjoint_and_exhaust_frozen_domain(self) -> None:
        intervals = []
        for shard in range(3):
            planned = runner.commands(
                self.manifest, "production", shard, Path("runner"), Path("out"), 8, COMMIT
            )
            self.assertEqual([item["samples"] for item in planned], [200000000, 200000000])
            self.assertEqual(planned[0]["counter_interval"], planned[1]["counter_interval"])
            intervals.append(planned[0]["counter_interval"])
        self.assertEqual(
            intervals,
            [[15551000000, 15751000000], [15751000000, 15951000000], [15951000000, 16151000000]],
        )

    def test_variance_pilot_withholds_means_and_freezes_first_mahalanobis_grid_hit(self) -> None:
        def campaign(offset: float):
            output = {}
            for n, scale in ((305, 1.0), (325, 0.98)):
                values = [offset + scale * (index - 9.5) * 0.004 for index in range(20)]
                output[n] = {"runs": [{
                    "metadata": {"samples": 20000, "batches": 20, "counter_first": 15550000000},
                    "delta_m": values,
                }]}
            return output

        first = scorer.pilot_result(self.manifest, campaign(0.0))
        shifted = scorer.pilot_result(self.manifest, campaign(1000.0))
        self.assertEqual(first["recommended_samples_per_design"], shifted["recommended_samples_per_design"])
        self.assertAlmostEqual(
            first["variance_estimates"]["305"]["centered_sampling_se"],
            shifted["variance_estimates"]["305"]["centered_sampling_se"],
            places=9,
        )
        self.assertEqual(
            shifted["variance_estimates"]["305"]["observed_target_mean"],
            "withheld_by_protocol",
        )
        distances = [
            row["conditional_h4_only_vs_equal_A12_mahalanobis"]
            for row in first["production_grid_projection"]
        ]
        chosen_index = [row["samples_per_design"] for row in first["production_grid_projection"]].index(
            first["recommended_samples_per_design"]
        )
        self.assertLess(distances[chosen_index - 1], 3.0)
        self.assertGreaterEqual(distances[chosen_index], 3.0)

    def test_two_column_score_recovers_a4_and_a12_not_repeated_h4_vote(self) -> None:
        a4, a12 = 0.7, 0.12
        campaigns = {}
        counter_base = self.manifest["production"]["counter_first"]
        for design in self.manifest["designs"]:
            n = int(design["N"])
            mean = n ** (-13.0 / 8.0) * (
                a4 * float(Fraction(design["delta_cos4"]))
                + a12 * float(Fraction(design["delta_cos12"]))
            )
            runs = []
            for shard in range(3):
                # Identical centered patterns make the exact campaign mean equal
                # to the requested two-column alternative while retaining >0 SE.
                values = [mean + (index - 49.5) * 1e-7 for index in range(100)]
                first = counter_base + shard * 200000000
                runs.append({
                    "metadata": {
                        "N": n,
                        "samples": 200000000,
                        "batches": 100,
                        "counter_first": first,
                        "counter_last": first + 200000000,
                        "commit": COMMIT,
                    },
                    "delta_m": values,
                })
            campaigns[n] = {"runs": runs}
        result = scorer.final_result(self.manifest, campaigns)
        fit = result["two_column_h4_h12"]
        self.assertAlmostEqual(fit["A4"], a4, places=12)
        self.assertAlmostEqual(fit["A12"], a12, places=12)
        self.assertEqual(fit["fit_df"], 0)
        self.assertIn("saturated", fit["fit_warning"])


if __name__ == "__main__":
    unittest.main()
