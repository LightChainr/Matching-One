from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from rng_domain_policy import derive_size_seed  # noqa: E402


class P321FreshMultiscaleDecisionFreezeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.freeze = json.loads(
            (ROOT / "analysis/p321_fresh_multiscale_decision_freeze.json")
            .read_text(encoding="utf-8")
        )

    def test_domains_are_exact_fresh_and_reproducible(self) -> None:
        rng = self.freeze["rng"]
        intervals = []
        for row in rng["domains"]:
            self.assertEqual(
                row["effective_seed"],
                derive_size_seed(rng["base_seed"], rng["experiment_tag"], row["N"]),
            )
            self.assertEqual(
                row["replica_counter_last_exclusive"] - row["replica_counter_first"],
                100000,
            )
            intervals.append(
                (row["replica_counter_first"], row["replica_counter_last_exclusive"])
            )
        self.assertEqual(len({row["effective_seed"] for row in rng["domains"]}), 3)
        for index, first in enumerate(intervals):
            for second in intervals[index + 1 :]:
                self.assertTrue(first[1] <= second[0] or second[1] <= first[0])
        pilot_seeds = {row["seed"] for row in rng["pilot_domains_excluded"]}
        self.assertFalse(pilot_seeds & {row["effective_seed"] for row in rng["domains"]})

    def test_sample_count_and_decisions_are_frozen(self) -> None:
        design = self.freeze["fresh_design"]
        self.assertEqual(design["samples_per_shape"], 100000)
        self.assertEqual(design["batches"], 50)
        self.assertEqual(design["samples_per_batch"], 2000)
        self.assertFalse(self.freeze["frozen_models"]["free_exponent_fit"])
        self.assertFalse(self.freeze["frozen_models"]["E4_curve_refit"])
        self.assertEqual(
            self.freeze["decision"]["conditional_E4"]["critical_chi_square"],
            11.34486673014437,
        )
        self.assertFalse(self.freeze["execution"]["authorized_by_this_file"])
        self.assertTrue(self.freeze["execution"]["no_interim_scoring"])

    def test_power_is_based_on_bias_corrected_pilot_noncentrality(self) -> None:
        pilot = self.freeze["pilot"]["frozen_E4_residual"]
        power = self.freeze["power_freeze"]
        lambda20 = max(pilot["chi_square"] - pilot["degrees_of_freedom"], 0)
        self.assertAlmostEqual(power["lambda_20k"], lambda20)
        self.assertAlmostEqual(power["lambda_100k"], 5 * lambda20)
        self.assertAlmostEqual(
            power["expected_chi_square_100k"],
            pilot["degrees_of_freedom"] + 5 * lambda20,
        )
        self.assertGreater(power["power_100k_at_alpha_0.01"], 0.94)
        self.assertLess(power["power_60k_at_alpha_0.01"], 0.73)


if __name__ == "__main__":
    unittest.main()
