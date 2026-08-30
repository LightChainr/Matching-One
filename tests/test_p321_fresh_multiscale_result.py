from __future__ import annotations

import hashlib
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results/local-20260830/P321-equal-area-fresh-100k"
FREEZE = json.loads(
    (ROOT / "analysis/p321_fresh_multiscale_decision_freeze.json").read_text(
        encoding="utf-8"
    )
)
MANIFEST = json.loads(
    (ROOT / "analysis/p321_fresh_multiscale_result_manifest.json").read_text(
        encoding="utf-8"
    )
)
SCORE = json.loads((RESULT / "multiscale_score.json").read_text(encoding="utf-8"))


class P321FreshMultiscaleResultTests(unittest.TestCase):
    def test_campaigns_match_every_frozen_domain(self) -> None:
        frozen = {row["N"]: row for row in FREEZE["rng"]["domains"]}
        for n in FREEZE["fresh_design"]["N_order"]:
            campaign = json.loads((RESULT / f"N{n}/campaign.json").read_text())
            domain = frozen[n]
            self.assertEqual(campaign["samples_per_shape"], 100000)
            self.assertEqual(campaign["batches"], 50)
            self.assertEqual(campaign["seed"], domain["effective_seed"])
            self.assertEqual(
                campaign["replica_counter_first"], domain["replica_counter_first"]
            )
            self.assertEqual(
                campaign["replica_counter_last_exclusive"],
                domain["replica_counter_last_exclusive"],
            )
            self.assertEqual(campaign["git_commit"], MANIFEST["freeze_commit"])

    def test_full_covariance_and_square_gates_pass(self) -> None:
        self.assertEqual([row["N"] for row in SCORE["campaigns"]], [144, 576, 1296])
        for campaign in SCORE["campaigns"]:
            self.assertTrue(campaign["square_histograms_byte_identical"])
            self.assertTrue(campaign["square_moments_byte_identical"])
            self.assertEqual(len(campaign["root_covariance"]), 5)
            self.assertTrue(all(len(row) == 5 for row in campaign["root_covariance"]))

    def test_frozen_decisions_are_applied_without_refit(self) -> None:
        fit = SCORE["scale_fit"]
        e4 = fit["conditional_thermal_Q4_E4_score"]
        self.assertFalse(fit["free_exponent_fit"])
        self.assertLess(
            fit["fit_chi_square"],
            FREEZE["decision"]["scale_fit"]["critical_chi_square"],
        )
        self.assertLess(
            e4["chi_square"],
            FREEZE["decision"]["conditional_E4"]["critical_chi_square"],
        )
        self.assertEqual(e4["rho_order"], FREEZE["fresh_design"]["primary_rho"])
        self.assertTrue(e4["endpoint_rho_9_is_diagnostic_only"])
        self.assertEqual(MANIFEST["decision"]["scale_fit"]["decision"], "pass")
        self.assertEqual(
            MANIFEST["decision"]["conditional_E4"]["decision"], "do_not_reject"
        )

    def test_package_checksums_verify(self) -> None:
        for line in (RESULT / "checksums.sha256").read_text().splitlines():
            expected, relative = line.split("  ", 1)
            observed = hashlib.sha256((RESULT / relative).read_bytes()).hexdigest()
            self.assertEqual(observed, expected, relative)


if __name__ == "__main__":
    unittest.main()
