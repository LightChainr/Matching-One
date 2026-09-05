#!/usr/bin/env python3
"""Lock the P3 manuscript evidence assembly to the committed ladder artifacts."""

from __future__ import annotations

import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import p3_manuscript_evidence_table as evidence  # noqa: E402


class P3EvidenceTable(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = evidence.assemble()
        cls.committed = json.loads(evidence.DEFAULT_OUTPUT.read_text(encoding="utf-8"))

    def test_the_committed_artifact_reproduces_exactly(self) -> None:
        """Stops us believing a table generated from a ladder artifact that has since moved."""
        self.assertEqual(json.loads(json.dumps(self.payload, sort_keys=True)), self.committed)

    def test_the_rendered_tables_do_not_drift_from_the_artifact(self) -> None:
        """Stops us believing a number in tables.md that the artifact no longer contains."""
        self.assertEqual(evidence.render(self.committed),
                         evidence.DEFAULT_MARKDOWN.read_text(encoding="utf-8"))

    def test_the_nominated_denominator_is_the_weakest_rung(self) -> None:
        """Stops us believing the denominator was weak by bad luck rather than by construction.

        The whole argument of the draft's section 3 is that r=1 was nominated as
        denominator while being the rung furthest from being resolved.  If some
        future response makes r=1 the strongest rung, that section is about a
        situation that no longer exists and must be rewritten, not re-run.
        """
        rows = {row["rung"]: row["sigma_from_zero"] for row in self.payload["response"]}
        self.assertEqual(min(rows, key=lambda rung: rows[rung]), 1)
        self.assertLess(rows[1], 4.0)
        self.assertGreater(rows[4], 20.0)

    def test_the_projective_statistic_is_fieller_squared_on_the_real_covariance(self) -> None:
        """Stops us believing a verdict change was caused by changing the statistic.

        Section 4.4 attributes every verdict change to the third rung.  That
        attribution is only honest if the projective statistic, restricted to the
        two rungs the frozen test used, is the frozen test.  The wrong number
        this stops us believing is a projective sigma that differs from Fieller
        for any reason other than the added rung.
        """
        control = self.payload["two_entry_control"]
        self.assertEqual(control["entry"], "r4_over_r1")
        self.assertEqual(len(control["rows"]), 8)
        self.assertLess(control["largest_relative_deviation"], 1e-12)
        for row in control["rows"]:
            self.assertAlmostEqual(row["projective_statistic_two_entries"],
                                   row["fieller_z_squared"],
                                   delta=1e-9 * max(abs(row["fieller_z_squared"]), 1.0))

    def test_exactly_two_verdicts_flip_and_both_flip_toward_exclusion(self) -> None:
        """Stops us believing the reanalysis is a uniform tightening of every verdict.

        A correction that pushed every competitor the same way would be the
        signature of an error, not of a fix.  This pins the actual shape: two
        competitors move from compatible to excluded, and no competitor moves
        from excluded to compatible.
        """
        changed = [row for row in self.payload["verdicts"] if row["verdict_changed"]]
        self.assertEqual(sorted(row["competitor"] for row in changed),
                         ["plain_area_scaling", "q4_jordan_weight4"])
        for row in changed:
            self.assertEqual(row["frozen_verdict"], "compatible")
            self.assertEqual(row["projective_verdict"], "excluded")
            self.assertGreater(row["projective_sigma_three_rungs"], 6.9)

    def test_at_least_one_exclusion_weakens_under_the_correction(self) -> None:
        """Stops us believing the third rung can only ever add significance.

        Section 4.4 rests on the correction moving numbers in both directions.
        If every projective sigma exceeded its frozen sigma, that sentence would
        be false and the draft would be claiming an even-handedness it does not
        have.
        """
        weakened = [row for row in self.payload["verdicts"]
                    if math.isfinite(row["projective_sigma_three_rungs"])
                    and row["projective_sigma_three_rungs"] < row["fieller_sigma_one_entry"]]
        self.assertTrue(weakened)
        self.assertIn("no_modulus_dependence", [row["competitor"] for row in weakened])

    def test_the_measured_curvature_is_negative_and_no_competitor_can_be(self) -> None:
        """Stops us believing a competitor might reach the measured concavity.

        Section 4.5's class-level claim is that the sign is unreachable, not that
        it is unlikely.  The wrong number this stops us believing is a predicted
        curvature that is negative for some competitor, which would turn the
        claim from a class statement into an amplitude comparison.
        """
        curvature = self.payload["curvature"]
        self.assertLess(curvature["measured"]["value"], 0.0)
        self.assertLess(curvature["measured"]["z"], -3.0)
        low, high = curvature["measured"]["z_over_admissible_correlations"]
        self.assertLess(high, 0.0, "the sign must not flip anywhere in the admissible range")
        self.assertLess(low, high)
        self.assertTrue(curvature["no_competitor_predicts_a_negative_curvature"])
        self.assertTrue(curvature["at_least_one_competitor_predicts_a_positive_curvature"])

    def test_the_curvature_weights_are_the_second_divided_difference(self) -> None:
        """Stops us believing a reweighted contrast that is not the divided difference.

        The functional's whole value is that it is exactly 1 on r squared and
        exactly 0 on every line.  A weight vector that has drifted would still
        produce a plausible z, and the class-level claim would silently stop
        holding.
        """
        weights = self.payload["curvature"]["weights"]
        rungs = (1.0, 2.0, 4.0)
        self.assertAlmostEqual(sum(w * r * r for w, r in zip(weights, rungs)), 1.0, places=12)
        self.assertAlmostEqual(sum(w for w in weights), 0.0, places=12)
        self.assertAlmostEqual(sum(w * r for w, r in zip(weights, rungs)), 0.0, places=12)

    def test_every_competitor_but_one_needs_a_spin8_amplitude_above_the_assumed_bound(self) -> None:
        """Stops us believing the dropped rung is consistent with the bound that dropped it.

        Section 4.6's fork exists only because the requirement exceeds 1 for
        every competitor except one, and because that one is excluded without
        the disputed rung.  A requirement table in which several competitors sat
        below 1 would dissolve the fork.
        """
        rows = {row["competitor"]: row["required_abs_A8_over_A4_to_reach_r2"]
                for row in self.payload["spin8"]["rows"]}
        below_one = [name for name, value in rows.items() if value < 1.0]
        self.assertEqual(below_one, ["no_modulus_dependence"])
        verdicts = {row["competitor"]: row for row in self.payload["verdicts"]}
        self.assertEqual(verdicts["no_modulus_dependence"]["projective_verdict"], "excluded")
        self.assertTrue(
            verdicts["no_modulus_dependence"]["verdict_survives_the_missing_covariance"])
        self.assertGreater(min(value for name, value in rows.items() if name not in below_one),
                           3.0)
        # the exact (u+1) solution, not the leading-order gap over lambda, which
        # would put the top of this column at 785 rather than 17.5
        self.assertLess(max(rows.values()), 20.0)

    def test_exactly_one_verdict_is_undetermined_by_the_missing_covariance(self) -> None:
        """Stops us believing a stable exclusion that the absent cov(r2,r4) could undo.

        Section 5.2 reports one undetermined verdict.  If a second competitor
        became unstable, the draft would be presenting as settled an exclusion
        that a covariance entry nobody stored could reverse.
        """
        unstable = [row["competitor"] for row in self.payload["verdicts"]
                    if not row["verdict_survives_the_missing_covariance"]]
        self.assertEqual(unstable, ["bare_aspect_ratio"])
        self.assertFalse(self.payload["covariance_is_complete"])
        low, high = self.payload["admissible_corr_r2_r4"]
        self.assertLess(-1.0, low)
        self.assertLess(high, 1.0)

    def test_no_competitor_survives_that_the_draft_reports_as_excluded(self) -> None:
        """Stops us believing the draft's seven-of-eight headline after a competitor moves."""
        excluded = [row["competitor"] for row in self.payload["verdicts"]
                    if row["projective_verdict"] == "excluded"]
        self.assertEqual(len(excluded), 7)
        self.assertNotIn("bare_aspect_ratio", excluded)


if __name__ == "__main__":
    unittest.main()
