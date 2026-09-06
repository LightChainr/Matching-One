#!/usr/bin/env python3
"""Lock the #582 shape-flow reading to its committed artifact and its controls."""

from __future__ import annotations

import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_wasserstein_shape_flow as flow  # noqa: E402


class ShapeFlowArtifact(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = json.loads(flow.DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        cls.primary = cls.payload["scored"][cls.payload["primary_weighting"]]

    def test_the_primary_reading_is_the_spin_four_corrected_one(self) -> None:
        """Stops the naive orientation average becoming the headline.

        The equal-weight residue alternates in sign along a lineage and has
        opposite sign between the two Gaussian lineages, which is a pattern
        easily mistaken for a shape flow that differs between lineages.  The
        wrong number this stops us reporting is a direction read off the
        uncorrected law.
        """
        self.assertEqual(self.payload["primary_weighting"], "spin0")
        self.assertIn("equal", self.payload["scored"])

    def test_the_equal_weight_residue_alternates_along_every_lineage(self) -> None:
        """Stops us believing the correction was unnecessary.

        If the naive net cos4theta were small everywhere, the whole second
        reconstruction would be dead weight and the note's account of why it
        exists would be false.
        """
        sizes = self.payload["sizes"]
        for lineage_sizes in self.payload["lineages"].values():
            residues = [sizes[str(size)]["equal_weight_net_cos4theta"]
                        for size in lineage_sizes]
            for left, right in zip(residues, residues[1:]):
                self.assertLess(left * right, 0.0,
                                f"residues do not alternate: {residues}")

    def test_the_same_production_null_control_is_consistent_with_zero(self) -> None:
        """Stops us believing a chi-square of 10^5 that the pipeline itself produced.

        Two disjoint halves of one production differ in nothing.  If this control
        showed shape flow, every transition number in the artifact would be the
        reconstruction, the binomial window, the bisection tolerance or the
        jackknife rather than the physics.
        """
        controls = self.payload["same_production_null_control"]
        self.assertGreaterEqual(len(controls), 3)
        for row in controls:
            self.assertLess(row["equivalent_sigma"], 3.0, row["n"])
            self.assertGreater(row["p_value"], 0.01, row["n"])

    def test_the_null_control_displacement_is_orders_below_a_transition(self) -> None:
        """Stops a null control that passes only because it had nothing to detect.

        A half-to-half displacement comparable to a real transition would mean
        the control has power; one a thousandth the size means it does not, and
        saying so is the difference between a control and a decoration.
        """
        largest_null = max(row["displacement_norm"]
                           for row in self.payload["same_production_null_control"])
        smallest_real = min(
            math.sqrt(sum(x * x for x in row["affine_only"]["displacement"]))
            for row in self.primary["transitions"])
        self.assertLess(largest_null * 100.0, smallest_real)

    def test_no_transition_is_consistent_with_pure_center_and_width(self) -> None:
        """Stops W0 being reported as a live option.

        The whole reading rests on the affine tangent failing at every
        transition.  One transition consistent with affine would make the
        headline false.
        """
        for row in self.primary["transitions"]:
            self.assertGreater(row["affine_only"]["statistic"], 1e4, row["label"])
            self.assertEqual(row["affine_only"]["degrees_of_freedom"], 7)

    def test_a_frozen_generator_beats_the_best_of_two_hundred_random_ones(self) -> None:
        """The control that makes the held-out number a finding rather than arithmetic.

        A direction fitted on four transitions removing 99 percent of the fifth's
        chi-square means nothing unless a direction that was never fitted removes
        much less.  The wrong number this stops us believing is a transfer that
        any smooth direction would have achieved.
        """
        held = {(row["held_out"], row["rank"]): row for row in self.primary["held_out"]}
        controls = {row["label"]: row for row in self.primary["random_direction_control"]}
        for row in self.primary["transitions"]:
            label = row["label"]
            trained = held[(label, 1)]
            removed = 1.0 - (trained["statistic_with_frozen_generators"]
                             / trained["statistic_affine_only"])
            self.assertGreater(removed, 0.98, label)
            self.assertLess(controls[label]["median_random_fraction_removed"], 0.10, label)
            self.assertLess(controls[label]["best_random_fraction_removed"], 0.40, label)
            self.assertGreater(removed, 3.0 * controls[label]["best_random_fraction_removed"],
                               label)

    def test_the_held_out_transition_is_never_in_its_own_training_set(self) -> None:
        """Stops a generator that was fitted on the transition it then explains."""
        for row in self.primary["held_out"]:
            self.assertNotIn(row["held_out"], row["trained_on"])
            self.assertEqual(len(row["trained_on"]),
                             len(self.primary["transitions"]) - 1)

    def test_rank_one_does_not_close_the_flow(self) -> None:
        """Stops W1 being claimed from a 99 percent reduction.

        Ninety-nine percent of 400,000 is still 4,000 on six degrees of freedom.
        The note's fifth outcome -- a dominant direction plus a resolved
        remainder -- is only honest if this stays true.
        """
        for row in self.primary["held_out"]:
            if row["rank"] != 1:
                continue
            self.assertGreater(row["statistic_with_frozen_generators"],
                               row["degrees_of_freedom_with_frozen_generators"] * 5.0,
                               row["held_out"])

    def test_the_two_lineage_groups_are_further_apart_than_they_are_wide(self) -> None:
        """Stops the structured remainder being smoothed into one direction.

        Within {85->170, 170->425} and within {65->130, 145->290} the shape
        directions agree to a few degrees; between the groups they sit tens of
        degrees apart.  The wrong reading this stops is 'one direction, plus
        noise'.
        """
        angles = {tuple(sorted(row["pair"])): row["acute_angle_degrees"]
                  for row in self.primary["pairwise_angles"]}
        within = max(angles[tuple(sorted(pair))] for pair in
                     (("85->170", "170->425"), ("65->130", "145->290")))
        between = min(angles[tuple(sorted(pair))] for pair in
                      (("65->130", "85->170"), ("145->290", "85->170"),
                       ("65->130", "170->425"), ("145->290", "170->425")))
        self.assertLess(within, 20.0)
        self.assertGreater(between, 2.0 * within)

    def test_the_answer_does_not_depend_on_the_frozen_level_grid(self) -> None:
        """Stops a rank read off a window that was chosen to produce it.

        #582 warns against picking u_min by whichever value minimises the rank.
        The wrong number is a transfer fraction that a narrower or wider grid
        would not reproduce.
        """
        rows = self.payload["level_grid_robustness"]
        self.assertGreaterEqual(len(rows), 4)
        for row in rows:
            low, high = row["rank1_held_out_fraction_removed_range"]
            self.assertGreater(low, 0.98, row["grid"])
            self.assertLessEqual(high, 1.0, row["grid"])

    def test_the_cross_size_covariance_is_dropped_on_evidence(self) -> None:
        """Stops an assumed-zero cross term.

        Two lineages reuse a seed across sizes.  The cross term is measured
        batch by batch rather than argued away, and the wrong number this stops
        us believing is a Var(Q_M - Q_N) that silently omitted a real coupling.
        """
        for row in self.primary["transitions"]:
            coupling = row["coupling"]
            self.assertGreater(coupling["batches_compared"], 50)
            self.assertTrue(coupling["cross_term_droppable"], row["label"])
            self.assertLess(coupling["largest_absolute_correlation"],
                            3.0 * coupling["correlation_expected_from_noise_alone"])

    def test_the_claim_boundary_refuses_a_mechanism_and_an_exponent(self) -> None:
        """Stops a direction being reported as a field, or three sizes as a power law."""
        joined = " ".join(self.payload["not_established"]).lower()
        self.assertIn("mechanism", joined)
        self.assertIn("exponent", joined)
        self.assertIn("spin 8", joined)


if __name__ == "__main__":
    unittest.main()
