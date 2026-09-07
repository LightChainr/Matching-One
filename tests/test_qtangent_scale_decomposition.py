#!/usr/bin/env python3
"""Lock #581's exact Q-tangent gate: the identity, the split, the scales, the spectrum."""

from __future__ import annotations

import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import qtangent_scale_decomposition as qt  # noqa: E402


class ExactGate(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.table = qt.enumerate_torus(2)
        cls.gate = qt.exact_identity_gate(cls.table)

    def test_the_duality_identity_holds_on_every_configuration(self) -> None:
        """Stops us splitting the Q score with a map that is not the duality map.

        ``T - T* = X`` is what makes ``T = V/2 + B_even/2 + X/2`` a decomposition
        at all.  Naive bit-complement -- the wrong dual, which this repository
        has already had to rule out once -- fails it on most configurations, and
        every scale profile built on top would then be a profile of nothing.
        """
        self.assertEqual(self.gate["failures"]["t_minus_t_dual_is_x"], 0)
        self.assertEqual(self.gate["configurations"], 256)

    def test_the_two_algebraic_consequences_also_hold(self) -> None:
        """Stops the two derived identities being quoted without ever being checked."""
        self.assertEqual(self.gate["failures"]["t_is_half_v_plus_half_b"], 0)
        self.assertEqual(self.gate["failures"]["b_minus_b_dual_is_two_x"], 0)

    def test_the_image_rank_is_the_largest_component_rank(self) -> None:
        """Stops us reading X off one component when two of them wrap.

        ``X = r - 1`` needs the rank of the *image* in H1(T^2), the span of every
        component's winding lattice.  Two disjoint non-contractible cycles on a
        torus are parallel, so the span equals the largest component's rank --
        but that is an argument, and the wrong number it would hide is an X that
        is too small on a configuration with two wrapping clusters.
        """
        self.assertEqual(self.gate["failures"]["span_rank_differs_from_max"], 0)

    def test_all_three_homology_ranks_are_realised(self) -> None:
        """Stops a gate that passed because the torus was too small to have topology."""
        self.assertEqual(self.gate["x_values_realised"], [-1, 0, 1])

    def test_a_wrong_dual_map_breaks_the_identity(self) -> None:
        """Shows the gate has power rather than being satisfied by construction.

        Naive bit-complement occupies the vacant primal edges instead of the
        identified dual edges.  If the identity survived that too, passing it
        would say nothing about the transport being right.
        """
        pairs = qt.square_bond_pairs(2)
        bonds = len(pairs)
        failures = 0
        for mask in range(1 << bonds):
            primal = qt.configuration_statistics(2, mask, pairs)
            complement = qt.configuration_statistics(2, ~mask & ((1 << bonds) - 1), pairs)
            left = (2 * primal["components"] + primal["open_edges"]
                    - 2 * complement["components"] - complement["open_edges"])
            if left != 2 * (primal["span_rank"] - 1):
                failures += 1
        self.assertGreater(failures, 0)


class Split(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.table = qt.enumerate_torus(2)
        cls.rows = qt.covariance_split(cls.table)

    def test_every_observable_splits_exactly(self) -> None:
        """Stops a Betti piece and a topological piece that do not add up.

        If the two do not sum to Cov(O,T), then asking which of them a response
        loads on is asking about two things that are not a decomposition of the
        score.
        """
        self.assertEqual(len(self.rows), len(qt.OBSERVABLES))
        for row in self.rows:
            self.assertTrue(row["split_is_exact"], row["observable"])
            self.assertEqual(Fraction(row["cov_with_q_score"]),
                             Fraction(row["betti_even_piece"])
                             + Fraction(row["ambient_homology_piece"]))

    def test_the_open_edge_count_is_purely_topological_on_this_torus(self) -> None:
        """Stops a split that reports a Betti piece where an exact zero belongs.

        ``|A|`` is duality-odd on a self-dual torus, so its duality-even Betti
        piece is exactly zero and its whole covariance with the score is
        ambient.  A nonzero value here would mean the even projection is not
        even.
        """
        row = next(r for r in self.rows if r["observable"] == "open_edges")
        self.assertEqual(Fraction(row["betti_even_piece"]), 0)
        self.assertEqual(row["topological_fraction"], 1.0)

    def test_the_two_pieces_can_carry_opposite_signs(self) -> None:
        """Stops us reading a topological fraction as if it were a proportion.

        For the component count the two pieces have opposite signs and largely
        cancel, so the fraction leaves [0,1].  A reader who took it for a share
        of variance would misread the whole table.
        """
        row = next(r for r in self.rows if r["observable"] == "components")
        self.assertLess(row["topological_fraction"], 0.0)


class Scales(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.table = qt.enumerate_torus(2)
        cls.decomposition = qt.scale_decomposition(
            cls.table, ("wrap_either", "twice_b_even", "x_source"))

    def test_the_filtration_is_strictly_nested_and_ends_at_every_bond(self) -> None:
        """Stops a 'filtration' whose levels are not a filtration.

        The telescoping identity needs F_(j-1) contained in F_j and F_last equal
        to everything; a level ordering that repeated or skipped bonds would
        still produce plausible-looking Gamma_j that summed to the wrong total.
        """
        for length in (2, 3):
            levels = qt.spatial_filtration_levels(length)
            self.assertEqual(levels[-1], list(range(2 * length * length)))
            for earlier, later in zip(levels, levels[1:]):
                self.assertLess(set(earlier), set(later))

    def test_every_increment_matrix_is_positive_semidefinite(self) -> None:
        """Stops a cross entry being read at a scale where the matrix is not PSD.

        Gamma_j[O, B_even] only has an unambiguous scale-resolved meaning because
        the whole matrix is PSD; a negative eigenvalue would mean the increment
        is not a covariance and the cross term is not bounded by the diagonals.
        """
        for level in self.decomposition["levels"]:
            self.assertTrue(level["is_positive_semidefinite"], level["level"])

    def test_the_scales_telescope_to_the_total_covariance(self) -> None:
        """Stops a scale profile that is not a decomposition of the total.

        The wrong number is a per-scale cross-covariance whose sum differs from
        Cov(O, B_even), which would mean the profile and the total describe
        different objects.
        """
        self.assertTrue(self.decomposition["telescopes_to_the_total_covariance"])


class BooleanDegree(unittest.TestCase):
    def test_the_cross_spectrum_sums_to_the_covariance(self) -> None:
        """Stops a degree profile whose total is not the covariance it resolves.

        This is the rho = 1 end of the noise semigroup, where the operator is the
        identity, so the sum over degrees must reproduce the plain covariance
        exactly.  Any discrepancy is a bug in the transform or in the
        normalisation, not a physical statement.
        """
        table = qt.enumerate_torus(2)
        for left, right in (("wrap_either", "twice_b_even"), ("wrap_either", "x_source")):
            report = qt.walsh_cross_spectrum(table, left, right)
            self.assertTrue(report["agrees_with_the_direct_covariance"], (left, right))
            self.assertEqual(Fraction(report["sum_over_degrees"]),
                             Fraction(report["direct_covariance"]))

    def test_the_degree_zero_term_is_excluded(self) -> None:
        """Stops the mean sneaking into a covariance.

        The empty subset carries E[O]E[Y]; including it would turn every
        cross-spectrum into a second moment and every reported scale profile
        into a profile of the means.
        """
        table = qt.enumerate_torus(2)
        report = qt.walsh_cross_spectrum(table, "wrap_either", "x_source")
        self.assertNotIn("0", report["by_degree"])


class DegreeOneParity(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = qt.degree_one_parity(qt.enumerate_torus(2))

    def test_the_betti_even_piece_has_no_degree_one_weight_at_all(self) -> None:
        """The cleanest separation the two typed pieces have, and it is exact.

        The dual transport sends a single-bond character to minus the character
        of the crossing bond, so the duality-even part keeps the antisymmetric
        degree-one combination, which vanishes when the two crossing bonds carry
        equal weight.  The wrong number this stops us believing is a degree-one
        Betti loading -- which would mean the degree-one part of a response's
        Q-score coupling is no longer purely ambient, and every reading built on
        that would need a fit instead of a parity argument.
        """
        piece = self.report["pieces"]["twice_b_even"]
        self.assertTrue(piece["all_degree_one_coefficients_vanish"])
        self.assertEqual(piece["distinct_values"], ["0"])

    def test_the_ambient_piece_does_carry_degree_one_weight(self) -> None:
        """Stops a vanishing that came from the transform rather than from parity.

        If both pieces had zero degree-one weight the result would be a bug in
        the Walsh normalisation, not a separation.
        """
        piece = self.report["pieces"]["x_source"]
        self.assertFalse(piece["all_degree_one_coefficients_vanish"])
        self.assertEqual(len(piece["distinct_values"]), 1,
                         "every bond must carry the same weight by symmetry")

    def test_the_symmetry_the_argument_needs_is_measured_not_assumed(self) -> None:
        """Stops the parity explanation resting on an unchecked lattice symmetry.

        The degree-one vanishing needs f^({i}) = f^({sigma(i)}) for every bond.
        That is what makes the antisymmetric combination zero, and it is the one
        step of the mechanism that could fail on a less symmetric lattice.
        """
        self.assertTrue(self.report["crossing_map_is_a_bijection"])
        for name, piece in self.report["pieces"].items():
            self.assertTrue(piece["invariant_under_the_crossing_map"], name)

    def test_the_mechanism_is_not_claimed_for_general_side_length(self) -> None:
        """Stops two enumerations being reported as a theorem."""
        self.assertIn("not proved here", self.report["mechanism"])


class ClaimBoundary(unittest.TestCase):
    def test_the_artifact_refuses_a_field_identification(self) -> None:
        """Stops B_even being cited as the energy field, or X as a defect theory."""
        payload = json.loads(qt.DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        joined = " ".join(payload["not_established"]).lower()
        self.assertIn("energy field", joined)
        self.assertIn("defect", joined)
        self.assertIn("h4", joined)


if __name__ == "__main__":
    unittest.main()
