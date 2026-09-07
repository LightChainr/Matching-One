"""Tests for the #598 reflection-parity selection rule.

Each test names the wrong number it would stop us believing.
"""

from __future__ import annotations

import math
import sys
import unittest
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p398_intervention_transport import Generator, dense_expm, matvec  # noqa: E402
from p398_reflection_parity import (  # noqa: E402
    D0_NAMES,
    EPSILONS,
    READOUTS,
    REPORTED_R_POSITIVE,
    SECTOR_TOLERANCE,
    Propagator,
    apply_permutation,
    burnside_check,
    canonical_labels,
    catalan,
    coarsest_lumping,
    combine_rows,
    conjugate,
    decide,
    equivariance,
    exact_rows,
    is_involution,
    lumping_comparison,
    move_indicator,
    observable_partition,
    orbit_count,
    orbit_labels,
    parity_defect,
    parity_split,
    reflection,
    reflection_survey,
    rows_equal,
    selection_rule,
    state_permutation,
    tilted_rates,
    transpose_rows,
)


class DihedralAction(unittest.TestCase):
    def test_reflections_are_involutions_on_states(self) -> None:
        """Stops us believing a parity decomposition built on a non-involution.

        ``H_even`` and ``H_odd`` are defined by ``R K R``, and the +-1
        eigenspaces only exist if ``R^2 = I`` on the state space.  If the induced
        state permutation were not an involution the split would still compute a
        number, and that number would mean nothing.
        """

        for width in (4, 5, 6, 7):
            generator = Generator(width)
            for k in range(width):
                pi = state_permutation(generator, reflection(width, k))
                self.assertTrue(is_involution(pi))
                self.assertEqual(sorted(pi), list(range(generator.size)))

    def test_every_reflection_fixes_the_same_number_of_states(self) -> None:
        """Stops us believing the orbit formula identifies which reflection R is.

        The headline arithmetic -- ``r_positive = (Catalan + central binomial)/2``
        -- holds for *every* reflection, so quoting it as evidence that a
        particular reflection is the surviving symmetry would be an overclaim.
        The identification has to come from the dictionary, and this test pins
        the fact that makes that necessary.
        """

        for width in (4, 5, 6, 7, 8):
            generator = Generator(width)
            counts = set()
            for k in range(width):
                pi = state_permutation(generator, reflection(width, k))
                counts.add(sum(1 for i in range(generator.size) if pi[i] == i))
            self.assertEqual(counts, {math.comb(width, width // 2)})


class ReflectionIdentification(unittest.TestCase):
    def test_the_declared_dictionary_names_exactly_one_reflection(self) -> None:
        """Stops us believing an R that was picked to make the theorem work.

        If several reflections preserved the declared readouts, choosing among
        them after seeing the response table would be fitting the symmetry to
        the answer.  The dictionary has to select one on its own, and it does --
        ``wrap`` is the only D0 readout that is not fully dihedral-invariant, and
        it fixes the reflection that swaps points 0 and w-1.
        """

        for width in (4, 5, 6, 7, 8):
            survey = reflection_survey(Generator(width))
            self.assertEqual(survey["reflections_preserving_D0"], [width - 1])
            self.assertEqual(survey["task_compatible_k"], width - 1)

    def test_burnside_matches_a_direct_orbit_enumeration(self) -> None:
        """Stops us believing an orbit count computed by a formula we never checked.

        ``r_positive`` is compared against ``(|X| + |Fix|)/2``.  Independent
        means: enumerate the orbits by union-find on the permutation and count
        them directly.  A wrong count here would make the whole quotient reading
        agree with #580 by accident.
        """

        for width in (4, 5, 6, 7):
            generator = Generator(width)
            pi = state_permutation(generator, reflection(width, width - 1))
            direct = len({frozenset((i, pi[i])) for i in range(generator.size)})
            self.assertEqual(orbit_count(pi), direct)
            check = burnside_check(generator, width - 1)
            self.assertEqual(check["orbits_counted"], REPORTED_R_POSITIVE[width])
            self.assertTrue(check["formula_matches_count"])
            self.assertEqual(
                check["r_reflect_formula"],
                (catalan(width) + math.comb(width, width // 2)) // 2,
            )


class Lumping(unittest.TestCase):
    def test_the_coarsest_lumping_is_the_orbit_partition_not_merely_its_size(self) -> None:
        """Stops us believing 750 blocks and 750 orbits are the same 750 blocks.

        Cardinality agreement is what #598 reports; it is much weaker than what
        Gate 1 asks.  Two partitions can have equal block counts and share no
        block at all, and the whole symmetry-quotient reading rests on them being
        the same partition.
        """

        for width in (4, 5, 6):
            generator = Generator(width)
            comparison = lumping_comparison(generator, width - 1)["D0"]
            self.assertTrue(comparison["coarsest_equals_orbit_partition"])
            self.assertEqual(comparison["coarsest_blocks"], REPORTED_R_POSITIVE[width])
            self.assertTrue(comparison["orbit_partition_is_admissible"])

    def test_a_non_invariant_readout_destroys_the_quotient(self) -> None:
        """Stops us believing a quotient that survives any dictionary we like.

        At odd width ``halves_linked`` is not R-invariant, so demanding the
        lumping carry it must collapse the partition to the identity.  If it did
        not, the quotient would not actually be the task's symmetry quotient and
        the parity story would be decoration.
        """

        generator = Generator(5)
        every = tuple(name for name, _ in READOUTS)
        initial = observable_partition(generator, every)
        joins = {move: (1.0 if move[0] == "join" else 0.0) for move in generator.moves}
        detaches = {move: (1.0 if move[0] == "detach" else 0.0) for move in generator.moves}
        coarsest = coarsest_lumping(generator, (joins, detaches), initial)
        self.assertEqual(len(set(coarsest)), generator.size)

    def test_refinement_stops_at_a_partition_that_is_actually_stable(self) -> None:
        """Stops us believing a lumping the refinement quit on too early.

        Re-running the refinement from its own output must be a fixed point.  A
        premature stop would report a coarser partition than exists and inflate
        the agreement with the orbit count.
        """

        generator = Generator(5)
        joins = {move: (1.0 if move[0] == "join" else 0.0) for move in generator.moves}
        detaches = {move: (1.0 if move[0] == "detach" else 0.0) for move in generator.moves}
        first = coarsest_lumping(generator, (joins, detaches), observable_partition(generator, D0_NAMES))
        again = coarsest_lumping(generator, (joins, detaches), first)
        self.assertEqual(canonical_labels(first), canonical_labels(again))
        orbits = orbit_labels(generator, 4)
        self.assertEqual(canonical_labels(first), canonical_labels(orbits))


class ExactSplit(unittest.TestCase):
    def test_equivariance_holds_and_the_move_map_explains_why(self) -> None:
        """Stops us believing R G R = G as a coincidence of two big matrices.

        The operator identity is checked, but so is its reason: R sends the join
        at p to the join at (k-1)-p and the detach at p to the detach at k-p.  If
        the move map were wrong, uniform sums would still be invariant and a
        localized tilt might be too -- and the asymmetry the theorem needs would
        not exist.
        """

        for width in (4, 5, 6):
            report = equivariance(Generator(width), width - 1)
            self.assertTrue(report["R_J_R_equals_J"])
            self.assertTrue(report["R_D_R_equals_D"])
            self.assertTrue(report["R_G0_R_equals_G0"])
            self.assertTrue(report["all_moves_map_as_predicted"])

    def test_the_split_is_exact_and_reassembles(self) -> None:
        """Stops us believing an H_odd that is not the odd part of H_single.

        Every second-order slope in the artifact is a statement about *this*
        decomposition.  If the halves did not sum back to the original tilt in
        exact rational arithmetic, the ladder would be measuring a different
        perturbation than the one #580 applied.
        """

        for width in (4, 5, 6, 7):
            split = parity_split(Generator(width), width - 1)
            self.assertTrue(split["even_plus_odd_is_H_single"])
            self.assertTrue(split["R_Heven_R_equals_plus_Heven"])
            self.assertTrue(split["R_Hodd_R_equals_minus_Hodd"])
            self.assertTrue(split["mirrored_single_is_the_partner_join"])
            self.assertEqual(split["mirror_partner_point"], width - 2)
            self.assertGreater(split["odd_share_of_norm"], 0.5)

    def test_conjugation_is_an_involution_that_preserves_the_norm(self) -> None:
        """Stops us believing a conjugation that silently transposes or drops entries.

        ``R K R`` applied twice must return K exactly.  A transposed or
        index-shifted conjugation would still produce an "even" and an "odd"
        piece that sum correctly, and both parity relations would then be
        checking the wrong operator.
        """

        generator = Generator(5)
        pi = state_permutation(generator, reflection(5, 4))
        rows = exact_rows(generator, move_indicator(generator, "join", 0))
        self.assertTrue(rows_equal(conjugate(conjugate(rows, pi), pi), rows))
        doubled = combine_rows(rows, rows, +1)
        self.assertTrue(rows_equal(doubled, rows))
        cancelled = combine_rows(rows, rows, -1)
        self.assertTrue(all(not row for row in cancelled))


class Propagation(unittest.TestCase):
    def test_uniformization_matches_a_dense_matrix_exponential(self) -> None:
        """Stops us believing an integrand computed with the wrong propagator.

        Every vanishing number in the selection rule is a product of two
        propagated vectors.  Independent means: at width 4 the generator is 14x14,
        so a dense expm is available and the two must agree.  A propagator that
        was subtly wrong could make integrands vanish for the wrong reason.
        """

        generator = Generator(4)
        rows = generator.rows({move: 1.0 for move in generator.moves})
        dense = [[0.0] * generator.size for _ in range(generator.size)]
        for source, row in enumerate(rows):
            for column, value in row:
                dense[source][column] = value
        propagator = Propagator(rows, generator.size)
        vector = [float(index % 3) - 1.0 for index in range(generator.size)]
        for t in (0.25, 1.0, 2.0):
            exact = matvec(
                [[(j, v) for j, v in enumerate(r)] for r in dense_expm(dense, t)], vector
            )
            got = propagator.apply(7, vector, t)
            scale = math.sqrt(sum(v * v for v in exact))
            defect = math.sqrt(sum((a - b) ** 2 for a, b in zip(exact, got)))
            self.assertLess(defect / scale, 1e-12)

    def test_transpose_rows_is_the_actual_transpose(self) -> None:
        """Stops us believing an evolved source measure that never left the row side.

        The integrand pairs a forward-evolved readout with a backward-evolved
        source.  If the "transpose" were the same operator, both sides would stay
        even for a trivial reason and the selection rule would pass vacuously.
        """

        generator = Generator(4)
        rows = generator.rows({move: 1.0 for move in generator.moves})
        transposed = transpose_rows(rows, generator.size)
        left = [float(i) for i in range(generator.size)]
        right = [float(generator.size - i) for i in range(generator.size)]
        a = sum(left[i] * v for i, v in enumerate(matvec(rows, right)))
        b = sum(right[i] * v for i, v in enumerate(matvec(transposed, left)))
        self.assertLess(abs(a - b), 1e-9 * (abs(a) + 1))


class SelectionRule(unittest.TestCase):
    def test_parity_defect_separates_the_two_eigenspaces(self) -> None:
        """Stops us believing a readout classified into the wrong parity sector.

        The whole verdict splits pairs into "protected" and "exposed" by this
        number.  Misclassifying an odd readout as even would move a nonzero
        integrand into the protected set and turn a passing theorem into a
        failing one, or the reverse.
        """

        generator = Generator(5)
        pi = state_permutation(generator, reflection(5, 4))
        even = [1.0] * generator.size
        odd = [0.0] * generator.size
        for i in range(generator.size):
            if pi[i] != i:
                odd[i] = 1.0 if i < pi[i] else -1.0
        self.assertLess(parity_defect(even, pi, +1), 1e-15)
        self.assertLess(parity_defect(odd, pi, -1), 1e-15)
        self.assertGreater(parity_defect(odd, pi, +1), 0.5)

    def test_an_odd_readout_gives_a_visibly_nonzero_integrand(self) -> None:
        """Stops us believing a selection rule that cannot fail.

        At width 5 ``halves_linked`` is not R-even, so the theorem *requires* a
        nonzero first-order response there.  A run in which every integrand
        vanished would not be a stronger result -- it would mean the integrand is
        being computed as zero for an unrelated reason.
        """

        generator = Generator(5)
        split = parity_split(generator, 4)
        rule = selection_rule(generator, 4, split["_odd"])
        self.assertEqual(rule["odd_containing_readouts"], ["halves_linked"])
        self.assertTrue(rule["protected_all_vanish"])
        self.assertLess(rule["worst_protected_relative_integrand"], SECTOR_TOLERANCE)
        self.assertTrue(rule["exposed_channel_is_visible"])
        self.assertGreater(rule["largest_exposed_relative_integrand"], 0.01)

    def test_a_vacuous_run_is_reported_as_vacuous(self) -> None:
        """Stops us believing a pass from widths where nothing could have failed.

        At even widths every declared readout is R-even, so the exposed set is
        empty.  If only such widths were run the verdict must say the test was
        vacuous rather than claim the rule was confirmed.
        """

        base = {
            "orbit_count": {"matches_reported_r_positive": True},
            "lumping": {"D0": {"coarsest_equals_orbit_partition": True}},
            "equivariance": {"R_G0_R_equals_G0": True, "all_moves_map_as_predicted": True},
            "parity_split": {
                "R_Hodd_R_equals_minus_Hodd": True,
                "R_Heven_R_equals_plus_Heven": True,
                "even_plus_odd_is_H_single": True,
            },
            "selection_rule": {"protected_all_vanish": True, "exposed_channel_is_visible": False},
            "second_order": {
                "odd_vs_none": {"consistent_with_second_order": True},
                "even_vs_none": {"consistent_with_first_order": True},
            },
        }
        self.assertIn("VACUOUS", decide([base])["verdict"])
        exposed = json_copy(base)
        exposed["selection_rule"]["exposed_channel_is_visible"] = True
        self.assertIn("BRANCH_A", decide([exposed])["verdict"])
        broken = json_copy(exposed)
        broken["selection_rule"]["protected_all_vanish"] = False
        self.assertIn("SELECTION_RULE_FAILS", decide([broken])["verdict"])


def json_copy(value):
    if isinstance(value, dict):
        return {key: json_copy(item) for key, item in value.items()}
    return value


class SecondOrder(unittest.TestCase):
    def test_the_three_tilts_stay_generators_and_reassemble(self) -> None:
        """Stops us believing a slope measured on an operator that is not a generator.

        Uniformization requires non-negative rates.  A tilt that left the cone
        would still return numbers, and the ``eps^2`` fit would be reporting the
        decay of an unphysical operator rather than a selection rule.
        """

        generator = Generator(6)
        for epsilon in EPSILONS:
            single = tilted_rates(generator, 5, "single", epsilon)
            even = tilted_rates(generator, 5, "even", epsilon)
            odd = tilted_rates(generator, 5, "odd", epsilon)
            base = tilted_rates(generator, 5, "none", 0.0)
            for move in generator.moves:
                self.assertGreaterEqual(single[move], 0.0)
                self.assertGreaterEqual(odd[move], 0.0)
                self.assertAlmostEqual(
                    single[move] - base[move],
                    (even[move] - base[move]) + (odd[move] - base[move]),
                    places=12,
                )


class DeclaredManifest(unittest.TestCase):
    def test_the_reported_positive_counts_match_the_orbit_formula(self) -> None:
        """Stops us believing a match against numbers we quietly regenerated.

        ``REPORTED_R_POSITIVE`` is #580's published table.  Checking the formula
        against a freshly recomputed lumping instead would let both drift
        together and still agree.
        """

        self.assertEqual(REPORTED_R_POSITIVE, {4: 10, 5: 26, 6: 76, 7: 232, 8: 750})
        for width, reported in REPORTED_R_POSITIVE.items():
            self.assertEqual(
                (catalan(width) + math.comb(width, width // 2)) // 2, reported
            )
        self.assertEqual(catalan(8), 1430)
        self.assertTrue(all(epsilon > 0 for epsilon in EPSILONS))
        self.assertEqual(set(D0_NAMES), {"blocks", "singletons", "wrap"})


if __name__ == "__main__":
    unittest.main()
