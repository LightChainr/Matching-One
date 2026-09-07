#!/usr/bin/env python3
"""Lock #580's exact intervention-transport control on the P398 state space."""

from __future__ import annotations

import math
import sys
import unittest
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import p398_intervention_transport as it  # noqa: E402


class GeneratorFamily(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.generator = it.Generator(4)

    def test_every_declared_intervention_is_a_markov_generator(self) -> None:
        """Stops us reading transport errors for an object that is not a generator.

        Row sums, non-negative off-diagonal rates across the whole declared eta
        range, exact affineness in eta and irreducibility are what make
        ``exp(t G_eta)`` a probability semigroup.  If any of them failed, every
        ``E_transported`` downstream would be an error about nothing.
        """

        for intervention in it.INTERVENTIONS:
            gate = it.generator_gate(self.generator, intervention)
            self.assertEqual(gate["row_sum_failures"], 0, intervention)
            self.assertEqual(gate["negative_offdiagonal_entries"], 0, intervention)
            self.assertEqual(gate["affine_in_eta_failures"], 0, intervention)
            self.assertTrue(gate["strongly_connected"], intervention)
            self.assertTrue(gate["passed"], intervention)

    def test_moves_that_do_not_change_the_state_are_dropped(self) -> None:
        """Stops a phantom exit rate inflating every relaxation time.

        Detaching an already-singleton point, or joining two points already in
        one block, is a no-op.  Recording it would add the same rate to an
        off-diagonal entry and to the exit-rate diagonal.  The two cancel in the
        generator but not in ``Lambda``, so the uniformization horizon and every
        fitted timescale would silently shift.
        """

        discrete = tuple(range(4))
        index = self.generator.states.index(discrete)
        self.assertTrue(
            all(self.generator.target[("detach", point)][index] is None for point in range(4)),
            "detaching an all-singleton state must be a no-op",
        )
        rows = self.generator.rows(self.generator.baseline_rates())
        for row in rows:
            self.assertAlmostEqual(sum(value for _, value in row), 0.0, places=12)

    def test_the_two_in_pencil_interventions_span_the_baseline_pencil(self) -> None:
        """Stops us calling an in-pencil control an independent intervention.

        ``G_0 = J + D`` and ``H = J - D`` span the same plane as ``{J, D}``, so
        the ``uniform_detach_only`` direction has to be an exact linear
        combination of the other two.  The out-of-pencil control must not be.
        """

        size = self.generator.size
        base = self.generator.rows(self.generator.baseline_rates())
        pencil = [base, self.generator.tangent_rows("uniform_join_minus_detach")]
        detach_only = self.generator.tangent_rows("uniform_detach_only")
        single = self.generator.tangent_rows("single_point_join")

        def dense(rows):
            matrix = [[0.0] * size for _ in range(size)]
            for source, row in enumerate(rows):
                for destination, value in row:
                    matrix[source][destination] = value
            return [value for row in matrix for value in row]

        columns = it.orthonormalize([dense(rows) for rows in pencil])
        for label, candidate, expected_inside in (
            ("uniform_detach_only", detach_only, True),
            ("single_point_join", single, False),
        ):
            vector = dense(candidate)
            residual = list(vector)
            for column in columns:
                overlap = it._dot(column, residual)
                for position in range(len(residual)):
                    residual[position] -= overlap * column[position]
            relative = it._norm(residual) / it._norm(vector)
            if expected_inside:
                self.assertLess(relative, 1e-12, label)
            else:
                self.assertGreater(relative, 1e-3, label)


class ExactControls(unittest.TestCase):
    def test_the_response_machinery_survives_four_independent_checks(self) -> None:
        """Stops a truncated series or an unconverged eigenvector being read as physics.

        Probability leaking, a Poisson horizon cut too early, a stationary power
        iteration stopped short, or the affine-transport claim being a coding
        error rather than a theorem would each produce plausible-looking but
        wrong transport errors.
        """

        generator = it.Generator(4)
        for intervention in it.INTERVENTIONS:
            control = it.exact_controls(generator, it.LAGS, intervention)
            self.assertLess(
                control["constant_function_is_conserved_max_drift"], 1e-10, intervention
            )
            self.assertLess(
                control["uniformization_matches_dense_taylor_max_drift"], 1e-9, intervention
            )
            self.assertLess(
                control["power_iteration_matches_exact_rational_stationary_max_drift"],
                1e-9,
                intervention,
            )
            self.assertLess(
                control["full_span_transport_is_exact_max_drift"], 1e-9, intervention
            )
            self.assertTrue(control["passed"], intervention)

    def test_galerkin_compression_of_an_affine_family_is_affine(self) -> None:
        """Stops us reporting a Contract II 'result' that is a theorem, not a finding.

        For a frozen span, ``Phi^T G_eta Phi = A_0 + eta B`` identically.  A
        transport gap can therefore only come from ``A_eta`` being allowed to
        see the target -- which is exactly what separates ``M_fixed_span`` from
        ``M_transport`` here.  If this drifted, the ladder would be measuring
        arithmetic noise.
        """

        generator = it.Generator(4)
        basis = it.deterministic_span(generator.size, 5)
        a0 = it.reduced_operator(basis, generator.rows(generator.baseline_rates()))
        for intervention in it.INTERVENTIONS:
            tangent = it.reduced_operator(basis, generator.tangent_rows(intervention))
            for eta in it.ETAS:
                galerkin = it.reduced_operator(
                    basis, generator.rows(generator.rates(intervention, eta))
                )
                drift = max(
                    abs(galerkin[i][j] - (a0[i][j] + eta * tangent[i][j]))
                    for i in range(len(basis))
                    for j in range(len(basis))
                )
                self.assertLess(drift, 1e-11, (intervention, eta))


class NumericalPrimitives(unittest.TestCase):
    def test_the_krylov_basis_does_not_overflow_at_the_ranks_we_use(self) -> None:
        """Stops a NaN span being scored as a perfect or a hopeless model.

        Raw Krylov powers grow like ``||G||^k``.  Before this was normalized the
        rank-12 span at width 4 came back with NaN columns, and the reported
        errors for that rank were NaN rather than a number anyone would notice.
        """

        generator = it.Generator(5)
        observables = [
            [function(state) for state in generator.states]
            for _, function in it.PRIMARY_READOUTS
        ]
        seeds = [[1.0] * generator.size] + observables
        basis = it.krylov_span(
            generator.rows(generator.baseline_rates()), seeds, 16
        )
        for column in basis:
            self.assertTrue(all(math.isfinite(value) for value in column))
            self.assertAlmostEqual(it._norm(column), 1.0, places=10)
        for i, left in enumerate(basis):
            for right in basis[i + 1 :]:
                self.assertLess(abs(it._dot(left, right)), 1e-9)

    def test_span_drift_separates_a_moved_span_from_a_frozen_one(self) -> None:
        """Stops a zero refit/fixed-span gap being read as a stable state.

        A refit that returns the same subspace produces an identically zero gap
        for a trivial reason.  Without this measurement that tautology and a
        genuine intervention-stable state look the same in the error table.
        """

        left = it.deterministic_span(20, 4)
        self.assertLess(it.span_drift(left, left), 1e-12)
        right = it.deterministic_span(20, 4, seed=7)
        self.assertGreater(it.span_drift(left, right), 1e-3)

    def test_the_dense_exponential_matches_a_closed_form(self) -> None:
        """Stops a mis-scaled squaring loop turning every prediction into a rescaled one.

        ``[[0,1],[0,0]]`` exponentiates to ``[[1,t],[0,1]]`` exactly, and a
        rotation generator to a rotation: both catch a wrong number of squarings,
        which would show up as a uniformly wrong relaxation rate rather than as
        an obvious failure.
        """

        nilpotent = it.dense_expm([[0.0, 1.0], [0.0, 0.0]], 2.5)
        self.assertAlmostEqual(nilpotent[0][1], 2.5, places=12)
        self.assertAlmostEqual(nilpotent[0][0], 1.0, places=12)
        rotation = it.dense_expm([[0.0, -1.0], [1.0, 0.0]], 1.3)
        self.assertAlmostEqual(rotation[0][0], math.cos(1.3), places=11)
        self.assertAlmostEqual(rotation[1][0], math.sin(1.3), places=11)

    def test_the_jacobi_eigensolver_reproduces_its_input(self) -> None:
        """Stops a wrong pseudo-inverse making the realizable estimator look hopeless.

        The generator fit is a foil for the transported model.  If the Gram
        truncation were wrong, the foil would fail for a numerical reason and the
        transported model would look better than it is.
        """

        matrix = [[4.0, 1.0, 0.5], [1.0, 3.0, -0.25], [0.5, -0.25, 2.0]]
        values, columns = it.symmetric_eigen(matrix)
        rebuilt = [[0.0] * 3 for _ in range(3)]
        for value, column in zip(values, columns):
            for i in range(3):
                for j in range(3):
                    rebuilt[i][j] += value * column[i] * column[j]
        for i in range(3):
            for j in range(3):
                self.assertAlmostEqual(rebuilt[i][j], matrix[i][j], places=10)

    def test_the_contrast_tensor_removes_a_constant_offset(self) -> None:
        """Stops the score being dominated by a part every model gets for free.

        The raw response carries the observable's mean and the common long-time
        limit.  Both are reproduced by any span containing the constant
        function, so scoring them would compress ``M_refit``, ``M_fixed_span``,
        ``M_transport`` and even the random control toward each other.
        """

        tensor = [[[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]]
        shifted = [[[value + 10.0 for value in row] for row in block] for block in tensor]
        self.assertEqual(it.contrast_tensor(tensor, 2), it.contrast_tensor(shifted, 2))


class ModelLadder(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = it.assemble(
            widths=(4,),
            families=("krylov", "random_control"),
            ranks=(4, 6),
            etas=it.PRIMARY_ETAS,
            lags=(0.5, 1.0, 2.0),
            interventions=("uniform_join_minus_detach", "single_point_join"),
            control_width=4,
        )

    def test_the_in_span_oracle_is_never_beaten_in_the_metric_it_minimizes(self) -> None:
        """Stops a bookkeeping slip inverting the ladder that carries the reading.

        ``E_fixed_span`` is the orthogonal projection, so in the state-space
        metric it is an exact infimum over every reduced trajectory the frozen
        span admits -- transported, fitted or time-dependent.  If a model in the
        same span ever scored below it there, the two numbers would not be
        measuring the same span and #580's 'span survives but the generator law
        is incomplete' branch could never be reached honestly.

        The same inequality does *not* hold on the source-sampled response,
        because four sources cannot see every direction the projection
        optimizes.  Asserting it there would be asserting something false, which
        is why both metrics are carried.
        """

        for block in self.result["blocks"]:
            for entry in block["spans"]:
                if entry.get("status") != "scored":
                    continue
                for key, row in entry["etas"].items():
                    where = (
                        block["intervention"],
                        entry["family"],
                        entry["requested_rank"],
                        key,
                    )
                    for label in ("transported", "fixed_span_generator_fit"):
                        self.assertLessEqual(
                            row["fixed_span"]["state_space_primary"],
                            row[label]["state_space_primary"] + 1e-9,
                            where + (label,),
                        )

    def test_the_random_span_control_is_far_worse_than_the_declared_span(self) -> None:
        """Stops a loose metric being read as evidence that a state transports.

        If an information-free span reached the same error, the pass threshold
        would be measuring the observables' trivial structure rather than the
        span.  The declared control margin exists to make that failure visible.
        """

        for block in self.result["blocks"]:
            krylov = it._find_span(block, "krylov", 6)
            control = it._find_span(block, "random_control", 6)
            self.assertIsNotNone(krylov)
            self.assertIsNotNone(control)
            best = it._worst_primary(krylov, "transported", "contrast_primary")
            worst = it._worst_primary(control, "transported", "contrast_primary")
            self.assertGreater(worst, it.CONTROL_MARGIN * best, block["intervention"])

    def test_the_decision_records_which_side_of_the_pencil_each_number_came_from(self) -> None:
        """Stops a transport pass being read as a physical state without its attribution.

        The declared intervention lies inside the same operator pencil that
        builds the state space, so a pass there is compatible with an entirely
        algebraic explanation.  The artifact has to carry the out-of-pencil
        comparison next to it or the verdict is unreadable.
        """

        pencil = self.result["decision"]["pencil_attribution"]
        self.assertIn("reading", pencil)
        self.assertIsNotNone(pencil["worst_in_pencil_transport_error"])
        self.assertIsNotNone(pencil["worst_out_of_pencil_transport_error"])
        for reading in self.result["decision"]["per_family"]:
            self.assertIn("intervention", reading)
            self.assertIn("intervention_is_inside_the_baseline_pencil", reading)


class FlagInvariance(unittest.TestCase):
    def test_the_pencil_interventions_keep_the_readouts_inside_the_baseline_flag(self) -> None:
        """Stops a transport pass being credited to physics when it is algebra.

        ``G_0 = J + D`` and the declared intervention ``J - D`` span the same
        plane, so one step of that intervention lands the declared readouts back
        inside the Krylov flag the baseline already builds -- a frozen span then
        transports it for a reason that has nothing to do with the state being
        physical.  A spatially inhomogeneous tilt does leave the flag.  Without
        this number the two cases are indistinguishable in the error table.
        """

        generator = it.Generator(5)
        seeds = [[1.0] * generator.size] + [
            [function(state) for state in generator.states]
            for _, function in it.PRIMARY_READOUTS
        ]
        rows = generator.rows(generator.baseline_rates())
        inside = it.flag_invariance(
            rows, generator.tangent_rows("uniform_join_minus_detach"), seeds, (4,)
        )
        detach = it.flag_invariance(
            rows, generator.tangent_rows("uniform_detach_only"), seeds, (4,)
        )
        outside = it.flag_invariance(
            rows, generator.tangent_rows("single_point_join"), seeds, (4,)
        )
        self.assertTrue(inside["declared_readouts_stay_inside_the_baseline_flag"])
        self.assertTrue(detach["declared_readouts_stay_inside_the_baseline_flag"])
        self.assertFalse(outside["declared_readouts_stay_inside_the_baseline_flag"])
        self.assertGreater(outside["at_the_declared_readouts"], 0.1)


class DeclaredManifest(unittest.TestCase):
    def test_the_held_out_readouts_never_touch_the_span(self) -> None:
        """Stops a held-out prediction that quietly used the held-out functions.

        The span is seeded by the constant function and the three primary
        readouts only.  If the two blocks overlapped, ``contrast_held_out``
        would be an in-sample number wearing an out-of-sample label.
        """

        primary = {name for name, _ in it.PRIMARY_READOUTS}
        held_out = {name for name, _ in it.HELD_OUT_READOUTS}
        self.assertEqual(primary & held_out, set())
        self.assertTrue(held_out)

    def test_the_pencil_labels_match_the_intervention_table(self) -> None:
        """Stops the attribution reading being computed from a stale label.

        ``IN_PENCIL`` is what the decision uses to split the interventions.  If
        it disagreed with ``INTERVENTIONS`` the artifact could report a pencil
        contrast that compares an intervention with itself.
        """

        self.assertEqual(set(it.INTERVENTIONS), set(it.IN_PENCIL))
        self.assertTrue(it.IN_PENCIL["uniform_join_minus_detach"])
        self.assertFalse(it.IN_PENCIL["single_point_join"])


if __name__ == "__main__":
    unittest.main()


class RankNotions(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.generator = it.Generator(5)
        cls.observables = {
            name: [function(state) for state in cls.generator.states]
            for name, function in list(it.PRIMARY_READOUTS) + list(it.HELD_OUT_READOUTS)
        }

    def test_the_lumped_chain_reproduces_the_full_one(self) -> None:
        """Stops a refinement bug being reported as a small positive realization.

        The block count is only a certified positive/Markov realization dimension
        if the aggregated chain actually reproduces the full chain's evolution of
        a function constant on blocks.  A signature bug would over-merge and make
        ``r_positive`` look smaller than it is -- the direction that most
        flatters the conclusion this artifact draws.
        """

        rates = self.generator.baseline_rates()
        for name in ("blocks", "covering_depth"):
            colours = [(value,) for value in self.observables[name]]
            result = it.exact_lumping(
                self.generator, colours, [rates], verify_with=self.observables[name]
            )
            self.assertLess(result["verified_max_response_drift"], 1e-10, name)
            self.assertLessEqual(result["blocks"], self.generator.size)

    def test_a_finer_colouring_can_only_give_a_finer_lumping(self) -> None:
        """Stops the dictionary ladder reporting a positive dimension that shrinks.

        ``D0`` is a subset of ``D2``, so the coarsest lumping keeping ``D2``
        refines the one keeping ``D0``.  If the block count ever fell as the
        dictionary grew, the refinement would not be computing what the artifact
        says it computes.
        """

        rates = self.generator.baseline_rates()
        counts = []
        for members in (
            it.DICTIONARIES["D0_additive_local_counts"],
            it.DICTIONARIES["D2_plus_nonlocal_topology"],
        ):
            colours = [
                tuple(self.observables[member][state] for member in members)
                for state in range(self.generator.size)
            ]
            counts.append(it.exact_lumping(self.generator, colours, [rates])["blocks"])
        self.assertLessEqual(counts[0], counts[1])

    def test_refining_against_the_intervention_cannot_coarsen_the_lumping(self) -> None:
        """Stops an intervention-stable positive realization being claimed for free.

        Asking the same partition to lump the baseline generator *and* the
        intervention tangent is a strictly stronger requirement, so it can only
        split blocks.  If it ever merged them, the "positive realization survives
        the whole affine family" line would be meaningless.
        """

        members = it.DICTIONARIES["D0_additive_local_counts"]
        colours = [
            tuple(self.observables[member][state] for member in members)
            for state in range(self.generator.size)
        ]
        rates = self.generator.baseline_rates()
        alone = it.exact_lumping(self.generator, colours, [rates])["blocks"]
        for intervention in it.INTERVENTIONS:
            joint = it.exact_lumping(
                self.generator,
                colours,
                [rates, self.generator.coefficients(intervention)],
            )["blocks"]
            self.assertGreaterEqual(joint, alone, intervention)

    def test_the_modular_rank_agrees_with_a_second_prime_and_beats_floating_point(self) -> None:
        """Stops ``r_linear`` being a number produced by a pivot tolerance.

        Repeated application of the generator collapses onto the dominant
        direction, so a float elimination loses dimensions: at width 5 it
        returned 26 against a true 42, and that error runs in the direction that
        would make the linear rank look close to the transported rank and hide
        the whole separation.  Two independent primes agreeing is the cheap
        check that the modular rank is the rational one.
        """

        vectors = [
            self.observables[name]
            for name in it.DICTIONARIES["D0_additive_local_counts"]
        ]
        first = it.observable_reachable_dimension(self.generator, vectors, 200)
        second = it.observable_reachable_dimension(
            self.generator, vectors, 200, modulus=2147483629
        )
        self.assertEqual(first["dimension"], second["dimension"])
        self.assertFalse(first["limited_by_budget"])
        self.assertLess(first["dimension"], self.generator.size)
        self.assertGreater(first["dimension"], 3 * max(it.RANKS[:3]))


    def test_the_linear_rank_never_exceeds_the_positive_one(self) -> None:
        """Stops the three rank notions being reported in an impossible order.

        The lumped chain is itself a realization, and every function in the
        observable Krylov space is constant on its blocks, so the linear rank
        cannot exceed the block count.  If the artifact ever printed the reverse,
        one of the two would be measuring something other than what its name
        says, and the whole `r_linear` / `r_positive` / `r_transport` separation
        would be unreadable.
        """

        rates = self.generator.baseline_rates()
        for members in it.DICTIONARIES.values():
            vectors = [self.observables[member] for member in members]
            colours = [
                tuple(self.observables[member][state] for member in members)
                for state in range(self.generator.size)
            ]
            linear = it.observable_reachable_dimension(
                self.generator, vectors, self.generator.size
            )["dimension"]
            positive = it.exact_lumping(self.generator, colours, [rates])["blocks"]
            self.assertLessEqual(linear, positive, members)


class ReadoutClassification(unittest.TestCase):
    def test_an_unrepresented_readout_gets_no_transport_verdict(self) -> None:
        """Stops a dictionary failure being counted as a transport result.

        A readout the frozen span never represented at ``eta = 0`` can show a
        tiny finite-``eta`` excess purely because both numbers are large and
        close.  Binning it as 'transports' would be the most flattering possible
        misreading of the whole artifact.
        """

        classification = it.classify_readouts(
            ("good", "bad", "faint"),
            ("good",),
            [[[1.0, 1.0, 1e-9], [2.0, 2.0, 2e-9]]],
            {"good": 0.01, "bad": 0.60, "faint": 0.01},
            {"good": 0.02, "bad": 0.61, "faint": 0.01},
            {"good": 0.005, "bad": 0.004, "faint": 0.0},
        )
        rows = classification["by_readout"]
        self.assertEqual(rows["good"]["bin"], "represented_and_transports")
        self.assertEqual(
            rows["bad"]["bin"],
            "unrepresented_at_baseline_transport_not_identifiable",
        )
        self.assertEqual(rows["faint"]["bin"], "weakly_identifiable_baseline_signal")

    def test_the_balanced_score_ignores_a_rescaling(self) -> None:
        """Stops a coordinate scale deciding whether the state transports.

        The declared pooled Frobenius norm is magnitude-weighted: doubling one
        readout changes the verdict it carries.  The balanced score gives each
        readout its own denominator, so it must be invariant under rescaling any
        single readout -- which the pooled one is not.
        """

        truth = [[[1.0, 10.0], [2.0, 20.0]]]
        prediction = [[[1.1, 12.0], [2.2, 24.0]]]
        scaled_truth = [[[1.0, 1000.0], [2.0, 2000.0]]]
        scaled_prediction = [[[1.1, 1200.0], [2.2, 2400.0]]]
        self.assertAlmostEqual(
            it.balanced_error(truth, prediction, (0, 1)),
            it.balanced_error(scaled_truth, scaled_prediction, (0, 1)),
            places=12,
        )
        self.assertNotAlmostEqual(
            it.relative_error(truth, prediction, (0, 1)),
            it.balanced_error(truth, prediction, (0, 1)),
            places=3,
        )
