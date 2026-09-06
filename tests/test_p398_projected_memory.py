"""Tests for the exact projected-memory machinery of #588 Phase A.

The failure mode this file is built against is a kernel routine that returns
something smooth, plausible and wrong.  ``K(t)`` has no independent measurement
to be checked against in the wild, so every test here names a specific wrong
number it would stop us believing.
"""

from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from scripts.p398_intervention_transport import (  # noqa: E402
    Generator,
    RANKS,
    _dot,
    _norm,
    _prefix,
    build_span,
    krylov_span,
    matvec,
    orthonormalize,
)
from scripts.p398_projected_memory import (  # noqa: E402
    FROZEN_RANK,
    KERNEL_GRID,
    KERNEL_STEP,
    LAGS,
    MEMORY_ETAS,
    MEMORY_INTERVENTIONS,
    INTERVENTIONS,
    adapted_blocks,
    block_hankel,
    dense_kernel,
    dense_rows,
    declared_observables,
    frobenius,
    frozen_span,
    hankel_subspace,
    integrate_gle,
    markov_only_trajectory,
    memory_kernel,
    project_out,
    relative_kernel_distance,
    resolved_operator,
    resolvent_identity,
    singular_spectrum,
    subspace_residual,
    synthetic_control,
    unresolved_columns,
    unresolved_forcing,
)


class ProjectionPrimitives(unittest.TestCase):
    def setUp(self) -> None:
        self.generator = Generator(4)
        self.rates = self.generator.baseline_rates()
        self.rows = self.generator.rows(self.rates)
        self.rate = self.generator.exit_rate(self.rates)
        self.basis = frozen_span(self.generator, self.rows)

    def test_the_complement_is_orthogonal_to_the_span(self) -> None:
        """Stops us believing a memory kernel that is really span leakage.

        ``Q`` is the only thing standing between "this is the discarded
        dynamics" and "this is the resolved dynamics counted twice".  A single
        classical Gram-Schmidt pass loses orthogonality at the condition number
        of the basis, and the leak would appear as a perfectly smooth kernel.
        """

        vector = [float(index % 7) - 3.0 for index in range(self.generator.size)]
        residual = project_out(self.basis, vector)
        for column in self.basis:
            self.assertLess(abs(_dot(column, residual)), 1e-14)

    def test_a_full_rank_projection_has_no_memory_at_all(self) -> None:
        """Stops us believing a nonzero kernel where the complement is empty.

        With ``P = I`` there is nothing unresolved, so ``K`` must be identically
        zero -- not small, zero.  Anything else is arithmetic in the kernel
        routine itself.
        """

        size = self.generator.size
        full = orthonormalize(
            [[1.0 if i == j else 0.0 for j in range(size)] for i in range(size)]
        )
        kernels = memory_kernel(full, self.rows, size, self.rate, (0.5, 1.0))
        for block in kernels:
            self.assertEqual(frobenius(block), 0.0)

    def test_a_generator_closed_span_has_no_memory(self) -> None:
        """Stops us believing memory that is really rank truncation.

        The saturated observable-Krylov space is exactly invariant under ``G``,
        so ``Q G Phi = 0`` and the kernel vanishes even though the span is a
        proper subspace.  A routine that mistook truncation for memory would
        return something of order ``||A||`` here.
        """

        size = self.generator.size
        observables, _, _ = declared_observables(self.generator)
        seeds = [[1.0] * size] + observables[:3]
        closed = krylov_span(self.rows, seeds, size)
        self.assertLess(len(closed), size)
        kernels = memory_kernel(closed, self.rows, size, self.rate, (0.5, 1.0, 2.0))
        reference = frobenius(resolved_operator(closed, self.rows))
        for block in kernels:
            self.assertLess(frobenius(block), 1e-8 * reference)


class BlockDecomposition(unittest.TestCase):
    def setUp(self) -> None:
        self.generator = Generator(4)
        self.rows = self.generator.rows(self.generator.baseline_rates())
        self.rate = self.generator.exit_rate(self.generator.baseline_rates())
        self.basis = frozen_span(self.generator, self.rows)

    def test_the_adapted_basis_reproduces_the_generator(self) -> None:
        """Stops us believing blocks read off a basis that is not a basis.

        If the complement missed a direction, ``[[A,B],[C,D]]`` would still be a
        well-formed matrix and would still produce a kernel; it would just be
        the kernel of a different operator.
        """

        size = self.generator.size
        resolved, out, into, unresolved = adapted_blocks(self.basis, self.rows, size)
        self.assertEqual(len(unresolved), size - len(self.basis))
        total = (
            sum(value * value for row in resolved for value in row)
            + sum(value * value for row in out for value in row)
            + sum(value * value for row in into for value in row)
            + sum(value * value for row in unresolved for value in row)
        )
        images = [matvec(self.rows, [1.0 if i == j else 0.0 for i in range(size)])
                  for j in range(size)]
        reference = sum(value * value for image in images for value in image)
        self.assertLess(abs(total - reference) / reference, 1e-12)

    def test_the_schur_complement_identity_holds(self) -> None:
        """Stops us believing a transposed B and C.

        Swapping ``B`` and ``C`` leaves every norm in this file unchanged and
        still yields a decaying kernel.  The resolvent identity does not
        survive it, and it involves no time integration at all, so a failure
        here cannot be blamed on a quadrature rule.
        """

        worst = resolvent_identity(self.basis, self.rows, self.generator.size, (0.5, 2.0))
        self.assertLess(worst, 1e-10)


class MemoryKernel(unittest.TestCase):
    def setUp(self) -> None:
        self.generator = Generator(4)
        self.rows = self.generator.rows(self.generator.baseline_rates())
        self.rate = self.generator.exit_rate(self.generator.baseline_rates())
        self.basis = frozen_span(self.generator, self.rows)

    def test_matrix_free_matches_an_explicitly_built_unresolved_block(self) -> None:
        """Stops us believing a kernel from a propagation that left range(Q).

        ``D`` is never formed in the production path -- it is 1424 x 1424 at
        width 8 -- so the only way to know the constrained propagation is
        ``exp(tau D)`` is to build ``D`` at a width where that is possible and
        compare.
        """

        lags = (0.25, 1.0, 2.5)
        _, out, into, unresolved = adapted_blocks(self.basis, self.rows, self.generator.size)
        expected = dense_kernel(out, unresolved, into, lags)
        measured = memory_kernel(
            self.basis, self.rows, self.generator.size, self.rate, lags
        )
        self.assertLess(relative_kernel_distance(measured, expected), 1e-9)

    def test_a_planted_block_is_recovered_exactly(self) -> None:
        """Stops us believing a routine that silently needs a generator.

        The planted matrix has no vanishing row sums and the uniformization rate
        it is given is not an exit rate of anything.  If either had been assumed
        rather than used, this control would be the one that shows it.
        """

        self.assertLess(synthetic_control(), 1e-10)

    def test_the_kernel_at_zero_lag_is_the_product_of_the_couplings(self) -> None:
        """Stops us believing an off-by-one in the propagated ladder.

        ``K(0) = B C`` with no propagator in between.  A ladder that recorded
        one step early or late would be smooth, plausible, and wrong by exactly
        one grid step -- which is invisible in every norm but this one.
        """

        size = self.generator.size
        kernels = memory_kernel(self.basis, self.rows, size, self.rate, (0.0, 0.25))
        columns = unresolved_columns(self.basis, self.rows)
        direct = [
            [
                _dot(self.basis[i], matvec(self.rows, columns[j]))
                for j in range(len(self.basis))
            ]
            for i in range(len(self.basis))
        ]
        difference = relative_kernel_distance([kernels[0]], [direct])
        self.assertLess(difference, 1e-12)


class UnresolvedForcing(unittest.TestCase):
    def setUp(self) -> None:
        self.generator = Generator(5)
        self.rows = self.generator.rows(self.generator.baseline_rates())
        self.rate = self.generator.exit_rate(self.generator.baseline_rates())
        self.basis = frozen_span(self.generator, self.rows)

    def test_a_span_seeding_readout_has_no_unresolved_initial_condition(self) -> None:
        """Stops us believing the two closure defects are the same defect.

        Every declared readout seeds the Krylov span, so ``Q f = 0`` and the
        inhomogeneous term is identically zero for it.  If that came back
        nonzero, the memory share and the unrepresented-readout share of the
        closure error would be mixing, and #588's whole question is which of the
        two the held-out failure is.
        """

        observables, names, _ = declared_observables(self.generator)
        forcing = unresolved_forcing(
            self.basis, self.rows, self.generator.size, self.rate, observables[:3], (0.5, 1.0)
        )
        for block in forcing:
            for vector in block:
                self.assertLess(_norm(vector), 1e-12)

    def test_a_held_out_readout_does_have_one(self) -> None:
        """Stops us believing a forcing term that is zero because of a bug.

        The previous test passes trivially if ``unresolved_forcing`` always
        returns zero.  A held-out readout is not in the span, so its term must
        be genuinely nonzero.
        """

        observables, _, held_out = declared_observables(self.generator)
        seeds = [observables[index] for index in held_out]
        forcing = unresolved_forcing(
            self.basis, self.rows, self.generator.size, self.rate, seeds, (0.5,)
        )
        self.assertGreater(max(_norm(vector) for vector in forcing[0]), 1e-3)


class GeneralizedLangevin(unittest.TestCase):
    def setUp(self) -> None:
        self.generator = Generator(4)
        self.rows = self.generator.rows(self.generator.baseline_rates())
        self.rate = self.generator.exit_rate(self.generator.baseline_rates())
        self.basis = frozen_span(self.generator, self.rows)
        self.resolved = resolved_operator(self.basis, self.rows)

    def _trajectories(self, substeps: int):
        from scripts.p398_intervention_transport import evolve_observables
        from scripts.p398_projected_memory import project_coordinates

        size = self.generator.size
        step = KERNEL_STEP / substeps
        grid = [step * index for index in range(8 * substeps + 1)]
        observables, names, held_out = declared_observables(self.generator)
        column = held_out[0]
        kernels = memory_kernel(self.basis, self.rows, size, self.rate, grid)
        forcing_blocks = unresolved_forcing(
            self.basis, self.rows, size, self.rate, [observables[column]], grid
        )
        evolved = evolve_observables(
            self.rows, size, self.rate, [observables[column]], grid
        )
        exact = [project_coordinates(self.basis, block[0]) for block in evolved]
        forcing = [block[0] for block in forcing_blocks]
        solved = integrate_gle(self.resolved, kernels, forcing, exact[0], step)
        reference = math.sqrt(sum(_dot(vector, vector) for vector in exact))
        error = math.sqrt(
            sum(
                sum((a - b) ** 2 for a, b in zip(left, right))
                for left, right in zip(solved, exact)
            )
        ) / reference
        return error, grid, exact

    def test_exact_memory_and_forcing_reproduce_the_exact_trajectory(self) -> None:
        """Stops us believing a kernel that is merely of the right size.

        The projected equation is an identity, not an approximation: with the
        exact ``K`` and the exact inhomogeneous term, the only thing left
        between the reduced solve and ``Phi^T exp(tG) f`` is the quadrature.
        A kernel that were, say, scaled by two would still decay and would still
        look like memory, but it would not close this.
        """

        error, _, _ = self._trajectories(4)
        self.assertLess(error, 5e-3)

    def test_the_remaining_error_is_quadrature(self) -> None:
        """Stops us believing a wrong kernel that happens to be close.

        The scheme is second order, so halving the step must divide the residual
        by about four.  A residual that stalls instead is a residual that is not
        quadrature, which is exactly what a slightly wrong ``B``, ``C`` or ``D``
        would produce.
        """

        coarse, _, _ = self._trajectories(2)
        fine, _, _ = self._trajectories(4)
        self.assertGreater(coarse / fine, 3.0)

    def test_dropping_both_terms_is_the_markov_closure(self) -> None:
        """Stops us believing a memory share inflated by the integrator.

        With a zero kernel and zero forcing the projected equation is just
        ``x' = A x``, so the integrator must reproduce ``exp(tA)``.  If it did
        not, every reported share of the closure error would carry the
        integrator's own discretization inside it.
        """

        from scripts.p398_intervention_transport import evolve_observables
        from scripts.p398_projected_memory import project_coordinates

        substeps = 8
        step = KERNEL_STEP / substeps
        grid = [step * index for index in range(4 * substeps + 1)]
        rank = len(self.basis)
        observables, _, _ = declared_observables(self.generator)
        start = project_coordinates(self.basis, observables[0])
        zero_kernels = [[[0.0] * rank for _ in range(rank)] for _ in grid]
        zero_forcing = [[0.0] * rank for _ in grid]
        integrated = integrate_gle(self.resolved, zero_kernels, zero_forcing, start, step)
        exponential = markov_only_trajectory(self.resolved, start, grid)
        worst = max(
            abs(a - b)
            for left, right in zip(integrated, exponential)
            for a, b in zip(left, right)
        )
        scale = max(abs(value) for vector in exponential for value in vector)
        self.assertLess(worst / scale, 1e-3)


class KernelStatistics(unittest.TestCase):
    def test_a_single_exponential_has_hankel_rank_one(self) -> None:
        """Stops us believing an effective pole count that counts noise.

        A one-pole kernel must produce a rank-one block Hankel matrix.  If the
        construction indexed the samples wrongly the matrix would still be
        square and would still have a spectrum -- just not this one.

        The bound is ``1e-7`` and not machine epsilon on purpose: the singular
        values come from an eigendecomposition of ``H^T H``, so the floor is the
        square root of epsilon.  That floor is why the declared numerical-rank
        tolerance is ``1e-6``.
        """

        kernels = [[[3.0 * math.exp(-1.7 * KERNEL_STEP * index)]] for index in range(17)]
        spectrum = singular_spectrum(block_hankel(kernels))
        self.assertGreater(spectrum[0], 0.0)
        self.assertLess(spectrum[1] / spectrum[0], 1e-7)

    def test_a_matrix_lies_in_its_own_leading_subspace(self) -> None:
        """Stops us believing a transport residual that is always small.

        ``subspace_residual`` is the A3 statistic for whether the baseline
        memory directions still carry the intervened kernel.  If it returned
        near zero for everything, an intervention that changed the memory
        spectrum in kind would read as one that transported.
        """

        kernels = [
            [[math.exp(-1.1 * KERNEL_STEP * index), 0.0],
             [0.0, math.exp(-2.9 * KERNEL_STEP * index)]]
            for index in range(17)
        ]
        matrix = block_hankel(kernels)
        subspace = hankel_subspace(matrix, 2)
        self.assertLess(subspace_residual(subspace, matrix), 1e-10)
        other = block_hankel(
            [
                [[math.exp(-0.2 * KERNEL_STEP * index), 0.0], [0.0, 0.0]]
                for index in range(17)
            ]
        )
        self.assertGreater(subspace_residual(hankel_subspace(other, 1), matrix), 1e-3)


class FrozenConfiguration(unittest.TestCase):
    def test_the_span_is_the_one_580_froze(self) -> None:
        """Stops us believing a memory measured on a different span.

        #588 is only meaningful on #580's basis.  A basis rebuilt with a
        different rank ladder, or built directly at rank 6 rather than prefixed
        from rank 12, would be a *different subspace* and would carry a
        different kernel -- and nothing downstream would notice.
        """

        for width in (4, 5, 6):
            generator = Generator(width)
            rows = generator.rows(generator.baseline_rates())
            observables, _, _ = declared_observables(generator)
            seeds = [[1.0] * generator.size] + observables[:3]
            reference = _prefix(
                build_span(
                    "krylov",
                    generator,
                    rows,
                    seeds,
                    min(max(RANKS), generator.size),
                    generator.exit_rate(generator.baseline_rates()),
                ),
                FROZEN_RANK,
            )
            measured = frozen_span(generator, rows)
            self.assertEqual(len(measured), len(reference))
            for left, right in zip(measured, reference):
                self.assertLess(max(abs(a - b) for a, b in zip(left, right)), 1e-14)


class DeclaredManifest(unittest.TestCase):
    def test_every_declared_lag_is_a_kernel_grid_point(self) -> None:
        """Stops us believing a kernel norm read off an interpolated lag.

        The reported lag-grid quantities are read from the uniform kernel
        samples.  A declared lag that was not on the grid would have to be
        interpolated, and the norm profile would then be partly a spline.
        """

        for lag in LAGS:
            self.assertTrue(
                any(abs(lag - point) < 1e-12 for point in KERNEL_GRID),
                f"declared lag {lag} is not a kernel grid point",
            )

    def test_the_baseline_is_in_the_declared_eta_ladder(self) -> None:
        """Stops us believing a transported number with no baseline.

        Every quantity here is an excess over ``eta = 0``.  A ladder without it
        would silently compare the first intervention against nothing.
        """

        self.assertIn(0.0, MEMORY_ETAS)

    def test_both_declared_interventions_exist_and_split_the_pencil(self) -> None:
        """Stops us believing an out-of-pencil control that is in the pencil.

        The whole force of A3 is that one of the two interventions is outside
        the operator pencil that built the span.  If both were inside, the
        transport statistic would be measuring an algebraic identity.
        """

        from scripts.p398_projected_memory import IN_PENCIL

        for name in MEMORY_INTERVENTIONS:
            self.assertIn(name, INTERVENTIONS)
        pencil = {IN_PENCIL[name] for name in MEMORY_INTERVENTIONS}
        self.assertEqual(pencil, {True, False})


if __name__ == "__main__":
    unittest.main()
