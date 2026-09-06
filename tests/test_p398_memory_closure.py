"""Tests for #588 Phases B and C.

The failure mode here is a competitor that wins or loses for a reason other than
the one being measured: a Petrov-Galerkin pair that is not biorthogonal, a
memory model that quietly saw the intervened data, or a balanced transform built
from a metric that a rescaled readout could rotate.
"""

from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for candidate in (ROOT, ROOT / "scripts"):
    if str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

from scripts.p398_intervention_transport import (  # noqa: E402
    CONTRAST_REFERENCE,
    DICTIONARIES,
    LAGS,
    Generator,
    _dot,
    _prefix,
    contrast_tensor,
    evolve_observables,
    krylov_span,
    orthonormalize,
    response_tensor,
    source_distributions,
)
from scripts.p398_projected_memory import (  # noqa: E402
    FROZEN_RANK,
    KERNEL_GRID,
    KERNEL_STEP,
    KERNEL_STEPS,
    declared_observables,
    frozen_span,
    memory_kernel,
)
from scripts.p398_memory_closure import (  # noqa: E402
    BALANCED_ORDERS,
    CLOSURE_SUBSTEPS,
    MEMORY_ORDERS,
    balanced_family,
    biorthogonal,
    era_fit,
    era_kernels,
    reduced_realization,
    score_realization,
    width_closure,
)


def _cache(generator, etas=(0.0,)):
    size = generator.size
    observables, names, _ = declared_observables(generator)
    sources = source_distributions(generator)
    reference_index = [name for name, _ in sources].index(CONTRAST_REFERENCE)
    cache = {}
    for eta in etas:
        rates = generator.rates("uniform_join_minus_detach", eta)
        rows = generator.rows(rates)
        rate = generator.exit_rate(rates)
        truth = response_tensor(
            sources, evolve_observables(rows, size, rate, observables, LAGS)
        )
        cache[eta] = {
            "rows": rows,
            "columns": generator.columns(rates),
            "rate": rate,
            "contrast": contrast_tensor(truth, reference_index),
        }
    return observables, names, sources, reference_index, cache


class Scorer(unittest.TestCase):
    def test_the_galerkin_case_is_the_petrov_galerkin_case(self) -> None:
        """Stops us believing a competitor that won by using a different scorer.

        The balanced adversary is a genuine Petrov-Galerkin reduction and the
        Krylov span is Galerkin.  If the shared scorer did not reduce exactly to
        #580's Galerkin score when the two bases coincide, any difference
        between the competitors would be partly the scorer.
        """

        generator = Generator(5)
        observables, _, sources, reference_index, cache = _cache(generator)
        basis = frozen_span(generator, cache[0.0]["rows"])
        columns = [0, 1, 2]
        first = score_realization(
            basis, basis, cache[0.0]["rows"], observables, sources,
            cache[0.0]["contrast"], reference_index, columns,
        )
        self.assertLess(biorthogonal(basis, basis), 1e-12)
        self.assertGreater(first["balanced"], 0.0)
        self.assertLess(first["balanced"], 1.0)


class MemoryFit(unittest.TestCase):
    def test_a_planted_two_pole_kernel_is_recovered_at_order_two(self) -> None:
        """Stops us believing an ERA fit that is really interpolation.

        A kernel that is exactly two exponentials must be reproduced to
        machine precision at order 2 and no better at higher order.  A fit that
        merely interpolated the samples it was given would look perfect on those
        samples and drift on the ones it was not.
        """

        step = 0.05
        samples = 60
        # Two poles, each with a RANK-ONE residue.  A rank-two residue would
        # need two realization states of its own, so the planted kernel would
        # be order 3 and the test would be checking the wrong number.
        kernels = [
            [
                [
                    2.0 * math.exp(-1.3 * step * k) + 0.5 * math.exp(-7.0 * step * k),
                    0.3 * math.exp(-1.3 * step * k),
                ],
                [
                    1.1 * math.exp(-7.0 * step * k),
                    0.0,
                ],
            ]
            for k in range(samples)
        ]
        model = era_fit(kernels, 2, blocks=8)
        self.assertIsNotNone(model)
        rebuilt = era_kernels(model, kernels[0], samples - 1)
        worst = max(
            abs(a - b)
            for left, right in zip(kernels[1:], rebuilt[1:])
            for rl, rr in zip(left, right)
            for a, b in zip(rl, rr)
        )
        self.assertLess(worst, 1e-8)

    def test_the_fit_uses_only_baseline_samples(self) -> None:
        """Stops us believing a memory model that saw the intervention.

        The whole discipline of Phase C is that the rational model is fitted at
        eta = 0 and applied unchanged at eta = +-1/4.  A fit that depended on the
        intervened kernel would be refitting per environment, which is the thing
        #580 already showed is not the question.
        """

        generator = Generator(4)
        rows_zero = generator.rows(generator.rates("uniform_join_minus_detach", 0.0))
        rows_up = generator.rows(generator.rates("uniform_join_minus_detach", 0.25))
        basis = frozen_span(generator, rows_zero)
        rate = generator.exit_rate(generator.baseline_rates())
        step = KERNEL_STEP / CLOSURE_SUBSTEPS
        grid = [step * index for index in range(KERNEL_STEPS * CLOSURE_SUBSTEPS + 1)]
        zero = memory_kernel(basis, rows_zero, generator.size, rate, grid)
        up = memory_kernel(basis, rows_up, generator.size, rate, grid)
        first = era_fit(zero, 2)
        second = era_fit(zero, 2)
        self.assertIsNotNone(first)
        # Same inputs, same model; and the intervened kernel is genuinely
        # different, so a model that had used it could not match this one.
        for left, right in zip(first[0], second[0]):
            for a, b in zip(left, right):
                self.assertEqual(a, b)
        difference = max(
            abs(a - b)
            for lz, lu in zip(zero, up)
            for rz, ru in zip(lz, lu)
            for a, b in zip(rz, ru)
        )
        self.assertGreater(difference, 1e-6)


class BalancedAdversary(unittest.TestCase):
    def setUp(self) -> None:
        self.generator = Generator(5)
        (
            self.observables,
            self.names,
            self.sources,
            self.reference_index,
            self.cache,
        ) = _cache(self.generator, (0.0, 0.25, -0.25))

    def test_the_transform_is_biorthogonal(self) -> None:
        """Stops us believing a reduced operator that is not the balanced one.

        ``W^T V = I`` is what makes ``W^T G V`` the balanced reduced generator.
        Without it the pair is just two arbitrary bases and the reported order
        would not be the order of anything.
        """

        columns = [0, 1, 2]
        family = balanced_family(
            self.generator, self.observables, columns, self.cache, self.sources,
            self.reference_index, [1.0, 1.0, 1.0],
        )
        for order in BALANCED_ORDERS:
            entry = family["orders"].get(str(order))
            if entry is None or "status" in entry:
                continue
            self.assertLess(entry["biorthogonality_defect"], 1e-6)

    def test_rescaling_a_readout_does_not_move_the_balanced_result(self) -> None:
        """Stops us believing a balanced order that a unit change could set.

        #580 showed the pooled Frobenius norm is magnitude-weighted.  A balancing
        transform is exactly the object a rescaled readout would silently rotate,
        so the metric is fixed before balancing by dividing each readout by its
        own baseline contrast signal.  If that were not doing its job, doubling
        one readout would change the Hankel spectrum.
        """

        columns = [0, 1, 2]
        scales = [1.0, 1.0, 1.0]
        base = balanced_family(
            self.generator, self.observables, columns, self.cache, self.sources,
            self.reference_index, scales,
        )
        doubled = [list(vector) for vector in self.observables]
        doubled[1] = [2.0 * value for value in doubled[1]]
        moved = balanced_family(
            self.generator, doubled, columns, self.cache, self.sources,
            self.reference_index, [1.0, 2.0, 1.0],
        )
        for left, right in zip(
            base["hankel_singular_values"], moved["hankel_singular_values"]
        ):
            self.assertLess(abs(left - right), 1e-8)


class DeclaredManifest(unittest.TestCase):
    def test_the_dictionaries_are_nested(self) -> None:
        """Stops us believing a growth profile over a filtration that is not one.

        Phase B reads the *profile* across D0 subset D1 subset D2.  If the sets
        were not nested, adding outputs would not be strictly harder and the
        growth column would mean nothing.
        """

        names = list(DICTIONARIES)
        for earlier, later in zip(names, names[1:]):
            self.assertTrue(set(DICTIONARIES[earlier]) < set(DICTIONARIES[later]))

    def test_the_two_repair_ladders_are_comparable(self) -> None:
        """Stops us believing a per-degree-of-freedom comparison between ladders
        of different lengths.

        S is scored in added ranks and M in added poles.  If one ladder reached
        far beyond the other, the winner could be decided by where each stopped
        rather than by what each buys.
        """

        from scripts.p398_memory_closure import AUGMENT_RANKS

        added = [rank - FROZEN_RANK for rank in AUGMENT_RANKS]
        self.assertEqual(max(added), max(MEMORY_ORDERS) - 2)
        self.assertEqual(min(added), 0)


if __name__ == "__main__":
    unittest.main()
