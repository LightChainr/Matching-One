#!/usr/bin/env python3
"""Exact intervention transport on the P398 noncrossing planar state space.

Issue #580 asks a question that no same-environment low-rank fit can answer:
can the state functions be frozen **once**, in the baseline environment, and
then predict a new interventional dynamics by changing only the generator?

P398 gives an unusually clean place to ask it because the intervention is
known exactly.  On the canonical noncrossing connectivity states of width
``w`` (issue #11's frozen state space) every boundary point carries two
competing moves, ``join_cyclic_adjacent`` and ``detach``.  Tilting their
rates by one scalar gives the declared generator family

    G_eta = (1 + eta) * sum(join moves)
          + (1 - eta) * sum(detach moves)
          - exit-rate diagonal
          = G_0 + eta * H,

so the intervention direction ``H`` is not estimated, it is written down.

That exactness cuts both ways, and the run says so out loud.  The baseline
generator is ``J + D`` and the declared intervention is ``J - D``, so the
intervened generator never leaves the two-dimensional operator pencil
``span{J, D}``.  Two further directions are therefore scored beside it: a
second uniform direction inside the same pencil, and a tilt of one boundary
point only, which is outside it.  Comparing the three is what separates
transport that reflects a state from transport that reflects the algebra the
state was built from.

Everything below is deterministic.  No Monte Carlo, no width increase and no
mark search is used, as #580 requires.

What is frozen and what is allowed to move
------------------------------------------

Frozen at ``eta = 0``, never relearned from an ``eta != 0`` outcome:

* the state functions ``Phi`` (an orthonormal basis of a declared subspace of
  functions on the microscopic state space);
* the readout functions and their coordinates ``Phi^T f``;
* the source distributions and their coordinates ``Phi^T mu``;
* the inner product used by the projector (plain counting measure);
* the reduced tangent ``B = Phi^T H Phi``.

Allowed to move with ``eta``, because the intervention physically moves them:

* the microscopic generator ``G_eta``;
* the stationary law ``pi_eta``, which the intervention really does move.  It
  is computed and cross-checked against an exact rational solve, but it is
  never used to choose a state function, so no result here depends on it.

The three models
----------------

``M_refit``       the same span construction redone from scratch at each eta,
                  then scored at its own best-possible in-span prediction.
``M_fixed_span``  ``Phi`` frozen at eta = 0, scored at *its* best-possible
                  in-span prediction at eta.
``M_transport``   ``Phi``, ``A_0`` and ``B`` all frozen at eta = 0, with
                  ``A_eta = A_0 + eta B`` and no target refit at all.

"Best-possible in-span prediction" is not a fit.  Every model in this family
reconstructs the evolved observable as ``Phi v(t)``, so replacing ``v`` by the
exact projected trajectory ``Phi^T exp(t G_eta) f`` is an infimum over *all*
reduced generators in the span -- time-dependent ones included.  Using it
removes every estimator choice from the two upper rungs of the ladder, which
matters because a badly conditioned snapshot fit can make a good span look bad.
One ordinary least-squares generator fit is still reported, as
``fixed_span_generator_fit``, to show what an estimator actually recovers.

The infimum is exact in the state-space metric, reported as
``state_space_primary``.  On the source-sampled response matrix it is not a
bound: four sources cannot see every direction the projection optimizes.  Both
channels are carried for every model so that neither claim has to be softened.

One exact fact worth stating before any number is read: for a generator that
is affine in eta and a frozen span, the Galerkin compression is *automatically*
affine, ``Phi^T G_eta Phi = A_0 + eta B``.  So Contract II adds no assumption
on top of Contract I as long as ``A_eta`` is obtained by projection; the
contracts only separate once ``A_eta`` is allowed to see the target.  The gap
between ``M_fixed_span`` and ``M_transport`` is therefore the gap between "the
frozen functions can still host *some* reduced dynamics" and "the frozen
functions host the reduced dynamics the intervention law predicts".

A second fact the run has to check rather than assume: a refit that returns
the same span is not a refit.  ``refit_span_drift_from_frozen`` is reported for
every rank, because a zero ``M_refit``/``M_fixed_span`` gap means one of two
opposite things -- a stable state, or a span that structurally cannot move.

Scientific boundary
-------------------

This is a falsification gate on finite low-rank state claims.  Surviving it is
not latent-state identification (Bing et al. arXiv:2312.03580), and it says
nothing about Virasoro/LCFT modules.  It is also not a statement about square
site percolation: the state space here is P398's, not Matching One's.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

try:  # pragma: no cover - import shim, exercised both ways in practice
    from scripts.noncrossing_connectivity_codec import noncrossing_states
    from scripts.planar_state_operations import (
        detach_rgs,
        join_cyclic_adjacent_rgs,
        rgs_to_blocks,
    )
except ModuleNotFoundError:  # Direct ``python scripts/...`` execution.
    from noncrossing_connectivity_codec import noncrossing_states
    from planar_state_operations import (
        detach_rgs,
        join_cyclic_adjacent_rgs,
        rgs_to_blocks,
    )


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "p398-intervention-transport" / "latest.json"
SCHEMA = "matching-one.p398-intervention-transport.v1"
ISSUE = 580

DEFAULT_WIDTHS: Tuple[int, ...] = (4, 5, 6, 7)
#: The declared intervention ladder.  +-1/4 is the pair #580 names; the two
#: smaller and two larger values exist only to separate "the frozen span is
#: wrong" from "finite eta is nonlinear", and are declared here before any
#: number is computed.
#: ``0.0`` is in the ladder on purpose: at eta = 0 the transported model *is*
#: the baseline Galerkin model, so its error there is the rank-truncation error
#: of the frozen span with no intervention at all.  Every finite-eta number has
#: to be read as an excess over it, or the truncation gets reported as a
#: transport failure.
ETAS: Tuple[float, ...] = (-0.5, -0.25, -0.125, 0.0, 0.125, 0.25, 0.5)
PRIMARY_ETAS: Tuple[float, ...] = (-0.25, 0.25)
#: Declared lags.  The full matrix over all of them is scored; no lag is picked
#: after the fact.
LAGS: Tuple[float, ...] = (0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0)
RANKS: Tuple[int, ...] = (3, 4, 6, 8, 12)
SPECTRAL_ITERATIONS = 400
#: Orthogonal iteration only needs re-orthogonalization often enough that the
#: columns do not collapse onto the dominant direction between sweeps; doing it
#: every step made the span construction the most expensive part of the run.
REORTHOGONALIZE_EVERY = 5
STATIONARY_ITERATIONS = 20000
STATIONARY_TOLERANCE = 1e-14
POISSON_TAIL = 1e-16


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


# --------------------------------------------------------------------------
# Microscopic observables on noncrossing connectivity states
# --------------------------------------------------------------------------


def observable_blocks(state: Sequence[int]) -> float:
    return float(max(state) + 1)


def observable_singletons(state: Sequence[int]) -> float:
    return float(sum(1 for block in rgs_to_blocks(state) if len(block) == 1))


def observable_wrap(state: Sequence[int]) -> float:
    return 1.0 if state[0] == state[-1] else 0.0


def observable_max_block(state: Sequence[int]) -> float:
    return float(max(len(block) for block in rgs_to_blocks(state)))


def observable_linked_pairs(state: Sequence[int]) -> float:
    return float(sum(len(block) * (len(block) - 1) // 2 for block in rgs_to_blocks(state)))


def observable_halves_linked(state: Sequence[int]) -> float:
    width = len(state)
    half = width // 2
    for block in rgs_to_blocks(state):
        if any(point < half for point in block) and any(point >= half for point in block):
            return 1.0
    return 0.0


def observable_covering_depth(state: Sequence[int]) -> float:
    blocks = [sorted(block) for block in rgs_to_blocks(state)]
    best = 0
    for point in range(len(state)):
        depth = sum(
            1
            for block in blocks
            if block[0] < point < block[-1] and point not in block
        )
        best = max(best, depth)
    return float(best)


def observable_boundary_span(state: Sequence[int]) -> float:
    return float(sum(max(block) - min(block) for block in rgs_to_blocks(state)))


#: The block used to seed the frozen span.  ``one`` is the trivial invariant
#: mode and is always physical; the other three are the coarsest exact
#: connectivity counts on the boundary.
PRIMARY_READOUTS: Tuple[Tuple[str, Callable[[Sequence[int]], float]], ...] = (
    ("blocks", observable_blocks),
    ("singletons", observable_singletons),
    ("wrap", observable_wrap),
)
#: Never used to build ``Phi``.  These exist so that every reported transport
#: number has a held-out counterpart.
HELD_OUT_READOUTS: Tuple[Tuple[str, Callable[[Sequence[int]], float]], ...] = (
    ("max_block", observable_max_block),
    ("linked_pairs", observable_linked_pairs),
    ("halves_linked", observable_halves_linked),
    ("covering_depth", observable_covering_depth),
    ("boundary_span", observable_boundary_span),
)


# --------------------------------------------------------------------------
# The exact generator family
# --------------------------------------------------------------------------


#: The intervention directions.  Every rate multiplier is ``1 + eta * c``, so
#: all three stay non-negative for ``|eta| <= 1`` and all three are exact.
#:
#: ``uniform_join_minus_detach`` is the family #580 declares.  The other two
#: exist to answer a question the declared family cannot: the baseline
#: generator is ``J + D`` and that intervention is ``J - D``, so the intervened
#: generator never leaves the two-dimensional operator pencil ``span{J, D}``
#: that the state space was built from.  ``uniform_detach_only`` is a second,
#: independent direction *inside* the same pencil; ``single_point_join`` tilts
#: one boundary point only and is therefore *outside* it.  If a frozen span
#: transports the first two and fails the third, transport success is a
#: statement about where the intervention sits relative to the algebra that
#: generated the state -- not about the state being physical.
INTERVENTIONS: Dict[str, Callable[[str, int, int], float]] = {
    "uniform_join_minus_detach": lambda operation, point, width: (
        1.0 if operation == "join" else -1.0
    ),
    "uniform_detach_only": lambda operation, point, width: (
        1.0 if operation == "detach" else 0.0
    ),
    "single_point_join": lambda operation, point, width: (
        1.0 if operation == "join" and point == 0 else 0.0
    ),
}
IN_PENCIL = {
    "uniform_join_minus_detach": True,
    "uniform_detach_only": True,
    "single_point_join": False,
}
OPERATIONS: Tuple[str, ...] = ("join", "detach")


class Generator:
    """The exact move targets behind every rate-tilted generator on this space.

    ``target[(operation, point)][source]`` is the state reached, or ``None``
    when the move does not change the state.  A self transition is a no-op for
    a continuous-time chain -- keeping it would add the same rate to an
    off-diagonal entry and to the exit rate, cancelling exactly -- so it is
    dropped rather than recorded.
    """

    def __init__(self, width: int) -> None:
        _require(type(width) is int and 1 <= width <= 8, "width must be in [1,8]")
        states = noncrossing_states(width)
        rank = {state: index for index, state in enumerate(states)}
        target: Dict[Tuple[str, int], List[Optional[int]]] = {}
        for point in range(width):
            for operation, implementation in (
                ("join", join_cyclic_adjacent_rgs),
                ("detach", detach_rgs),
            ):
                column: List[Optional[int]] = []
                for source, state in enumerate(states):
                    image = rank[implementation(state, point)]
                    column.append(None if image == source else image)
                target[(operation, point)] = column
        self.width = width
        self.states = states
        self.size = len(states)
        self.target = target
        self.moves: Tuple[Tuple[str, int], ...] = tuple(
            (operation, point) for point in range(width) for operation in OPERATIONS
        )

    # -- rate vectors ------------------------------------------------------

    def baseline_rates(self) -> Dict[Tuple[str, int], float]:
        return {move: 1.0 for move in self.moves}

    def rates(self, intervention: str, eta: float) -> Dict[Tuple[str, int], float]:
        coefficient = INTERVENTIONS[intervention]
        return {
            (operation, point): 1.0 + eta * coefficient(operation, point, self.width)
            for operation, point in self.moves
        }

    def coefficients(self, intervention: str) -> Dict[Tuple[str, int], float]:
        coefficient = INTERVENTIONS[intervention]
        return {
            (operation, point): coefficient(operation, point, self.width)
            for operation, point in self.moves
        }

    # -- exact structure ---------------------------------------------------

    def exact_entry(
        self,
        rates: Mapping[Tuple[str, int], Fraction],
        source: int,
        destination: int,
    ) -> Fraction:
        value = Fraction(0)
        for move, rate in rates.items():
            image = self.target[move][source]
            if image is None:
                continue
            if image == destination:
                value += rate
            if source == destination:
                value -= rate
        return value

    def exact_rates(self, intervention: str, eta: Fraction) -> Dict[Tuple[str, int], Fraction]:
        coefficient = INTERVENTIONS[intervention]
        return {
            (operation, point): Fraction(1) + eta * Fraction(
                coefficient(operation, point, self.width)
            ).limit_denominator(10**6)
            for operation, point in self.moves
        }

    def exit_rate(self, rates: Mapping[Tuple[str, int], float]) -> float:
        best = 0.0
        for source in range(self.size):
            total = 0.0
            for move, rate in rates.items():
                if self.target[move][source] is not None:
                    total += rate
            best = max(best, total)
        return best

    # -- float sparse forms ------------------------------------------------

    def rows(
        self, rates: Mapping[Tuple[str, int], float]
    ) -> List[List[Tuple[int, float]]]:
        """Row-oriented operator; ``matvec(rows, f)`` acts on functions."""

        out: List[List[Tuple[int, float]]] = []
        for source in range(self.size):
            entries: Dict[int, float] = {}
            diagonal = 0.0
            for move, rate in rates.items():
                image = self.target[move][source]
                if image is None or rate == 0.0:
                    continue
                entries[image] = entries.get(image, 0.0) + rate
                diagonal -= rate
            entries[source] = entries.get(source, 0.0) + diagonal
            out.append(sorted(entries.items()))
        return out

    def columns(
        self, rates: Mapping[Tuple[str, int], float]
    ) -> List[List[Tuple[int, float]]]:
        """Column-oriented operator; used only for the stationary law."""

        gathered: List[Dict[int, float]] = [dict() for _ in range(self.size)]
        for source, row in enumerate(self.rows(rates)):
            for destination, value in row:
                gathered[destination][source] = (
                    gathered[destination].get(source, 0.0) + value
                )
        return [sorted(column.items()) for column in gathered]

    def tangent_rows(self, intervention: str) -> List[List[Tuple[int, float]]]:
        """``H = dG/deta`` for the named intervention, known exactly."""

        return self.rows(self.coefficients(intervention))


def matvec(rows: Sequence[Sequence[Tuple[int, float]]], vector: Sequence[float]) -> List[float]:
    return [sum(value * vector[column] for column, value in row) for row in rows]


# --------------------------------------------------------------------------
# Exact gates on the generator family
# --------------------------------------------------------------------------


def generator_gate(generator: Generator, intervention: str) -> Dict[str, Any]:
    """Everything about the generator family that must hold before any fit is read.

    Row sums, positivity of the off-diagonal rates on the whole declared eta
    range, exact affineness in eta, and irreducibility.  If any of these fails
    the object being compressed is not a Markov generator and every error
    reported downstream is an error about nothing.
    """

    etas = (Fraction(-1), Fraction(-1, 4), Fraction(0), Fraction(1, 4), Fraction(1))
    row_sum_failures = 0
    negative_offdiagonal = 0
    affine_failures = 0
    baseline_rates = {move: Fraction(1) for move in generator.moves}
    tangent_rates = generator.exact_rates(intervention, Fraction(1))
    for eta in etas:
        rates = generator.exact_rates(intervention, eta)
        for source in range(generator.size):
            total = Fraction(0)
            for destination in range(generator.size):
                entry = generator.exact_entry(rates, source, destination)
                total += entry
                if destination != source and entry < 0:
                    negative_offdiagonal += 1
                baseline = generator.exact_entry(baseline_rates, source, destination)
                tangent = (
                    generator.exact_entry(tangent_rates, source, destination) - baseline
                )
                if entry != baseline + eta * tangent:
                    affine_failures += 1
            if total != 0:
                row_sum_failures += 1
    reachable = _strongly_connected(generator)
    return {
        "width": generator.width,
        "intervention": intervention,
        "states": generator.size,
        "checked_etas": [str(eta) for eta in etas],
        "row_sum_failures": row_sum_failures,
        "negative_offdiagonal_entries": negative_offdiagonal,
        "affine_in_eta_failures": affine_failures,
        "strongly_connected": reachable,
        "maximum_exit_rate_at_eta_0": generator.exit_rate(generator.baseline_rates()),
        "passed": (
            row_sum_failures == 0
            and negative_offdiagonal == 0
            and affine_failures == 0
            and reachable
        ),
    }


def _strongly_connected(generator: Generator) -> bool:
    forward: List[List[int]] = []
    for source in range(generator.size):
        neighbours = set()
        for move in generator.moves:
            image = generator.target[move][source]
            if image is not None:
                neighbours.add(image)
        forward.append(sorted(neighbours))
    backward: List[List[int]] = [[] for _ in range(generator.size)]
    for source, targets in enumerate(forward):
        for destination in targets:
            backward[destination].append(source)

    def reach(adjacency: Sequence[Sequence[int]]) -> int:
        seen = {0}
        stack = [0]
        while stack:
            node = stack.pop()
            for neighbour in adjacency[node]:
                if neighbour not in seen:
                    seen.add(neighbour)
                    stack.append(neighbour)
        return len(seen)

    return reach(forward) == generator.size and reach(backward) == generator.size


# --------------------------------------------------------------------------
# Dense linear algebra on the reduced (r x r) blocks
# --------------------------------------------------------------------------


def _dot(left: Sequence[float], right: Sequence[float]) -> float:
    return sum(a * b for a, b in zip(left, right))


def _norm(vector: Sequence[float]) -> float:
    return math.sqrt(_dot(vector, vector))


def orthonormalize(vectors: Sequence[Sequence[float]], tolerance: float = 1e-10) -> List[List[float]]:
    """Modified Gram-Schmidt with a rank-revealing drop."""

    basis: List[List[float]] = []
    for candidate in vectors:
        working = list(candidate)
        original = _norm(working)
        if original == 0.0:
            continue
        for _ in range(2):  # one reorthogonalization pass
            for column in basis:
                overlap = _dot(column, working)
                for index in range(len(working)):
                    working[index] -= overlap * column[index]
        length = _norm(working)
        if length <= tolerance * original:
            continue
        basis.append([value / length for value in working])
    return basis


def solve_dense(matrix: Sequence[Sequence[float]], rhs: Sequence[Sequence[float]]) -> List[List[float]]:
    """Gaussian elimination with partial pivoting, ``matrix @ X = rhs``."""

    size = len(matrix)
    augmented = [list(matrix[index]) + list(rhs[index]) for index in range(size)]
    width = len(augmented[0])
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(augmented[row][column]))
        _require(abs(augmented[pivot][column]) > 1e-13, "reduced least-squares system is singular")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        for index in range(column, width):
            augmented[column][index] /= scale
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor == 0.0:
                continue
            for index in range(column, width):
                augmented[row][index] -= factor * augmented[column][index]
    return [row[size:] for row in augmented]


def dense_expm(matrix: Sequence[Sequence[float]], time: float) -> List[List[float]]:
    """Scaling-and-squaring Taylor exponential for the small reduced blocks."""

    size = len(matrix)
    scale = max(
        (abs(value) for row in matrix for value in row),
        default=0.0,
    ) * abs(time)
    squarings = 0
    while scale > 0.5:
        scale /= 2.0
        squarings += 1
    factor = time / (2.0**squarings)
    scaled = [[value * factor for value in row] for row in matrix]
    result = [[1.0 if i == j else 0.0 for j in range(size)] for i in range(size)]
    term = [row[:] for row in result]
    for order in range(1, 19):
        term = [
            [sum(term[i][k] * scaled[k][j] for k in range(size)) / order for j in range(size)]
            for i in range(size)
        ]
        result = [[result[i][j] + term[i][j] for j in range(size)] for i in range(size)]
    for _ in range(squarings):
        result = [
            [sum(result[i][k] * result[k][j] for k in range(size)) for j in range(size)]
            for i in range(size)
        ]
    return result


# --------------------------------------------------------------------------
# Frozen spans
# --------------------------------------------------------------------------


def krylov_span(
    rows: Sequence[Sequence[Tuple[int, float]]],
    seeds: Sequence[Sequence[float]],
    rank: int,
) -> List[List[float]]:
    """Block Krylov space of the declared readouts under the declared generator.

    This is deliberately the ordinary Hankel/low-rank state: it is exactly the
    span a same-environment rank-``r`` fit of these observables would recover.

    Two details are load-bearing.  Each new frontier vector is normalized
    before the next application of ``G``, because the raw powers grow like
    ``||G||^k`` and reach NaN well before the ranks this run needs.  And the
    frontier keeps the declared seed order, so that a rank cutting through the
    middle of a Krylov level is still a *declared* span rather than an artifact
    of how the previous level happened to be orthogonalized.
    """

    basis = orthonormalize(list(seeds))[:rank]
    frontier = [list(seed) for seed in seeds]
    for _ in range(len(seeds) * rank + 1):
        if len(basis) >= rank:
            break
        stepped: List[List[float]] = []
        for vector in frontier:
            image = matvec(rows, vector)
            length = _norm(image)
            if length == 0.0:
                continue
            stepped.append([value / length for value in image])
        if not stepped:
            break
        grown = orthonormalize(list(basis) + stepped)[:rank]
        if len(grown) == len(basis):
            break
        basis = grown
        frontier = stepped
    return basis[:rank]


def spectral_span(
    rows: Sequence[Sequence[Tuple[int, float]]],
    size: int,
    rank: int,
    lazy_rate: float,
    iterations: int = SPECTRAL_ITERATIONS,
) -> List[List[float]]:
    """Dominant invariant subspace of the lazy uniformized chain.

    ``P = I + G/(2*Lambda)`` has spectrum ``1 + lambda/(2*Lambda)`` inside the
    unit disc with the slow modes on the outside, so subspace iteration on
    ``P`` converges to the slow invariant subspace of ``G`` -- the span that
    closes *exactly* at baseline, which is what makes any failure at
    ``eta != 0`` attributable to the intervention alone.
    """

    vectors = [
        [1.0 if index == column else 0.0 for index in range(size)] for column in range(rank)
    ]
    vectors[0] = [1.0] * size  # the trivial invariant mode, exactly
    basis = orthonormalize(vectors)
    for step in range(iterations):
        moved = []
        for vector in basis:
            image = matvec(rows, vector)
            stepped = [vector[i] + image[i] / lazy_rate for i in range(size)]
            length = _norm(stepped)
            moved.append(
                [value / length for value in stepped] if length > 0 else stepped
            )
        if step % REORTHOGONALIZE_EVERY == REORTHOGONALIZE_EVERY - 1:
            moved = orthonormalize(moved)
        basis = moved
        if len(basis) < rank:
            break
    basis = orthonormalize(basis)
    return basis[:rank]


def deterministic_span(size: int, rank: int, seed: int = 20260906) -> List[List[float]]:
    """A control span from a fixed linear congruential stream.

    It carries no information about the dynamics, so any transport quality it
    reaches is the quality that costs nothing.
    """

    state = seed
    vectors: List[List[float]] = []
    for _ in range(rank):
        vector = []
        for _ in range(size):
            state = (1103515245 * state + 12345) % (1 << 31)
            vector.append(state / float(1 << 30) - 1.0)
        vectors.append(vector)
    return orthonormalize(vectors)[:rank]


def project_coordinates(basis: Sequence[Sequence[float]], vector: Sequence[float]) -> List[float]:
    return [_dot(column, vector) for column in basis]


def reduced_operator(
    basis: Sequence[Sequence[float]],
    rows: Sequence[Sequence[Tuple[int, float]]],
) -> List[List[float]]:
    """``Phi^T G Phi`` with ``Phi`` orthonormal in the counting inner product."""

    images = [matvec(rows, column) for column in basis]
    return [
        [_dot(basis[i], images[j]) for j in range(len(basis))] for i in range(len(basis))
    ]


def closure_residual(
    basis: Sequence[Sequence[float]],
    rows: Sequence[Sequence[Tuple[int, float]]],
) -> Tuple[float, float, List[List[float]]]:
    """``(I - Phi Phi^T) G Phi`` -- its Frobenius norm, relative norm, columns."""

    residual_columns: List[List[float]] = []
    total = 0.0
    reference = 0.0
    for column in basis:
        image = matvec(rows, column)
        reference += _dot(image, image)
        residual = list(image)
        for other in basis:
            overlap = _dot(other, residual)
            for index in range(len(residual)):
                residual[index] -= overlap * other[index]
        total += _dot(residual, residual)
        residual_columns.append(residual)
    absolute = math.sqrt(total)
    relative = absolute / math.sqrt(reference) if reference > 0 else 0.0
    return absolute, relative, residual_columns


def leading_direction(columns: Sequence[Sequence[float]], iterations: int = 500) -> Tuple[List[float], float]:
    """Top left singular vector of the residual block, by power iteration on ``R R^T``."""

    if not columns:
        return [], 0.0
    size = len(columns[0])
    vector = [0.0] * size
    for column in columns:
        for index in range(size):
            vector[index] += column[index]
    length = _norm(vector)
    if length == 0.0:
        return [0.0] * size, 0.0
    vector = [value / length for value in vector]
    value = 0.0
    for _ in range(iterations):
        coefficients = [_dot(column, vector) for column in columns]
        image = [0.0] * size
        for coefficient, column in zip(coefficients, columns):
            if coefficient == 0.0:
                continue
            for index in range(size):
                image[index] += coefficient * column[index]
        length = _norm(image)
        if length == 0.0:
            return vector, 0.0
        vector = [entry / length for entry in image]
        value = math.sqrt(length)
    return vector, value


# --------------------------------------------------------------------------
# Exact microscopic response
# --------------------------------------------------------------------------


def poisson_horizon(rate_time: float) -> int:
    horizon = int(rate_time + 12.0 * math.sqrt(rate_time + 1.0) + 40.0)
    return max(horizon, 40)


def evolve_observables(
    rows: Sequence[Sequence[Tuple[int, float]]],
    size: int,
    uniform_rate: float,
    observables: Sequence[Sequence[float]],
    lags: Sequence[float],
) -> List[List[List[float]]]:
    """``exp(t G) f`` for every declared readout and lag, by uniformization.

    Uniformization is used rather than a Taylor series in ``G`` because every
    term is a Poisson probability times a non-negative power of a stochastic
    matrix: there is no cancellation, so the answer stays accurate at the large
    ``t`` the slow modes need.  Returned as ``[lag][observable][state]``.
    """

    horizon = poisson_horizon(uniform_rate * max(lags))
    means = [uniform_rate * lag for lag in lags]
    weights = [math.exp(-mean) for mean in means]
    accumulated = [
        [[0.0] * size for _ in observables] for _ in lags
    ]
    current = [list(observable) for observable in observables]
    for order in range(horizon + 1):
        for lag_index in range(len(lags)):
            weight = weights[lag_index]
            if weight <= POISSON_TAIL:
                continue
            block = accumulated[lag_index]
            for observable_index, vector in enumerate(current):
                row = block[observable_index]
                for position in range(size):
                    row[position] += weight * vector[position]
        for lag_index, mean in enumerate(means):
            weights[lag_index] *= mean / (order + 1)
        if order == horizon:
            break
        stepped = []
        for vector in current:
            image = matvec(rows, vector)
            stepped.append([vector[i] + image[i] / uniform_rate for i in range(size)])
        current = stepped
    return accumulated


def stationary_law(
    columns: Sequence[Sequence[Tuple[int, float]]],
    size: int,
    uniform_rate: float,
    iterations: int = STATIONARY_ITERATIONS,
) -> Tuple[List[float], float]:
    """Left null vector of ``G``, by power iteration on the lazy chain."""

    vector = [1.0 / size] * size
    for _ in range(iterations):
        image = [sum(value * vector[source] for source, value in column) for column in columns]
        stepped = [vector[i] + image[i] / uniform_rate for i in range(size)]
        total = sum(stepped)
        stepped = [value / total for value in stepped]
        shift = max(abs(stepped[i] - vector[i]) for i in range(size))
        vector = stepped
        if shift < STATIONARY_TOLERANCE:
            break
    residual = [sum(value * vector[source] for source, value in column) for column in columns]
    return vector, max(abs(value) for value in residual)


def exact_stationary(
    generator: Generator, rates: Mapping[Tuple[str, int], Fraction]
) -> List[Fraction]:
    """Exact rational stationary law, affordable only at the smallest widths."""

    size = generator.size
    matrix = [
        [generator.exact_entry(rates, source, target) for source in range(size)]
        for target in range(size)
    ]
    matrix[size - 1] = [Fraction(1)] * size
    rhs = [Fraction(0)] * size
    rhs[size - 1] = Fraction(1)
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if matrix[row][column] != 0), None
        )
        _require(pivot is not None, "exact stationary system is singular")
        matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
        rhs[column], rhs[pivot] = rhs[pivot], rhs[column]
        scale = matrix[column][column]
        matrix[column] = [value / scale for value in matrix[column]]
        rhs[column] /= scale
        for row in range(size):
            if row == column or matrix[row][column] == 0:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                matrix[row][index] - factor * matrix[column][index] for index in range(size)
            ]
            rhs[row] -= factor * rhs[column]
    return rhs


# --------------------------------------------------------------------------
# Frozen sources and the scored response tensor
# --------------------------------------------------------------------------


def source_distributions(generator: Generator) -> List[Tuple[str, List[float]]]:
    """The frozen initial laws.  ``uniform`` is the declared contrast reference."""

    size = generator.size
    rank = {state: index for index, state in enumerate(generator.states)}
    width = generator.width
    discrete = tuple(range(width))
    single = tuple(0 for _ in range(width))
    wrapped = tuple([0] + list(range(1, width - 1)) + [0]) if width >= 3 else single
    named: List[Tuple[str, List[float]]] = []
    for name, state in (
        ("delta_all_singletons", discrete),
        ("delta_single_block", single),
        ("delta_wrapped_pair", wrapped),
    ):
        vector = [0.0] * size
        vector[rank[state]] = 1.0
        named.append((name, vector))
    named.append(("uniform", [1.0 / size] * size))
    return named


CONTRAST_REFERENCE = "uniform"


def response_tensor(
    sources: Sequence[Tuple[str, Sequence[float]]],
    evolved: Sequence[Sequence[Sequence[float]]],
) -> List[List[List[float]]]:
    """``C(t)[source][readout] = <mu, exp(t G) f>`` over the declared lags."""

    return [
        [[_dot(source, vector) for vector in block] for _, source in sources]
        for block in evolved
    ]


def contrast_tensor(
    tensor: Sequence[Sequence[Sequence[float]]],
    reference_index: int,
) -> List[List[List[float]]]:
    """Each source response minus the declared reference source's response.

    The raw response carries a constant part -- the observable's mean and the
    common ``t -> infinity`` limit -- that every model with the constant
    function in its span reproduces for free.  Scoring that part would make
    each of the three models look equally good and would hide the entire
    question.  The contrast is what the dynamics actually has to get right.
    """

    return [
        [
            [row[readout] - block[reference_index][readout] for readout in range(len(row))]
            for index, row in enumerate(block)
            if index != reference_index
        ]
        for block in tensor
    ]


def relative_error(
    truth: Sequence[Sequence[Sequence[float]]],
    prediction: Sequence[Sequence[Sequence[float]]],
    columns: Optional[Sequence[int]] = None,
) -> float:
    numerator = 0.0
    denominator = 0.0
    for lag_index, block in enumerate(truth):
        for row_index, row in enumerate(block):
            indices = range(len(row)) if columns is None else columns
            for readout in indices:
                difference = row[readout] - prediction[lag_index][row_index][readout]
                numerator += difference * difference
                denominator += row[readout] * row[readout]
    if denominator == 0.0:
        return 0.0
    return math.sqrt(numerator / denominator)


def _multiply(
    left: Sequence[Sequence[float]], right: Sequence[Sequence[float]]
) -> List[List[float]]:
    size = len(left)
    return [
        [sum(left[i][k] * right[k][j] for k in range(size)) for j in range(size)]
        for i in range(size)
    ]


def propagator_ladder(
    reduced: Sequence[Sequence[float]], lags: Sequence[float]
) -> List[List[List[float]]]:
    """``exp(t A)`` for every declared lag, from one scaling-and-squaring.

    The declared lag grid is a set of integer multiples of its smallest entry,
    so every propagator is a power of the first one and the whole grid costs one
    matrix exponential plus a short product ladder.  Exponentiating each lag
    separately was the single largest cost in the run, and it recomputed the
    same squaring sequence seven times over.
    """

    step = min(lags)
    multiples = [lag / step for lag in lags]
    if step <= 0 or any(abs(value - round(value)) > 1e-9 for value in multiples):
        return [dense_expm(reduced, lag) for lag in lags]
    size = len(reduced)
    identity = [[1.0 if i == j else 0.0 for j in range(size)] for i in range(size)]
    powers: Dict[int, List[List[float]]] = {0: identity, 1: dense_expm(reduced, step)}

    def power(exponent: int) -> List[List[float]]:
        if exponent in powers:
            return powers[exponent]
        half = exponent // 2
        value = _multiply(power(half), power(exponent - half))
        powers[exponent] = value
        return value

    return [power(int(round(value))) for value in multiples]


def balanced_error(
    truth: Sequence[Sequence[Sequence[float]]],
    prediction: Sequence[Sequence[Sequence[float]]],
    columns: Sequence[int],
) -> float:
    """Root mean square of the *per readout* relative errors.

    Added after the first reading, at the owner's request on #580, and post-hoc
    by construction.  The declared pooled relative Frobenius norm is
    magnitude-weighted: at width 8 the declared block scores 0.068 while the
    declared ``wrap`` readout inside it is at 0.230, because ``blocks`` and
    ``singletons`` are ``O(w)`` and ``wrap`` is ``O(1)``.  Multiplying one
    observable by a constant must not decide whether a state transports.  Each
    readout carries its own denominator here, so it cannot.
    """

    values = [relative_error(truth, prediction, (column,)) for column in columns]
    if not values:
        return 0.0
    return math.sqrt(sum(value * value for value in values) / len(values))


def contrast_signal(
    truth: Sequence[Sequence[Sequence[float]]], column: int
) -> float:
    """Frobenius norm of one readout's contrast response -- its own denominator.

    A readout with almost no baseline contrast signal gets an unstable relative
    error, so it is classified as weakly identifiable rather than being given a
    near-zero denominator and a spectacular score in either direction.
    """

    total = 0.0
    for block in truth:
        for row in block:
            total += row[column] * row[column]
    return math.sqrt(total)


def reduced_trajectories(
    basis: Sequence[Sequence[float]],
    reduced: Sequence[Sequence[float]],
    observables: Sequence[Sequence[float]],
    lags: Sequence[float],
    propagators: Optional[Sequence[Sequence[Sequence[float]]]] = None,
) -> List[List[List[float]]]:
    """``exp(t A) Phi^T f`` -- one model's reduced trajectory per readout and lag."""

    coordinates = [project_coordinates(basis, vector) for vector in observables]
    rank = len(basis)
    if propagators is None:
        propagators = propagator_ladder(reduced, lags)
    out: List[List[List[float]]] = []
    for propagator in propagators:
        out.append(
            [
                [sum(propagator[i][j] * column[j] for j in range(rank)) for i in range(rank)]
                for column in coordinates
            ]
        )
    return out


def oracle_trajectories(
    basis: Sequence[Sequence[float]],
    evolved: Sequence[Sequence[Sequence[float]]],
) -> List[List[List[float]]]:
    """``Phi^T exp(t G_eta) f`` -- the best reduced trajectory the span admits.

    Every model in this family reconstructs the evolved observable as
    ``Phi v(t)``.  The orthogonal projection is the exact minimizer of
    ``||exp(t G_eta) f - Phi v(t)||`` over all ``v``, time-dependent ones
    included, so this is an infimum rather than a fit: it separates "the frozen
    functions are the wrong functions" from "the frozen functions are right but
    the transported generator law is wrong", which no single fitted ``A`` can.

    The infimum is in the state-space metric.  On the source-sampled response
    it is not a bound -- four sources cannot see every direction the projection
    optimizes -- which is why both metrics are reported.
    """

    return [[project_coordinates(basis, vector) for vector in block] for block in evolved]


def response_from_trajectories(
    source_coordinates: Sequence[Sequence[float]],
    trajectories: Sequence[Sequence[Sequence[float]]],
) -> List[List[List[float]]]:
    """``<mu, Phi v(t)> = (Phi^T mu) . v(t)`` over sources, readouts and lags."""

    return [
        [[_dot(source, vector) for vector in block] for source in source_coordinates]
        for block in trajectories
    ]


def reconstruction_relative(
    squared_norms: Sequence[Sequence[float]],
    projected: Sequence[Sequence[Sequence[float]]],
    trajectories: Sequence[Sequence[Sequence[float]]],
    columns: Sequence[int],
) -> float:
    """``||exp(t G) f - Phi v(t)|| / ||exp(t G) f||`` over the whole state space.

    Computed from ``||g||^2 - 2 v.c + |v|^2`` with ``c = Phi^T g``, which is
    exact for an orthonormal ``Phi`` and costs ``O(r)`` per readout and lag
    instead of ``O(n r)``: rebuilding ``Phi v`` for every model, rank, eta and
    width was the single largest cost in the run.

    This is the metric in which the in-span projection is provably minimal, so
    ``fixed_span <= transported`` is an inequality here rather than an empirical
    ordering.  On the source-sampled response it is not a bound.
    """

    numerator = 0.0
    denominator = 0.0
    for lag_index, block in enumerate(trajectories):
        for index in columns:
            coordinates = projected[lag_index][index]
            reduced = block[index]
            total = squared_norms[lag_index][index]
            for value, target in zip(reduced, coordinates):
                total += value * value - 2.0 * value * target
            numerator += max(total, 0.0)
            denominator += squared_norms[lag_index][index]
    if denominator == 0.0:
        return 0.0
    return math.sqrt(numerator / denominator)


def symmetric_eigen(matrix: Sequence[Sequence[float]], sweeps: int = 60) -> Tuple[List[float], List[List[float]]]:
    """Cyclic Jacobi eigendecomposition of a small symmetric matrix.

    Needed because the snapshot Gram matrix of decaying trajectories is badly
    conditioned -- every trajectory ends up near the constant function -- and a
    naive solve turns that into a reduced generator with a spurious growing
    mode.  Returns eigenvalues and the columns of the eigenvector matrix.
    """

    size = len(matrix)
    working = [list(row) for row in matrix]
    vectors = [[1.0 if i == j else 0.0 for j in range(size)] for i in range(size)]
    for _ in range(sweeps):
        off = math.sqrt(
            sum(working[i][j] ** 2 for i in range(size) for j in range(size) if i != j)
        )
        if off < 1e-15:
            break
        for p in range(size - 1):
            for q in range(p + 1, size):
                if abs(working[p][q]) < 1e-18:
                    continue
                theta = (working[q][q] - working[p][p]) / (2.0 * working[p][q])
                sign = 1.0 if theta >= 0 else -1.0
                t = sign / (abs(theta) + math.sqrt(theta * theta + 1.0))
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c
                for k in range(size):
                    akp = working[k][p]
                    akq = working[k][q]
                    working[k][p] = c * akp - s * akq
                    working[k][q] = s * akp + c * akq
                for k in range(size):
                    apk = working[p][k]
                    aqk = working[q][k]
                    working[p][k] = c * apk - s * aqk
                    working[q][k] = s * apk + c * aqk
                for k in range(size):
                    vkp = vectors[k][p]
                    vkq = vectors[k][q]
                    vectors[k][p] = c * vkp - s * vkq
                    vectors[k][q] = s * vkp + c * vkq
    values = [working[i][i] for i in range(size)]
    columns = [[vectors[row][column] for row in range(size)] for column in range(size)]
    return values, columns


def snapshot_moments(
    basis: Sequence[Sequence[float]],
    images: Sequence[Tuple[Sequence[float], Sequence[float]]],
) -> Tuple[List[List[float]], List[List[float]]]:
    """Gram and cross moments of the projected snapshot pairs.

    Split out of the fit so that the pseudo-inverse cut-off can be scanned
    without recomputing the microscopic matrix-vector products, which dominate
    the run at the larger widths.
    """

    rank = len(basis)
    gram = [[0.0] * rank for _ in range(rank)]
    cross = [[0.0] * rank for _ in range(rank)]
    for vector, image in images:
        coordinates = project_coordinates(basis, vector)
        moved = project_coordinates(basis, image)
        for i in range(rank):
            for j in range(rank):
                gram[i][j] += coordinates[i] * coordinates[j]
                cross[i][j] += moved[i] * coordinates[j]
    return gram, cross


def solve_reduced_fit(
    gram: Sequence[Sequence[float]],
    cross: Sequence[Sequence[float]],
    tolerance: float = 1e-8,
    eigen: Optional[Tuple[List[float], List[List[float]]]] = None,
) -> List[List[float]]:
    """One realizable reduced generator, least-squares fit to the eta target.

    ``M_fixed_span``'s projection says whether the frozen functions *could* host
    the intervened dynamics; this says whether an ordinary estimator actually
    recovers a generator from them.  The snapshots are the projected
    trajectories of the declared readouts, paired as
    ``(Phi^T g, Phi^T G_eta g)``, so an exactly invariant span returns the
    Galerkin operator.  The Gram inverse is truncated at ``tolerance`` times its
    top eigenvalue.
    """

    rank = len(gram)
    values, columns = symmetric_eigen(gram) if eigen is None else eigen
    top = max(values) if values else 0.0
    inverse = [[0.0] * rank for _ in range(rank)]
    for value, column in zip(values, columns):
        if value <= tolerance * top or value <= 0.0:
            continue
        for i in range(rank):
            for j in range(rank):
                inverse[i][j] += column[i] * column[j] / value
    return [
        [sum(cross[i][k] * inverse[k][j] for k in range(rank)) for j in range(rank)]
        for i in range(rank)
    ]


def span_drift(left: Sequence[Sequence[float]], right: Sequence[Sequence[float]]) -> float:
    """``||(I - P_left) Phi_right||_F / sqrt(rank)`` -- zero iff the spans agree.

    A refit that returns the same span is not a refit, and reading a
    ``M_refit``/``M_fixed_span`` gap that is identically zero as evidence for a
    stable state would be reading a tautology.
    """

    if not right:
        return 0.0
    total = 0.0
    for column in right:
        residual = list(column)
        for other in left:
            overlap = _dot(other, residual)
            for index in range(len(residual)):
                residual[index] -= overlap * other[index]
        total += _dot(residual, residual)
    return math.sqrt(total / len(right))


# --------------------------------------------------------------------------
# One (width, span family, rank) experiment
# --------------------------------------------------------------------------


def build_span(
    name: str,
    generator: Generator,
    rows: Sequence[Sequence[Tuple[int, float]]],
    seeds: Sequence[Sequence[float]],
    rank: int,
    exit_rate: float,
) -> List[List[float]]:
    if name == "krylov":
        return krylov_span(rows, seeds, rank)
    if name == "spectral":
        return spectral_span(rows, generator.size, rank, 2.0 * exit_rate)
    if name == "random_control":
        return deterministic_span(generator.size, rank)
    raise ValueError(f"unknown span family {name!r}")


def _best_generator_fit(
    basis: Sequence[Sequence[float]],
    images: Sequence[Tuple[Sequence[float], Sequence[float]]],
    observables: Sequence[Sequence[float]],
    source_coordinates: Sequence[Sequence[float]],
    lags: Sequence[float],
    truth_contrast: Sequence[Sequence[Sequence[float]]],
    score_columns: Sequence[int],
    reference_index: int,
) -> List[List[float]]:
    """The realizable estimator, given its best shot over the truncation level.

    The snapshot Gram matrix of decaying trajectories is badly conditioned, so
    the pseudo-inverse cut-off is not a neutral choice.  Scanning it and keeping
    the best is generous to the estimator on purpose: this number is a foil for
    the transported model, and a foil crippled by a tuning constant would prove
    nothing.
    """

    rank = len(basis)
    gram, cross = snapshot_moments(basis, images)
    eigen = symmetric_eigen(gram)
    best: Optional[List[List[float]]] = None
    best_error = float("inf")
    for tolerance in (1e-4, 1e-6, 1e-8, 1e-10, 1e-12):
        reduced = solve_reduced_fit(gram, cross, tolerance, eigen)
        trajectories = reduced_trajectories(basis, reduced, observables, lags)
        prediction = response_from_trajectories(source_coordinates, trajectories)
        error = relative_error(
            truth_contrast, contrast_tensor(prediction, reference_index), score_columns
        )
        if error == error and error < best_error:
            best_error = error
            best = reduced
    if best is None:  # pragma: no cover - every tolerance produced NaN
        return [[0.0] * rank for _ in range(rank)]
    return best


def flag_invariance(
    rows: Sequence[Sequence[Tuple[int, float]]],
    tangent: Sequence[Sequence[Tuple[int, float]]],
    seeds: Sequence[Sequence[float]],
    ranks: Sequence[int],
) -> Dict[str, Any]:
    """Does the intervention keep the observables inside the flag the baseline builds?

    ``K_r`` is the block Krylov flag of the declared readouts under ``G_0``.
    This measures ``||(I - P_{r+s}) H Phi_r|| / ||H Phi_r||`` with ``s`` the
    number of seeds -- that is, whether one step of the intervention lands
    inside the next Krylov level.

    It is the quantity that decides how a transport pass should be read.  If the
    intervention direction never leaves the flag the observables already
    generate, a frozen span will transport it for an algebraic reason, and
    calling that a physical state would be a category error.
    """

    seed_count = len(seeds)
    levels: Dict[str, float] = {}
    for rank in ranks:
        basis = krylov_span(rows, seeds, rank)
        if len(basis) < rank:
            continue
        outer = krylov_span(rows, seeds, rank + seed_count)
        numerator = 0.0
        denominator = 0.0
        for column in basis:
            image = matvec(tangent, column)
            denominator += _dot(image, image)
            residual = list(image)
            for other in outer:
                overlap = _dot(other, residual)
                for index in range(len(residual)):
                    residual[index] -= overlap * other[index]
            numerator += _dot(residual, residual)
        if denominator > 0:
            levels[str(rank)] = math.sqrt(numerator / denominator)
    declared = levels.get(str(seed_count))
    return {
        "measures": "||(I - P_{r+s}) H Phi_r|| / ||H Phi_r||, Phi_r the baseline Krylov flag",
        "by_rank": levels,
        "at_the_declared_readouts": declared,
        "declared_readouts_stay_inside_the_baseline_flag": (
            declared is not None and declared < 1e-10
        ),
    }


#: A readout whose baseline contrast signal is below this share of the largest
#: readout's has no stable relative denominator; it is classified as weakly
#: identifiable rather than scored.
WEAK_SIGNAL_FLOOR = 1e-3
#: Modular elimination costs ``r^2 n / 2``, so the whole rank is affordable up
#: to 429 states and only the largest width needs a floor instead of a value.
#: A cost control, not a modelling choice.
LINEAR_RANK_BUDGET_FULL_UP_TO = 429
LINEAR_RANK_BUDGET_ABOVE = 150


def linear_rank_budget(size: int) -> int:
    return size if size <= LINEAR_RANK_BUDGET_FULL_UP_TO else LINEAR_RANK_BUDGET_ABOVE
_LINEAR_RANK_CACHE: Dict[Tuple[int, str, int], Dict[str, Any]] = {}
#: The nested observable dictionaries #588 Phase B declares, built only from
#: readouts #580 had already declared.
DICTIONARIES: Dict[str, Tuple[str, ...]] = {
    "D0_additive_local_counts": ("blocks", "singletons", "wrap"),
    "D1_plus_size_and_extent": (
        "blocks",
        "singletons",
        "wrap",
        "max_block",
        "linked_pairs",
        "boundary_span",
    ),
    "D2_plus_nonlocal_topology": (
        "blocks",
        "singletons",
        "wrap",
        "max_block",
        "linked_pairs",
        "boundary_span",
        "halves_linked",
        "covering_depth",
    ),
}


def classify_readouts(
    readout_names: Sequence[str],
    declared: Sequence[str],
    truth_contrast_at_zero: Sequence[Sequence[Sequence[float]]],
    baseline_by_readout: Mapping[str, float],
    model_baseline_by_readout: Mapping[str, float],
    excess_by_readout: Mapping[str, float],
) -> Dict[str, Any]:
    """Separate baseline representability from intervention excess, per readout.

    A readout the frozen span cannot represent at ``eta = 0`` carries no
    transport verdict at all: its finite-eta error is dominated by a dictionary
    failure that was there before any intervention was applied.  Counting such a
    readout as a transport failure -- or, worse, as a success because its excess
    happens to be small -- is the mistake this classification exists to stop.

    Four bins, at the same 0.10 threshold the run already declared:
    represented and transports; represented but transport fails; unrepresented
    at baseline, so not identifiable; and weakly identifiable, where the
    readout's own baseline contrast signal is too small to be a denominator.
    """

    signals = {
        name: contrast_signal(truth_contrast_at_zero, index)
        for index, name in enumerate(readout_names)
    }
    strongest = max(signals.values()) if signals else 0.0
    rows: Dict[str, Any] = {}
    for name in readout_names:
        share = signals[name] / strongest if strongest > 0 else 0.0
        representability = baseline_by_readout.get(name)
        excess = excess_by_readout.get(name)
        if share < WEAK_SIGNAL_FLOOR:
            bin_name = "weakly_identifiable_baseline_signal"
        elif representability is None or representability > PASS_THRESHOLD:
            bin_name = "unrepresented_at_baseline_transport_not_identifiable"
        elif excess is not None and excess <= PASS_THRESHOLD:
            bin_name = "represented_and_transports"
        else:
            bin_name = "represented_but_transport_fails"
        rows[name] = {
            "in_the_declared_dictionary": name in declared,
            "R_baseline_representability": representability,
            "R_baseline_of_the_transported_model": model_baseline_by_readout.get(name),
            "T_intervention_excess": excess,
            "baseline_contrast_signal_share": share,
            "bin": bin_name,
        }
    tally: Dict[str, int] = {}
    for row in rows.values():
        tally[row["bin"]] = tally.get(row["bin"], 0) + 1
    return {
        "threshold": PASS_THRESHOLD,
        "weak_signal_floor": WEAK_SIGNAL_FLOOR,
        "R_is": "the best in-span reconstruction error for that readout at eta = 0",
        "T_is": "the transported model's excess over its own eta = 0 error",
        "by_readout": rows,
        "tally": tally,
    }


def exact_lumping(
    generator: Generator,
    colours: Sequence[Any],
    operators: Sequence[Mapping[Tuple[str, int], float]],
    verify_with: Optional[Sequence[float]] = None,
) -> Dict[str, Any]:
    """Coarsest exact strong lumping of the chain that keeps the colouring.

    Strong (ordinary) lumpability asks that every state of a block send the same
    total rate into every other block.  When it holds, the aggregated process is
    an exact Markov chain **for every initial distribution**, so its block count
    is a certified positive/Markov realization dimension for any observable
    constant on the blocks -- not a fitted rank.

    This is the number that separates ``r_linear`` from ``r_positive``.  A
    signed low-rank realization can be tiny while every genuine positive
    realization of the same responses is large, and reporting one dimension for
    both is exactly the conflation #580's frontier update warns about.

    Refining against several operators at once (the baseline generator and the
    intervention tangent) returns the coarsest lumping valid for the whole
    affine family, which is the positive analogue of an intervention-stable
    state.
    """

    size = generator.size
    block_of = {}
    for index, colour in enumerate(sorted(set(colours))):
        block_of[colour] = index
    partition = [block_of[colour] for colour in colours]
    edges: List[List[Tuple[int, Tuple[float, ...]]]] = []
    for source in range(size):
        gathered: Dict[int, List[float]] = {}
        for operator_index, rates in enumerate(operators):
            for move, rate in rates.items():
                image = generator.target[move][source]
                if image is None or rate == 0.0:
                    continue
                row = gathered.setdefault(image, [0.0] * len(operators))
                row[operator_index] += rate
        edges.append([(image, tuple(row)) for image, row in sorted(gathered.items())])

    rounds = 0
    while True:
        rounds += 1
        signatures: List[Tuple[Any, ...]] = []
        for source in range(size):
            totals: Dict[int, List[float]] = {}
            for image, rates in edges[source]:
                row = totals.setdefault(partition[image], [0.0] * len(operators))
                for index, rate in enumerate(rates):
                    row[index] += rate
            signatures.append(
                (partition[source],)
                + tuple(sorted((block, tuple(row)) for block, row in totals.items()))
            )
        ordering = {value: index for index, value in enumerate(sorted(set(signatures)))}
        refined = [ordering[signature] for signature in signatures]
        if len(set(refined)) == len(set(partition)):
            partition = refined
            break
        partition = refined
        if rounds > size:  # pragma: no cover - refinement terminates in <= n rounds
            break
    result = {
        "blocks": len(set(partition)),
        "states": size,
        "refinement_rounds": rounds,
        "compression": len(set(partition)) / size,
    }
    if verify_with is not None:
        result["verified_max_response_drift"] = verify_lumping(
            generator, partition, operators[0], verify_with, LAGS
        )
    return result


def verify_lumping(
    generator: Generator,
    partition: Sequence[int],
    rates: Mapping[Tuple[str, int], float],
    observable: Sequence[float],
    lags: Sequence[float],
) -> float:
    """Check the aggregated chain against the full one, on a function it should keep.

    Strong lumpability says ``exp(t G) f`` stays constant on blocks whenever
    ``f`` is.  So the lumped generator's own evolution has to reproduce the full
    chain's, state by state.  Without this the block count would be a number
    produced by a refinement loop rather than a certified positive realization,
    and a bug in the signature would show up as a spuriously small
    ``r_positive`` -- the direction that would most flatter the conclusion.
    """

    size = generator.size
    blocks = sorted(set(partition))
    index_of = {block: index for index, block in enumerate(blocks)}
    representative = {}
    for state, block in enumerate(partition):
        representative.setdefault(block, state)
    lumped_rows: List[List[Tuple[int, float]]] = []
    for block in blocks:
        source = representative[block]
        entries: Dict[int, float] = {}
        diagonal = 0.0
        for move, rate in rates.items():
            image = generator.target[move][source]
            if image is None or rate == 0.0:
                continue
            entries[index_of[partition[image]]] = (
                entries.get(index_of[partition[image]], 0.0) + rate
            )
            diagonal -= rate
        entries[index_of[block]] = entries.get(index_of[block], 0.0) + diagonal
        lumped_rows.append(sorted(entries.items()))
    lumped_observable = [observable[representative[block]] for block in blocks]
    exit_rate = max(
        sum(rate for move, rate in rates.items() if generator.target[move][state] is not None)
        for state in range(size)
    )
    full = evolve_observables(
        generator.rows(rates), size, exit_rate, [list(observable)], lags
    )
    small = evolve_observables(
        lumped_rows, len(blocks), exit_rate, [lumped_observable], lags
    )
    drift = 0.0
    for lag_index in range(len(lags)):
        for state in range(size):
            drift = max(
                drift,
                abs(
                    full[lag_index][0][state]
                    - small[lag_index][0][index_of[partition[state]]]
                ),
            )
    return drift


#: A large Mersenne prime.  The generator and every declared readout are
#: integer-valued, so linear independence can be decided by exact modular
#: elimination: no pivot tolerance, no conditioning, and a rank that is a
#: certified lower bound on the rational rank.
RANK_PRIME = 2147483647


def observable_reachable_dimension(
    generator: Generator,
    observables: Sequence[Sequence[float]],
    budget: int,
    modulus: int = RANK_PRIME,
) -> Dict[str, Any]:
    """``dim span{G_0^k f}`` -- the dimension an ordinary linear realization needs.

    This is the observability side of the minimal realization order.  It neither
    bounds nor is bounded by ``r_transport``: it is simply the other number, and
    reporting it beside the rank that actually transports is the point.

    Done by exact elimination over a large prime field rather than in floating
    point.  Repeated application of the generator collapses onto the dominant
    direction, so a float elimination silently loses dimensions -- at width 5 it
    returned 26 where the true answer is 42, and that error is in the direction
    that would make the linear rank look close to the transported rank.  Integer
    entries make the modular rank exact, and it is a certified lower bound on
    the rational rank (equal to it unless the prime divides a maximal minor).

    Block Arnoldi with deflation: each level applies the generator to the
    directions the previous level added, which is correct because the generator
    applied to older directions already lies inside the space.
    """

    size = generator.size
    limit = min(budget, size)
    rows: List[List[Tuple[int, int]]] = []
    for row in generator.rows(generator.baseline_rates()):
        rows.append([(column, int(round(value)) % modulus) for column, value in row])

    pivots: Dict[int, List[int]] = {}

    def insert(vector: List[int]) -> Optional[List[int]]:
        working = list(vector)
        for column, pivot in pivots.items():
            if working[column]:
                factor = working[column]
                working = [
                    (working[i] - factor * pivot[i]) % modulus for i in range(size)
                ]
        column = next((i for i in range(size) if working[i]), None)
        if column is None:
            return None
        inverse = pow(working[column], modulus - 2, modulus)
        working = [(value * inverse) % modulus for value in working]
        for other, pivot in list(pivots.items()):
            if pivot[column]:
                factor = pivot[column]
                pivots[other] = [
                    (pivot[i] - factor * working[i]) % modulus for i in range(size)
                ]
        pivots[column] = working
        return working

    frontier = [
        [int(round(value)) % modulus for value in vector] for vector in observables
    ]
    while frontier and len(pivots) < limit:
        added: List[List[int]] = []
        for vector in frontier:
            inserted = insert(vector)
            if inserted is not None:
                added.append(inserted)
            if len(pivots) >= limit:
                break
        if not added:
            break
        frontier = [
            [
                sum(value * vector[column] for column, value in row) % modulus
                for row in rows
            ]
            for vector in added
        ]
    dimension = len(pivots)
    return {
        "dimension": dimension,
        "arithmetic": f"exact modulo the prime {modulus}",
        "budget": limit,
        "reached_the_whole_state_space": dimension >= size,
        "limited_by_budget": dimension >= budget and budget < size,
        "is_a_certified_lower_bound_on_the_rational_rank": True,
    }


def _prefix(basis: Sequence[Sequence[float]], rank: int) -> List[List[float]]:
    """Leading ``rank`` columns.

    Both span constructions are ordered: block Krylov adds powers in the declared
    seed order and orthogonal iteration converges its leading ``k`` columns to
    the dominant ``k``-dimensional invariant subspace, so a prefix is the same
    span the construction would have produced had it been asked for rank ``k``.
    """

    return [list(column) for column in basis[:rank]]


def _rms(values: Sequence[float]) -> Optional[float]:
    values = [value for value in values if value is not None]
    if not values:
        return None
    return math.sqrt(sum(value * value for value in values) / len(values))


def dictionary_scores(entry: Mapping[str, Any]) -> Dict[str, Dict[str, Optional[float]]]:
    """Readout-balanced baseline representability and intervention excess, per dictionary.

    Balanced rather than pooled so that a dictionary's verdict cannot be carried
    by its largest-magnitude member, and split into ``R`` and ``T`` so that a
    dictionary the span never represented is not credited with transporting.
    """

    out: Dict[str, Dict[str, Optional[float]]] = {}
    zero = entry["etas"].get(_eta_key(0.0))
    if zero is None:
        return out
    baseline = zero["fixed_span"].get("contrast_by_readout", {})
    model_zero = zero["transported"].get("contrast_by_readout", {})
    for name, members in DICTIONARIES.items():
        representability = _rms([baseline.get(member) for member in members])
        at_zero = _rms([model_zero.get(member) for member in members])
        excess = None
        for eta in PRIMARY_ETAS:
            block = entry["etas"].get(_eta_key(eta))
            if block is None or at_zero is None:
                continue
            moved = _rms(
                [
                    block["transported"].get("contrast_by_readout", {}).get(member)
                    for member in members
                ]
            )
            if moved is None:
                continue
            gap = moved - at_zero
            excess = gap if excess is None else max(excess, gap)
        out[name] = {
            "R_baseline_representability": representability,
            "transported_at_eta_zero": at_zero,
            "T_intervention_excess": excess,
        }
    return out


def rank_notions(
    generator: Generator,
    intervention: str,
    observables: Mapping[str, Sequence[float]],
    spans: Sequence[Mapping[str, Any]],
    budget: Optional[int] = None,
    verify: bool = False,
) -> Dict[str, Any]:
    """``r_linear``, ``r_positive`` and ``r_transport``, reported as three numbers.

    #580's frontier update asks for exactly this separation, and P398 is one of
    the few places in the repository where all three are computable rather than
    fitted.  They are not the same question: an ordinary signed realization, a
    positive/Markov realization and a frozen realization that survives a change
    of generator can differ by two orders of magnitude, and calling any one of
    them "the state dimension" is the conflation the issue warns against.
    """

    baseline_rates = generator.baseline_rates()
    tangent_rates = generator.coefficients(intervention)
    if budget is None:
        budget = linear_rank_budget(generator.size)
    rows: Dict[str, Any] = {}
    for name, members in DICTIONARIES.items():
        vectors = [observables[member] for member in members]
        colours = [
            tuple(observables[member][state] for member in members)
            for state in range(generator.size)
        ]
        key = (generator.width, name, budget)
        if key not in _LINEAR_RANK_CACHE:
            _LINEAR_RANK_CACHE[key] = observable_reachable_dimension(
                generator, vectors, budget
            )
        linear = _LINEAR_RANK_CACHE[key]
        positive_baseline = exact_lumping(
            generator,
            colours,
            [baseline_rates],
            verify_with=vectors[0] if verify else None,
        )
        positive_family = exact_lumping(
            generator, colours, [baseline_rates, tangent_rates]
        )
        transport: Optional[int] = None
        for span in spans:
            if span.get("status") != "scored" or span["family"] != "krylov":
                continue
            scores = dictionary_scores(span).get(name)
            if scores is None:
                continue
            representability = scores["R_baseline_representability"]
            excess = scores["T_intervention_excess"]
            if (
                representability is not None
                and representability <= PASS_THRESHOLD
                and excess is not None
                and excess <= PASS_THRESHOLD
            ):
                rank = span["requested_rank"]
                transport = rank if transport is None else min(transport, rank)
        rows[name] = {
            "readouts": list(members),
            "r_linear": linear,
            "r_positive_baseline_lumping": positive_baseline,
            "r_positive_whole_affine_family": positive_family,
            "r_transport": transport,
            "r_transport_note": (
                "smallest declared Krylov rank whose readout-balanced baseline "
                "representability and intervention excess are both at or below "
                f"{PASS_THRESHOLD}; None means no declared rank reached it"
            ),
        }
    return rows


def width_experiment(
    generator: Generator,
    intervention: str,
    families: Sequence[str],
    ranks: Sequence[int],
    etas: Sequence[float],
    lags: Sequence[float],
) -> Dict[str, Any]:
    size = generator.size
    baseline_rates = generator.baseline_rates()
    baseline_exit = generator.exit_rate(baseline_rates)
    baseline_rows = generator.rows(baseline_rates)
    tangent = generator.tangent_rows(intervention)
    readouts = list(PRIMARY_READOUTS) + list(HELD_OUT_READOUTS)
    observables = [
        [function(state) for state in generator.states] for _, function in readouts
    ]
    primary_columns = list(range(len(PRIMARY_READOUTS)))
    held_out_columns = list(
        range(len(PRIMARY_READOUTS), len(PRIMARY_READOUTS) + len(HELD_OUT_READOUTS))
    )
    seeds = [[1.0] * size] + [observables[index] for index in primary_columns]
    sources = source_distributions(generator)
    reference_index = [name for name, _ in sources].index(CONTRAST_REFERENCE)
    maximum_rank = min(max(ranks), size)

    cache: Dict[float, Dict[str, Any]] = {}
    for eta in etas:
        rates = generator.rates(intervention, eta)
        rows = generator.rows(rates)
        rate = generator.exit_rate(rates)
        evolved = evolve_observables(rows, size, rate, observables, lags)
        truth = response_tensor(sources, evolved)
        snapshots: List[Tuple[Sequence[float], Sequence[float]]] = [
            (observables[index], matvec(rows, observables[index]))
            for index in primary_columns
        ]
        for block in evolved:
            for index in primary_columns:
                snapshots.append((block[index], matvec(rows, block[index])))
        cache[eta] = {
            "rows": rows,
            "exit_rate": rate,
            "evolved": evolved,
            "truth": truth,
            "contrast": contrast_tensor(truth, reference_index),
            "snapshots": snapshots,
            "squared_norms": [
                [_dot(vector, vector) for vector in block] for block in evolved
            ],
        }

    span_cache: Dict[Tuple[str, float], List[List[float]]] = {}
    for family in families:
        span_cache[(family, 0.0)] = build_span(
            family, generator, baseline_rows, seeds, maximum_rank, baseline_exit
        )
        for eta in etas:
            span_cache[(family, eta)] = build_span(
                family,
                generator,
                cache[eta]["rows"],
                seeds,
                maximum_rank,
                cache[eta]["exit_rate"],
            )

    results: List[Dict[str, Any]] = []
    for family in families:
        frozen_full = span_cache[(family, 0.0)]
        for rank in ranks:
            if rank > len(frozen_full):
                results.append(
                    {
                        "family": family,
                        "requested_rank": rank,
                        "achieved_rank": len(frozen_full),
                        "status": "span_degenerate_at_this_width",
                    }
                )
                continue
            basis = _prefix(frozen_full, rank)
            source_coordinates = [
                project_coordinates(basis, source) for _, source in sources
            ]
            baseline_absolute, baseline_relative, _ = closure_residual(basis, baseline_rows)
            a0 = reduced_operator(basis, baseline_rows)
            b_tangent = reduced_operator(basis, tangent)
            leak_absolute, leak_relative, leak_columns = closure_residual(basis, tangent)
            direction, leak_top = leading_direction(leak_columns)
            augmented = (
                orthonormalize(list(basis) + [direction]) if leak_top > 1e-12 else list(basis)
            )
            augmented_a0 = reduced_operator(augmented, baseline_rows)
            augmented_b = reduced_operator(augmented, tangent)
            _, augmented_relative, _ = closure_residual(augmented, tangent)

            per_eta: Dict[str, Any] = {}
            galerkin_drift = 0.0
            for eta in etas:
                rows = cache[eta]["rows"]
                evolved = cache[eta]["evolved"]
                truth = cache[eta]["truth"]
                truth_contrast = cache[eta]["contrast"]
                galerkin = reduced_operator(basis, rows)
                affine = [
                    [a0[i][j] + eta * b_tangent[i][j] for j in range(rank)]
                    for i in range(rank)
                ]
                galerkin_drift = max(
                    galerkin_drift,
                    max(
                        abs(galerkin[i][j] - affine[i][j])
                        for i in range(rank)
                        for j in range(rank)
                    ),
                )
                refit_basis = _prefix(span_cache[(family, eta)], rank)
                entry: Dict[str, Any] = {
                    "closure_relative_frozen_span": closure_residual(basis, rows)[1],
                    "closure_relative_refit_span": closure_residual(refit_basis, rows)[1],
                    "refit_span_drift_from_frozen": span_drift(basis, refit_basis),
                }
                models: Dict[str, Tuple[List[List[float]], List[List[List[float]]]]] = {}
                oracle_projection = oracle_trajectories(basis, evolved)
                models["fixed_span"] = (list(basis), oracle_projection)
                models["transported"] = (
                    list(basis),
                    reduced_trajectories(basis, affine, observables, lags),
                )
                models["fixed_span_generator_fit"] = (
                    list(basis),
                    reduced_trajectories(
                        basis,
                        _best_generator_fit(
                            basis,
                            cache[eta]["snapshots"],
                            observables,
                            source_coordinates,
                            lags,
                            truth_contrast,
                            list(primary_columns) + list(held_out_columns),
                            reference_index,
                        ),
                        observables,
                        lags,
                    ),
                )
                augmented_affine = [
                    [
                        augmented_a0[i][j] + eta * augmented_b[i][j]
                        for j in range(len(augmented))
                    ]
                    for i in range(len(augmented))
                ]
                oracle_augmented = oracle_trajectories(augmented, evolved)
                models["transported_augmented"] = (
                    list(augmented),
                    reduced_trajectories(augmented, augmented_affine, observables, lags),
                )
                oracle_refit: List[List[List[float]]] = []
                if len(refit_basis) == rank:
                    oracle_refit = oracle_trajectories(refit_basis, evolved)
                    models["refit"] = (list(refit_basis), oracle_refit)
                    models["refit_generator"] = (
                        list(refit_basis),
                        reduced_trajectories(
                            refit_basis,
                            reduced_operator(refit_basis, rows),
                            observables,
                            lags,
                        ),
                    )
                augmented_coordinates = [
                    project_coordinates(augmented, source) for _, source in sources
                ]
                refit_coordinates = (
                    [project_coordinates(refit_basis, source) for _, source in sources]
                    if len(refit_basis) == rank
                    else []
                )
                squared_norms = cache[eta]["squared_norms"]
                for label, (model_basis, trajectories) in models.items():
                    if label == "transported_augmented":
                        coordinates = augmented_coordinates
                        projected = oracle_augmented
                    elif label in ("refit", "refit_generator"):
                        coordinates = refit_coordinates
                        projected = oracle_refit
                    else:
                        coordinates = source_coordinates
                        projected = oracle_projection
                    prediction = response_from_trajectories(coordinates, trajectories)
                    prediction_contrast = contrast_tensor(prediction, reference_index)
                    entry[label] = {
                        "contrast_primary": relative_error(
                            truth_contrast, prediction_contrast, primary_columns
                        ),
                        "contrast_held_out": relative_error(
                            truth_contrast, prediction_contrast, held_out_columns
                        ),
                        "raw_primary": relative_error(truth, prediction, primary_columns),
                        "raw_held_out": relative_error(truth, prediction, held_out_columns),
                        "state_space_primary": reconstruction_relative(
                            squared_norms, projected, trajectories, primary_columns
                        ),
                        "state_space_held_out": reconstruction_relative(
                            squared_norms, projected, trajectories, held_out_columns
                        ),
                    }
                    entry[label]["contrast_balanced_primary"] = balanced_error(
                        truth_contrast, prediction_contrast, primary_columns
                    )
                    entry[label]["contrast_balanced_held_out"] = balanced_error(
                        truth_contrast, prediction_contrast, held_out_columns
                    )
                    if label in ("transported", "fixed_span") and (
                        eta in PRIMARY_ETAS or eta == 0.0
                    ):
                        entry[label]["contrast_by_readout"] = {
                            readouts[column][0]: relative_error(
                                truth_contrast, prediction_contrast, (column,)
                            )
                            for column in range(len(readouts))
                        }
                per_eta[_eta_key(eta)] = entry

            results.append(
                {
                    "family": family,
                    "requested_rank": rank,
                    "achieved_rank": rank,
                    "status": "scored",
                    "baseline": {
                        "closure_absolute": baseline_absolute,
                        "closure_relative": baseline_relative,
                        "galerkin_is_affine_max_drift": galerkin_drift,
                    },
                    "tangent_leakage": {
                        "absolute": leak_absolute,
                        "relative": leak_relative,
                        "leading_singular_value": leak_top,
                        "augmented_relative": augmented_relative,
                        "augmented_rank": len(augmented),
                    },
                    "etas": per_eta,
                }
            )
    named_observables = {
        name: observables[index] for index, (name, _) in enumerate(readouts)
    }
    for entry in results:
        if entry.get("status") != "scored":
            continue
        zero = entry["etas"].get(_eta_key(0.0))
        if zero is None:
            continue
        baseline_by_readout = zero["fixed_span"].get("contrast_by_readout", {})
        model_zero = zero["transported"].get("contrast_by_readout", {})
        excess_by_readout: Dict[str, float] = {}
        for name in named_observables:
            base = model_zero.get(name)
            if base is None:
                continue
            for eta in PRIMARY_ETAS:
                block = entry["etas"].get(_eta_key(eta))
                if block is None:
                    continue
                moved = block["transported"].get("contrast_by_readout", {}).get(name)
                if moved is None:
                    continue
                gap = moved - base
                excess_by_readout[name] = max(excess_by_readout.get(name, gap), gap)
        entry["readout_classification"] = classify_readouts(
            [name for name, _ in readouts],
            [name for name, _ in PRIMARY_READOUTS],
            cache[0.0]["contrast"],
            baseline_by_readout,
            model_zero,
            excess_by_readout,
        )
        entry["dictionary_scores"] = dictionary_scores(entry)

    return {
        "width": generator.width,
        "intervention": intervention,
        "intervention_is_inside_the_baseline_pencil": IN_PENCIL[intervention],
        "rank_notions": rank_notions(
            generator,
            intervention,
            named_observables,
            results,
            verify=intervention == VERDICT_INTERVENTION,
        ),
        "flag_invariance": flag_invariance(
            baseline_rows, tangent, seeds, sorted(set(ranks) | {len(seeds)})
        ),
        "states": size,
        "readouts": {
            "used_to_build_the_span": [name for name, _ in PRIMARY_READOUTS],
            "held_out": [name for name, _ in HELD_OUT_READOUTS],
        },
        "sources": [name for name, _ in sources],
        "contrast_reference_source": CONTRAST_REFERENCE,
        "spans": results,
    }


def _eta_key(eta: float) -> str:
    return f"{eta:+.4f}"


# --------------------------------------------------------------------------
# Exact controls
# --------------------------------------------------------------------------


def exact_controls(
    generator: Generator, lags: Sequence[float], intervention: str
) -> Dict[str, Any]:
    """Four independent ways the response machinery could be silently wrong.

    Probability could leak (the constant function would stop being conserved),
    the uniformization could be truncated too early (a dense Taylor exponential
    would disagree), the stationary power iteration could stop short (the exact
    rational null vector would disagree), and the whole affine-transport claim
    could be a coding error rather than a theorem (a full-rank span would not
    reproduce the intervened generator exactly).
    """

    size = generator.size
    eta = 0.25
    rates = generator.rates(intervention, eta)
    rows = generator.rows(rates)
    rate = generator.exit_rate(rates)
    constant = [1.0] * size
    conserved = evolve_observables(rows, size, rate, [constant], lags)
    conservation_drift = max(
        abs(value - 1.0) for block in conserved for vector in block for value in vector
    )

    dense = [[0.0] * size for _ in range(size)]
    for source, row in enumerate(rows):
        for destination, value in row:
            dense[source][destination] = value
    observable = [observable_blocks(state) for state in generator.states]
    dense_drift = 0.0
    for lag in lags:
        propagator = dense_expm(dense, lag)
        direct = [
            sum(propagator[i][j] * observable[j] for j in range(size)) for i in range(size)
        ]
        reference = evolve_observables(rows, size, rate, [observable], [lag])[0][0]
        dense_drift = max(
            dense_drift, max(abs(direct[i] - reference[i]) for i in range(size))
        )

    numeric, numeric_residual = stationary_law(
        generator.columns(rates), size, 2.0 * rate
    )
    exact = exact_stationary(
        generator, generator.exact_rates(intervention, Fraction(1, 4))
    )
    stationary_drift = max(abs(float(exact[i]) - numeric[i]) for i in range(size))

    full_basis = orthonormalize(
        [[1.0 if i == j else 0.0 for i in range(size)] for j in range(size)]
    )
    a0 = reduced_operator(full_basis, generator.rows(generator.baseline_rates()))
    b_tangent = reduced_operator(full_basis, generator.tangent_rows(intervention))
    galerkin = reduced_operator(full_basis, rows)
    full_rank_drift = max(
        abs(galerkin[i][j] - (a0[i][j] + eta * b_tangent[i][j]))
        for i in range(size)
        for j in range(size)
    )

    return {
        "width": generator.width,
        "intervention": intervention,
        "eta": eta,
        "constant_function_is_conserved_max_drift": conservation_drift,
        "uniformization_matches_dense_taylor_max_drift": dense_drift,
        "power_iteration_matches_exact_rational_stationary_max_drift": stationary_drift,
        "exact_stationary_left_null_residual": numeric_residual,
        "full_span_transport_is_exact_max_drift": full_rank_drift,
        "passed": (
            conservation_drift < 1e-10
            and dense_drift < 1e-9
            and stationary_drift < 1e-9
            and full_rank_drift < 1e-9
        ),
    }


# --------------------------------------------------------------------------
# The declared decision
# --------------------------------------------------------------------------

#: Declared before the run.  A relative error on the contrast tensor at or
#: below this counts as "predicts the intervention"; anything above does not.
PASS_THRESHOLD = 0.10
#: The random-span control must be at least this many times worse than the
#: transported model, otherwise the metric is too loose to decide anything.
CONTROL_MARGIN = 3.0
VERDICT_FAMILY = "krylov"
VERDICT_INTERVENTION = "uniform_join_minus_detach"
VERDICT_RANK = 6


def _find_span(width_block: Mapping[str, Any], family: str, rank: int) -> Optional[Mapping[str, Any]]:
    for entry in width_block["spans"]:
        if entry["family"] == family and entry["requested_rank"] == rank:
            return entry if entry.get("status") == "scored" else None
    return None


def _worst_primary(entry: Optional[Mapping[str, Any]], label: str, metric: str) -> Optional[float]:
    if entry is None:
        return None
    values = []
    for eta in PRIMARY_ETAS:
        block = entry["etas"].get(_eta_key(eta))
        if block is None or label not in block:
            return None
        values.append(block[label][metric])
    return max(values) if values else None


def _classify(
    refit: Optional[float],
    fixed: Optional[float],
    transported: Optional[float],
    held_out: Optional[float],
    augmented: Optional[float],
) -> str:
    """#580's decision table, with one branch the table did not anticipate.

    The table assumes the declared and the held-out readouts stand or fall
    together.  They need not: a frozen span built from three observables can
    predict those three across an intervention and still miss five microscopic
    readouts it never saw.  That case gets its own label rather than being
    swept into the refinement branch, which would credit an augmentation that
    fixes nothing.
    """

    passes = lambda value: value is not None and value <= PASS_THRESHOLD
    if passes(transported) and passes(held_out):
        return "COMMON_STATE_TRANSPORTS"
    if passes(transported):
        return "STATE_TRANSPORTS_ON_ITS_OWN_READOUTS_ONLY"
    if passes(augmented):
        return "ONE_FROZEN_REFINEMENT_CLOSES_THE_INTERVENTION"
    if passes(fixed):
        return "SPAN_SURVIVES_BUT_THE_GENERATOR_LAW_IS_INCOMPLETE"
    if passes(refit):
        return "LOW_RANK_IS_ENVIRONMENT_SPECIFIC"
    return "NO_SMALL_FROZEN_SPAN_TRANSPORTS"


def _family_reading(
    width_block: Mapping[str, Any], family: str, rank: int, control: Optional[Mapping[str, Any]]
) -> Dict[str, Any]:
    entry = _find_span(width_block, family, rank)
    if entry is None:
        return {"family": family, "rank": rank, "verdict": "NOT_SCORED"}
    refit = _worst_primary(entry, "refit", "contrast_primary")
    fixed = _worst_primary(entry, "fixed_span", "contrast_primary")
    transported = _worst_primary(entry, "transported", "contrast_primary")
    augmented = _worst_primary(entry, "transported_augmented", "contrast_primary")
    held_out = _worst_primary(entry, "transported", "contrast_held_out")
    fitted = _worst_primary(entry, "fixed_span_generator_fit", "contrast_primary")
    control_value = _worst_primary(control, "transported", "contrast_primary")
    drift = max(
        entry["etas"][_eta_key(eta)]["refit_span_drift_from_frozen"] for eta in PRIMARY_ETAS
    )
    flag = width_block.get("flag_invariance", {})
    by_readout: Dict[str, float] = {}
    for eta in PRIMARY_ETAS:
        block = entry["etas"][_eta_key(eta)]["transported"].get("contrast_by_readout", {})
        for name, value in block.items():
            by_readout[name] = max(by_readout.get(name, 0.0), value)
    declared = {name for name, _ in PRIMARY_READOUTS}
    worst_declared = (
        max((row for row in by_readout.items() if row[0] in declared), key=lambda row: row[1])
        if by_readout
        else None
    )
    worst_held_out = (
        max((row for row in by_readout.items() if row[0] not in declared), key=lambda row: row[1])
        if by_readout
        else None
    )
    scaling = {
        _eta_key(eta): entry["etas"][_eta_key(eta)]["transported"]["contrast_primary"]
        for eta in sorted(
            float(key) for key in (float(k) for k in entry["etas"])
        )
        if _eta_key(eta) in entry["etas"]
    }
    baseline_block = entry["etas"].get(_eta_key(0.0))
    at_zero = (
        baseline_block["transported"]["contrast_primary"] if baseline_block else None
    )
    at_zero_held_out = (
        baseline_block["transported"]["contrast_held_out"] if baseline_block else None
    )
    return {
        "E_transported_at_eta_zero": at_zero,
        "E_transported_excess_over_eta_zero": (
            None
            if at_zero is None or transported is None
            else transported - at_zero
        ),
        "E_transported_held_out_excess_over_eta_zero": (
            None
            if at_zero_held_out is None or held_out is None
            else held_out - at_zero_held_out
        ),
        "declared_readouts_stay_inside_the_baseline_flag": flag.get(
            "declared_readouts_stay_inside_the_baseline_flag"
        ),
        "intervention_leakage_out_of_the_baseline_flag": flag.get(
            "at_the_declared_readouts"
        ),
        "family": family,
        "rank": rank,
        "verdict": _classify(refit, fixed, transported, held_out, augmented),
        "E_refit": refit,
        "E_fixed_span": fixed,
        "E_transported": transported,
        "E_transported_held_out_readouts": held_out,
        "E_transported_with_one_frozen_refinement": augmented,
        "E_fixed_span_generator_fit": fitted,
        "E_random_span_control": control_value,
        "readout_bins": (entry.get("readout_classification") or {}).get("tally"),
        "readout_classification": (entry.get("readout_classification") or {}).get(
            "by_readout"
        ),
        "dictionary_scores": entry.get("dictionary_scores"),
        "rank_notions": width_block.get("rank_notions"),
        "E_transported_balanced_declared": _worst_primary(
            entry, "transported", "contrast_balanced_primary"
        ),
        "E_transported_balanced_held_out": _worst_primary(
            entry, "transported", "contrast_balanced_held_out"
        ),
        "E_transported_worst_declared_readout": worst_declared,
        "E_transported_worst_held_out_readout": worst_held_out,
        "verdict_by_worst_single_readout": _classify(
            None,
            None,
            worst_declared[1] if worst_declared else None,
            worst_held_out[1] if worst_held_out else None,
            None,
        ),
        "why_that_second_verdict_exists": (
            "the declared metric is a pooled relative Frobenius norm, so a "
            "large-magnitude readout can carry the score while a small one "
            "fails inside it; this re-reads the same numbers one readout at a "
            "time, at the same threshold, and is post-hoc by construction"
        ),
        "transported_error_by_eta": scaling,
        "E_transported_by_readout": (
            entry["etas"][_eta_key(PRIMARY_ETAS[-1])]["transported"].get(
                "contrast_by_readout"
            )
        ),
        "state_space_refit": _worst_primary(entry, "refit", "state_space_primary"),
        "state_space_fixed_span": _worst_primary(entry, "fixed_span", "state_space_primary"),
        "state_space_transported": _worst_primary(
            entry, "transported", "state_space_primary"
        ),
        "refit_span_drift_from_frozen": drift,
        "the_refit_span_actually_moved": drift > 1e-8,
        "baseline_closure_relative": entry["baseline"]["closure_relative"],
        "galerkin_is_affine_max_drift": entry["baseline"]["galerkin_is_affine_max_drift"],
        "tangent_leakage_relative": entry["tangent_leakage"]["relative"],
        "tangent_leakage_after_one_refinement": entry["tangent_leakage"]["augmented_relative"],
        "control_margin_satisfied": (
            control_value is not None
            and transported is not None
            and control_value >= CONTROL_MARGIN * max(transported, 1e-12)
        ),
    }


def decide(blocks: Sequence[Mapping[str, Any]]) -> Dict[str, Any]:
    largest = max(block["width"] for block in blocks)
    readings: List[Dict[str, Any]] = []
    for block in blocks:
        if block["width"] != largest:
            continue
        control = _find_span(block, "random_control", VERDICT_RANK)
        for family in ("krylov", "spectral"):
            reading = _family_reading(block, family, VERDICT_RANK, control)
            reading["intervention"] = block["intervention"]
            reading["intervention_is_inside_the_baseline_pencil"] = block[
                "intervention_is_inside_the_baseline_pencil"
            ]
            readings.append(reading)

    headline = next(
        (
            row
            for row in readings
            if row["family"] == VERDICT_FAMILY
            and row["intervention"] == VERDICT_INTERVENTION
        ),
        readings[0] if readings else {"verdict": "NOT_SCORED"},
    )

    inside = [
        row
        for row in readings
        if row["family"] == VERDICT_FAMILY
        and row.get("intervention_is_inside_the_baseline_pencil")
        and row.get("E_transported") is not None
    ]
    outside = [
        row
        for row in readings
        if row["family"] == VERDICT_FAMILY
        and row.get("intervention_is_inside_the_baseline_pencil") is False
        and row.get("E_transported") is not None
    ]
    pencil: Dict[str, Any] = {
        "what_this_asks": (
            "the baseline generator is J + D and the declared intervention is "
            "J - D, so that intervention never leaves the two-dimensional "
            "operator pencil span{J, D}; a spatially inhomogeneous tilt does. "
            "Comparing them separates transport that reflects a state from "
            "transport that reflects the algebra the state was built from."
        ),
        "worst_in_pencil_transport_error": (
            max(row["E_transported"] for row in inside) if inside else None
        ),
        "worst_out_of_pencil_transport_error": (
            max(row["E_transported"] for row in outside) if outside else None
        ),
        "span_moves_under_the_in_pencil_intervention": (
            any(row["the_refit_span_actually_moved"] for row in inside) if inside else None
        ),
        "span_moves_under_the_out_of_pencil_intervention": (
            any(row["the_refit_span_actually_moved"] for row in outside) if outside else None
        ),
    }
    if inside and outside:
        worst_in = max(row["E_transported"] for row in inside)
        worst_out = max(row["E_transported"] for row in outside)
        if worst_in <= PASS_THRESHOLD < worst_out:
            pencil["reading"] = "TRANSPORT_TRACKS_THE_PENCIL_NOT_THE_STATE"
        elif worst_in <= PASS_THRESHOLD and worst_out <= PASS_THRESHOLD:
            pencil["reading"] = "TRANSPORT_SURVIVES_AN_INTERVENTION_OUTSIDE_THE_PENCIL"
        else:
            pencil["reading"] = "NO_TRANSPORT_TO_ATTRIBUTE"
        pencil["out_over_in_ratio"] = worst_out / max(worst_in, 1e-12)

    return {
        "verdict": headline["verdict"],
        "scored_at": {
            "width": largest,
            "headline_family": VERDICT_FAMILY,
            "headline_intervention": VERDICT_INTERVENTION,
            "rank": VERDICT_RANK,
            "etas": list(PRIMARY_ETAS),
            "metric": "relative Frobenius error on the source-contrast response tensor",
        },
        "per_family": readings,
        "pencil_attribution": pencil,
        "pass_threshold": PASS_THRESHOLD,
        "control_margin": CONTROL_MARGIN,
        "boundary": (
            "A transport pass is a falsification gate on finite low-rank state "
            "claims, not latent-state identification, and it is a statement "
            "about the P398 noncrossing process, not about square-site "
            "Matching One."
        ),
    }


def assemble(
    widths: Sequence[int] = DEFAULT_WIDTHS,
    families: Sequence[str] = ("krylov", "spectral", "random_control"),
    ranks: Sequence[int] = RANKS,
    etas: Sequence[float] = ETAS,
    lags: Sequence[float] = LAGS,
    interventions: Sequence[str] = tuple(INTERVENTIONS),
    control_width: int = 4,
) -> Dict[str, Any]:
    blocks: List[Dict[str, Any]] = []
    gates: List[Dict[str, Any]] = []
    controls: List[Dict[str, Any]] = []
    for width in widths:
        generator = Generator(width)
        for intervention in interventions:
            if width <= 5:
                gates.append(generator_gate(generator, intervention))
            if width == control_width:
                controls.append(exact_controls(generator, lags, intervention))
            blocks.append(
                width_experiment(generator, intervention, families, ranks, etas, lags)
            )
    return {
        "schema": SCHEMA,
        "issue": ISSUE,
        "generated_by": "scripts/p398_intervention_transport.py",
        "frozen_manifest": {
            "state_space": "canonical noncrossing connectivity states, issue #11 codec",
            "baseline_generator": (
                "G_0 = sum over boundary points of (join_cyclic_adjacent + detach) "
                "minus the exit-rate diagonal"
            ),
            "intervention_family": (
                "every rate multiplier is 1 + eta * c(operation, point), exact"
            ),
            "interventions": {
                name: {
                    "inside_the_baseline_pencil": IN_PENCIL[name],
                    "declared_by_issue_580": name == "uniform_join_minus_detach",
                }
                for name in interventions
            },
            "inner_product": (
                "counting measure on the microscopic state space, eta-independent"
            ),
            "frozen_at_eta_zero": [
                "the state functions Phi",
                "the readout functions and their coordinates Phi^T f",
                "the source distributions and their coordinates Phi^T mu",
                "the reduced tangent B = Phi^T H Phi",
            ],
            "allowed_to_move_with_eta": [
                "the microscopic generator G_eta",
                "the stationary law pi_eta (reported, never used to choose functions)",
            ],
            "span_families": list(families),
            "ranks": list(ranks),
            "etas": list(etas),
            "lags": list(lags),
            "declared_before_the_run": [
                "the eta ladder",
                "the lag grid",
                "the readout split into span-building and held-out blocks",
                "the contrast reference source",
                "the pass threshold and the control margin",
                "the three intervention directions and which sit inside the pencil",
            ],
        },
        "exact_generator_gates": gates,
        "exact_controls": controls,
        "blocks": blocks,
        "decision": decide(blocks),
        "scientific_boundary": (
            "Deterministic P398 calculation only.  No Monte Carlo, no width "
            "increase and no mark search.  Intervention invariance is a "
            "falsification gate, not latent-state identification."
        ),
    }


def render(result: Mapping[str, Any]) -> str:
    lines: List[str] = []
    decision = result["decision"]
    scored = decision.get("scored_at", {})
    lines.append(f"issue #{result['issue']}  verdict: {decision['verdict']}")
    lines.append(
        f"  scored at width {scored.get('width')}, rank {scored.get('rank')}, "
        f"eta = {scored.get('etas')}, intervention {scored.get('headline_intervention')}"
    )
    pencil = decision.get("pencil_attribution", {})
    lines.append(
        "  pencil attribution: {reading}  in={inside}  out={outside}".format(
            reading=pencil.get("reading"),
            inside=_format(pencil.get("worst_in_pencil_transport_error")),
            outside=_format(pencil.get("worst_out_of_pencil_transport_error")),
        )
    )
    for reading in decision.get("per_family", []):
        lines.append(
            f"  [{reading.get('intervention')} / {reading['family']}] {reading['verdict']}"
        )
        for key in (
            "E_refit",
            "E_fixed_span",
            "E_transported",
            "E_transported_held_out_readouts",
            "E_transported_with_one_frozen_refinement",
            "E_fixed_span_generator_fit",
            "E_random_span_control",
            "E_transported_at_eta_zero",
            "E_transported_excess_over_eta_zero",
            "E_transported_held_out_excess_over_eta_zero",
            "E_transported_balanced_declared",
            "E_transported_balanced_held_out",
        ):
            lines.append(f"    {key:44} {_format(reading.get(key))}")
        lines.append(
            f"    {'refit_span_drift_from_frozen':44} "
            f"{_format(reading.get('refit_span_drift_from_frozen'))} "
            f"(moved={reading.get('the_refit_span_actually_moved')})"
        )
        bins = reading.get("readout_bins") or {}
        if bins:
            lines.append(
                f"    {'readout bins':44} "
                + "  ".join(f"{key}={value}" for key, value in sorted(bins.items()))
            )
        notions = reading.get("rank_notions") or {}
        for name, row in sorted(notions.items()):
            lines.append(
                f"    {name:44} r_linear="
                f"{row['r_linear']['dimension']}"
                f"{'+' if row['r_linear']['limited_by_budget'] else ''}"
                f"{' (=n)' if row['r_linear']['reached_the_whole_state_space'] else ''}  "
                f"r_positive={row['r_positive_baseline_lumping']['blocks']}"
                f"/{row['r_positive_whole_affine_family']['blocks']}  "
                f"r_transport={row['r_transport']}"
            )
        worst = reading.get("E_transported_worst_declared_readout")
        held = reading.get("E_transported_worst_held_out_readout")
        if worst or held:
            lines.append(
                f"    {'worst single readout':44} declared "
                f"{worst[0] if worst else '-'}={_format(worst[1] if worst else None)}  "
                f"held-out {held[0] if held else '-'}={_format(held[1] if held else None)}"
                f"  -> {reading.get('verdict_by_worst_single_readout')}"
            )
        scaling = reading.get("transported_error_by_eta") or {}
        if scaling:
            lines.append(
                "    {label:44} ".format(label="transported error by eta")
                + "  ".join(f"{key}:{value:.4f}" for key, value in sorted(scaling.items()))
            )
        lines.append(
            f"    {'intervention_leakage_out_of_the_flag':44} "
            f"{_format(reading.get('intervention_leakage_out_of_the_baseline_flag'))} "
            f"(inside={reading.get('declared_readouts_stay_inside_the_baseline_flag')})"
        )
        lines.append(
            f"    {'tangent_leakage_relative':44} "
            f"{_format(reading.get('tangent_leakage_relative'))} -> "
            f"{_format(reading.get('tangent_leakage_after_one_refinement'))} "
            "after one frozen refinement"
        )
    for gate in result["exact_generator_gates"]:
        lines.append(
            f"  generator gate width {gate['width']} {gate['intervention']}: "
            f"passed={gate['passed']} rows={gate['row_sum_failures']} "
            f"negatives={gate['negative_offdiagonal_entries']} "
            f"affine={gate['affine_in_eta_failures']} connected={gate['strongly_connected']}"
        )
    for control in result["exact_controls"]:
        lines.append(
            f"  exact controls width {control['width']} {control['intervention']}: "
            f"passed={control['passed']}"
        )
    for block in result["blocks"]:
        lines.append(
            f"  width {block['width']} ({block['states']} states), "
            f"{block['intervention']} (in pencil="
            f"{block['intervention_is_inside_the_baseline_pencil']}), eta = +0.25"
        )
        for entry in block["spans"]:
            if entry.get("status") != "scored":
                lines.append(
                    f"    {entry['family']:14} r={entry['requested_rank']} {entry['status']}"
                )
                continue
            row = entry["etas"].get(_eta_key(0.25))
            if row is None:
                continue
            lines.append(
                "    {family:14} r={rank:<3} closure={closure:.2e} leak={leak:.3f} "
                "drift={drift:.1e} refit={refit} fixed={fixed:.4f} "
                "transported={transported:.4f} augmented={augmented:.4f}".format(
                    family=entry["family"],
                    rank=entry["requested_rank"],
                    closure=entry["baseline"]["closure_relative"],
                    leak=entry["tangent_leakage"]["relative"],
                    drift=row["refit_span_drift_from_frozen"],
                    refit=(
                        f"{row['refit']['contrast_primary']:.4f}" if "refit" in row else "------"
                    ),
                    fixed=row["fixed_span"]["contrast_primary"],
                    transported=row["transported"]["contrast_primary"],
                    augmented=row["transported_augmented"]["contrast_primary"],
                )
            )
    return "\n".join(lines)


def _format(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:.6f}" if abs(value) >= 1e-4 else f"{value:.3e}"
    return str(value)


def rounded(value: Any, digits: int = 10) -> Any:
    """Round every float in the artifact to ``digits`` significant figures.

    The run is deterministic, so the full double expansion of every entry is
    reproducible from the script; carrying seventeen digits for a hundred
    thousand entries only makes the artifact large enough to be awkward to
    read and to review.
    """

    if isinstance(value, float):
        if value == 0.0 or not math.isfinite(value):
            return value
        exponent = math.floor(math.log10(abs(value)))
        return round(value, digits - 1 - exponent)
    if isinstance(value, dict):
        return {key: rounded(item, digits) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [rounded(item, digits) for item in value]
    return value


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--widths", type=int, nargs="+", default=list(DEFAULT_WIDTHS))
    parser.add_argument("--ranks", type=int, nargs="+", default=list(RANKS))
    parser.add_argument(
        "--families", nargs="+", default=["krylov", "spectral", "random_control"]
    )
    parser.add_argument("--interventions", nargs="+", default=list(INTERVENTIONS))
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args(argv)
    result = assemble(
        widths=tuple(args.widths),
        families=tuple(args.families),
        ranks=tuple(args.ranks),
        interventions=tuple(args.interventions),
        control_width=min(args.widths),
    )
    print(render(result))
    if not args.no_write:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(rounded(result), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
