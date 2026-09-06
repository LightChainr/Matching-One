#!/usr/bin/env python3
"""#588 Phase A: the exact Mori-Zwanzig memory of #580's frozen span.

#580 left one fork open that its own language cannot close.  A rank-6
observable-Krylov span frozen at ``eta = 0`` transports a declared intervention
for ``+0.006`` excess relative error from width 4 to width 8, and yet fails
badly on five readouts it was not built from -- already at ``eta = 0``.  Two
completely different things produce that pattern:

    (i)  the dictionary is missing an instantaneous state coordinate, or
    (ii) an exact Markov process, projected onto an observable-relative
         subspace, generates *memory*.

They call for opposite repairs.  Under (i) the right move is more state; under
(ii) more state is the wrong object entirely and any "fiber label" found
downstream may be history wearing a label.

P398 is unusually well suited to settle this because the microscopic generator
is known exactly, so the memory kernel is *computed*, not estimated.  With
``P = Phi Phi^T`` the frozen projection, ``Q = I - P``, and ``G`` acting on
functions,

    A = Phi^T G Phi        B y = Phi^T G y        (y in range Q)
    C a = Q G Phi a        D   = Q G Q

    d x_R/dt = A x_R(t) + int_0^t K(t-s) x_R(s) ds + B exp(tD) x_U(0)
    K(tau)   = B exp(tau D) C

``K`` is exactly the dynamics the Markov closure discarded, and the third term
is exactly the part of the readout that the span does not see at ``t = 0``.

Three implementation points are load-bearing.

``D`` is never formed.  At width 8 the unresolved block is 1424 x 1424 and a
dense exponential of it is out of reach in pure Python; but ``K`` is only
``r x r``, so the columns ``c_j = Q G phi_j`` are propagated as full-space
vectors constrained to ``range(Q)`` and ``K`` is read off by two inner
products.  ``exp(tau D) v`` is taken by uniformization on ``I + D/Lambda``
rather than a Taylor series in ``D``, for the same reason #580 uniformizes:
no cancellation at the large ``tau`` the slow modes need.

The block formula is *verified against the propagator*, not trusted from
notation.  The repository stores observables in one convention and the
Mori-Zwanzig literature in the other, and a transposed ``B``/``C`` would
produce a plausible, wrong, and perfectly smooth kernel.

And the span is #580's, unchanged.  Choosing a basis that makes the memory
small would answer a different question.

Boundary: memory here is a property of a chosen projection and dictionary.  It
is not a physical field, a low-order kernel is not evidence for a Jordan block,
and P398 is a calibration model -- nothing here transports to square-site
Matching One without a declared map between microscopic state spaces.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

try:  # pragma: no cover - import shim, exercised both ways in practice
    from scripts.p398_intervention_transport import (
        CONTRAST_REFERENCE,
        DICTIONARIES,
        HELD_OUT_READOUTS,
        INTERVENTIONS,
        IN_PENCIL,
        LAGS,
        POISSON_TAIL,
        PRIMARY_ETAS,
        PRIMARY_READOUTS,
        RANKS,
        Generator,
        _dot,
        _norm,
        _prefix,
        _require,
        balanced_error,
        contrast_tensor,
        dense_expm,
        evolve_observables,
        krylov_span,
        matvec,
        orthonormalize,
        poisson_horizon,
        project_coordinates,
        propagator_ladder,
        reduced_operator,
        relative_error,
        response_tensor,
        rounded,
        solve_dense,
        source_distributions,
        symmetric_eigen,
    )
except ModuleNotFoundError:  # Direct ``python scripts/...`` execution.
    from p398_intervention_transport import (  # type: ignore[no-redef]
        CONTRAST_REFERENCE,
        DICTIONARIES,
        HELD_OUT_READOUTS,
        INTERVENTIONS,
        IN_PENCIL,
        LAGS,
        POISSON_TAIL,
        PRIMARY_ETAS,
        PRIMARY_READOUTS,
        RANKS,
        Generator,
        _dot,
        _norm,
        _prefix,
        _require,
        balanced_error,
        contrast_tensor,
        dense_expm,
        evolve_observables,
        krylov_span,
        matvec,
        orthonormalize,
        poisson_horizon,
        project_coordinates,
        propagator_ladder,
        reduced_operator,
        relative_error,
        response_tensor,
        rounded,
        solve_dense,
        source_distributions,
        symmetric_eigen,
    )


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "p398-projected-memory" / "latest.json"
SCHEMA = "matching-one.p398-projected-memory.v1"
ISSUE = 588

DEFAULT_WIDTHS: Tuple[int, ...] = (4, 5, 6, 7, 8)

#: #580's verdict configuration, quoted rather than re-chosen.  The span is
#: built at ``max(RANKS)`` and then prefixed, exactly as #580 built it, so the
#: rank-6 basis here is the same basis to the last digit.
FROZEN_FAMILY = "krylov"
FROZEN_RANK = 6
FROZEN_DICTIONARY = "D0_additive_local_counts"

#: Baseline plus #580's two primary interventions.  The wider ladder is not
#: repeated: the memory question is about the structure at the declared
#: operating point, and every extra ``eta`` is a full propagation of the
#: unresolved block.
MEMORY_ETAS: Tuple[float, ...] = (-0.25, 0.0, 0.25)
MEMORY_INTERVENTIONS: Tuple[str, ...] = ("uniform_join_minus_detach", "single_point_join")

#: ``K`` is sampled on a uniform grid so that a block Hankel matrix is defined.
#: #580's declared lag grid is a subset of it, and every reported lag-grid
#: quantity is read off these samples rather than recomputed.
KERNEL_STEP = 0.25
KERNEL_STEPS = 16  # t = 0, 0.25, ..., 4.0
KERNEL_GRID: Tuple[float, ...] = tuple(
    KERNEL_STEP * index for index in range(KERNEL_STEPS + 1)
)

#: Declared before any kernel is looked at.  The numerical-rank tolerance is
#: ``1e-6`` rather than something tighter because the singular values come from
#: an eigendecomposition of ``H^T H``, which caps their accuracy near the square
#: root of machine epsilon: an exactly one-pole kernel returns a second singular
#: value around ``1e-8``, so a tolerance at that level would count arithmetic as
#: a pole.  The energy levels are the honest "effective pole count" and are what
#: the reading rests on, because a kernel generated by a 1424-dimensional
#: unresolved block has full numerical rank unless it is genuinely degenerate.
HANKEL_TOLERANCE = 1e-6
HANKEL_ENERGY_LEVELS: Tuple[float, ...] = (0.99, 0.999)

#: Re-projection cadence for the constrained propagation.  ``Q G v`` lands in
#: ``range(Q)`` exactly, but the ``v`` it is added to carries drift, so the
#: iterate is cleaned periodically rather than never.
REPROJECT_EVERY = 8

#: Grid for the closing control that integrates the generalized Langevin
#: equation.  Run at two step sizes so the residual can be attributed to
#: quadrature rather than to the kernel.
#: Two step sizes wherever the width allows it, so the residual can be shown to
#: fall like ``h^2``; a single coarser pass at the two largest widths, where a
#: fine grid over the unresolved block is the dominant cost of the whole run and
#: the convergence has already been established below.
GLE_SUBSTEPS: Tuple[int, ...] = (4, 8)
GLE_LARGE_SUBSTEPS: Tuple[int, ...] = (4,)
GLE_FULL_LADDER_MAX_WIDTH = 6
GLE_MAX_WIDTH = 8
GLE_ETAS: Tuple[float, ...] = (0.0, 0.25)
GLE_LARGE_ETAS: Tuple[float, ...] = (0.0,)

#: Widths at which the dense resolvent/Schur identity is checked.
RESOLVENT_MAX_WIDTH = 5
RESOLVENT_SHIFTS: Tuple[float, ...] = (0.5, 1.0, 2.0)


# --------------------------------------------------------------------------
# Projection primitives
# --------------------------------------------------------------------------


def project_out(basis: Sequence[Sequence[float]], vector: Sequence[float]) -> List[float]:
    """``Q v = (I - Phi Phi^T) v``, with one repeat of classical Gram-Schmidt.

    The repeat is not decoration.  A single classical pass loses orthogonality
    at roughly the condition number of the basis, and the entire measurement
    here is the size of what is left *outside* the span: an error in ``Q`` shows
    up as a spurious kernel, in the direction that would manufacture memory.
    """

    out = list(vector)
    for _ in range(2):
        coefficients = [_dot(column, out) for column in basis]
        for column, coefficient in zip(basis, coefficients):
            if coefficient == 0.0:
                continue
            for index in range(len(out)):
                out[index] -= coefficient * column[index]
    return out


def unresolved_columns(
    basis: Sequence[Sequence[float]],
    rows: Sequence[Sequence[Tuple[int, float]]],
) -> List[List[float]]:
    """``C = Q G Phi`` -- one full-space column per resolved coordinate."""

    return [project_out(basis, matvec(rows, column)) for column in basis]


def coupling_rank(
    basis: Sequence[Sequence[float]],
    rows: Sequence[Sequence[Tuple[int, float]]],
) -> Dict[str, Any]:
    """Rank of ``C``, and how many of its columns vanish identically.

    This is the structural cap on the kernel: ``rank K(tau) <= rank C`` for
    every ``tau``, so if ``rank C`` is small and width-independent then the
    kernel's *instantaneous* rank is too, for reasons that are arithmetic about
    the Krylov construction rather than anything about the process.

    Reporting it is what stops a small leading Hankel order from being read as a
    discovery.  A rank-6 Krylov prefix over a 3-dimensional seed span contains
    its own first level entirely, so ``Q G`` annihilates the seed directions and
    only the frontier survives -- and the frontier is the rank.
    """

    columns = unresolved_columns(basis, rows)
    norms = [_norm(column) for column in columns]
    largest = max(norms) if norms else 0.0
    gram = [[_dot(left, right) for right in columns] for left in columns]
    spectrum = singular_spectrum(gram)
    top = spectrum[0] if spectrum else 0.0
    rank = sum(1 for value in spectrum if top > 0 and value > 1e-10 * top)
    return {
        "rank": rank,
        "vanishing_columns": sum(1 for value in norms if largest > 0 and value <= 1e-12 * largest),
        "column_norms": norms,
    }


def resolved_operator(
    basis: Sequence[Sequence[float]],
    rows: Sequence[Sequence[Tuple[int, float]]],
) -> List[List[float]]:
    """``A = Phi^T G Phi``."""

    return reduced_operator(basis, rows)


def apply_B(
    basis: Sequence[Sequence[float]],
    rows: Sequence[Sequence[Tuple[int, float]]],
    vector: Sequence[float],
) -> List[float]:
    """``B y = Phi^T G y`` for ``y`` in the unresolved block."""

    image = matvec(rows, vector)
    return [_dot(column, image) for column in basis]


# --------------------------------------------------------------------------
# exp(tau D) without ever forming D
# --------------------------------------------------------------------------


def propagate_unresolved(
    basis: Sequence[Sequence[float]],
    rows: Sequence[Sequence[Tuple[int, float]]],
    size: int,
    uniform_rate: float,
    seeds: Sequence[Sequence[float]],
    lags: Sequence[float],
) -> List[List[List[float]]]:
    """``exp(tau D) v`` for every seed and lag, by uniformization on ``Q``.

    ``D = Q G Q`` is 1424 x 1424 at width 8 and is never built.  The seeds
    already lie in ``range(Q)``, ``Q G v`` lands there exactly, and the
    uniformized step ``v -> v + (Q G v)/Lambda`` therefore preserves the block.
    All lags share one sweep, as in #580's observable evolution: the Poisson
    weights differ, the powers do not.

    Returned as ``[lag][seed][state]``.
    """

    horizon = poisson_horizon(uniform_rate * max(lags))
    means = [uniform_rate * lag for lag in lags]
    weights = [math.exp(-mean) for mean in means]
    accumulated = [[[0.0] * size for _ in seeds] for _ in lags]
    current = [list(seed) for seed in seeds]
    for order in range(horizon + 1):
        for lag_index in range(len(lags)):
            weight = weights[lag_index]
            if weight <= POISSON_TAIL:
                continue
            block = accumulated[lag_index]
            for seed_index, vector in enumerate(current):
                row = block[seed_index]
                for position in range(size):
                    row[position] += weight * vector[position]
        for lag_index, mean in enumerate(means):
            weights[lag_index] *= mean / (order + 1)
        if order == horizon:
            break
        stepped: List[List[float]] = []
        clean = (order % REPROJECT_EVERY) == REPROJECT_EVERY - 1
        for vector in current:
            image = project_out(basis, matvec(rows, vector))
            moved = [vector[i] + image[i] / uniform_rate for i in range(size)]
            stepped.append(project_out(basis, moved) if clean else moved)
        current = stepped
    return accumulated


def memory_kernel(
    basis: Sequence[Sequence[float]],
    rows: Sequence[Sequence[Tuple[int, float]]],
    size: int,
    uniform_rate: float,
    lags: Sequence[float],
) -> List[List[List[float]]]:
    """``K(tau) = B exp(tau D) C`` on the given lag grid, as ``[lag][i][j]``."""

    columns = unresolved_columns(basis, rows)
    propagated = propagate_unresolved(basis, rows, size, uniform_rate, columns, lags)
    out: List[List[List[float]]] = []
    for block in propagated:
        kernel = [[0.0] * len(basis) for _ in basis]
        for j, vector in enumerate(block):
            coordinates = apply_B(basis, rows, vector)
            for i in range(len(basis)):
                kernel[i][j] = coordinates[i]
        out.append(kernel)
    return out


def unresolved_forcing(
    basis: Sequence[Sequence[float]],
    rows: Sequence[Sequence[Tuple[int, float]]],
    size: int,
    uniform_rate: float,
    observables: Sequence[Sequence[float]],
    lags: Sequence[float],
) -> List[List[List[float]]]:
    """``B exp(tau D) Q f`` -- the unresolved initial-condition term.

    For any readout that seeded the Krylov span this is identically zero,
    because ``Q f = 0``.  For a held-out readout it is not, and it is a
    *different* failure from memory: the span does not represent the readout at
    ``t = 0`` at all.  Keeping the two terms apart is the whole point of
    reporting this separately.

    Returned as ``[lag][observable][resolved coordinate]``.
    """

    seeds = [project_out(basis, vector) for vector in observables]
    propagated = propagate_unresolved(basis, rows, size, uniform_rate, seeds, lags)
    return [
        [apply_B(basis, rows, vector) for vector in block] for block in propagated
    ]


# --------------------------------------------------------------------------
# Kernel statistics
# --------------------------------------------------------------------------


def frobenius(matrix: Sequence[Sequence[float]]) -> float:
    return math.sqrt(sum(value * value for row in matrix for value in row))


def _trapezoid(values: Sequence[float], step: float) -> float:
    if len(values) < 2:
        return 0.0
    total = 0.5 * (values[0] + values[-1])
    total += sum(values[1:-1])
    return total * step


def kernel_statistics(
    kernels: Sequence[Sequence[Sequence[float]]],
    grid: Sequence[float],
) -> Dict[str, Any]:
    """Norm profile, integrated weight, decay time and tail mass of ``K``.

    ``M_memory`` has units of ``1/time`` once integrated, the same units as
    ``||A||``, so the ratio of the two is dimensionless and is the number that
    can be compared across widths and interventions.
    """

    norms = [frobenius(kernel) for kernel in kernels]
    step = grid[1] - grid[0]
    weight = _trapezoid(norms, step)
    first_moment = _trapezoid([t * value for t, value in zip(grid, norms)], step)
    half = len(grid) // 2
    tail = _trapezoid(norms[half:], step)
    return {
        "norm_profile": list(norms),
        "norm_at_zero": norms[0],
        "integrated_weight": weight,
        "decay_time": (first_moment / weight) if weight > 0 else None,
        "tail_mass_fraction": (tail / weight) if weight > 0 else None,
        "half_grid_time": grid[half],
    }


def block_hankel(
    kernels: Sequence[Sequence[Sequence[float]]],
) -> List[List[float]]:
    """``H[p][q] = K(t_{p+q+1})`` from the uniform samples, excluding ``t=0``.

    The Hankel matrix of a kernel is what carries its pole count: a kernel that
    is a sum of ``m`` exponentials has Hankel rank ``m`` times the block size,
    whatever its amplitude.  Dropping ``t = 0`` keeps the first sample of the
    ladder at one full step, so every block is a genuine propagated sample.
    """

    samples = list(kernels[1:])
    blocks = len(samples) // 2
    rank = len(kernels[0])
    matrix = [[0.0] * (blocks * rank) for _ in range(blocks * rank)]
    for p in range(blocks):
        for q in range(blocks):
            block = samples[p + q]
            for i in range(rank):
                for j in range(rank):
                    matrix[p * rank + i][q * rank + j] = block[i][j]
    return matrix


def singular_spectrum(matrix: Sequence[Sequence[float]]) -> List[float]:
    """Singular values from the eigenvalues of ``H^T H``.

    The squared formulation halves the achievable precision, which is why the
    declared numerical tolerance is ``1e-8`` and not machine epsilon; the
    effective counts that carry the scientific reading are energy based and are
    insensitive to that.
    """

    size = len(matrix[0])
    gram = [
        [
            sum(matrix[k][i] * matrix[k][j] for k in range(len(matrix)))
            for j in range(size)
        ]
        for i in range(size)
    ]
    values, _ = symmetric_eigen(gram)
    return sorted((math.sqrt(max(value, 0.0)) for value in values), reverse=True)


def hankel_orders(spectrum: Sequence[float]) -> Dict[str, Any]:
    if not spectrum or spectrum[0] == 0.0:
        return {"numerical_rank": 0, "effective_orders": {}, "spectrum": []}
    largest = spectrum[0]
    numerical = sum(1 for value in spectrum if value > HANKEL_TOLERANCE * largest)
    total = sum(value * value for value in spectrum)
    effective: Dict[str, int] = {}
    for level in HANKEL_ENERGY_LEVELS:
        running = 0.0
        count = 0
        for value in spectrum:
            running += value * value
            count += 1
            if running >= level * total:
                break
        effective[f"{level:g}"] = count
    return {
        "numerical_rank": numerical,
        "effective_orders": effective,
        "spectrum": [value / largest for value in spectrum[:16]],
    }


def hankel_subspace(matrix: Sequence[Sequence[float]], order: int) -> List[List[float]]:
    """Leading ``order`` left singular directions of a block Hankel matrix.

    Obtained by orthonormalizing the images of the leading right singular
    vectors, which avoids a second eigensolve and is exact for the directions
    that matter here.
    """

    size = len(matrix[0])
    gram = [
        [
            sum(matrix[k][i] * matrix[k][j] for k in range(len(matrix)))
            for j in range(size)
        ]
        for i in range(size)
    ]
    values, vectors = symmetric_eigen(gram)
    ordering = sorted(range(len(values)), key=lambda index: -values[index])
    images: List[List[float]] = []
    for index in ordering[:order]:
        column = vectors[index]
        image = [sum(matrix[r][c] * column[c] for c in range(size)) for r in range(len(matrix))]
        if _norm(image) > 0:
            images.append(image)
    return orthonormalize(images)


def subspace_residual(
    basis: Sequence[Sequence[float]], matrix: Sequence[Sequence[float]]
) -> Optional[float]:
    """``||(I - U U^T) H||_F / ||H||_F`` -- does the memory structure transport?

    The useful question in A3 is not whether ``K_eta`` changes under an
    intervention; it must.  It is whether the *baseline* memory directions still
    carry it, or whether the memory spectrum changes in kind.
    """

    total = 0.0
    reference = 0.0
    height = len(matrix)
    width = len(matrix[0])
    for column_index in range(width):
        column = [matrix[row][column_index] for row in range(height)]
        reference += _dot(column, column)
        residual = list(column)
        for direction in basis:
            coefficient = _dot(direction, residual)
            for row in range(height):
                residual[row] -= coefficient * direction[row]
        total += _dot(residual, residual)
    if reference == 0.0:
        return None
    return math.sqrt(total / reference)


def relative_kernel_distance(
    left: Sequence[Sequence[Sequence[float]]],
    right: Sequence[Sequence[Sequence[float]]],
) -> Optional[float]:
    numerator = 0.0
    denominator = 0.0
    for block_left, block_right in zip(left, right):
        for row_left, row_right in zip(block_left, block_right):
            for a, b in zip(row_left, row_right):
                numerator += (a - b) * (a - b)
                denominator += b * b
    if denominator == 0.0:
        return None
    return math.sqrt(numerator / denominator)


# --------------------------------------------------------------------------
# The generalized Langevin equation, integrated
# --------------------------------------------------------------------------


def integrate_gle(
    resolved: Sequence[Sequence[float]],
    kernels: Sequence[Sequence[Sequence[float]]],
    forcing: Sequence[Sequence[float]],
    initial: Sequence[float],
    step: float,
) -> List[List[float]]:
    """Trapezoidal solve of the projected equation on a uniform grid.

    Both the ODE and the convolution use the trapezoid rule, so the scheme is
    second order and the residual against the exact resolved trajectory should
    fall by four when the step is halved.  That convergence is the point: it is
    what lets the remaining error be attributed to quadrature rather than to a
    wrong ``B``, ``C`` or ``D``.

    The last convolution node depends on the unknown, which is why ``K(0)``
    appears in the left-hand matrix rather than in the right-hand side.
    """

    rank = len(resolved)
    steps = len(kernels) - 1
    identity = [[1.0 if i == j else 0.0 for j in range(rank)] for i in range(rank)]
    system = [
        [
            identity[i][j]
            - 0.5 * step * resolved[i][j]
            - 0.25 * step * step * kernels[0][i][j]
            for j in range(rank)
        ]
        for i in range(rank)
    ]
    inverse = solve_dense(system, identity)

    def apply(matrix: Sequence[Sequence[float]], vector: Sequence[float]) -> List[float]:
        return [sum(matrix[i][j] * vector[j] for j in range(rank)) for i in range(rank)]

    trajectory = [list(initial)]
    for index in range(steps):
        current = trajectory[index]
        history_now = [0.0] * rank
        history_next = [0.0] * rank
        for node in range(index + 1):
            weight_now = 0.5 if node in (0, index) else 1.0
            weight_next = 0.5 if node == 0 else 1.0
            state = trajectory[node]
            kernel_now = kernels[index - node]
            kernel_next = kernels[index + 1 - node]
            for i in range(rank):
                row_now = kernel_now[i]
                row_next = kernel_next[i]
                accumulate_now = 0.0
                accumulate_next = 0.0
                for j in range(rank):
                    accumulate_now += row_now[j] * state[j]
                    accumulate_next += row_next[j] * state[j]
                history_now[i] += weight_now * accumulate_now
                history_next[i] += weight_next * accumulate_next
        history_now = [step * value for value in history_now]
        history_next = [step * value for value in history_next]
        drift = apply(resolved, current)
        right = [
            current[i]
            + 0.5 * step * (drift[i] + history_now[i] + forcing[index][i])
            + 0.5 * step * (history_next[i] + forcing[index + 1][i])
            for i in range(rank)
        ]
        trajectory.append(apply(inverse, right))
    return trajectory


def markov_only_trajectory(
    resolved: Sequence[Sequence[float]],
    initial: Sequence[float],
    grid: Sequence[float],
) -> List[List[float]]:
    """``exp(t A) x_R(0)`` -- the same closure with the memory simply dropped."""

    rank = len(resolved)
    propagators = [dense_expm(resolved, time) for time in grid]
    return [
        [sum(propagator[i][j] * initial[j] for j in range(rank)) for i in range(rank)]
        for propagator in propagators
    ]


# --------------------------------------------------------------------------
# A1 controls
# --------------------------------------------------------------------------


def dense_rows(matrix: Sequence[Sequence[float]]) -> List[List[Tuple[int, float]]]:
    return [
        [(index, value) for index, value in enumerate(row) if value != 0.0]
        for row in matrix
    ]


def adapted_blocks(
    basis: Sequence[Sequence[float]],
    rows: Sequence[Sequence[Tuple[int, float]]],
    size: int,
) -> Tuple[List[List[float]], List[List[float]], List[List[float]], List[List[float]]]:
    """Dense ``A,B,C,D`` in an orthonormal basis adapted to ``(P,Q)``.

    Only used at the control widths.  Its whole purpose is to compute the same
    kernel by a completely different route -- an explicit ``exp(tau D)`` of an
    explicitly built ``D`` -- so that a transposed ``B``/``C``, which would
    still produce a smooth and plausible kernel, cannot survive.
    """

    complement: List[List[float]] = []
    for index in range(size):
        candidate = [1.0 if position == index else 0.0 for position in range(size)]
        residual = project_out(basis, candidate)
        for direction in complement:
            coefficient = _dot(direction, residual)
            for position in range(size):
                residual[position] -= coefficient * direction[position]
        length = _norm(residual)
        if length > 1e-8:
            complement.append([value / length for value in residual])
        if len(complement) == size - len(basis):
            break
    images_resolved = [matvec(rows, column) for column in basis]
    images_unresolved = [matvec(rows, column) for column in complement]
    resolved = [
        [_dot(basis[i], images_resolved[j]) for j in range(len(basis))]
        for i in range(len(basis))
    ]
    coupling_out = [
        [_dot(basis[i], images_unresolved[j]) for j in range(len(complement))]
        for i in range(len(basis))
    ]
    coupling_in = [
        [_dot(complement[i], images_resolved[j]) for j in range(len(basis))]
        for i in range(len(complement))
    ]
    unresolved = [
        [_dot(complement[i], images_unresolved[j]) for j in range(len(complement))]
        for i in range(len(complement))
    ]
    return resolved, coupling_out, coupling_in, unresolved


def dense_kernel(
    coupling_out: Sequence[Sequence[float]],
    unresolved: Sequence[Sequence[float]],
    coupling_in: Sequence[Sequence[float]],
    lags: Sequence[float],
) -> List[List[List[float]]]:
    out: List[List[List[float]]] = []
    rank = len(coupling_out)
    hidden = len(unresolved)
    for time in lags:
        propagator = dense_expm(unresolved, time)
        middle = [
            [sum(propagator[i][k] * coupling_in[k][j] for k in range(hidden)) for j in range(rank)]
            for i in range(hidden)
        ]
        out.append(
            [
                [sum(coupling_out[i][k] * middle[k][j] for k in range(hidden)) for j in range(rank)]
                for i in range(rank)
            ]
        )
    return out


def resolvent_identity(
    basis: Sequence[Sequence[float]],
    rows: Sequence[Sequence[Tuple[int, float]]],
    size: int,
    shifts: Sequence[float],
) -> float:
    """Schur-complement check ``Phi^T (sI-G)^{-1} Phi = (sI - A - K^(s))^{-1}``.

    A time-free verification of the whole block decomposition at once, with no
    quadrature anywhere in it: if ``A``, ``B``, ``C`` or ``D`` were extracted
    wrongly, or if the resolved/unresolved split were transposed, this identity
    would fail even though every individual matrix would still look reasonable.
    """

    rank = len(basis)
    resolved, coupling_out, coupling_in, unresolved = adapted_blocks(basis, rows, size)
    hidden = len(unresolved)
    dense = [[0.0] * size for _ in range(size)]
    for source, row in enumerate(rows):
        for destination, value in row:
            dense[source][destination] = value
    worst = 0.0
    for shift in shifts:
        shifted = [
            [(shift if i == j else 0.0) - dense[i][j] for j in range(size)]
            for i in range(size)
        ]
        columns = [[column[index] for column in basis] for index in range(size)]
        solved = solve_dense(shifted, columns)
        direct = [
            [
                sum(basis[i][state] * solved[state][j] for state in range(size))
                for j in range(rank)
            ]
            for i in range(rank)
        ]
        shifted_hidden = [
            [(shift if i == j else 0.0) - unresolved[i][j] for j in range(hidden)]
            for i in range(hidden)
        ]
        inner = solve_dense(shifted_hidden, coupling_in)
        kernel_hat = [
            [sum(coupling_out[i][k] * inner[k][j] for k in range(hidden)) for j in range(rank)]
            for i in range(rank)
        ]
        schur = [
            [
                (shift if i == j else 0.0) - resolved[i][j] - kernel_hat[i][j]
                for j in range(rank)
            ]
            for i in range(rank)
        ]
        identity = [[1.0 if i == j else 0.0 for j in range(rank)] for i in range(rank)]
        predicted = solve_dense(schur, identity)
        numerator = 0.0
        denominator = 0.0
        for i in range(rank):
            for j in range(rank):
                numerator += (direct[i][j] - predicted[i][j]) ** 2
                denominator += direct[i][j] ** 2
        if denominator > 0:
            worst = max(worst, math.sqrt(numerator / denominator))
    return worst


def exact_controls(width: int = 4) -> Dict[str, Any]:
    """A1.  Every one of these has a known answer before it is run.

    A kernel routine that silently returns something smooth and plausible is
    the failure mode here, so the controls are chosen to be ones a plausible
    wrong implementation would fail: a projection with no complement, a span
    that is exactly closed under the generator, a planted block with a known
    answer, and two independent recomputations of the same kernel.
    """

    generator = Generator(width)
    size = generator.size
    rows = generator.rows(generator.baseline_rates())
    rate = generator.exit_rate(generator.baseline_rates())
    lags = (0.5, 1.0, 2.0)

    full_basis = orthonormalize(
        [[1.0 if i == j else 0.0 for j in range(size)] for i in range(size)]
    )
    full_kernel = memory_kernel(full_basis, rows, size, rate, lags)
    full_rank_norm = max(frobenius(block) for block in full_kernel)

    observables, _, _ = declared_observables(generator)
    seeds = [[1.0] * size] + observables[:len(PRIMARY_READOUTS)]
    closed_basis = krylov_span(rows, seeds, size)
    closed_kernel = memory_kernel(closed_basis, rows, size, rate, lags)
    closed_norm = max(frobenius(block) for block in closed_kernel)
    closed_reference = frobenius(resolved_operator(closed_basis, rows))

    planted = synthetic_control()

    frozen = frozen_span(generator, rows)
    dense_blocks = adapted_blocks(frozen, rows, size)
    dense = dense_kernel(dense_blocks[1], dense_blocks[3], dense_blocks[2], lags)
    measured = memory_kernel(frozen, rows, size, rate, lags)
    agreement = relative_kernel_distance(measured, dense)

    columns = unresolved_columns(frozen, rows)
    propagated = propagate_unresolved(frozen, rows, size, rate, columns, lags)
    # Measured against the seed, not against the propagated vector.  ``exp(tau D) c``
    # has decayed by orders of magnitude at the far end of the grid, so dividing by
    # its own norm reports ordinary roundoff as a large leakage and would hide a real
    # one behind a threshold loosened to accommodate it.
    scales = [max(_norm(column), 1e-300) for column in columns]
    leakage = 0.0
    for block in propagated:
        for index, vector in enumerate(block):
            resolved_part = math.sqrt(
                sum(_dot(column, vector) ** 2 for column in frozen)
            )
            leakage = max(leakage, resolved_part / scales[index])

    return {
        "width": width,
        "states": size,
        "full_rank_projection_kernel_norm": full_rank_norm,
        "generator_closed_span_dimension": len(closed_basis),
        "generator_closed_span_kernel_norm": closed_norm,
        "generator_closed_span_reference_norm": closed_reference,
        "synthetic_planted_block_relative_error": planted,
        "frozen_span_dense_versus_matrix_free": agreement,
        "unresolved_propagation_leakage_into_span": leakage,
        "resolvent_schur_identity": {
            str(control_width): resolvent_identity(
                frozen_span(
                    Generator(control_width),
                    Generator(control_width).rows(
                        Generator(control_width).baseline_rates()
                    ),
                ),
                Generator(control_width).rows(Generator(control_width).baseline_rates()),
                Generator(control_width).size,
                RESOLVENT_SHIFTS,
            )
            for control_width in range(4, RESOLVENT_MAX_WIDTH + 1)
        },
        "passed": (
            full_rank_norm < 1e-10
            and closed_norm < 1e-8 * max(closed_reference, 1.0)
            and (planted is not None and planted < 1e-10)
            and (agreement is not None and agreement < 1e-9)
            and leakage < 1e-10
        ),
    }


def synthetic_control() -> Optional[float]:
    """A planted ``(A,B,C,D)`` whose kernel is known in closed form.

    Deliberately not a generator: if the matrix-free path secretly relied on
    row sums vanishing, or on the uniformization rate being a genuine exit
    rate, this control would expose it.
    """

    size = 9
    rank = 3
    state = 20260906
    matrix: List[List[float]] = []
    for _ in range(size):
        row = []
        for _ in range(size):
            state = (1103515245 * state + 12345) % (1 << 31)
            row.append(state / float(1 << 30) - 1.0)
        matrix.append(row)
    for index in range(size):
        matrix[index][index] -= 4.0
    rows = dense_rows(matrix)
    basis = [[1.0 if i == j else 0.0 for j in range(size)] for i in range(rank)]
    lags = (0.25, 1.0, 2.0)
    coupling_out = [[matrix[i][rank + j] for j in range(size - rank)] for i in range(rank)]
    coupling_in = [[matrix[rank + i][j] for j in range(rank)] for i in range(size - rank)]
    unresolved = [
        [matrix[rank + i][rank + j] for j in range(size - rank)] for i in range(size - rank)
    ]
    expected = dense_kernel(coupling_out, unresolved, coupling_in, lags)
    measured = memory_kernel(basis, rows, size, 12.0, lags)
    return relative_kernel_distance(measured, expected)


# --------------------------------------------------------------------------
# The frozen #580 configuration
# --------------------------------------------------------------------------


def declared_observables(
    generator: Generator,
) -> Tuple[List[List[float]], List[str], List[int]]:
    readouts = list(PRIMARY_READOUTS) + list(HELD_OUT_READOUTS)
    observables = [
        [function(state) for state in generator.states] for _, function in readouts
    ]
    names = [name for name, _ in readouts]
    held_out = list(range(len(PRIMARY_READOUTS), len(readouts)))
    return observables, names, held_out


def frozen_span(
    generator: Generator, baseline_rows: Sequence[Sequence[Tuple[int, float]]]
) -> List[List[float]]:
    """#580's rank-6 baseline Krylov span, rebuilt by the identical recipe.

    Built at ``max(RANKS)`` and then prefixed, because that is how #580 built
    it: block Krylov adds powers in declared seed order, so the prefix of the
    rank-12 construction and a direct rank-6 construction are the same span,
    and quoting the same code path removes any doubt about which one this is.
    """

    size = generator.size
    observables, _, _ = declared_observables(generator)
    seeds = [[1.0] * size] + observables[: len(PRIMARY_READOUTS)]
    full = krylov_span(baseline_rows, seeds, min(max(RANKS), size))
    return _prefix(full, min(FROZEN_RANK, len(full)))


def _grid_index(time: float) -> int:
    index = int(round(time / KERNEL_STEP))
    _require(abs(index * KERNEL_STEP - time) < 1e-12, "declared lag is not on the kernel grid")
    return index


def _response(
    source_coordinates: Sequence[Sequence[float]],
    trajectories: Sequence[Sequence[Sequence[float]]],
) -> List[List[List[float]]]:
    """``<Phi^T mu, v(t)>`` over sources and readouts, as ``[lag][source][readout]``."""

    return [
        [[_dot(source, vector) for vector in block] for source in source_coordinates]
        for block in trajectories
    ]


def term_decomposition(
    resolved: Sequence[Sequence[float]],
    kernels: Sequence[Sequence[Sequence[float]]],
    coordinates: Sequence[Sequence[float]],
    forcing: Optional[Sequence[Sequence[float]]],
    grid: Sequence[float],
) -> Dict[str, Optional[float]]:
    """Size of each term of the projected equation, for one readout.

    The kernel belongs to the span, not to a readout; what differs between
    readouts is which trajectory the kernel is convolved against and how much
    of the readout the span fails to see at ``t = 0``.  So the readout-resolved
    object is this decomposition, and it is what separates "the span generates
    memory for this readout" from "the span does not contain this readout".

    ``coordinates[k]`` is the exact ``Phi^T exp(t_k G) f``, so the convolution
    is quadrature on the true trajectory rather than on a simulated one.
    """

    rank = len(resolved)
    step = grid[1] - grid[0]
    drift_norms: List[float] = []
    memory_norms: List[float] = []
    forcing_norms: List[float] = []
    for index in range(len(grid)):
        state = coordinates[index]
        drift = [sum(resolved[i][j] * state[j] for j in range(rank)) for i in range(rank)]
        drift_norms.append(_norm(drift))
        memory = [0.0] * rank
        for node in range(index + 1):
            weight = 0.5 if node in (0, index) and index > 0 else (1.0 if index > 0 else 0.0)
            if weight == 0.0:
                continue
            kernel = kernels[index - node]
            past = coordinates[node]
            for i in range(rank):
                row = kernel[i]
                memory[i] += weight * step * sum(row[j] * past[j] for j in range(rank))
        memory_norms.append(_norm(memory))
        if forcing is not None:
            forcing_norms.append(_norm(forcing[index]))

    def _rms(values: Sequence[float]) -> float:
        return math.sqrt(sum(value * value for value in values) / len(values))

    drift = _rms(drift_norms)
    return {
        "drift_rms": drift,
        "memory_rms": _rms(memory_norms),
        "forcing_rms": _rms(forcing_norms) if forcing is not None else 0.0,
        "memory_over_drift": (_rms(memory_norms) / drift) if drift > 0 else None,
        "forcing_over_drift": (
            (_rms(forcing_norms) / drift) if (forcing is not None and drift > 0) else 0.0
        ),
    }


def width_memory(
    generator: Generator,
    intervention: str,
    etas: Sequence[float],
) -> Dict[str, Any]:
    """A2/A3 for one width and one intervention on the frozen span."""

    size = generator.size
    baseline_rates = generator.baseline_rates()
    baseline_rows = generator.rows(baseline_rates)
    basis = frozen_span(generator, baseline_rows)
    rank = len(basis)
    observables, names, held_out_columns = declared_observables(generator)
    declared_columns = list(range(len(PRIMARY_READOUTS)))
    sources = source_distributions(generator)
    reference_index = [name for name, _ in sources].index(CONTRAST_REFERENCE)
    source_coordinates = [project_coordinates(basis, source) for _, source in sources]
    lag_indices = [_grid_index(lag) for lag in LAGS]

    representation = []
    for name, vector in zip(names, observables):
        length = _norm(vector)
        residual = _norm(project_out(basis, vector))
        representation.append(
            {
                "readout": name,
                "declared": name in DICTIONARIES[FROZEN_DICTIONARY],
                "unresolved_initial_fraction": (residual / length) if length > 0 else None,
            }
        )

    # The two declared interventions are not the same size.  ``uniform_join
    # _minus_detach`` tilts all ``2w`` moves; ``single_point_join`` tilts one.
    # At equal eta the in-pencil perturbation is therefore several times larger
    # in operator norm, and a raw comparison of how far each moves the kernel
    # would report that size difference as a pencil effect.  The tangent norms
    # are recorded so the comparison can be made per unit of perturbation.
    tangent_norm = _sparse_frobenius(generator.tangent_rows(intervention))
    generator_norm = _sparse_frobenius(baseline_rows)
    coupling = coupling_rank(basis, baseline_rows)

    baseline_hankel: Optional[List[List[float]]] = None
    baseline_kernels: Optional[List[List[List[float]]]] = None
    baseline_subspace: Optional[List[List[float]]] = None
    entries: Dict[str, Any] = {}
    _require(0.0 in etas, "the baseline is the reference for every transported quantity")
    # The baseline is computed first because every transported quantity is
    # measured against it; iterating the declared ladder in its natural order
    # would silently compare the first eta against nothing.
    ordered = [0.0] + [value for value in etas if value != 0.0]
    for eta in ordered:
        rates = generator.rates(intervention, eta)
        rows = generator.rows(rates)
        rate = generator.exit_rate(rates)
        resolved = resolved_operator(basis, rows)
        drift_norm = frobenius(resolved)

        kernels = memory_kernel(basis, rows, size, rate, KERNEL_GRID)
        statistics = kernel_statistics(kernels, KERNEL_GRID)
        hankel = block_hankel(kernels)
        spectrum = singular_spectrum(hankel)
        orders = hankel_orders(spectrum)

        evolved = evolve_observables(rows, size, rate, observables, KERNEL_GRID)
        coordinates = [
            [project_coordinates(basis, block[index]) for block in evolved]
            for index in range(len(observables))
        ]
        held_out_seeds = [observables[index] for index in held_out_columns]
        forcing_blocks = (
            unresolved_forcing(basis, rows, size, rate, held_out_seeds, KERNEL_GRID)
            if held_out_seeds
            else []
        )
        forcing_by_readout: Dict[int, List[List[float]]] = {}
        for position, column in enumerate(held_out_columns):
            forcing_by_readout[column] = [block[position] for block in forcing_blocks]

        terms = []
        for column, name in enumerate(names):
            terms.append(
                {
                    "readout": name,
                    "declared": name in DICTIONARIES[FROZEN_DICTIONARY],
                    **term_decomposition(
                        resolved,
                        kernels,
                        coordinates[column],
                        forcing_by_readout.get(column),
                        KERNEL_GRID,
                    ),
                }
            )

        truth = response_tensor(sources, [evolved[index] for index in lag_indices])
        truth_contrast = contrast_tensor(truth, reference_index)
        propagators = propagator_ladder(resolved, LAGS)
        initial = [project_coordinates(basis, vector) for vector in observables]
        markov = [
            [
                [
                    sum(propagator[i][j] * start[j] for j in range(rank))
                    for i in range(rank)
                ]
                for start in initial
            ]
            for propagator in propagators
        ]
        prediction = _response(source_coordinates, markov)
        prediction_contrast = contrast_tensor(prediction, reference_index)

        entry: Dict[str, Any] = {
            "eta": eta,
            "resolved_operator_norm": drift_norm,
            "kernel": statistics,
            "kernel_at_declared_lags": {
                f"{lag:g}": frobenius(kernels[index])
                for lag, index in zip(LAGS, lag_indices)
            },
            "hankel": orders,
            "memory_over_markov": (
                statistics["integrated_weight"] / drift_norm if drift_norm > 0 else None
            ),
            "terms": terms,
            "markov_only_response_error": {
                "declared_pooled": relative_error(
                    truth_contrast, prediction_contrast, declared_columns
                ),
                "declared_balanced": balanced_error(
                    truth_contrast, prediction_contrast, declared_columns
                ),
                "held_out_pooled": relative_error(
                    truth_contrast, prediction_contrast, held_out_columns
                ),
                "held_out_balanced": balanced_error(
                    truth_contrast, prediction_contrast, held_out_columns
                ),
            },
        }
        if eta == 0.0:
            baseline_kernels = kernels
            baseline_hankel = hankel
            order = max(orders["effective_orders"].get("0.999", rank), rank)
            baseline_subspace = hankel_subspace(hankel, min(order, len(hankel)))
        else:
            if baseline_kernels is not None:
                distance = relative_kernel_distance(kernels, baseline_kernels)
                entry["relative_distance_from_baseline_kernel"] = distance
                scale = abs(eta) * tangent_norm / generator_norm
                entry["kernel_distance_per_unit_perturbation"] = (
                    (distance / scale) if (distance is not None and scale > 0) else None
                )
            if baseline_subspace is not None:
                entry["baseline_memory_subspace_residual"] = subspace_residual(
                    baseline_subspace, hankel
                )
        entries[f"{eta:+.3f}"] = entry

    return {
        "width": generator.width,
        "states": size,
        "intervention": intervention,
        "in_pencil": IN_PENCIL[intervention],
        "rank": rank,
        "coupling_rank": coupling,
        "tangent_frobenius_norm": tangent_norm,
        "baseline_generator_frobenius_norm": generator_norm,
        "representation": representation,
        "etas": entries,
    }


def _sparse_frobenius(rows: Sequence[Sequence[Tuple[int, float]]]) -> float:
    return math.sqrt(sum(value * value for row in rows for _, value in row))


def gle_control(
    generator: Generator,
    intervention: str,
    eta: float,
    substeps: int,
) -> Dict[str, Any]:
    """The closing control: does exact memory actually remove the closure error?

    The Markov closure drops the convolution and the unresolved initial term.
    Putting both back must reproduce ``Phi^T exp(tG) f`` exactly, so what is
    left is quadrature.  Running two step sizes and watching the residual fall
    like ``h^2`` is what distinguishes "the kernel is right and the integrator
    is finite" from "the kernel is wrong".
    """

    size = generator.size
    baseline_rows = generator.rows(generator.baseline_rates())
    basis = frozen_span(generator, baseline_rows)
    rank = len(basis)
    observables, names, held_out_columns = declared_observables(generator)
    rates = generator.rates(intervention, eta)
    rows = generator.rows(rates)
    rate = generator.exit_rate(rates)
    resolved = resolved_operator(basis, rows)

    step = KERNEL_STEP / substeps
    grid = [step * index for index in range(KERNEL_STEPS * substeps + 1)]
    kernels = memory_kernel(basis, rows, size, rate, grid)
    forcing_blocks = unresolved_forcing(basis, rows, size, rate, observables, grid)
    evolved = evolve_observables(rows, size, rate, observables, grid)

    zero_kernels = [
        [[0.0] * rank for _ in range(rank)] for _ in range(len(grid))
    ]
    zero_forcing = [[0.0] * rank for _ in range(len(grid))]

    def _relative(
        trajectory: Sequence[Sequence[float]],
        exact: Sequence[Sequence[float]],
        reference: float,
    ) -> float:
        return (
            math.sqrt(
                sum(
                    sum((a - b) ** 2 for a, b in zip(left, right))
                    for left, right in zip(trajectory, exact)
                )
            )
            / reference
        )

    per_readout: List[Dict[str, Any]] = []
    for column, name in enumerate(names):
        exact = [project_coordinates(basis, block[column]) for block in evolved]
        forcing = [block[column] for block in forcing_blocks]
        initial = exact[0]
        reference = math.sqrt(sum(_dot(vector, vector) for vector in exact))
        if reference == 0.0:
            continue
        # Four closures on one trajectory.  The Markov error is the sum of two
        # completely different defects -- discarded history and a readout the
        # span does not contain at t=0 -- and #588's whole question is which one
        # it is.  Turning each term on alone attributes it exactly, with no
        # appeal to a norm ratio.
        # Every one of the four uses the *same* integrator, so the differences
        # between them are the terms and nothing else.  Scoring the Markov
        # closure with ``exp(tA)`` instead would fold the trapezoid rule's own
        # discretization into the memory share -- at h = 1/16 that was already
        # a few percent, in the direction that flatters memory.  The exponential
        # is kept as a separate check that the integrator reproduces it.
        errors = {
            "markov_exponential": _relative(
                markov_only_trajectory(resolved, initial, grid), exact, reference
            ),
            "markov_only": _relative(
                integrate_gle(resolved, zero_kernels, zero_forcing, initial, step),
                exact,
                reference,
            ),
            "plus_memory": _relative(
                integrate_gle(resolved, kernels, zero_forcing, initial, step),
                exact,
                reference,
            ),
            "plus_forcing": _relative(
                integrate_gle(resolved, zero_kernels, forcing, initial, step),
                exact,
                reference,
            ),
            "plus_both": _relative(
                integrate_gle(resolved, kernels, forcing, initial, step),
                exact,
                reference,
            ),
        }
        per_readout.append(
            {
                "readout": name,
                "declared": name in DICTIONARIES[FROZEN_DICTIONARY],
                **errors,
            }
        )
    return {
        "width": generator.width,
        "intervention": intervention,
        "eta": eta,
        "substeps": substeps,
        "step": step,
        "declared": _attribution(per_readout, True),
        "held_out": _attribution(per_readout, False),
        "per_readout": per_readout,
    }


def _attribution(rows: Sequence[Mapping[str, Any]], declared: bool) -> Dict[str, Any]:
    """Pooled closure error of one readout block, and what each term removes.

    Pooling declared and held-out readouts together would be meaningless here:
    for any readout that seeded the span ``Q f = 0``, so the forcing term is
    identically zero and its 'contribution' is exactly nothing.  A maximum over
    all eight readouts would report that structural zero as evidence.
    """

    selected = [row for row in rows if bool(row["declared"]) is declared]
    if not selected:
        return {}

    def pooled(key: str) -> float:
        return math.sqrt(sum(row[key] ** 2 for row in selected) / len(selected))

    markov = pooled("markov_only")
    memory = pooled("plus_memory")
    forcing = pooled("plus_forcing")
    both = pooled("plus_both")
    return {
        "readouts": [row["readout"] for row in selected],
        "markov_exponential": pooled("markov_exponential"),
        "markov_only": markov,
        "plus_memory_only": memory,
        "plus_forcing_only": forcing,
        "plus_both": both,
        "memory_share": ((markov - memory) / markov) if markov > 0 else None,
        "forcing_share": ((markov - forcing) / markov) if markov > 0 else None,
        "residual_share": (both / markov) if markov > 0 else None,
    }


# --------------------------------------------------------------------------
# The decision, declared before the numbers exist
# --------------------------------------------------------------------------

#: Thresholds fixed in ``notes/analysis-plan-20260906.md`` before this script
#: was run.  They are coarse on purpose: this is a four-way classification, not
#: a measurement, and a rule that needs a third digit to fire is a rule that was
#: chosen after seeing the data.
WEAK_MEMORY = 0.10
STRONG_MEMORY = 0.25
LOW_ORDER_MULTIPLE = 2.0
ORDER_GROWTH_STEP = 4
LARGE_HELD_OUT_ERROR = 0.25
#: A term "dominates" the closure error when turning it on alone removes at
#: least half of it.  The natural half, fixed before widths 6-8 were computed.
DOMINANT_SHARE = 0.5
VERDICT_INTERVENTION = "uniform_join_minus_detach"
VERDICT_ETA = "+0.000"


def _mean(values: Sequence[Optional[float]]) -> Optional[float]:
    present = [value for value in values if value is not None]
    if not present:
        return None
    return sum(present) / len(present)


def _block_terms(entry: Mapping[str, Any], declared: bool) -> Dict[str, Optional[float]]:
    rows = [row for row in entry["terms"] if bool(row["declared"]) is declared]
    return {
        "memory_over_drift": _mean([row["memory_over_drift"] for row in rows]),
        "forcing_over_drift": _mean([row["forcing_over_drift"] for row in rows]),
    }


def decide(
    blocks: Sequence[Mapping[str, Any]],
    attribution: Sequence[Mapping[str, Any]],
) -> Dict[str, Any]:
    """The four-way reading, from the attribution first and the norms second.

    The primary discriminator is the closure-error attribution, not a ratio of
    norms.  ``||int K x|| / ||A x||`` is small here at every width, but ``A``
    carries the fast modes and that ratio therefore understates how much of the
    *error* memory accounts for -- turning each term on alone answers the same
    question without a normalization convention in it.

    Provenance, per GOVERNANCE 2C: the norm-ratio thresholds were fixed in
    ``notes/analysis-plan-20260906.md`` before any kernel existed.  The
    attribution and its 0.5 shares were added after the width 4 and 5 controls
    were read and before widths 6, 7 and 8 were computed; 0.5 is the natural
    half, not a tuned value.
    """

    chosen = sorted(
        (block for block in blocks if block["intervention"] == VERDICT_INTERVENTION),
        key=lambda block: block["width"],
    )
    if not chosen:
        return {"verdict": "NO_DECLARED_INTERVENTION_BLOCK"}
    baselines = [block["etas"][VERDICT_ETA] for block in chosen]
    declared = [_block_terms(entry, True) for entry in baselines]
    held_out = [_block_terms(entry, False) for entry in baselines]
    orders = [entry["hankel"]["effective_orders"].get("0.999") for entry in baselines]
    tails = [entry["kernel"]["tail_mass_fraction"] for entry in baselines]
    decays = [entry["kernel"]["decay_time"] for entry in baselines]
    rank = chosen[-1]["rank"]

    finest: Dict[int, Mapping[str, Any]] = {}
    for row in attribution:
        if row["eta"] != 0.0 or row["intervention"] != VERDICT_INTERVENTION:
            continue
        current = finest.get(row["width"])
        if current is None or row["substeps"] > current["substeps"]:
            finest[row["width"]] = row
    ladder = [finest[width] for width in sorted(finest)]
    memory_shares = [row["held_out"].get("memory_share") for row in ladder]
    forcing_shares = [row["held_out"].get("forcing_share") for row in ladder]
    residual_shares = [row["held_out"].get("residual_share") for row in ladder]
    declared_memory_shares = [row["declared"].get("memory_share") for row in ladder]

    numerical = [entry["hankel"]["numerical_rank"] for entry in baselines]
    couplings = [block["coupling_rank"]["rank"] for block in chosen]
    last_order = orders[-1]
    first_order = orders[0]
    low_order = last_order is not None and last_order <= LOW_ORDER_MULTIPLE * rank
    # Two order notions, and they do not agree.  The energy-effective order is
    # how many poles reproduce the kernel to a stated accuracy; the numerical
    # rank is how many it has at all.  #588's table asks whether "memory
    # rank/tail grows with width" as if that were one question.  It is two, and
    # reporting only the one that happens to saturate would be choosing the
    # answer.
    grew_energy = (
        last_order is not None
        and first_order is not None
        and last_order - first_order >= ORDER_GROWTH_STEP
    )
    grew_numerical = numerical[-1] - numerical[0] >= ORDER_GROWTH_STEP
    grew = grew_energy and grew_numerical
    tail_grew = (
        tails[0] is not None and tails[-1] is not None and tails[-1] > tails[0] + 0.10
    )
    if grew_energy and grew_numerical:
        complexity = "MEMORY_ORDER_GROWS_WITH_WIDTH"
    elif not grew_energy and not grew_numerical:
        complexity = "MEMORY_ORDER_SATURATES"
    else:
        complexity = "BOUNDED_MEMORY_ORDER_IS_ACCURACY_DEPENDENT"
    memory_share = memory_shares[-1] if memory_shares else None
    forcing_share = forcing_shares[-1] if forcing_shares else None

    if grew or tail_grew:
        verdict = "MEMORY_COMPLEXITY_GROWS_WITH_WIDTH"
    elif memory_share is None or forcing_share is None:
        verdict = "MEMORY_STRUCTURE_UNRESOLVED_AT_PHASE_A"
    elif memory_share >= DOMINANT_SHARE and forcing_share < DOMINANT_SHARE and low_order:
        verdict = "COMPACT_DICTIONARY_RELATIVE_MEMORY"
    elif forcing_share >= DOMINANT_SHARE and memory_share < DOMINANT_SHARE:
        verdict = "PROJECTION_OR_READOUT_FAILURE_NOT_HIDDEN_HISTORY"
    elif memory_share < DOMINANT_SHARE and forcing_share < DOMINANT_SHARE:
        # Not a branch #588's table has.  Both terms are individually
        # sub-dominant and together they close the error, which is what an
        # interacting pair looks like: neither repair alone is the object.
        verdict = "MEMORY_AND_UNREPRESENTED_READOUT_ARE_JOINTLY_REQUIRED"
    else:
        verdict = "MEMORY_AND_READOUT_BOTH_DOMINANT_CHECK_ATTRIBUTION"

    transport: Dict[str, List[float]] = {}
    normalized: Dict[str, List[Optional[float]]] = {}
    for block in blocks:
        for key, entry in block["etas"].items():
            if key == VERDICT_ETA:
                continue
            residual = entry.get("baseline_memory_subspace_residual")
            if residual is None:
                continue
            transport.setdefault(block["intervention"], []).append(residual)
            normalized.setdefault(block["intervention"], []).append(
                entry.get("kernel_distance_per_unit_perturbation")
            )
    return {
        "verdict": verdict,
        "declared_intervention": VERDICT_INTERVENTION,
        "widths": [block["width"] for block in chosen],
        "declared_memory_over_drift": [row["memory_over_drift"] for row in declared],
        "held_out_memory_over_drift": [row["memory_over_drift"] for row in held_out],
        "held_out_forcing_over_drift": [row["forcing_over_drift"] for row in held_out],
        "effective_memory_order_p999": orders,
        "numerical_memory_order": numerical,
        "coupling_rank_by_width": couplings,
        "memory_complexity": complexity,
        "memory_complexity_note": (
            "rank K(tau) <= rank C at every tau, and rank C is fixed by the "
            "Krylov construction rather than by the width, so a small leading "
            "order is not by itself a finding about the process."
        ),
        "memory_decay_time": decays,
        "tail_mass_fraction": tails,
        "attribution_widths": sorted(finest),
        "held_out_memory_share": memory_shares,
        "held_out_forcing_share": forcing_shares,
        "held_out_residual_share": residual_shares,
        "declared_memory_share": declared_memory_shares,
        "worst_held_out_balanced_markov_error": max(
            entry["markov_only_response_error"]["held_out_balanced"] for entry in baselines
        ),
        "memory_subspace_transport_residual": {
            name: max(values) for name, values in transport.items()
        },
        "kernel_distance_per_unit_perturbation": {
            name: max(value for value in values if value is not None)
            for name, values in normalized.items()
            if any(value is not None for value in values)
        },
        "thresholds": {
            "weak_memory": WEAK_MEMORY,
            "strong_memory": STRONG_MEMORY,
            "low_order_multiple": LOW_ORDER_MULTIPLE,
            "order_growth_step": ORDER_GROWTH_STEP,
            "large_held_out_error": LARGE_HELD_OUT_ERROR,
            "dominant_share": DOMINANT_SHARE,
        },
    }


# --------------------------------------------------------------------------
# Assembly and rendering
# --------------------------------------------------------------------------


def assemble(
    widths: Sequence[int],
    interventions: Sequence[str],
    etas: Sequence[float],
    run_gle: bool,
) -> Dict[str, Any]:
    controls = exact_controls()
    blocks: List[Dict[str, Any]] = []
    for width in widths:
        generator = Generator(width)
        for intervention in interventions:
            blocks.append(width_memory(generator, intervention, etas))
    gle: List[Dict[str, Any]] = []
    if run_gle:
        for width in widths:
            if width > GLE_MAX_WIDTH:
                continue
            generator = Generator(width)
            large = width > GLE_FULL_LADDER_MAX_WIDTH
            ladder = GLE_LARGE_SUBSTEPS if large else GLE_SUBSTEPS
            eta_ladder = GLE_LARGE_ETAS if large else GLE_ETAS
            for eta in eta_ladder:
                for substeps in ladder:
                    gle.append(
                        gle_control(generator, VERDICT_INTERVENTION, eta, substeps)
                    )
    return {
        "schema": SCHEMA,
        "issue": ISSUE,
        "phase": "A",
        "frozen_from": {
            "issue": 580,
            "family": FROZEN_FAMILY,
            "rank": FROZEN_RANK,
            "dictionary": FROZEN_DICTIONARY,
            "note": "notes/p398-intervention-transport-20260906.md",
        },
        "declared": {
            "widths": list(widths),
            "interventions": list(interventions),
            "etas": list(etas),
            "kernel_grid": list(KERNEL_GRID),
            "declared_lags": list(LAGS),
            "hankel_tolerance": HANKEL_TOLERANCE,
            "hankel_energy_levels": list(HANKEL_ENERGY_LEVELS),
            "gle_substeps": list(GLE_SUBSTEPS),
        },
        "exact_controls": controls,
        "blocks": blocks,
        "gle_control": gle,
        "decision": decide(blocks, gle),
        "claim_boundary": (
            "Mori-Zwanzig memory of one declared projection on one exactly known "
            "finite process.  Not a physical field, not evidence for a Jordan "
            "block, and not transportable to square-site Matching One without a "
            "declared map between microscopic state spaces."
        ),
    }


def render(result: Mapping[str, Any]) -> str:
    lines: List[str] = []
    controls = result["exact_controls"]
    lines.append(
        f"#588 Phase A -- exact projected memory of #580's rank-"
        f"{result['frozen_from']['rank']} span"
    )
    lines.append("")
    lines.append(f"exact controls (width {controls['width']}): passed={controls['passed']}")
    lines.append(
        f"  full-rank projection kernel      {controls['full_rank_projection_kernel_norm']:.3e}"
    )
    lines.append(
        f"  generator-closed span (dim "
        f"{controls['generator_closed_span_dimension']:>3}) {controls['generator_closed_span_kernel_norm']:.3e}"
    )
    lines.append(
        f"  synthetic planted block          {controls['synthetic_planted_block_relative_error']:.3e}"
    )
    lines.append(
        f"  dense D versus matrix-free       {controls['frozen_span_dense_versus_matrix_free']:.3e}"
    )
    lines.append(
        f"  propagation leakage into span    {controls['unresolved_propagation_leakage_into_span']:.3e}"
    )
    for width, value in sorted(controls["resolvent_schur_identity"].items()):
        lines.append(f"  resolvent Schur identity w={width}   {value:.3e}")
    lines.append("")

    for block in result["blocks"]:
        pencil = "in pencil" if block["in_pencil"] else "OUT of pencil"
        lines.append(
            f"width {block['width']:>2} ({block['states']:>4} states)  "
            f"{block['intervention']}  [{pencil}]  rank {block['rank']}  "
            f"rank(C) = {block['coupling_rank']['rank']}  "
            f"||H||/||G_0|| = "
            f"{block['tangent_frobenius_norm'] / block['baseline_generator_frobenius_norm']:.4f}"
        )
        header = (
            "   eta      ||A||    |K(0)|   int|K|  mem/mkv   decay   tail   "
            "ord99 ord999   decl_mem  held_mem  held_frc   E_decl  E_held"
        )
        lines.append(header)
        for key in sorted(block["etas"]):
            entry = block["etas"][key]
            kernel = entry["kernel"]
            declared = _block_terms(entry, True)
            held_out = _block_terms(entry, False)
            errors = entry["markov_only_response_error"]
            lines.append(
                f"  {key:>7}  {entry['resolved_operator_norm']:7.3f}"
                f"  {kernel['norm_at_zero']:7.3f}"
                f"  {kernel['integrated_weight']:7.3f}"
                f"  {_format(entry['memory_over_markov'])}"
                f"  {_format(kernel['decay_time'])}"
                f"  {_format(kernel['tail_mass_fraction'])}"
                f"  {entry['hankel']['effective_orders'].get('0.99', 0):>5}"
                f"  {entry['hankel']['effective_orders'].get('0.999', 0):>6}"
                f"   {_format(declared['memory_over_drift'])}"
                f"  {_format(held_out['memory_over_drift'])}"
                f"  {_format(held_out['forcing_over_drift'])}"
                f"  {errors['declared_balanced']:6.3f}"
                f"  {errors['held_out_balanced']:6.3f}"
            )
            residual = entry.get("baseline_memory_subspace_residual")
            if residual is not None:
                lines.append(
                    f"           baseline memory subspace residual "
                    f"{residual:.4f}   kernel distance "
                    f"{_format(entry.get('relative_distance_from_baseline_kernel'))}"
                    f"   per unit perturbation "
                    f"{_format(entry.get('kernel_distance_per_unit_perturbation'))}"
                )
        lines.append("")

    if result["gle_control"]:
        lines.append("closing control -- which term is the closure error?")
        lines.append(
            "   width  eta    h        block       markov    +memory   +forcing      +both"
            "   mem%  frc%"
        )
        for row in result["gle_control"]:
            for label, key in (("declared", "declared"), ("held-out", "held_out")):
                block = row[key]
                if not block:
                    continue
                lines.append(
                    f"   {row['width']:>5}  {row['eta']:+.2f}  {row['step']:.5f}"
                    f"  {label:>8}  {block['markov_only']:9.3e}"
                    f"  {block['plus_memory_only']:9.3e}"
                    f"  {block['plus_forcing_only']:9.3e}"
                    f"  {block['plus_both']:9.3e}"
                    f"  {_percent(block['memory_share'])}"
                    f"  {_percent(block['forcing_share'])}"
                )
        lines.append("")

    decision = result["decision"]
    lines.append(f"decision: {decision['verdict']}")
    lines.append(f"  widths                    {decision['widths']}")
    lines.append(
        f"  declared memory/drift     "
        f"{[None if v is None else round(v, 4) for v in decision['declared_memory_over_drift']]}"
    )
    lines.append(
        f"  held-out memory/drift     "
        f"{[None if v is None else round(v, 4) for v in decision['held_out_memory_over_drift']]}"
    )
    lines.append(
        f"  held-out forcing/drift    "
        f"{[None if v is None else round(v, 4) for v in decision['held_out_forcing_over_drift']]}"
    )
    lines.append(f"  effective order (99.9%)   {decision['effective_memory_order_p999']}")
    lines.append(f"  numerical order           {decision.get('numerical_memory_order')}")
    lines.append(f"  rank(C) by width          {decision.get('coupling_rank_by_width')}")
    lines.append(f"  memory complexity         {decision.get('memory_complexity')}")
    lines.append(f"  attribution widths        {decision.get('attribution_widths')}")
    lines.append(
        f"  held-out memory share     "
        f"{[None if v is None else round(v, 3) for v in decision.get('held_out_memory_share', [])]}"
    )
    lines.append(
        f"  held-out forcing share    "
        f"{[None if v is None else round(v, 3) for v in decision.get('held_out_forcing_share', [])]}"
    )
    lines.append(
        f"  held-out residual share   "
        f"{[None if v is None else round(v, 3) for v in decision.get('held_out_residual_share', [])]}"
    )
    lines.append(
        f"  tail mass fraction        "
        f"{[None if v is None else round(v, 4) for v in decision['tail_mass_fraction']]}"
    )
    for name, value in sorted(decision["memory_subspace_transport_residual"].items()):
        lines.append(f"  memory subspace residual  {name}: {value:.4f}")
    for name, value in sorted(
        decision.get("kernel_distance_per_unit_perturbation", {}).items()
    ):
        lines.append(f"  kernel move per unit |H|  {name}: {value:.4f}")
    lines.append("")
    lines.append(result["claim_boundary"])
    return "\n".join(lines)


def _percent(value: Optional[float]) -> str:
    if value is None:
        return "   --"
    return f"{100.0 * value:5.1f}"


def _format(value: Optional[float]) -> str:
    if value is None:
        return "     --"
    return f"{value:7.4f}"


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--widths", type=int, nargs="+", default=list(DEFAULT_WIDTHS))
    parser.add_argument(
        "--interventions", nargs="+", default=list(MEMORY_INTERVENTIONS)
    )
    parser.add_argument("--etas", type=float, nargs="+", default=list(MEMORY_ETAS))
    parser.add_argument("--no-gle", action="store_true")
    parser.add_argument("--controls-only", action="store_true")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args(list(argv) if argv is not None else None)

    for name in args.interventions:
        _require(name in INTERVENTIONS, f"unknown intervention {name!r}")

    if args.controls_only:
        controls = exact_controls()
        print(json.dumps(rounded(controls), indent=2, sort_keys=True))
        return 0 if controls["passed"] else 1

    result = assemble(args.widths, args.interventions, args.etas, not args.no_gle)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(rounded(result), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    if args.json:
        print(json.dumps(rounded(result), indent=2, sort_keys=True))
    else:
        print(render(result))
    return 0 if result["exact_controls"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
