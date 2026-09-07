#!/usr/bin/env python3
"""#588 Phases B and C: dictionary growth, and memory against more state.

Phase A established what the closure error *is*: on the readouts that seed the
span it is 96% projected memory, and on the held-out readouts it is memory and
an unrepresented readout in roughly equal parts.  It also established the caveat
that keeps that from being a discovery -- ``rank C = 3`` uniformly, because a
rank-6 Krylov prefix over a 3-dimensional seed span contains its own first
level, so ``rank K(tau) <= 3`` for structural reasons.

Two questions follow, and this script answers both on the same object.

**Phase B.**  Not "what is the state dimension" but "how does the description
grow when the dictionary is refined".  For the predeclared filtration
D0 subset D1 subset D2, five numbers are reported *separately* per dictionary --
``r_linear``, ``E_markov``, ``M_memory``, ``r_memory``, ``r_transport`` -- because
collapsing them into one number is exactly the conflation #580 warned about.
Dictionary membership was closed as of #580; no observable is added after a
kernel has been looked at.

**Phase C.**  Two repairs for the held-out failure, scored per degree of freedom
rather than on raw fit:

    S   add the smallest baseline-only Krylov directions (more state);
    M   keep the small state and approximate K(t) with a fixed low-order
        rational model fitted to baseline data only (more memory);
    B   the owner's adversary -- finite-horizon balanced truncation on the
        zero-sum contrast subspace, transform frozen at eta = 0.

The balanced competitor matters because P398 is explicitly an input/output
problem, and #580 already showed that spectral slowness is not the same thing as
controllability plus observability: the dominant invariant subspace fails at
every rank while an observable-Krylov span works.  If a balanced state of the
same order also transports, then #580's Krylov choice was not special and there
is a genuine compact input/output module; if the balanced order jumps when the
nonlocal readouts enter, the dictionary complexity is real and instantaneous.

Two metric disciplines are carried throughout, because #580 showed a pooled
Frobenius norm is magnitude-weighted and can rotate a verdict: every score is
reported both pooled and readout-balanced, and the balancing transform is built
from a metric fixed before it is computed rather than from whatever scale the
readouts happen to carry.

Boundary: everything here is relative to a declared projection, dictionary and
source set on one exactly known finite process.  None of it is a percolation
statement.
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
        LAGS,
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
        observable_reachable_dimension,
        orthonormalize,
        project_coordinates,
        propagator_ladder,
        relative_error,
        response_tensor,
        rounded,
        solve_dense,
        source_distributions,
        symmetric_eigen,
    )
    from scripts.p398_projected_memory import (
        FROZEN_RANK,
        KERNEL_GRID,
        KERNEL_STEP,
        KERNEL_STEPS,
        block_hankel,
        coupling_rank,
        declared_observables,
        frobenius,
        hankel_orders,
        integrate_gle,
        kernel_statistics,
        memory_kernel,
        project_out,
        singular_spectrum,
        unresolved_forcing,
    )
except ModuleNotFoundError:  # Direct ``python scripts/...`` execution.
    from p398_intervention_transport import (  # type: ignore[no-redef]
        CONTRAST_REFERENCE,
        DICTIONARIES,
        HELD_OUT_READOUTS,
        LAGS,
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
        observable_reachable_dimension,
        orthonormalize,
        project_coordinates,
        propagator_ladder,
        relative_error,
        response_tensor,
        rounded,
        solve_dense,
        source_distributions,
        symmetric_eigen,
    )
    from p398_projected_memory import (  # type: ignore[no-redef]
        FROZEN_RANK,
        KERNEL_GRID,
        KERNEL_STEP,
        KERNEL_STEPS,
        block_hankel,
        coupling_rank,
        declared_observables,
        frobenius,
        hankel_orders,
        integrate_gle,
        kernel_statistics,
        memory_kernel,
        project_out,
        singular_spectrum,
        unresolved_forcing,
    )


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "p398-memory-closure" / "latest.json"
SCHEMA = "matching-one.p398-memory-closure.v1"
ISSUE = 588

DEFAULT_WIDTHS: Tuple[int, ...] = (4, 5, 6, 7, 8)
CLOSURE_ETAS: Tuple[float, ...] = (-0.25, 0.0, 0.25)
INTERVENTION = "uniform_join_minus_detach"

#: Phase B ladder.  The same ranks #580 declared, so a dictionary's r_transport
#: here is comparable to #580's without re-deriving anything.
DICTIONARY_RANKS: Tuple[int, ...] = RANKS
PASS_THRESHOLD = 0.10

#: Phase C.  S adds baseline-only Krylov directions on top of #580's rank 6;
#: M keeps rank 6 and adds poles.  The two ladders are matched in added degrees
#: of freedom so that "error per degree of freedom" is a fair comparison and not
#: an artifact of one ladder being longer.
AUGMENT_RANKS: Tuple[int, ...] = (6, 7, 8, 9, 10, 12)
MEMORY_ORDERS: Tuple[int, ...] = (1, 2, 3, 4, 6, 8)

#: The rational fit is done on a fine uniform grid so the resulting model can be
#: used by the same integrator the exact kernel is used by.  Ten Hankel blocks
#: reach t = 1.25 at this step, roughly nine memory decay times.
CLOSURE_SUBSTEPS = 4
ERA_BLOCKS = 10
CLOSURE_MAX_WIDTH = 6

#: Balanced truncation quadrature.  Finite horizon, on the zero-sum contrast
#: subspace; the stationary mode never enters because the sources are contrasts
#: against the declared reference and therefore sum to zero, which annihilates
#: the constant function exactly.  No infinite-horizon Lyapunov solve is used.
BALANCED_HORIZON_NODES = KERNEL_GRID
BALANCED_ORDERS: Tuple[int, ...] = (3, 4, 6, 8, 12)


# --------------------------------------------------------------------------
# One scorer for every realization in this file
# --------------------------------------------------------------------------


def reduced_realization(
    left: Sequence[Sequence[float]],
    right: Sequence[Sequence[float]],
    rows: Sequence[Sequence[Tuple[int, float]]],
) -> List[List[float]]:
    """``W^T G V``.  For an orthonormal Galerkin span ``W = V`` and this is Galerkin."""

    images = [matvec(rows, column) for column in right]
    return [[_dot(left[i], images[j]) for j in range(len(right))] for i in range(len(left))]


def score_realization(
    left: Sequence[Sequence[float]],
    right: Sequence[Sequence[float]],
    rows: Sequence[Sequence[Tuple[int, float]]],
    observables: Sequence[Sequence[float]],
    sources: Sequence[Tuple[str, Sequence[float]]],
    truth_contrast: Sequence[Sequence[Sequence[float]]],
    reference_index: int,
    columns: Sequence[int],
    lags: Sequence[float] = LAGS,
) -> Dict[str, float]:
    """Frozen-model response error, pooled and readout-balanced.

    A biorthogonal pair is accepted rather than only an orthonormal span,
    because the balanced adversary is a genuine Petrov-Galerkin reduction and
    scoring it by orthogonally projecting its span would hand it a handicap the
    method does not have.  With ``W = V`` orthonormal this reduces exactly to
    #580's Galerkin score.
    """

    rank = len(right)
    operator = reduced_realization(left, right, rows)
    starts = [[_dot(row, vector) for row in left] for vector in observables]
    source_coordinates = [
        [_dot(column, source) for column in right] for _, source in sources
    ]
    propagators = propagator_ladder(operator, lags)
    prediction = []
    for propagator in propagators:
        block = []
        for source in source_coordinates:
            row = []
            for start in starts:
                moved = [
                    sum(propagator[i][j] * start[j] for j in range(rank))
                    for i in range(rank)
                ]
                row.append(_dot(source, moved))
            block.append(row)
        prediction.append(block)
    prediction_contrast = contrast_tensor(prediction, reference_index)
    return {
        "pooled": relative_error(truth_contrast, prediction_contrast, columns),
        "balanced": balanced_error(truth_contrast, prediction_contrast, columns),
    }


def biorthogonal(
    left: Sequence[Sequence[float]], right: Sequence[Sequence[float]]
) -> float:
    """``max |W^T V - I|`` -- a Petrov-Galerkin pair that is not biorthogonal is not one."""

    rank = len(left)
    worst = 0.0
    for i in range(rank):
        for j in range(rank):
            value = _dot(left[i], right[j])
            worst = max(worst, abs(value - (1.0 if i == j else 0.0)))
    return worst


# --------------------------------------------------------------------------
# Phase B -- the growth profile under dictionary refinement
# --------------------------------------------------------------------------


def dictionary_profile(
    generator: Generator,
    name: str,
    members: Sequence[str],
    observables: Sequence[Sequence[float]],
    names: Sequence[str],
    baseline_rows: Sequence[Sequence[Tuple[int, float]]],
    cache: Mapping[float, Mapping[str, Any]],
    sources: Sequence[Tuple[str, Sequence[float]]],
    reference_index: int,
    budget: int,
) -> Dict[str, Any]:
    """``r_linear / E_markov / M_memory / r_memory / r_transport``, reported apart.

    Adding outputs can only make the declared prediction problem harder, so the
    object of interest is the *profile* across D0 subset D1 subset D2, not any
    one entry.  The five numbers are kept separate because they are five
    different questions: what an ordinary linear realization needs, how well a
    frozen one does, how much memory the projection generates, how many poles
    that memory has, and how small a frozen span still transports.
    """

    size = generator.size
    columns = [names.index(member) for member in members]
    member_observables = [observables[index] for index in columns]
    seeds = [[1.0] * size] + member_observables
    seed_span = len(orthonormalize([list(seed) for seed in seeds]))
    linear = observable_reachable_dimension(generator, member_observables, budget)
    full = krylov_span(baseline_rows, seeds, min(max(DICTIONARY_RANKS), size))

    ladder: List[Dict[str, Any]] = []
    transport_rank: Optional[int] = None
    for rank in DICTIONARY_RANKS:
        if rank > len(full):
            ladder.append({"rank": rank, "status": "span_degenerate_at_this_width"})
            continue
        basis = _prefix(full, rank)
        scores = {}
        for eta, entry in sorted(cache.items()):
            scores[f"{eta:+.3f}"] = score_realization(
                basis,
                basis,
                entry["rows"],
                observables,
                sources,
                entry["contrast"],
                reference_index,
                columns,
            )
        worst = max(
            scores[f"{eta:+.3f}"]["balanced"] for eta in PRIMARY_ETAS if f"{eta:+.3f}" in scores
        )
        if transport_rank is None and worst <= PASS_THRESHOLD:
            transport_rank = rank
        ladder.append({"rank": rank, "E_markov": scores, "worst_primary_balanced": worst})

    memory: Dict[str, Any] = {}
    wanted = {FROZEN_RANK}
    if transport_rank is not None:
        wanted.add(transport_rank)
    for rank in sorted(wanted):
        if rank > len(full):
            continue
        basis = _prefix(full, rank)
        zero = cache[0.0]
        kernels = memory_kernel(basis, zero["rows"], size, zero["rate"], KERNEL_GRID)
        statistics = kernel_statistics(kernels, KERNEL_GRID)
        orders = hankel_orders(singular_spectrum(block_hankel(kernels)))
        drift = frobenius(reduced_realization(basis, basis, zero["rows"]))
        memory[str(rank)] = {
            "coupling_rank": coupling_rank(basis, baseline_rows),
            "M_memory": statistics["integrated_weight"],
            "M_memory_over_markov": statistics["integrated_weight"] / drift if drift else None,
            "decay_time": statistics["decay_time"],
            "tail_mass_fraction": statistics["tail_mass_fraction"],
            "r_memory_numerical": orders["numerical_rank"],
            "r_memory_effective": orders["effective_orders"],
        }

    return {
        "dictionary": name,
        "members": list(members),
        "seed_span_dimension": seed_span,
        "r_linear": linear,
        "achieved_span": len(full),
        "ladder": ladder,
        "r_transport": transport_rank,
        "memory": memory,
    }


# --------------------------------------------------------------------------
# Phase C1 -- S, more state
# --------------------------------------------------------------------------


def augmentation_ladder(
    generator: Generator,
    observables: Sequence[Sequence[float]],
    baseline_rows: Sequence[Sequence[Tuple[int, float]]],
    cache: Mapping[float, Mapping[str, Any]],
    sources: Sequence[Tuple[str, Sequence[float]]],
    reference_index: int,
    declared_columns: Sequence[int],
    held_out_columns: Sequence[int],
) -> List[Dict[str, Any]]:
    """More Markov state, added the only way that keeps it prospective.

    The extra directions are the *baseline* Krylov continuation of the same
    declared seeds.  Choosing directions that help the held-out readouts would
    be choosing the repair after seeing what it has to repair, and would make
    the comparison with the memory closure meaningless -- that side is fitted to
    baseline data only.
    """

    size = generator.size
    seeds = [[1.0] * size] + [observables[index] for index in declared_columns]
    full = krylov_span(baseline_rows, seeds, min(max(AUGMENT_RANKS), size))
    out: List[Dict[str, Any]] = []
    for rank in AUGMENT_RANKS:
        if rank > len(full):
            out.append({"rank": rank, "status": "span_degenerate_at_this_width"})
            continue
        basis = _prefix(full, rank)
        entry: Dict[str, Any] = {
            "rank": rank,
            "added_degrees_of_freedom": rank - FROZEN_RANK,
            "coupling_rank": coupling_rank(basis, baseline_rows)["rank"],
        }
        for eta, block in sorted(cache.items()):
            entry[f"{eta:+.3f}"] = {
                "declared": score_realization(
                    basis, basis, block["rows"], observables, sources,
                    block["contrast"], reference_index, declared_columns,
                ),
                "held_out": score_realization(
                    basis, basis, block["rows"], observables, sources,
                    block["contrast"], reference_index, held_out_columns,
                ),
            }
        out.append(entry)
    return out


# --------------------------------------------------------------------------
# Phase C2 -- M, more memory
# --------------------------------------------------------------------------


def era_fit(
    kernels: Sequence[Sequence[Sequence[float]]],
    order: int,
    blocks: int = ERA_BLOCKS,
) -> Optional[Tuple[List[List[float]], List[List[float]], List[List[float]]]]:
    """Ho-Kalman / ERA realization ``K_k = C A^(k-1) B`` from baseline samples only.

    A discrete-time realization on the sampling grid, not a continuous-time pole
    fit, because the model has to be usable by the same convolution integrator
    the exact kernel is used by -- converting through a matrix logarithm would
    put a second approximation between the fit and the comparison, in the
    direction that flatters more state over more memory.
    """

    rank = len(kernels[0])
    samples = list(kernels[1:])
    if len(samples) < 2 * blocks:
        return None
    first = [
        [samples[p + q][i][j] for q in range(blocks) for j in range(rank)]
        for p in range(blocks)
        for i in range(rank)
    ]
    shifted = [
        [samples[p + q + 1][i][j] for q in range(blocks) for j in range(rank)]
        for p in range(blocks)
        for i in range(rank)
    ]
    size = len(first[0])
    gram = [
        [sum(first[k][i] * first[k][j] for k in range(len(first))) for j in range(size)]
        for i in range(size)
    ]
    values, vectors = symmetric_eigen(gram)
    ordering = sorted(range(len(values)), key=lambda index: -values[index])
    keep = [index for index in ordering[:order] if values[index] > 0.0]
    if not keep:
        return None
    singular = [math.sqrt(values[index]) for index in keep]
    right = [vectors[index] for index in keep]
    left = []
    for position, index in enumerate(keep):
        column = right[position]
        image = [
            sum(first[r][c] * column[c] for c in range(size)) for r in range(len(first))
        ]
        left.append([value / singular[position] for value in image])

    kept = len(keep)
    middle = [
        [
            sum(
                left[i][r] * sum(shifted[r][c] * right[j][c] for c in range(size))
                for r in range(len(shifted))
            )
            / math.sqrt(singular[i] * singular[j])
            for j in range(kept)
        ]
        for i in range(kept)
    ]
    inflow = [
        [math.sqrt(singular[i]) * right[i][j] for j in range(rank)] for i in range(kept)
    ]
    outflow = [
        [math.sqrt(singular[j]) * left[j][i] for j in range(kept)] for i in range(rank)
    ]
    return middle, inflow, outflow


def era_kernels(
    model: Tuple[List[List[float]], List[List[float]], List[List[float]]],
    zero_lag: Sequence[Sequence[float]],
    steps: int,
) -> List[List[List[float]]]:
    """Evaluate the fitted realization on the whole grid, keeping the exact ``K(0)``.

    ``K(0) = B C`` is available exactly and costs nothing, so there is no reason
    to make the fit carry it; every fitted order then differs from the exact
    kernel only in the propagated part.
    """

    middle, inflow, outflow = model
    order = len(middle)
    rank = len(outflow)
    out = [[list(row) for row in zero_lag]]
    power = [[1.0 if i == j else 0.0 for j in range(order)] for i in range(order)]
    for _ in range(steps):
        block = [
            [
                sum(
                    outflow[i][a] * sum(power[a][b] * inflow[b][j] for b in range(order))
                    for a in range(order)
                )
                for j in range(rank)
            ]
            for i in range(rank)
        ]
        out.append(block)
        power = [
            [sum(power[i][k] * middle[k][j] for k in range(order)) for j in range(order)]
            for i in range(order)
        ]
    return out


def _response_from_trajectory(
    source_coordinates: Sequence[Sequence[float]],
    trajectories: Sequence[Sequence[Sequence[float]]],
    lag_indices: Sequence[int],
) -> List[List[List[float]]]:
    return [
        [
            [_dot(source, trajectories[readout][index]) for readout in range(len(trajectories))]
            for source in source_coordinates
        ]
        for index in lag_indices
    ]


def memory_closure_ladder(
    generator: Generator,
    observables: Sequence[Sequence[float]],
    baseline_rows: Sequence[Sequence[Tuple[int, float]]],
    cache: Mapping[float, Mapping[str, Any]],
    sources: Sequence[Tuple[str, Sequence[float]]],
    reference_index: int,
    declared_columns: Sequence[int],
    held_out_columns: Sequence[int],
) -> Dict[str, Any]:
    """More memory instead of more state, at the same frozen rank 6.

    The rational model is fitted to the ``eta = 0`` kernel and then frozen --
    poles and residues both -- and applied unchanged at ``eta = +-1/4``.  That
    is the same prospective discipline the state ladder is held to.

    Three variants are scored because they answer different questions.  The
    information-matched one carries no forcing term at all, exactly like the
    state ladder, and is the fair comparison.  The exact-forcing variant is the
    ceiling: it says what a memory closure could reach if the unresolved initial
    condition were modelled too, which Phase A showed is a separate ~43% of the
    held-out error.  The exact-kernel variant is the ceiling for any pole count.
    """

    size = generator.size
    step = KERNEL_STEP / CLOSURE_SUBSTEPS
    grid = [step * index for index in range(KERNEL_STEPS * CLOSURE_SUBSTEPS + 1)]
    lag_indices = [int(round(lag / step)) for lag in LAGS]
    seeds = [[1.0] * size] + [observables[index] for index in declared_columns]
    full = krylov_span(baseline_rows, seeds, min(max(RANKS), size))
    basis = _prefix(full, min(FROZEN_RANK, len(full)))
    rank = len(basis)
    source_coordinates = [project_coordinates(basis, source) for _, source in sources]

    zero = cache[0.0]
    baseline_kernels = memory_kernel(basis, zero["rows"], size, zero["rate"], grid)
    models = {}
    for order in MEMORY_ORDERS:
        fitted = era_fit(baseline_kernels, order)
        if fitted is None:
            continue
        models[order] = era_kernels(fitted, baseline_kernels[0], len(grid) - 1)
    fit_quality = {
        str(order): _kernel_fit_error(baseline_kernels, kernels)
        for order, kernels in models.items()
    }

    zero_kernels = [[[0.0] * rank for _ in range(rank)] for _ in grid]
    results: Dict[str, Any] = {"rank": rank, "baseline_fit_relative_error": fit_quality}
    for eta, block in sorted(cache.items()):
        forcing_blocks = unresolved_forcing(
            basis, block["rows"], size, block["rate"], observables, grid
        )
        exact_forcing = [
            [row[column] for row in forcing_blocks] for column in range(len(observables))
        ]
        zero_forcing = [[0.0] * rank for _ in grid]
        resolved = reduced_realization(basis, basis, block["rows"])
        starts = [project_coordinates(basis, vector) for vector in observables]
        variants: Dict[str, Any] = {}
        candidates: List[Tuple[str, Sequence[Sequence[Sequence[float]]], bool]] = [
            ("markov_only", zero_kernels, False),
            ("exact_kernel", memory_kernel(basis, block["rows"], size, block["rate"], grid), False),
            ("exact_kernel_with_forcing", None, True),
        ]
        exact_here = candidates[1][1]
        candidates[2] = ("exact_kernel_with_forcing", exact_here, True)
        for order, kernels in sorted(models.items()):
            candidates.append((f"poles_{order}", kernels, False))
            candidates.append((f"poles_{order}_with_forcing", kernels, True))
        for label, kernels, with_forcing in candidates:
            trajectories = []
            for column in range(len(observables)):
                forcing = exact_forcing[column] if with_forcing else zero_forcing
                trajectories.append(
                    integrate_gle(resolved, kernels, forcing, starts[column], step)
                )
            prediction = _response_from_trajectory(
                source_coordinates, trajectories, lag_indices
            )
            prediction_contrast = contrast_tensor(prediction, reference_index)
            variants[label] = {
                "declared": {
                    "pooled": relative_error(
                        block["contrast"], prediction_contrast, declared_columns
                    ),
                    "balanced": balanced_error(
                        block["contrast"], prediction_contrast, declared_columns
                    ),
                },
                "held_out": {
                    "pooled": relative_error(
                        block["contrast"], prediction_contrast, held_out_columns
                    ),
                    "balanced": balanced_error(
                        block["contrast"], prediction_contrast, held_out_columns
                    ),
                },
            }
        results[f"{eta:+.3f}"] = variants
    return results


def _kernel_fit_error(
    exact: Sequence[Sequence[Sequence[float]]],
    fitted: Sequence[Sequence[Sequence[float]]],
) -> Optional[float]:
    numerator = 0.0
    denominator = 0.0
    for left, right in zip(exact, fitted):
        for row_left, row_right in zip(left, right):
            for a, b in zip(row_left, row_right):
                numerator += (a - b) * (a - b)
                denominator += a * a
    if denominator == 0.0:
        return None
    return math.sqrt(numerator / denominator)


# --------------------------------------------------------------------------
# Phase C3 -- the balanced adversary
# --------------------------------------------------------------------------


def _long_horizon_growth(operator: Sequence[Sequence[float]]) -> float:
    """``||exp(2 T A_r)||`` at twice the scoring horizon.

    A Petrov-Galerkin reduction of a generator is not guaranteed stable, and an
    unstable reduced model can still score well on a finite lag grid before it
    diverges.  Reporting the growth one horizon past the grid is what stops a
    good score from being an artifact of where the scoring stopped: a value near
    or below 1 means the model decays like the process it replaces.
    """

    horizon = 2.0 * max(LAGS)
    propagator = dense_expm(operator, horizon)
    return math.sqrt(sum(value * value for row in propagator for value in row))


def _rebiorthogonalize(
    left: Sequence[Sequence[float]], right: Sequence[Sequence[float]]
) -> List[List[float]]:
    """Restore ``W^T V = I`` exactly, without changing either span.

    The pair comes out biorthogonal in exact arithmetic, but the transform
    carries a factor ``sigma^{-1/2}`` and the smallest retained Hankel value is
    six orders below the largest, so the defect reaches ``7e-3`` at order 12 --
    enough that ``W^T G V`` would no longer be the balanced reduced generator.
    Replacing ``W`` by ``S^{-T} W`` with ``S = W^T V`` fixes it inside the same
    two spans, so it corrects the coordinates and not the method.
    """

    rank = len(left)
    gram = [[_dot(left[i], right[j]) for j in range(rank)] for i in range(rank)]
    identity = [[1.0 if i == j else 0.0 for j in range(rank)] for i in range(rank)]
    inverse = solve_dense(gram, identity)
    size = len(left[0])
    return [
        [sum(inverse[j][i] * left[j][state] for j in range(rank)) for state in range(size)]
        for i in range(rank)
    ]


def balanced_family(
    generator: Generator,
    observables: Sequence[Sequence[float]],
    columns: Sequence[int],
    cache: Mapping[float, Mapping[str, Any]],
    sources: Sequence[Tuple[str, Sequence[float]]],
    reference_index: int,
    readout_scales: Sequence[float],
) -> Dict[str, Any]:
    """Finite-horizon balanced truncation on the zero-sum contrast subspace.

    The stationary mode is removed by construction rather than by regularization:
    the sources enter as contrasts against the declared reference, so they sum to
    zero, and a zero-sum functional annihilates the constant function exactly.
    Nothing marginally stable is ever fed to a Lyapunov solve, and the horizon is
    the declared lag grid rather than infinity.

    The output metric is fixed before balancing.  #580 showed a pooled Frobenius
    norm is magnitude-weighted, and a balancing transform is exactly the kind of
    object that a rescaled readout would silently rotate, so each readout is
    divided by its own baseline contrast signal before the Gramians are built.
    """

    size = generator.size
    zero = cache[0.0]
    step = KERNEL_GRID[1] - KERNEL_GRID[0]
    weights = [
        math.sqrt(step * (0.5 if index in (0, len(KERNEL_GRID) - 1) else 1.0))
        for index in range(len(KERNEL_GRID))
    ]
    scaled = [
        [value / scale for value in observables[column]]
        for column, scale in zip(columns, readout_scales)
    ]
    forward = evolve_observables(zero["rows"], size, zero["rate"], scaled, KERNEL_GRID)
    reach: List[List[float]] = []
    for index, block in enumerate(forward):
        for vector in block:
            reach.append([weights[index] * value for value in vector])

    reference = sources[reference_index][1]
    contrasts = [
        [source[state] - reference[state] for state in range(size)]
        for name, source in sources
        if name != sources[reference_index][0]
    ]
    backward = evolve_observables(
        zero["columns"], size, zero["rate"], contrasts, KERNEL_GRID
    )
    observe: List[List[float]] = []
    for index, block in enumerate(backward):
        for vector in block:
            observe.append([weights[index] * value for value in vector])

    cross = [[_dot(left, right) for right in reach] for left in observe]
    height = len(cross)
    gram = [
        [sum(cross[i][k] * cross[j][k] for k in range(len(reach))) for j in range(height)]
        for i in range(height)
    ]
    values, vectors = symmetric_eigen(gram)
    ordering = sorted(range(len(values)), key=lambda index: -values[index])
    spectrum = [math.sqrt(max(values[index], 0.0)) for index in ordering]
    largest = spectrum[0] if spectrum else 0.0

    results: Dict[str, Any] = {
        "hankel_singular_values": [
            value / largest for value in spectrum[: min(16, len(spectrum))]
        ],
        "orders": {},
    }
    for order in BALANCED_ORDERS:
        # The Hankel values come from an eigendecomposition of a Gram matrix, so
        # they carry absolute error of order eps * sigma_max^2: below sigma =
        # 1e-6 sigma_max the retained direction is arithmetic, not a Hankel
        # direction, and sigma^{-1/2} then amplifies it into the transform.  An
        # order that cannot be filled is reported as rank deficient rather than
        # filled with noise -- how many input/output coordinates the task
        # actually exposes is the question, so refusing to invent one is part of
        # the answer.
        keep = [index for index in ordering[:order] if values[index] > 1e-12 * max(values)]
        if len(keep) < order:
            results["orders"][str(order)] = {"status": "gramian_rank_deficient"}
            continue
        left_basis: List[List[float]] = []
        right_basis: List[List[float]] = []
        for index in keep:
            sigma = math.sqrt(math.sqrt(max(values[index], 0.0)))
            if sigma == 0.0:
                continue
            u = vectors[index]
            left_basis.append(
                [
                    sum(u[c] * observe[c][state] for c in range(height)) / sigma
                    for state in range(size)
                ]
            )
            v = [
                sum(u[c] * cross[c][k] for c in range(height))
                / math.sqrt(max(values[index], 1e-300))
                for k in range(len(reach))
            ]
            right_basis.append(
                [
                    sum(v[k] * reach[k][state] for k in range(len(reach))) / sigma
                    for state in range(size)
                ]
            )
        if len(left_basis) < order:
            results["orders"][str(order)] = {"status": "gramian_rank_deficient"}
            continue
        raw_defect = biorthogonal(left_basis, right_basis)
        left_basis = _rebiorthogonalize(left_basis, right_basis)
        entry: Dict[str, Any] = {
            "biorthogonality_defect_before_correction": raw_defect,
            "biorthogonality_defect": biorthogonal(left_basis, right_basis),
            "long_horizon_growth": _long_horizon_growth(
                reduced_realization(left_basis, right_basis, cache[0.0]["rows"])
            ),
        }
        for eta, block in sorted(cache.items()):
            entry[f"{eta:+.3f}"] = score_realization(
                left_basis, right_basis, block["rows"], observables, sources,
                block["contrast"], reference_index, columns,
            )
        results["orders"][str(order)] = entry
    return results


# --------------------------------------------------------------------------
# Assembly
# --------------------------------------------------------------------------


def _contrast_scale(
    contrast: Sequence[Sequence[Sequence[float]]], column: int
) -> float:
    total = 0.0
    for block in contrast:
        for row in block:
            total += row[column] * row[column]
    return math.sqrt(total) if total > 0 else 1.0


def width_closure(generator: Generator, run_memory: bool) -> Dict[str, Any]:
    size = generator.size
    baseline_rates = generator.baseline_rates()
    baseline_rows = generator.rows(baseline_rates)
    observables, names, held_out_columns = declared_observables(generator)
    declared_columns = list(range(len(PRIMARY_READOUTS)))
    sources = source_distributions(generator)
    reference_index = [name for name, _ in sources].index(CONTRAST_REFERENCE)

    cache: Dict[float, Dict[str, Any]] = {}
    for eta in CLOSURE_ETAS:
        rates = generator.rates(INTERVENTION, eta)
        rows = generator.rows(rates)
        rate = generator.exit_rate(rates)
        evolved = evolve_observables(rows, size, rate, observables, LAGS)
        truth = response_tensor(sources, evolved)
        cache[eta] = {
            "rows": rows,
            "columns": generator.columns(rates),
            "rate": rate,
            "contrast": contrast_tensor(truth, reference_index),
        }

    budget = size if size <= 429 else 150
    profiles = [
        dictionary_profile(
            generator, name, members, observables, names, baseline_rows,
            cache, sources, reference_index, budget,
        )
        for name, members in DICTIONARIES.items()
    ]

    scales_by_name = {
        name: _contrast_scale(cache[0.0]["contrast"], index)
        for index, name in enumerate(names)
    }
    balanced = {}
    for name, members in DICTIONARIES.items():
        columns = [names.index(member) for member in members]
        balanced[name] = balanced_family(
            generator, observables, columns, cache, sources, reference_index,
            [scales_by_name[member] for member in members],
        )

    block: Dict[str, Any] = {
        "width": generator.width,
        "states": size,
        "intervention": INTERVENTION,
        "phase_b_dictionaries": profiles,
        "phase_c_state_augmentation": augmentation_ladder(
            generator, observables, baseline_rows, cache, sources,
            reference_index, declared_columns, held_out_columns,
        ),
        "phase_c_balanced": balanced,
    }
    if run_memory and generator.width <= CLOSURE_MAX_WIDTH:
        block["phase_c_memory_closure"] = memory_closure_ladder(
            generator, observables, baseline_rows, cache, sources,
            reference_index, declared_columns, held_out_columns,
        )
    return block


def decide(blocks: Sequence[Mapping[str, Any]]) -> Dict[str, Any]:
    """Three readings, kept apart because they are three questions.

    Declared before the numbers: the dictionary reading is about whether the
    description grows when outputs are added; the competitor reading is about
    whether #580's Krylov choice was special; the repair reading is about which
    of more state or more memory buys more per degree of freedom.  A single
    verdict string would have to pick one of the three to be about.
    """

    ordered = sorted(blocks, key=lambda block: block["width"])
    profile: Dict[str, Any] = {}
    for name in DICTIONARIES:
        row = []
        for block in ordered:
            entry = next(
                item for item in block["phase_b_dictionaries"] if item["dictionary"] == name
            )
            memory = entry["memory"].get(str(FROZEN_RANK), {})
            row.append(
                {
                    "width": block["width"],
                    "r_linear": entry["r_linear"].get("dimension"),
                    "r_transport": entry["r_transport"],
                    "M_memory": memory.get("M_memory"),
                    "r_memory_numerical": memory.get("r_memory_numerical"),
                    "r_memory_effective_999": (
                        memory.get("r_memory_effective", {}) or {}
                    ).get("0.999"),
                    "coupling_rank": (memory.get("coupling_rank") or {}).get("rank"),
                }
            )
        profile[name] = row

    dictionary_growth = {}
    for name, row in profile.items():
        transports = [item["r_transport"] for item in row]
        dictionary_growth[name] = {
            "r_transport_by_width": transports,
            "r_linear_by_width": [item["r_linear"] for item in row],
            "r_memory_numerical_by_width": [item["r_memory_numerical"] for item in row],
        }

    balanced_orders: Dict[str, Any] = {}
    for name in DICTIONARIES:
        best: Dict[str, Any] = {}
        for block in ordered:
            entry = block["phase_c_balanced"][name]["orders"]
            smallest = None
            for order in sorted(int(key) for key in entry):
                row = entry[str(order)]
                if "status" in row:
                    continue
                worst = max(
                    row[f"{eta:+.3f}"]["balanced"]
                    for eta in PRIMARY_ETAS
                    if f"{eta:+.3f}" in row
                )
                if worst <= PASS_THRESHOLD:
                    smallest = order
                    break
            best[str(block["width"])] = smallest
        balanced_orders[name] = best

    return {
        "dictionary_profile": profile,
        "dictionary_growth": dictionary_growth,
        "balanced_transport_order": balanced_orders,
        "pass_threshold": PASS_THRESHOLD,
    }


def assemble(widths: Sequence[int], run_memory: bool) -> Dict[str, Any]:
    blocks = [width_closure(Generator(width), run_memory) for width in widths]
    return {
        "schema": SCHEMA,
        "issue": ISSUE,
        "phase": "B+C",
        "frozen_from": {
            "issue": 580,
            "rank": FROZEN_RANK,
            "note": "notes/p398-projected-memory-20260906.md",
        },
        "declared": {
            "widths": list(widths),
            "intervention": INTERVENTION,
            "etas": list(CLOSURE_ETAS),
            "dictionary_ranks": list(DICTIONARY_RANKS),
            "augmentation_ranks": list(AUGMENT_RANKS),
            "memory_orders": list(MEMORY_ORDERS),
            "balanced_orders": list(BALANCED_ORDERS),
            "pass_threshold": PASS_THRESHOLD,
            "closure_substeps": CLOSURE_SUBSTEPS,
            "era_blocks": ERA_BLOCKS,
        },
        "blocks": blocks,
        "decision": decide(blocks),
        "claim_boundary": (
            "Dictionary growth and repair comparison for one declared projection "
            "and source set on one exactly known finite process.  Not a "
            "percolation statement, and not transportable to square-site "
            "Matching One without a declared map between microscopic state spaces."
        ),
    }


def render(result: Mapping[str, Any]) -> str:
    lines: List[str] = ["#588 Phases B and C -- dictionary growth, and memory versus state", ""]
    for block in result["blocks"]:
        lines.append(f"width {block['width']:>2} ({block['states']:>4} states)")
        lines.append("  Phase B -- per dictionary, at the frozen rank 6")
        lines.append(
            "    dictionary                 seeds  r_linear  r_transport  rank(C)"
            "   M_mem  decay   r_mem(num/99.9)"
        )
        for entry in block["phase_b_dictionaries"]:
            memory = entry["memory"].get(str(FROZEN_RANK), {})
            effective = (memory.get("r_memory_effective") or {}).get("0.999")
            lines.append(
                f"    {entry['dictionary']:<26} {entry['seed_span_dimension']:>5}"
                f"  {_linear(entry['r_linear']):>8}"
                f"  {str(entry['r_transport']):>11}"
                f"  {str((memory.get('coupling_rank') or {}).get('rank')):>7}"
                f"  {_num(memory.get('M_memory'))}"
                f"  {_num(memory.get('decay_time'))}"
                f"   {memory.get('r_memory_numerical')}/{effective}"
            )
        lines.append("  Phase C1 -- more state (held-out balanced, declared intervention)")
        lines.append("    rank  +dof   eta=0    eta=+1/4   eta=-1/4   declared(0)")
        for entry in block["phase_c_state_augmentation"]:
            if "status" in entry:
                continue
            lines.append(
                f"    {entry['rank']:>4}  {entry['added_degrees_of_freedom']:>4}"
                f"  {entry['+0.000']['held_out']['balanced']:7.4f}"
                f"   {entry['+0.250']['held_out']['balanced']:8.4f}"
                f"   {entry['-0.250']['held_out']['balanced']:8.4f}"
                f"   {entry['+0.000']['declared']['balanced']:11.4f}"
            )
        memory_block = block.get("phase_c_memory_closure")
        if memory_block:
            lines.append("  Phase C2 -- more memory at rank 6 (eta = 0)")
            lines.append("    variant                        declared   held-out   kernel fit")
            zero = memory_block["+0.000"]
            for label in sorted(zero):
                fit = memory_block["baseline_fit_relative_error"].get(
                    label.replace("poles_", "").replace("_with_forcing", "")
                )
                lines.append(
                    f"    {label:<28}  {zero[label]['declared']['balanced']:8.4f}"
                    f"   {zero[label]['held_out']['balanced']:8.4f}"
                    f"   {_num(fit)}"
                )
        lines.append("  Phase C3 -- balanced adversary (frozen at eta=0, transported)")
        for name, family in block["phase_c_balanced"].items():
            spectrum = " ".join(f"{value:.3f}" for value in family["hankel_singular_values"][:6])
            lines.append(f"    {name:<26} hankel sigma: {spectrum}")
            row = []
            for order in sorted(family["orders"], key=int):
                entry = family["orders"][order]
                if "status" in entry:
                    row.append(f"r{order}:--")
                    continue
                worst = max(
                    entry[f"{eta:+.3f}"]["balanced"]
                    for eta in PRIMARY_ETAS
                    if f"{eta:+.3f}" in entry
                )
                row.append(
                    f"r{order}:{entry['+0.000']['balanced']:.4f}/{worst:.4f}"
                    f"(x{entry['long_horizon_growth']:.2f})"
                )
            lines.append("                               " + "  ".join(row))
        lines.append("")

    decision = result["decision"]
    lines.append("decision")
    for name, row in decision["dictionary_growth"].items():
        lines.append(f"  {name}")
        lines.append(f"    r_linear      {row['r_linear_by_width']}")
        lines.append(f"    r_transport   {row['r_transport_by_width']}")
        lines.append(f"    r_memory(num) {row['r_memory_numerical_by_width']}")
        lines.append(
            f"    balanced order {list(decision['balanced_transport_order'][name].values())}"
        )
    lines.append("")
    lines.append(result["claim_boundary"])
    return "\n".join(lines)


def _linear(payload: Mapping[str, Any]) -> str:
    value = payload.get("dimension")
    if payload.get("reached_the_whole_state_space"):
        return f"{value}="
    if payload.get("limited_by_budget"):
        return f">={value}"
    return str(value)


def _num(value: Optional[float]) -> str:
    if value is None:
        return "     --"
    return f"{value:7.4f}"


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--widths", type=int, nargs="+", default=list(DEFAULT_WIDTHS))
    parser.add_argument("--no-memory-closure", action="store_true")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args(list(argv) if argv is not None else None)
    result = assemble(args.widths, not args.no_memory_closure)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(rounded(result), indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    if args.json:
        print(json.dumps(rounded(result), indent=2, sort_keys=True))
    else:
        print(render(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
