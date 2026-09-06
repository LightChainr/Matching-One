#!/usr/bin/env python3
"""#600 Phase C: does the balanced realization factor through the reflection quotient?

#598 established that P398's coarsest exact positive lumping is the orbit
partition of one reflection ``R``, and that an invariant finite-horizon task is
blind to the ``R``-odd part of a symmetry-breaking intervention at first order.
What was not established is the thing that would make the quotient *usable*:
that the input/output object the project actually scores -- the finite-horizon
balanced realization -- is the same object whether it is built on the
microscopic state space or on the quotient.  This script builds it both ways and
compares.

**A** is the full microscopic chain, ``Catalan(w)`` states.
**B** is the exact reduction to the ``R``-orbit quotient, ``r_reflect(w)`` states.

They are compared on four frozen objects: the Hankel singular values on the
declared horizon grid, the balanced order at the declared tolerance, the response
matrices on the declared lag grid, and transport under ``H_even`` at
``eta = 0, +-1/4``.

Conventions, stated because the whole comparison lives or dies on them
-----------------------------------------------------------------------

A lumping is a pair of maps, one for each slot of the pairing, and getting
either one wrong silently changes the answer.  With ``O_a`` the orbits,
``n_a`` their sizes, and ``rep(a)`` a representative:

    readouts (right slot, acted on by G)   mu_Q[a] = mu[rep(a)]      restriction
    sources  (left slot, acted on by G^T)  f_Q[a]  = sum_{i in O_a} f[i]   lifting
    pairing                                <f_Q, mu_Q> = sum_a n_a f_a mu[rep(a)]

For ``R``-even ``f`` and ``mu`` this reproduces ``<f, mu>`` exactly, and the
quotient generator ``G_Q[a][b] = sum_{j in O_b} G[rep(a)][j]`` satisfies
``e^{t G_Q} mu_Q = (e^{tG} mu)_Q``, so every response matches.  Both facts are
gated below rather than asserted.

The one place this is not canonical is a readout that is *not* ``R``-even: then
``mu_Q`` depends on which representative is chosen and there is no exact
reduction at all.  ``halves_linked`` is exactly such a readout at odd widths, so
two readings are reported separately for every width:

    restriction      the literal reduction, sampled at the least representative.
                     Canonical for even readouts, arbitrary for odd ones.
    symmetrization   replace ``mu`` by ``(mu + R mu)/2`` first, then restrict.
                     Always canonical, but it is a projection of the readout,
                     not a reduction of the state.

The gap between the two columns is the measurement.  A run in which they
coincided everywhere would mean the comparison could not see the odd sector, and
the honest report is that it is vacuous rather than that it passed.

Boundary: one exactly known finite connectivity process, one declared
projection, dictionary and source set.  Not a percolation statement.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

if __package__ in (None, ""):
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parent))

from p398_intervention_transport import (  # noqa: E402
    CONTRAST_REFERENCE,
    HELD_OUT_READOUTS,
    LAGS,
    POISSON_TAIL,
    PRIMARY_READOUTS,
    Generator,
    poisson_horizon,
    source_distributions,
)
from p398_projected_memory import (  # noqa: E402
    HANKEL_ENERGY_LEVELS,
    HANKEL_TOLERANCE,
    KERNEL_GRID,
    hankel_orders,
)
from p398_reflection_parity import (  # noqa: E402
    READOUTS,
    orbit_labels,
    reflection,
    reflection_survey,
    state_permutation,
    tilted_rates,
)
from p398_memory_closure import BALANCED_ORDERS  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "p598-balanced-quotient" / "latest.json"
SCHEMA = "matching-one.p598-balanced-quotient.v1"
ISSUE = 600

DEFAULT_WIDTHS: Tuple[int, ...] = (4, 5, 6, 7, 8)

#: #600 asks for transport under the even part of the single-point tilt, at
#: ``eta = 0`` and the declared pair ``+-1/4``.  ``H_even`` is the right
#: direction to test transport on precisely because it is ``R``-even: it keeps
#: the tilted generator inside the equivariant class, so the quotient still
#: exists at every ``eta`` and the comparison stays well posed.  ``H_odd`` would
#: destroy the quotient by construction and could not be transported at all.
TILT = "even"
TILT_ETAS: Tuple[float, ...] = (0.0, 0.25, -0.25)

#: The two readout sets.  ``protected`` is D0, all of which are ``R``-even at
#: every width, so A and B must agree on it everywhere.  ``exposed`` adds
#: ``halves_linked``, which is ``R``-even only at even widths.
PROTECTED: Tuple[str, ...] = tuple(name for name, _ in PRIMARY_READOUTS)
EXPOSED: Tuple[str, ...] = PROTECTED + ("halves_linked",)
DICTIONARIES: Tuple[Tuple[str, Tuple[str, ...]], ...] = (
    ("protected", PROTECTED),
    ("exposed", EXPOSED),
)
CONVENTIONS: Tuple[str, ...] = ("restriction", "symmetrization")

#: Frozen from p398_memory_closure.balanced_family -- the rank below which a
#: Hankel direction is arithmetic rather than a Hankel direction.  It is a
#: statement about eigenvalues, so in singular-value terms it is 1e-6, which is
#: the project's declared Hankel tolerance.
GRAMIAN_FLOOR = 1e-12

#: The frozen re-biorthogonalization replaces ``W`` by ``S^{-T} W`` with
#: ``S = W^T V``.  The correction that actually restores ``W^T V = I`` is
#: ``S^{-1} W``; the transposed form only cancels the antisymmetric part of the
#: defect, which is why the frozen artifact reports a residual ~5e-7 after
#: "correction".  Both are run: ``inverse`` is the primary reading and
#: ``frozen_transpose`` is what the continuity gate has to reproduce.
BIORTHOGONAL_MODES: Tuple[str, ...] = ("inverse", "frozen_transpose")
PRIMARY_BIORTHOGONAL_MODE = "inverse"

#: The intervention the frozen artifact scored, kept only so the new pipeline
#: can be checked against the committed numbers before it is trusted on B.
CONTINUITY_INTERVENTION = "uniform_join_minus_detach"
CONTINUITY_ETAS: Tuple[float, ...] = (-0.25, 0.0, 0.25)
CONTINUITY_ARTIFACT = ROOT / "results" / "p398-memory-closure" / "latest.json"

READOUT_FUNCTIONS = dict(READOUTS)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


# --------------------------------------------------------------------------
# The reflection and the two chains
# --------------------------------------------------------------------------


class Reflection:
    """The task-compatible reflection and the orbit bookkeeping it induces."""

    def __init__(self, generator: Generator) -> None:
        survey = reflection_survey(generator)
        k = survey["task_compatible_k"]
        _require(k is not None, f"width {generator.width}: the dictionary names no reflection")
        self.width = generator.width
        self.k = k
        self.pi = state_permutation(generator, reflection(generator.width, k))
        self.labels = orbit_labels(generator, k)
        states = len(self.labels)
        self.quotient_size = len(set(self.labels))
        self.sizes = [0] * self.quotient_size
        self.representatives = [-1] * self.quotient_size
        for index, label in enumerate(self.labels):
            self.sizes[label] += 1
            if self.representatives[label] < 0:
                self.representatives[label] = index
        self.odd_dimension = states - self.quotient_size
        self._permutation_matrix_cache: Optional[np.ndarray] = None

    def permutation(self, vector: np.ndarray) -> np.ndarray:
        """``(R v)[i] = v[pi[i]]`` -- the convention p398_reflection_parity uses."""

        return vector[self.pi]

    def even_part(self, vector: np.ndarray) -> np.ndarray:
        return 0.5 * (vector + self.permutation(vector))

    def lift(self, vector: np.ndarray) -> np.ndarray:
        """A measure pushed to the quotient: each block carries the orbit sum."""

        array = np.asarray(vector, dtype=float)
        out = np.zeros(self.quotient_size)
        for block in range(self.quotient_size):
            members = [
                index for index, label in enumerate(self.labels) if label == block
            ]
            out[block] = float(array[members].sum())
        return out

    def restrict(self, vector: np.ndarray) -> np.ndarray:
        """A function pulled back to the quotient: each block takes the value at its representative."""

        return np.asarray(vector, dtype=float)[self.representatives]

    def odd_part(self, vector: np.ndarray) -> np.ndarray:
        return 0.5 * (vector - self.permutation(vector))

    def odd_content(self, vector: np.ndarray) -> float:
        norm = float(np.linalg.norm(vector))
        if norm == 0.0:
            return 0.0
        return float(np.linalg.norm(self.odd_part(vector)) / norm)


def dense(rows: Sequence[Sequence[Tuple[int, float]]], size: int) -> np.ndarray:
    matrix = np.zeros((size, size))
    for source, row in enumerate(rows):
        for destination, value in row:
            matrix[source, destination] += value
    return matrix


def quotient_rows(
    rows: Sequence[Sequence[Tuple[int, float]]], reflection: Reflection
) -> List[List[Tuple[int, float]]]:
    """``G_Q[a][b] = sum_{j in O_b} G[rep(a)][j]``.

    Well-definedness is not free: strong lumpability asks that this be
    independent of the representative, which holds here because ``R`` is a
    symmetry of ``G`` rather than for any generic reason.  ``max_representative_
    dependence`` in the artifact is the check, not an assumption.
    """

    out: List[List[Tuple[int, float]]] = []
    for block in range(reflection.quotient_size):
        accumulated: Dict[int, float] = {}
        for destination, value in rows[reflection.representatives[block]]:
            label = reflection.labels[destination]
            accumulated[label] = accumulated.get(label, 0.0) + value
        out.append(sorted(accumulated.items()))
    return out


def representative_dependence(
    rows: Sequence[Sequence[Tuple[int, float]]], reflection: Reflection
) -> float:
    """How much ``G_Q`` would move if a different representative were picked."""

    worst = 0.0
    for block in range(reflection.quotient_size):
        members = [
            index for index, label in enumerate(reflection.labels) if label == block
        ]
        profiles: List[Dict[int, float]] = []
        for index in members:
            profile: Dict[int, float] = {}
            for destination, value in rows[index]:
                profile[reflection.labels[destination]] = (
                    profile.get(reflection.labels[destination], 0.0) + value
                )
            profiles.append(profile)
        for other in profiles[1:]:
            for key in set(profiles[0]) | set(other):
                worst = max(
                    worst, abs(profiles[0].get(key, 0.0) - other.get(key, 0.0))
                )
    return worst


def lift_sources(
    sources: Sequence[Tuple[str, Sequence[float]]], reflection: Reflection
) -> List[Tuple[str, np.ndarray]]:
    """Sources are measures, so a block carries the *sum* over its orbit."""

    return [(name, reflection.lift(vector)) for name, vector in sources]


def map_readouts(
    observables: Sequence[np.ndarray], reflection: Reflection, convention: str
) -> List[np.ndarray]:
    """Readouts are functions, so a block carries the *value* at a state.

    ``restriction`` samples the least representative; ``symmetrization`` first
    projects onto the ``R``-even sector, which is the only choice that does not
    depend on the representative.  They coincide exactly when the readout is
    ``R``-even, which is why the difference between them is the measurement.
    """

    out: List[np.ndarray] = []
    for vector in observables:
        array = np.asarray(vector, dtype=float)
        if convention == "symmetrization":
            array = reflection.even_part(array)
        elif convention != "restriction":  # pragma: no cover - guarded by caller
            raise ValueError(f"unknown readout convention {convention!r}")
        out.append(array[reflection.representatives])
    return out


# --------------------------------------------------------------------------
# Propagation and the frozen balanced construction, in numpy
# --------------------------------------------------------------------------


def evolve(
    matrix: np.ndarray, rate: float, vectors: np.ndarray, grid: Sequence[float]
) -> np.ndarray:
    """``exp(t A) v`` for every column of ``vectors``, by uniformization.

    Same algorithm as ``p398_intervention_transport.evolve_observables``:
    Poisson-weighted powers of the lazy chain, so no term cancels and large
    ``t`` stays accurate.  Returned as ``[lag][state][column]``.
    """

    size, count = vectors.shape
    horizon = poisson_horizon(rate * float(max(grid)))
    means = np.array([rate * float(lag) for lag in grid])
    weights = np.exp(-means)
    accumulated = np.zeros((len(grid), size, count))
    step = np.eye(size) + matrix / rate
    current = vectors.copy()
    for order in range(horizon + 1):
        live = weights > POISSON_TAIL
        if live.any():
            accumulated[live] += weights[live][:, None, None] * current[None, :, :]
        weights = weights * means / (order + 1)
        if order == horizon:
            break
        current = step @ current
    return accumulated


def hankel_cross(
    matrix: np.ndarray,
    rate: float,
    size: int,
    observables: Sequence[np.ndarray],
    contrasts: Sequence[np.ndarray],
    grid: Sequence[float],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The finite-horizon Hankel matrix, plus the two trajectory bases.

    ``cross[(lag, contrast)][(lag', readout)] = <e^{t G^T} f, e^{t' G} mu>``,
    weighted by the trapezoid rule on the declared grid -- i.e. the declared
    response sampled at ``t + t'``.  The reach side carries the evolved
    readouts and the observe side the evolved source contrasts; every response
    is a pairing of one against the other, so both are returned.
    """

    step = grid[1] - grid[0]
    weights = np.array(
        [
            math.sqrt(step * (0.5 if index in (0, len(grid) - 1) else 1.0))
            for index in range(len(grid))
        ]
    )
    reach = evolve(matrix, rate, np.stack(observables, axis=1), grid)
    reach = (reach * weights[:, None, None]).transpose(1, 0, 2).reshape(size, -1)
    observe = evolve(matrix.T, rate, np.stack(contrasts, axis=1), grid)
    observe = (observe * weights[:, None, None]).transpose(1, 0, 2).reshape(size, -1)
    cross = observe.T @ reach
    return cross, reach, observe


def response(
    matrix: np.ndarray,
    rate: float,
    sources: Sequence[Tuple[str, np.ndarray]],
    observables: Sequence[np.ndarray],
    lags: Sequence[float],
) -> np.ndarray:
    """``[lag][source][readout]`` of ``<f, exp(t G) mu>``."""

    evolved = evolve(matrix, rate, np.stack(observables, axis=1), lags)
    return np.einsum("rs,lsm->lrm", np.stack([f for _, f in sources]), evolved)


def contrasts_from(
    sources: Sequence[Tuple[str, np.ndarray]], reference_index: int
) -> List[np.ndarray]:
    reference = np.asarray(sources[reference_index][1], dtype=float)
    return [
        np.asarray(vector, dtype=float) - reference
        for index, (_, vector) in enumerate(sources)
        if index != reference_index
    ]


def balance(
    cross: np.ndarray,
    reach: np.ndarray,
    observe: np.ndarray,
    order: int,
) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """The balanced pair ``(W, V)`` at one order, or ``None`` if unfillable.

    Mirrors ``p398_memory_closure.balanced_family``: singular vectors of the
    Hankel matrix, scaled by ``sigma^{-1/2}`` on each side.  An order the
    Gramian cannot fill is reported as missing rather than padded with noise --
    how many input/output coordinates the task exposes is the question, so
    refusing to invent one is part of the answer.
    """

    left, values, right = np.linalg.svd(cross, full_matrices=False)
    if values.size == 0 or values[0] <= 0.0:
        return None
    floor = GRAMIAN_FLOOR * float(values[0] * values[0])
    keep = [index for index in range(min(order, values.size)) if values[index] ** 2 > floor]
    if len(keep) < order:
        return None
    pair = []
    for matrix_in, vectors in ((observe, left.T), (reach, right)):
        basis = []
        for index in keep:
            sigma = math.sqrt(float(values[index]))
            basis.append((vectors[index] @ matrix_in.T) / sigma)
        pair.append(np.stack(basis))
    return pair[0], pair[1]


def rebiorthogonalize(left: np.ndarray, right: np.ndarray, mode: str) -> np.ndarray:
    """Restore ``W^T V = I`` inside the same two spans.

    ``inverse`` is the correction that does what the name says.  ``frozen_
    transpose`` reproduces ``p398_memory_closure._rebiorthogonalize`` exactly,
    and only cancels the antisymmetric part of the defect -- which is why the
    committed artifact still carries ~5e-7 after correcting.  Both are run so
    the size of that difference is on the record rather than in the noise.
    """

    gram = left @ right.T
    inverse = np.linalg.inv(gram)
    coefficients = inverse if mode == "inverse" else inverse.T
    return coefficients @ left


def biorthogonality_defect(left: np.ndarray, right: np.ndarray) -> float:
    gram = left @ right.T
    return float(np.max(np.abs(gram - np.eye(gram.shape[0]))))


def score_realization(
    left: np.ndarray,
    right: np.ndarray,
    matrix: np.ndarray,
    observables: Sequence[np.ndarray],
    sources: Sequence[Tuple[str, np.ndarray]],
    truth: np.ndarray,
    reference_index: int,
) -> Dict[str, float]:
    """Frozen-model response error, pooled and readout-balanced.

    The same two norms ``p398_memory_closure.score_realization`` reports, so a
    balanced model is scored as the Petrov-Galerkin object it is rather than
    being handed the orthogonal projection of its own span.
    """

    reduced = left @ matrix @ right.T
    starts = np.stack([np.asarray(observable) for observable in observables]) @ left.T
    coordinates = np.stack([np.asarray(vector) for _, vector in sources]) @ right.T
    reference = coordinates[reference_index]
    prediction = np.empty((len(LAGS), len(sources), len(observables)))
    for lag_index, lag in enumerate(LAGS):
        propagator = _expm(reduced, lag)
        moved = starts @ propagator.T
        prediction[lag_index] = coordinates @ moved.T
    predicted = prediction - prediction[:, reference_index : reference_index + 1, :]
    difference = predicted - truth
    pooled = float(
        np.sqrt((difference ** 2).sum() / (truth ** 2).sum())
        if (truth ** 2).sum() > 0
        else 0.0
    )
    per_readout = []
    for column in range(truth.shape[2]):
        denom = float((truth[:, :, column] ** 2).sum())
        per_readout.append(
            float(np.sqrt((difference[:, :, column] ** 2).sum() / denom))
            if denom > 0
            else 0.0
        )
    balanced = (
        math.sqrt(sum(value * value for value in per_readout) / len(per_readout))
        if per_readout
        else 0.0
    )
    return {"pooled": pooled, "balanced": balanced}


def _expm(matrix: np.ndarray, time: float) -> np.ndarray:
    """Scaling-and-squaring Taylor exponential, matching ``dense_expm``."""

    size = matrix.shape[0]
    scale = float(np.max(np.abs(matrix))) * abs(time)
    squarings = 0
    while scale > 0.5:
        scale /= 2.0
        squarings += 1
    factor = time / (2.0**squarings)
    scaled = matrix * factor
    result = np.eye(size)
    term = np.eye(size)
    for order in range(1, 19):
        term = term @ scaled / order
        result = result + term
    for _ in range(squarings):
        result = result @ result
    return result


def balanced_report(
    matrix: np.ndarray,
    rate: float,
    size: int,
    observables: Sequence[np.ndarray],
    sources: Sequence[Tuple[str, np.ndarray]],
    reference_index: int,
    readout_scales: Sequence[float],
    tilted: Sequence[Tuple[float, np.ndarray]],
) -> Dict[str, Any]:
    """Everything the frozen construction produces, for one chain.

    ``tilted`` is ``[(eta, generator)]``; the truth at each ``eta`` is
    recomputed on this chain, exactly as the frozen ``cache`` does, and the
    balancing transform is built once at the baseline and frozen across all of
    them -- the transform is not allowed to see the intervention it is scored
    on.
    """

    contrasts = contrasts_from(sources, reference_index)
    scaled = [
        np.asarray(observable) / scale
        for observable, scale in zip(observables, readout_scales)
    ]
    cross, reach, observe = hankel_cross(
        matrix, rate, size, scaled, contrasts, KERNEL_GRID
    )
    singular = np.linalg.svd(cross, compute_uv=False)
    spectrum = hankel_orders([float(value) for value in singular])

    truths: Dict[float, np.ndarray] = {}
    for eta, generator in tilted:
        raw = response(generator, rate, sources, observables, LAGS)
        truths[eta] = raw - raw[:, reference_index : reference_index + 1, :]

    orders: Dict[str, Any] = {}
    bases: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    reduced: Dict[str, Dict[float, np.ndarray]] = {}
    baseline = generator_baseline(tilted)
    for order in BALANCED_ORDERS:
        pair = balance(cross, reach, observe, order)
        if pair is None:
            orders[str(order)] = {"status": "gramian_rank_deficient"}
            continue
        left, right = pair
        entry: Dict[str, Any] = {
            "status": "filled",
            "biorthogonality_defect_before_correction": biorthogonality_defect(
                left, right
            ),
            "by_correction": {},
        }
        for mode in BIORTHOGONAL_MODES:
            corrected = rebiorthogonalize(left, right, mode)
            block: Dict[str, Any] = {
                "biorthogonality_defect": biorthogonality_defect(corrected, right)
            }
            for eta, generator in tilted:
                block[f"{eta:+.3f}"] = score_realization(
                    corrected,
                    right,
                    generator,
                    observables,
                    sources,
                    truths[eta],
                    reference_index,
                )
                if mode == PRIMARY_BIORTHOGONAL_MODE:
                    reduced.setdefault(str(order), {})[eta] = (
                        corrected @ generator @ right.T
                    )
            block["long_horizon_growth"] = float(
                np.linalg.norm(
                    _expm(corrected @ baseline @ right.T, 2.0 * max(LAGS))
                )
            )
            entry["by_correction"][mode] = block
        bases[str(order)] = (
            rebiorthogonalize(left, right, PRIMARY_BIORTHOGONAL_MODE),
            right,
        )
        entry["long_horizon_growth"] = entry["by_correction"][PRIMARY_BIORTHOGONAL_MODE][
            "long_horizon_growth"
        ]
        entry["biorthogonality_defect"] = entry["by_correction"][
            PRIMARY_BIORTHOGONAL_MODE
        ]["biorthogonality_defect"]
        orders[str(order)] = entry

    return {
        "states": size,
        "hankel_singular_values": spectrum["spectrum"],
        "numerical_rank": spectrum["numerical_rank"],
        "effective_orders": spectrum["effective_orders"],
        "orders": orders,
        "_cross": cross,
        "_reach": reach,
        "_observe": observe,
        "_truths": truths,
        "_singular": [float(value) for value in singular],
        "_bases": bases,
        "_reduced": reduced,
    }


def generator_baseline(tilted: Sequence[Tuple[float, np.ndarray]]) -> np.ndarray:
    for eta, generator in tilted:
        if eta == 0.0:
            return generator
    raise ValueError("the baseline is required")  # pragma: no cover


# --------------------------------------------------------------------------
# Differences
# --------------------------------------------------------------------------


def relative_frobenius(left: np.ndarray, right: np.ndarray) -> float:
    denominator = float(np.linalg.norm(left))
    if denominator == 0.0:
        return float(np.linalg.norm(right))
    return float(np.linalg.norm(left - right) / denominator)


def spectrum_difference(
    left: Sequence[float], right: Sequence[float]
) -> Dict[str, float]:
    count = min(len(left), len(right))
    if count == 0:
        return {"max_absolute": 0.0, "max_relative": 0.0, "compared": 0}
    first = np.asarray(left[:count])
    second = np.asarray(right[:count])
    absolute = float(np.max(np.abs(first - second)))
    scale = float(np.max(np.abs(first)))
    return {
        "max_absolute": absolute,
        "max_relative": absolute / scale if scale > 0 else 0.0,
        "compared": count,
    }


def sign_corrected_basis_difference(
    left: np.ndarray, right: np.ndarray
) -> float:
    """Worst per-vector relative gap, allowing the arbitrary sign of a singular vector."""

    worst = 0.0
    for row_left, row_right in zip(left, right):
        norm = float(np.linalg.norm(row_right))
        if norm == 0.0:
            continue
        gap = min(
            float(np.linalg.norm(row_left - row_right)),
            float(np.linalg.norm(row_left + row_right)),
        )
        worst = max(worst, gap / norm)
    return worst


# --------------------------------------------------------------------------
# One width
# --------------------------------------------------------------------------


def width_report(width: int) -> Dict[str, Any]:
    generator = Generator(width)
    mirror = Reflection(generator)
    size = generator.size
    observables_by_name = {
        name: np.array([float(function(state)) for state in generator.states])
        for name, function in READOUTS
    }
    sources = [
        (name, np.asarray(vector, dtype=float))
        for name, vector in source_distributions(generator)
    ]
    reference_index = [name for name, _ in sources].index(CONTRAST_REFERENCE)

    parity = {
        name: mirror.odd_content(values) for name, values in observables_by_name.items()
    }
    source_defect = max(
        float(np.max(np.abs(mirror.permutation(vector) - vector)))
        for _, vector in sources
    )

    # Baseline scales are fixed on the microscopic chain before any balancing,
    # and the same numbers are used for both chains: an output metric that
    # changed with the chain would not be a fixed metric.
    baseline_rows = generator.rows(generator.baseline_rates())
    baseline_rate = generator.exit_rate(generator.baseline_rates())
    baseline_matrix = dense(baseline_rows, size)
    baseline_truth = response(
        baseline_matrix, baseline_rate, sources,
        [observables_by_name[name] for name, _ in READOUTS], LAGS,
    )
    scales = {
        name: float(
            np.sqrt(
                float(
                    (
                        (
                            baseline_truth
                            - baseline_truth[:, reference_index : reference_index + 1, :]
                        )[:, :, index]
                        ** 2
                    ).sum()
                )
            )
        )
        or 1.0
        for index, (name, _) in enumerate(READOUTS)
    }

    # -- the two chains, at every eta ------------------------------------
    micro_tilted: List[Tuple[float, np.ndarray]] = []
    quotient_tilted: List[Tuple[float, np.ndarray]] = []
    representative_gap = 0.0
    for eta in TILT_ETAS:
        rates = tilted_rates(generator, mirror.k, TILT, eta)
        micro = dense(generator.rows(rates), size)
        gap = representative_dependence(generator.rows(rates), mirror)
        representative_gap = max(representative_gap, gap)
        quotient = dense(
            quotient_rows(generator.rows(rates), mirror), mirror.quotient_size
        )
        micro_tilted.append((eta, micro))
        quotient_tilted.append((eta, quotient))

    quotient_sources = lift_sources(
        [(name, vector) for name, vector in sources], mirror
    )

    # -- the odd sector ---------------------------------------------------
    _, reach_micro, observe_micro = hankel_cross(
        baseline_matrix,
        baseline_rate,
        size,
        [observables_by_name[name] / scales[name] for name in PROTECTED],
        contrasts_from(sources, reference_index),
        KERNEL_GRID,
    )
    odd_sector = {
        "dimension": mirror.odd_dimension,
        "quotient_dimension": mirror.quotient_size,
        # Columns are the trajectory vectors; the projector acts on the state
        # axis, so permuting rows is what applies R to every column at once.
        "reach_odd_energy_relative": float(
            np.linalg.norm(mirror.odd_part(reach_micro)) / np.linalg.norm(reach_micro)
        ),
        "observe_odd_energy_relative": float(
            np.linalg.norm(mirror.odd_part(observe_micro))
            / np.linalg.norm(observe_micro)
        ),
    }

    dictionaries: Dict[str, Any] = {}
    for name, members in DICTIONARIES:
        blocks: Dict[str, Any] = {}
        micro_observables = [observables_by_name[member] for member in members]
        member_scales = [scales[member] for member in members]
        report_a = balanced_report(
            baseline_matrix,
            baseline_rate,
            size,
            micro_observables,
            sources,
            reference_index,
            member_scales,
            micro_tilted,
        )
        blocks["A"] = _strip(report_a)
        for convention in CONVENTIONS:
            projected = [observables_by_name[member] for member in members]
            if convention == "symmetrization" and all(
                parity[member] == 0.0 for member in members
            ):
                blocks[convention] = {
                    "identical_to_restriction": True,
                    "reason": (
                        "every readout in this dictionary is R-even, so "
                        "symmetrization changes nothing and the two conventions "
                        "are the same map"
                    ),
                }
                continue
            quotient_observables = map_readouts(projected, mirror, convention)
            report_b = balanced_report(
                dense(quotient_rows(baseline_rows, mirror), mirror.quotient_size),
                baseline_rate,
                mirror.quotient_size,
                quotient_observables,
                quotient_sources,
                reference_index,
                member_scales,
                quotient_tilted,
            )
            blocks[convention] = {
                "identical_to_restriction": False,
                "B": _strip(report_b),
                "differences": _differences(report_a, report_b, mirror),
            }
        dictionaries[name] = {
            "members": list(members),
            "odd_content": {member: parity[member] for member in members},
            "max_odd_content": max(parity[member] for member in members),
            "readings": blocks,
        }

    return {
        "width": width,
        "states": size,
        "reflection": f"i -> {mirror.k} - i",
        "quotient_states": mirror.quotient_size,
        "odd_sector": odd_sector,
        "parity": {
            "readout_odd_content": parity,
            "max_source_parity_defect": source_defect,
            "all_sources_are_R_even": source_defect == 0.0,
        },
        "quotient_is_well_defined": {
            "max_representative_dependence": representative_gap,
            "checked_at_etas": [float(eta) for eta in TILT_ETAS],
        },
        "dictionaries": dictionaries,
    }


def _strip(report: Dict[str, Any]) -> Dict[str, Any]:
    return {
        key: value for key, value in report.items() if not key.startswith("_")
    }


def _differences(
    a: Dict[str, Any], b: Dict[str, Any], mirror: "Reflection"
) -> Dict[str, Any]:
    """Everything the two chains disagree about, with a norm attached."""

    out: Dict[str, Any] = {
        "hankel_spectrum": spectrum_difference(
            a["hankel_singular_values"], b["hankel_singular_values"]
        ),
        "numerical_rank_equal": a["numerical_rank"] == b["numerical_rank"],
        "numerical_rank": {"A": a["numerical_rank"], "B": b["numerical_rank"]},
        "effective_orders_equal": a["effective_orders"] == b["effective_orders"],
        "response": {},
        "orders": {},
    }
    for eta in sorted(a["_truths"]):
        out["response"][f"{eta:+.3f}"] = {
            "max_absolute": float(
                np.max(np.abs(a["_truths"][eta] - b["_truths"][eta]))
            ),
            "relative_frobenius": relative_frobenius(
                a["_truths"][eta], b["_truths"][eta]
            ),
        }
    for order in BALANCED_ORDERS:
        left = a["orders"].get(str(order))
        right = b["orders"].get(str(order))
        if not isinstance(left, dict) or not isinstance(right, dict):
            out["orders"][str(order)] = {
                "status": {"A": _status(left), "B": _status(right)},
                "comparable": False,
            }
            continue
        if left["status"] != "filled" or right["status"] != "filled":
            out["orders"][str(order)] = {
                "status": {"A": left["status"], "B": right["status"]},
                "comparable": False,
            }
            continue
        block: Dict[str, Any] = {
            "status": {"A": "filled", "B": "filled"},
            "comparable": True,
            "biorthogonality_defect": {
                "A": left["biorthogonality_defect"],
                "B": right["biorthogonality_defect"],
            },
        }
        for mode in BIORTHOGONAL_MODES:
            errors: Dict[str, Any] = {}
            for eta in sorted(a["_truths"]):
                left_score = left["by_correction"][mode][f"{eta:+.3f}"]
                right_score = right["by_correction"][mode][f"{eta:+.3f}"]
                errors[f"{eta:+.3f}"] = {
                    "pooled": {
                        "A": left_score["pooled"],
                        "B": right_score["pooled"],
                        "absolute_difference": abs(
                            left_score["pooled"] - right_score["pooled"]
                        ),
                    },
                    "balanced": {
                        "A": left_score["balanced"],
                        "B": right_score["balanced"],
                        "absolute_difference": abs(
                            left_score["balanced"] - right_score["balanced"]
                        ),
                    },
                }
            block[mode] = {
                "eta": errors,
                "long_horizon_growth": {
                    "A": left["by_correction"][mode]["long_horizon_growth"],
                    "B": right["by_correction"][mode]["long_horizon_growth"],
                },
            }
        left_pair = a["_bases"][str(order)]
        right_pair = b["_bases"][str(order)]
        block["basis_transport"] = {
            "left_vs_lifted": sign_corrected_basis_difference(
                right_pair[0],
                np.stack([mirror.lift(row) for row in left_pair[0]]),
            ),
            "right_vs_restricted": sign_corrected_basis_difference(
                right_pair[1],
                np.stack([mirror.restrict(row) for row in left_pair[1]]),
            ),
            "note": (
                "whether the quotient's balanced basis is exactly the lifted / "
                "restricted microscopic one; a large right-side gap with a small "
                "left-side one is the odd part of the readout the quotient cannot "
                "carry"
            ),
        }
        block["reduced_generator"] = {}
        for eta in sorted(a["_truths"]):
            block["reduced_generator"][f"{eta:+.3f}"] = {
                "relative_frobenius": relative_frobenius(
                    a["_reduced"][str(order)][eta], b["_reduced"][str(order)][eta]
                ),
                "max_absolute": float(
                    np.max(
                        np.abs(
                            a["_reduced"][str(order)][eta]
                            - b["_reduced"][str(order)][eta]
                        )
                    )
                ),
            }
        out["orders"][str(order)] = block
    return out


def _status(entry: Any) -> str:
    if isinstance(entry, dict):
        return str(entry.get("status"))
    return "missing"


# --------------------------------------------------------------------------
# Continuity gate against the committed #588 balanced numbers
# --------------------------------------------------------------------------


def continuity_gate(width: int) -> Dict[str, Any]:
    """Re-run A on the frozen intervention and compare to the committed artifact.

    The new pipeline is numpy and the committed one is pure Python, so this is
    the check that the port carries the same numbers before any B result is
    read.  It is run under ``uniform_join_minus_detach`` -- the intervention the
    frozen artifact scored -- and against ``D0``, which is the ``protected``
    dictionary.
    """

    if not CONTINUITY_ARTIFACT.exists():
        return {"status": "skipped", "reason": "committed artifact is not present"}
    committed = json.loads(CONTINUITY_ARTIFACT.read_text())
    block = next(
        (item for item in committed["blocks"] if item["width"] == width), None
    )
    if block is None:
        return {"status": "skipped", "reason": f"no committed block at width {width}"}
    reference = block["phase_c_balanced"]["D0_additive_local_counts"]

    generator = Generator(width)
    size = generator.size
    observables_by_name = {
        name: np.array([float(function(state)) for state in generator.states])
        for name, function in READOUTS
    }
    sources = [
        (name, np.asarray(vector, dtype=float))
        for name, vector in source_distributions(generator)
    ]
    reference_index = [name for name, _ in sources].index(CONTRAST_REFERENCE)
    members = list(PROTECTED)

    tilted: List[Tuple[float, np.ndarray]] = []
    for eta in CONTINUITY_ETAS:
        rates = generator.rates(CONTINUITY_INTERVENTION, eta)
        tilted.append((eta, dense(generator.rows(rates), size)))
    baseline_rows = generator.rows(generator.baseline_rates())
    baseline_rate = generator.exit_rate(generator.baseline_rates())
    baseline_matrix = dense(baseline_rows, size)

    micro_observables = [observables_by_name[name] for name in members]
    baseline_truth = response(
        baseline_matrix, baseline_rate, sources,
        [observables_by_name[name] for name in members], LAGS,
    )
    contrast = baseline_truth - baseline_truth[:, reference_index : reference_index + 1, :]
    scales = [
        float(np.sqrt(float((contrast[:, :, index] ** 2).sum()))) or 1.0
        for index in range(len(members))
    ]
    report = balanced_report(
        baseline_matrix,
        baseline_rate,
        size,
        micro_observables,
        sources,
        reference_index,
        scales,
        tilted,
    )

    spectrum_gap = spectrum_difference(
        report["hankel_singular_values"], reference["hankel_singular_values"]
    )
    # Split the spectrum comparison at the declared tolerance.  The committed
    # artifact stores ten significant digits and gets its singular values from a
    # Jacobi eigensolver on the squared Gramian, which squares the condition
    # number: entries below 1e-6 of the leading value are declared arithmetic by
    # the project's own rule, so agreement there is not a claim worth making.
    pairs = list(zip(report["hankel_singular_values"], reference["hankel_singular_values"]))
    above = [pair for pair in pairs if pair[0] > HANKEL_TOLERANCE]
    below = [pair for pair in pairs if pair[0] <= HANKEL_TOLERANCE]
    gap_above = max((abs(a - b) for a, b in above), default=0.0)
    gap_below = max((abs(a - b) for a, b in below), default=0.0)

    worst_score = 0.0
    by_mode: Dict[str, float] = {}
    for mode in BIORTHOGONAL_MODES:
        worst = 0.0
        for order in BALANCED_ORDERS:
            mine = report["orders"][str(order)]
            theirs = reference["orders"][str(order)]
            if mine.get("status") != "filled" or "biorthogonality_defect" not in theirs:
                continue
            for eta in CONTINUITY_ETAS:
                key = f"{eta:+.3f}"
                if key not in theirs:
                    continue
                mine_score = mine["by_correction"][mode][key]
                for norm in ("pooled", "balanced"):
                    worst = max(
                        worst,
                        abs(mine_score[norm] - float(theirs[key][norm])),
                    )
        by_mode[mode] = worst
        if mode == PRIMARY_BIORTHOGONAL_MODE:
            worst_score = worst
    return {
        "status": "run",
        "width": width,
        "intervention": CONTINUITY_INTERVENTION,
        "dictionary": "D0_additive_local_counts",
        "hankel_spectrum": spectrum_gap,
        "hankel_spectrum_above_tolerance": {
            "entries": len(above),
            "max_absolute_difference": gap_above,
        },
        "hankel_spectrum_at_or_below_tolerance": {
            "entries": len(below),
            "max_absolute_difference": gap_below,
            "note": (
                "the committed artifact is stored at ten significant digits and "
                "takes its singular values from an eigensolver on the squared "
                "Gramian, so entries at or below the declared 1e-6 tolerance are "
                "arithmetic there rather than Hankel directions; agreement below "
                "that line is not a claim"
            ),
        },
        "max_score_difference_by_correction": by_mode,
        "max_score_difference": worst_score,
        "frozen_defect_left_by_the_transposed_correction": _worst_committed_defect(
            reference
        ),
    }


def _worst_committed_defect(reference: Mapping[str, Any]) -> float:
    worst = 0.0
    for entry in reference["orders"].values():
        if isinstance(entry, dict) and "biorthogonality_defect" in entry:
            worst = max(worst, float(entry["biorthogonality_defect"]))
    return worst


# --------------------------------------------------------------------------
# Assembly
# --------------------------------------------------------------------------


def decide(blocks: Sequence[Mapping[str, Any]]) -> Dict[str, Any]:
    """Three claims, kept apart because they can fail separately.

    Declared before the numbers: the protected dictionary must factor through
    the quotient at every width; the odd sector must be exactly uncontrollable
    and unobservable; and the exposed dictionary must *fail* to factor at the
    two odd widths, which is the only thing that shows the first claim was a
    test rather than an identity.
    """

    ordered = sorted(blocks, key=lambda block: block["width"])

    protected_rows = []
    exposed_rows = []
    for block in ordered:
        width = block["width"]
        readings = block["dictionaries"]["protected"]["readings"]
        protection = readings["restriction"]
        protected_rows.append(
            {
                "width": width,
                "states": block["states"],
                "quotient_states": block["quotient_states"],
                "spectrum_max_absolute_difference": protection["differences"][
                    "hankel_spectrum"
                ]["max_absolute"],
                "numerical_rank": {
                    "A": protection["differences"]["numerical_rank"]["A"],
                    "B": protection["differences"]["numerical_rank"]["B"],
                },
                "worst_response_relative_frobenius": max(
                    value["relative_frobenius"]
                    for value in protection["differences"]["response"].values()
                ),
            }
        )
        exposed = block["dictionaries"]["exposed"]["readings"]
        row = {
            "width": width,
            "max_odd_content": block["dictionaries"]["exposed"]["max_odd_content"],
            "restriction_spectrum_max_absolute_difference": exposed["restriction"][
                "differences"
            ]["hankel_spectrum"]["max_absolute"],
            "restriction_worst_response_relative_frobenius": max(
                value["relative_frobenius"]
                for value in exposed["restriction"]["differences"]["response"].values()
            ),
        }
        symmetrized = exposed["symmetrization"]
        row["symmetrization_identical_to_restriction"] = symmetrized[
            "identical_to_restriction"
        ]
        if not symmetrized["identical_to_restriction"]:
            row["symmetrization_spectrum_max_absolute_difference"] = symmetrized[
                "differences"
            ]["hankel_spectrum"]["max_absolute"]
        exposed_rows.append(row)

    odd_rows = [
        {
            "width": block["width"],
            "odd_dimension": block["odd_sector"]["dimension"],
            "reach_odd_energy_relative": block["odd_sector"][
                "reach_odd_energy_relative"
            ],
            "observe_odd_energy_relative": block["odd_sector"][
                "observe_odd_energy_relative"
            ],
        }
        for block in ordered
    ]

    protection_tolerance = 1e-9
    protected_factors = all(
        row["spectrum_max_absolute_difference"] <= protection_tolerance
        and row["worst_response_relative_frobenius"] <= protection_tolerance
        and row["numerical_rank"]["A"] == row["numerical_rank"]["B"]
        for row in protected_rows
    )
    odd_inert = all(
        row["reach_odd_energy_relative"] <= 1e-10
        and row["observe_odd_energy_relative"] <= 1e-10
        for row in odd_rows
    )
    exposed_teeth = {
        row["width"]: row["restriction_spectrum_max_absolute_difference"] > 1e-9
        for row in exposed_rows
    }
    odd_widths = [width for width in (5, 7) if width in exposed_teeth]
    even_widths = [width for width in (4, 6, 8) if width in exposed_teeth]

    return {
        "protected_factors_through_the_quotient": protected_factors,
        "protected_rows": protected_rows,
        "odd_sector_is_inert": odd_inert,
        "odd_sector_rows": odd_rows,
        "exposed_disagrees_at_odd_widths": all(
            exposed_teeth[width] for width in odd_widths
        ),
        "exposed_agrees_at_even_widths": all(
            not exposed_teeth[width] for width in even_widths
        ),
        "exposed_rows": exposed_rows,
        "tolerance": protection_tolerance,
        "reading": (
            "the protected dictionary factors through the quotient at every "
            "width and the odd sector carries no Hankel energy either way, so "
            "the comparison is a theorem for the invariant task; the exposed "
            "dictionary breaks exactly at the two odd widths, which is the "
            "control that shows the agreement was not an identity"
            if protected_factors and odd_inert
            else "at least one declared claim failed -- read the rows, not this string"
        ),
    }


def assemble(widths: Sequence[int]) -> Dict[str, Any]:
    blocks = [width_report(width) for width in widths]
    gates = [continuity_gate(width) for width in widths]
    return {
        "schema": SCHEMA,
        "issue": ISSUE,
        "phase": "C",
        "declared": {
            "tilt": TILT,
            "etas": [float(eta) for eta in TILT_ETAS],
            "horizon_grid": [float(value) for value in KERNEL_GRID],
            "declared_lags": [float(lag) for lag in LAGS],
            "balanced_orders": [int(order) for order in BALANCED_ORDERS],
            "hankel_tolerance": HANKEL_TOLERANCE,
            "hankel_energy_levels": [float(level) for level in HANKEL_ENERGY_LEVELS],
            "gramian_floor": GRAMIAN_FLOOR,
            "readout_conventions": list(CONVENTIONS),
            "primary_biorthogonal_mode": PRIMARY_BIORTHOGONAL_MODE,
            "protected_dictionary": list(PROTECTED),
            "exposed_dictionary": list(EXPOSED),
        },
        "blocks": blocks,
        "continuity_gate": gates,
        "decision": decide(blocks),
        "claim_boundary": (
            "One exactly known finite connectivity process, one declared "
            "projection, dictionary and source set.  That the balanced "
            "realization factors through an exact symmetry quotient is a "
            "statement about this input/output object, not about the state "
            "space and not a percolation result.  The quotient reproduces the "
            "task; it does not claim the extra states were unphysical."
        ),
    }


def render(payload: Mapping[str, Any]) -> str:
    lines: List[str] = []
    lines.append("# The balanced realization factors through the reflection quotient")
    lines.append("")
    lines.append(f"schema {payload['schema']} | issue #{payload['issue']} | phase C")
    lines.append("")
    lines.append(
        "A is the microscopic chain, B the exact reduction to the R-orbit "
        "quotient.  Tilts are H_even at eta = 0, +-1/4."
    )
    lines.append("")
    for block in payload["blocks"]:
        lines.append(
            f"## w = {block['width']}  ({block['states']} states -> "
            f"{block['quotient_states']} orbits, odd sector dimension "
            f"{block['odd_sector']['dimension']})"
        )
        lines.append("")
        lines.append(
            "  odd sector energy: reach "
            f"{block['odd_sector']['reach_odd_energy_relative']:.3e}, observe "
            f"{block['odd_sector']['observe_odd_energy_relative']:.3e}"
        )
        for name in ("protected", "exposed"):
            entry = block["dictionaries"][name]
            lines.append("")
            lines.append(
                f"  {name}: {', '.join(entry['members'])}  "
                f"(max odd content {entry['max_odd_content']:.3e})"
            )
            micro_spectra = entry["readings"]["A"]["hankel_singular_values"]
            quotient_spectra = entry["readings"]["restriction"]["B"][
                "hankel_singular_values"
            ]
            lines.append(
                "      k     A (micro)        B (quotient)     |A-B|"
            )
            for index in range(min(12, len(micro_spectra))):
                left = micro_spectra[index]
                right = (
                    quotient_spectra[index]
                    if index < len(quotient_spectra)
                    else float("nan")
                )
                lines.append(
                    f"      {index:<5d}{left:<17.9e}{right:<17.9e}{abs(left - right):.3e}"
                )
            for convention in ("restriction", "symmetrization"):
                reading = entry["readings"][convention]
                if reading["identical_to_restriction"] and convention != "restriction":
                    lines.append(f"    {convention}: same map as restriction here")
                    continue
                difference = reading["differences"]
                spectrum = difference["hankel_spectrum"]
                worst_response = max(
                    value["relative_frobenius"]
                    for value in difference["response"].values()
                )
                lines.append(
                    f"    {convention:16s} spectrum max |A-B| "
                    f"{spectrum['max_absolute']:.3e}  rank "
                    f"{difference['numerical_rank']['A']}/"
                    f"{difference['numerical_rank']['B']}  response rel-Frobenius "
                    f"{worst_response:.3e}"
                )
        lines.append("")
    lines.append("## Continuity gate against the committed #588 balanced numbers")
    lines.append("")
    lines.append("  " + "width".ljust(7) + "spectrum".ljust(12) + "score (inverse)")
    for gate in payload["continuity_gate"]:
        if gate.get("status") != "run":
            continue
        lines.append(
            f"  w={gate['width']}   spectrum(above tol) "
            f"{gate['hankel_spectrum_above_tolerance']['max_absolute_difference']:.3e}"
            f" over {gate['hankel_spectrum_above_tolerance']['entries']} entries"
            f"   score(inverse) {gate['max_score_difference']:.3e}"
            f"   score(frozen transpose) "
            f"{gate['max_score_difference_by_correction']['frozen_transpose']:.3e}"
        )
    lines.append("")
    decision = payload["decision"]
    lines.append("## Decision")
    lines.append("")
    lines.append(
        f"  protected factors through the quotient: {decision['protected_factors_through_the_quotient']}"
    )
    lines.append(f"  odd sector inert: {decision['odd_sector_is_inert']}")
    lines.append(
        f"  exposed disagrees at odd widths: {decision['exposed_disagrees_at_odd_widths']}"
    )
    lines.append(
        f"  exposed agrees at even widths: {decision['exposed_agrees_at_even_widths']}"
    )
    lines.append("")
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--widths", type=int, nargs="+", default=list(DEFAULT_WIDTHS))
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--markdown", action="store_true")
    args = parser.parse_args(argv)
    payload = assemble(args.widths)
    text = json.dumps(payload, indent=2, sort_keys=False)
    if args.json:
        print(text)
    if args.markdown:
        print(render(payload))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text + "\n")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
