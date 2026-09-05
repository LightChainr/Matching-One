#!/usr/bin/env python3
"""Test a measured response vector against a model's ray or subspace, without
ever forming a ratio.

Why this module exists.  The N=580 modulus ladder measured a three-entry
response `(A4(i), A4(2i), A4(4i))` and then scored competitors on the ratio
`A4(4i)/A4(i)`.  `A4(i)` came in 3.6 sigma from zero, so the ratio's sampling
distribution was nothing like normal and its standard error at the observed
point understated the spread out where the large predictions sit.  Every
verdict computed that way was wrong; three of eight flipped when the test was
moved onto the linear contrast `Y - R0 X` instead (Fieller).

Fieller is the two-entry, one-ray case of a general fact: **a model that
predicts only proportions predicts a ray, and a model that predicts a
low-dimensional family predicts a subspace.**  Neither is a statement about any
single coordinate, so neither should be tested by dividing one coordinate by
another.  Dividing chooses a denominator, and the choice is arbitrary; when the
chosen denominator is poorly resolved the test degrades for a reason that has
nothing to do with the physics.

What to compute instead.  For a model ray `v` and a measured `y` with covariance
`S`, the statistic is the covariance-weighted distance from `y` to the line
through `v`:

    D(v) = min_a (y - a v)^T S^+ (y - a v),

referred to chi-square with `rank(S) - 1` degrees of freedom.  For a model whose
image is an `r`-dimensional subspace `V` the same expression minimises over all
of `V` and carries `rank(S) - r` degrees of freedom.  There is no denominator,
so there is nothing to be badly conditioned except `S` itself -- and `S`'s
condition number is reported rather than hidden, because that *is* the honest
statement of how much the geometry can support.

This also buys something the ratio test could not have.  A systematic that
enters along a known direction -- the spin-8 leakage in the ladder's r=2 rung,
say -- can be carried as an extra basis vector of the model subspace instead of
by discarding the contaminated entry.  A nuisance direction costs one degree of
freedom; discarding an entry costs the whole entry.

Precision.  The linear algebra runs in mpmath at a caller-chosen precision, not
in float64.  The subject of the module is conditioning, and a routine that
silently lost three digits inverting an ill-conditioned covariance would be
answering a question about itself rather than about the data.
"""

from __future__ import annotations

from typing import Any, Mapping, Sequence

from mpmath import mp


DEFAULT_DPS = 50
# Eigenvalues below this fraction of the largest are treated as null directions
# of the covariance rather than as very small variances.  A measured covariance
# with a genuine null direction is a covariance whose data cannot speak in that
# direction at all; inverting it there manufactures certainty.
DEFAULT_RELATIVE_TOLERANCE = "1e-12"


def _as_matrix(rows: Sequence[Sequence[Any]]):
    size = len(rows)
    if any(len(row) != size for row in rows):
        raise ValueError("covariance must be square")
    return mp.matrix([[mp.mpf(str(value)) for value in row] for row in rows])


def _as_vector(values: Sequence[Any]):
    return mp.matrix([mp.mpf(str(value)) for value in values])


def spectral_pseudo_inverse(covariance: Sequence[Sequence[Any]],
                            relative_tolerance: str = DEFAULT_RELATIVE_TOLERANCE):
    """Moore-Penrose inverse of a symmetric covariance, with its rank and condition.

    Returns ``(pinv, rank, condition_number, eigenvalues)``.  Directions whose
    eigenvalue falls below ``relative_tolerance`` times the largest are dropped:
    a covariance estimated from a finite jackknife can be singular by
    construction (100 batches cannot support more than 99 directions), and
    inverting such a direction invents information.
    """
    matrix = _as_matrix(covariance)
    size = matrix.rows
    for i in range(size):
        for j in range(i + 1, size):
            if abs(matrix[i, j] - matrix[j, i]) > mp.mpf("1e-30") * (
                    abs(matrix[i, j]) + abs(matrix[j, i]) + mp.mpf(1)):
                raise ValueError(f"covariance is not symmetric at ({i}, {j})")
            matrix[j, i] = matrix[i, j]
    values, vectors = mp.eigsy(matrix)
    eigenvalues = [values[i] for i in range(size)]
    largest = max(eigenvalues)
    if largest <= 0:
        raise ValueError("covariance has no positive direction")
    floor = mp.mpf(relative_tolerance) * largest
    kept = [i for i in range(size) if eigenvalues[i] > floor]
    if not kept:
        raise ValueError("covariance has no direction above the tolerance")
    pinv = mp.zeros(size, size)
    for index in kept:
        inverse = 1 / eigenvalues[index]
        for i in range(size):
            for j in range(size):
                pinv[i, j] += inverse * vectors[i, index] * vectors[j, index]
    smallest_kept = min(eigenvalues[i] for i in kept)
    return pinv, len(kept), largest / smallest_kept, eigenvalues


def subspace_residual(observed: Sequence[Any],
                      covariance: Sequence[Sequence[Any]],
                      basis: Sequence[Sequence[Any]],
                      relative_tolerance: str = DEFAULT_RELATIVE_TOLERANCE,
                      dps: int = DEFAULT_DPS) -> dict[str, Any]:
    """Covariance-weighted distance from ``observed`` to the span of ``basis``.

    ``basis`` is a sequence of direction vectors.  One direction is a ray -- a
    model that predicts proportions only.  Two or more span the image of a model
    with that many free amplitudes, or a model plus a known nuisance direction.
    """
    with mp.workdps(dps):
        y = _as_vector(observed)
        size = y.rows
        directions = [_as_vector(row) for row in basis]
        if not directions:
            raise ValueError("need at least one model direction")
        for direction in directions:
            if direction.rows != size:
                raise ValueError("model direction has the wrong length")
        pinv, rank, condition, _ = spectral_pseudo_inverse(covariance, relative_tolerance)
        if pinv.rows != size:
            raise ValueError("covariance and observation disagree on dimension")

        # Normal equations for min_a (y - Va)^T S^+ (y - Va): (V^T S^+ V) a = V^T S^+ y.
        width = len(directions)
        gram = mp.zeros(width, width)
        rhs = mp.zeros(width, 1)
        weighted = [pinv * direction for direction in directions]
        for i in range(width):
            for j in range(width):
                gram[i, j] = (directions[i].T * weighted[j])[0]
            rhs[i] = (directions[i].T * (pinv * y))[0]
        model_rank = _rank_of(gram, relative_tolerance)
        if model_rank < width:
            raise ValueError(
                "model directions are linearly dependent after weighting; drop one"
            )
        amplitudes = mp.lu_solve(gram, rhs)
        residual = y - sum((amplitudes[i] * directions[i] for i in range(width)),
                           mp.zeros(size, 1))
        statistic = (residual.T * (pinv * residual))[0]
        degrees = rank - width
        if degrees < 0:
            raise ValueError("model has more directions than the covariance has rank")
        return {
            "statistic": float(statistic),
            "degrees_of_freedom": degrees,
            "amplitudes": [float(amplitudes[i]) for i in range(width)],
            "covariance_rank": rank,
            "covariance_condition_number": float(condition),
            "p_value": chi_square_upper_tail(float(statistic), degrees),
            "equivalent_sigma": _equivalent_sigma(float(statistic), degrees),
        }


def ray_residual(observed: Sequence[Any],
                 covariance: Sequence[Sequence[Any]],
                 direction: Sequence[Any],
                 relative_tolerance: str = DEFAULT_RELATIVE_TOLERANCE,
                 dps: int = DEFAULT_DPS) -> dict[str, Any]:
    """The one-direction case: a model that predicts proportions and nothing else.

    For a two-entry observation this is exactly Fieller's test on the ratio of
    the two entries, squared -- see the tests.  For three or more entries it is
    the thing Fieller cannot do: it uses every entry at once, and never has to
    nominate one of them as a denominator.
    """
    return subspace_residual(observed, covariance, [direction], relative_tolerance, dps)


def _rank_of(matrix, relative_tolerance: str) -> int:
    values, _ = mp.eigsy(matrix)
    magnitudes = [abs(values[i]) for i in range(matrix.rows)]
    largest = max(magnitudes)
    if largest == 0:
        return 0
    floor = mp.mpf(relative_tolerance) * largest
    return sum(1 for value in magnitudes if value > floor)


def chi_square_upper_tail(statistic: float, degrees: int) -> float:
    """P(chi^2_degrees >= statistic), via the regularised incomplete gamma."""
    if degrees <= 0:
        return 1.0 if statistic <= 0 else 0.0
    if statistic <= 0:
        return 1.0
    with mp.workdps(40):
        return float(mp.gammainc(mp.mpf(degrees) / 2, mp.mpf(statistic) / 2,
                                 mp.inf, regularized=True))


def _equivalent_sigma(statistic: float, degrees: int) -> float:
    """The two-sided normal deviate with the same tail probability.

    Reported so that a many-degree-of-freedom result can be compared with the
    one-degree-of-freedom numbers the project already quotes, without anyone
    having to take a square root that is only correct when ``degrees == 1``.
    """
    tail = chi_square_upper_tail(statistic, degrees)
    if tail <= 0.0:
        return float("inf")
    if tail >= 1.0:
        return 0.0
    with mp.workdps(40):
        return float(-mp.sqrt(2) * mp.erfinv(mp.mpf(tail) - 1))
