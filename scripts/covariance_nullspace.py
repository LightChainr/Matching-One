#!/usr/bin/env python3
"""Shared generalized chi-square diagnostics with covariance-nullspace checks.

The active-subspace quadratic form is only interpretable after checking that
the standardized residual is compatible with every discarded eigendirection.
Estimated covariance matrices retain the active statistic but report cutoff
sensitivity. A caller that declares an exact structural covariance may ask for
fail-closed behavior instead.
"""

from __future__ import annotations

from typing import Callable, Sequence

import mpmath as mp


DEFAULT_SENSITIVITY_RELATIVE_CUTOFFS = (
    "1e-14",
    "1e-12",
    "1e-10",
    "1e-8",
    "1e-6",
)


class CovarianceNullspaceViolation(ValueError):
    """Raised when a structural covariance constraint is violated."""


def _number(value: object) -> mp.mpf:
    if isinstance(value, mp.mpf):
        return value
    return mp.mpf(value)


def _survival(chi_square: mp.mpf, degrees: int) -> mp.mpf:
    return mp.gammainc(
        mp.mpf(degrees) / 2, chi_square / 2, mp.inf
    ) / mp.gamma(mp.mpf(degrees) / 2)


def covariance_spectral_diagnostics(
    residual: Sequence[object],
    covariance: Sequence[Sequence[object]],
    relative_cutoff: object = "1e-10",
    *,
    nullspace_policy: str = "estimated",
    null_projection_tolerance: object = "1e-12",
    standardize: bool = True,
    sensitivity_relative_cutoffs: Sequence[object] = (
        DEFAULT_SENSITIVITY_RELATIVE_CUTOFFS
    ),
) -> dict[str, object]:
    """Return an active-subspace score plus explicit discarded-mode checks.

    ``nullspace_policy="estimated"`` is appropriate for jackknife or sample
    covariance matrices. It reports an incompatibility and the statistic's
    cutoff sensitivity without silently treating the active-subspace survival
    probability as a passing score. ``"structural"`` is reserved for a contract
    that declares every below-cutoff eigendirection to be an exact deterministic
    constraint; it raises when a discarded projection exceeds
    ``null_projection_tolerance``.
    """
    if nullspace_policy not in {"estimated", "structural"}:
        raise ValueError("nullspace_policy must be 'estimated' or 'structural'")
    if not residual:
        raise ValueError("residual vector is empty")
    dimension = len(residual)
    if len(covariance) != dimension or any(
        len(row) != dimension for row in covariance
    ):
        raise ValueError("residual covariance shape differs from residual vector")

    values = [_number(value) for value in residual]
    matrix = [[_number(value) for value in row] for row in covariance]
    if standardize:
        scales = [mp.sqrt(matrix[index][index]) for index in range(dimension)]
        if any(scale <= 0 for scale in scales):
            raise ValueError("residual covariance has nonpositive diagonal")
        correlation = mp.matrix(
            [
                [
                    matrix[left][right] / (scales[left] * scales[right])
                    for right in range(dimension)
                ]
                for left in range(dimension)
            ]
        )
        standardized = mp.matrix(
            [value / scale for value, scale in zip(values, scales)]
        )
        spectral_matrix = correlation
        spectral_vector = standardized
        spectral_basis = "correlation_standardized"
        nullspace_projection_basis = "component_standardized_correlation_basis"
    else:
        spectral_matrix = mp.matrix(matrix)
        spectral_vector = mp.matrix(values)
        spectral_basis = "raw_covariance"
        standardized = []
        nullspace_projection_basis = "spectral_scale_normalized_raw_covariance_basis"
    eigenvalues, eigenvectors = mp.eigsy(spectral_matrix)
    spectrum = [eigenvalues[index] for index in range(dimension)]
    largest = max(spectrum)
    if largest <= 0:
        raise ValueError("residual covariance has no positive eigenvalue")

    relative = _number(relative_cutoff)
    if relative <= 0:
        raise ValueError("relative eigenvalue cutoff must be positive")
    cutoff = largest * relative
    if min(spectrum) < -cutoff:
        raise ValueError("residual covariance spectrum is materially indefinite")

    tolerance = _number(null_projection_tolerance)
    if tolerance < 0:
        raise ValueError("null projection tolerance must be nonnegative")
    projections = [
        (eigenvectors[:, index].T * spectral_vector)[0]
        for index in range(dimension)
    ]
    compatibility_projections = (
        projections
        if standardize
        else [projection / mp.sqrt(largest) for projection in projections]
    )

    def at_cutoff(candidate: mp.mpf) -> dict[str, object]:
        absolute = largest * candidate
        active = [
            index for index, eigenvalue in enumerate(spectrum)
            if eigenvalue > absolute
        ]
        discarded = [index for index in range(dimension) if index not in active]
        discarded_l2 = mp.sqrt(
            mp.fsum(projections[index] ** 2 for index in discarded)
        )
        discarded_max = max(
            (abs(projections[index]) for index in discarded), default=mp.mpf(0)
        )
        compatibility_l2 = mp.sqrt(
            mp.fsum(compatibility_projections[index] ** 2 for index in discarded)
        )
        compatibility_max = max(
            (abs(compatibility_projections[index]) for index in discarded),
            default=mp.mpf(0),
        )
        compatible = compatibility_max <= tolerance
        if active:
            chi_square = mp.fsum(
                projections[index] ** 2 / spectrum[index] for index in active
            )
            survival = _survival(chi_square, len(active))
            condition = largest / min(spectrum[index] for index in active)
        else:
            chi_square = None
            survival = None
            condition = None
        return {
            "relative_eigenvalue_cutoff": candidate,
            "absolute_eigenvalue_cutoff": absolute,
            "numerical_rank": len(active),
            "chi_square": chi_square,
            "degrees_of_freedom": len(active),
            "chi_square_survival": survival,
            "active_condition_number": condition,
            "discarded_residual_projection_l2": discarded_l2,
            "max_abs_discarded_residual_projection": discarded_max,
            "discarded_nullspace_projection_l2": compatibility_l2,
            "max_abs_discarded_nullspace_projection": compatibility_max,
            "nullspace_compatible": compatible,
            "active_indices": active,
            "discarded_indices": discarded,
        }

    primary = at_cutoff(relative)
    if not primary["active_indices"]:
        raise ValueError("residual covariance has zero numerical rank")

    discarded_directions = [
        {
            "eigen_index": index,
            "spectral_eigenvalue": spectrum[index],
            "correlation_eigenvalue": spectrum[index] if standardize else None,
            "residual_projection": projections[index],
            "abs_residual_projection": abs(projections[index]),
            "nullspace_compatibility_projection": compatibility_projections[index],
            "abs_nullspace_compatibility_projection": abs(
                compatibility_projections[index]
            ),
            "standardized_residual_projection": (
                projections[index] if standardize else None
            ),
            "abs_standardized_residual_projection": (
                abs(projections[index]) if standardize else None
            ),
            "spectral_scale_normalized_residual_projection": (
                compatibility_projections[index] if not standardize else None
            ),
            "abs_spectral_scale_normalized_residual_projection": (
                abs(compatibility_projections[index]) if not standardize else None
            ),
        }
        for index in primary["discarded_indices"]
    ]
    compatible = bool(primary["nullspace_compatible"])
    if nullspace_policy == "structural" and not compatible:
        maximum = primary["max_abs_discarded_nullspace_projection"]
        raise CovarianceNullspaceViolation(
            "residual violates a structural covariance null "
            f"constraint: max discarded projection {mp.nstr(maximum, 12)} "
            f"> tolerance {mp.nstr(tolerance, 12)}"
        )

    if not primary["discarded_indices"]:
        status = "full_rank"
        interpretation = "ordinary_full_rank_generalized_chi_square"
    elif nullspace_policy == "structural":
        status = "structural_null_compatible"
        interpretation = "structural_constraints_satisfied"
    elif compatible:
        status = "estimated_near_null_compatible"
        interpretation = "active_subspace_statistic_null_compatible"
    else:
        status = "estimated_near_null_incompatibility"
        interpretation = "active_subspace_statistic_requires_cutoff_sensitivity_review"

    sensitivity_values = sorted(
        {
            _number(candidate)
            for candidate in (*sensitivity_relative_cutoffs, relative)
            if _number(candidate) > 0
        }
    )
    sensitivity = [at_cutoff(candidate) for candidate in sensitivity_values]

    return {
        "chi_square": primary["chi_square"],
        "degrees_of_freedom": primary["degrees_of_freedom"],
        "chi_square_survival": primary["chi_square_survival"],
        "numerical_rank": primary["numerical_rank"],
        "relative_eigenvalue_cutoff": relative,
        "absolute_eigenvalue_cutoff": cutoff,
        "spectral_basis": spectral_basis,
        "nullspace_projection_basis": nullspace_projection_basis,
        "active_condition_number": primary["active_condition_number"],
        "spectral_eigenvalues": spectrum,
        "correlation_eigenvalues": spectrum if standardize else [],
        "component_standardized_residuals": list(standardized),
        "eigenbasis_residual_projections": projections,
        "eigenbasis_standardized_residuals": projections if standardize else [],
        "discarded_eigendirections": discarded_directions,
        "discarded_residual_projection_l2": primary[
            "discarded_residual_projection_l2"
        ],
        "max_abs_discarded_residual_projection": primary[
            "max_abs_discarded_residual_projection"
        ],
        "discarded_nullspace_projection_l2": primary[
            "discarded_nullspace_projection_l2"
        ],
        "max_abs_discarded_nullspace_projection": primary[
            "max_abs_discarded_nullspace_projection"
        ],
        "discarded_standardized_residual_l2": (
            primary["discarded_residual_projection_l2"] if standardize else None
        ),
        "max_abs_discarded_standardized_residual": (
            primary["max_abs_discarded_residual_projection"] if standardize else None
        ),
        "null_projection_tolerance": tolerance,
        "covariance_nullspace_policy": nullspace_policy,
        "nullspace_status": status,
        "nullspace_compatible": compatible,
        "chi_square_interpretation": interpretation,
        "cutoff_sensitivity": sensitivity,
    }


def serialize_diagnostics(
    value: object, number_formatter: Callable[[mp.mpf], object]
) -> object:
    """Recursively convert mpmath numbers while preserving schema values."""
    if isinstance(value, dict):
        return {
            key: serialize_diagnostics(item, number_formatter)
            for key, item in value.items()
        }
    if isinstance(value, (list, tuple)):
        return [serialize_diagnostics(item, number_formatter) for item in value]
    if isinstance(value, mp.mpf):
        return number_formatter(value)
    return value
