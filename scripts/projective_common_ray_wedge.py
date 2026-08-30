#!/usr/bin/env python3
"""Exact aligned-batch wedge and delete-one jackknife control for Issue 439."""

from __future__ import annotations

import argparse
from fractions import Fraction
import itertools
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence, Union


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "analysis" / "projective_common_ray_wedge_certificate.json"
SCHEMA = "matching-one/projective-common-ray-wedge/v1"
ExactInput = Union[int, str, Fraction]


def exact_fraction(value: ExactInput, *, field: str) -> Fraction:
    if isinstance(value, bool) or isinstance(value, float):
        raise TypeError(f"{field} must be exact; floats and booleans are forbidden")
    try:
        return value if isinstance(value, Fraction) else Fraction(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"invalid exact value for {field}") from exc


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}" if value.denominator != 1 else str(value.numerator)


def determinant(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    size = len(matrix)
    if any(len(row) != size for row in matrix):
        raise ValueError("determinant requires a square matrix")
    if size == 0:
        return Fraction(1)
    if size == 1:
        return matrix[0][0]
    return sum(
        (-1) ** column
        * matrix[0][column]
        * determinant(
            [
                [matrix[row][other] for other in range(size) if other != column]
                for row in range(1, size)
            ]
        )
        for column in range(size)
    )


def matrix_rank(matrix: Sequence[Sequence[Fraction]]) -> int:
    work = [list(row) for row in matrix]
    if not work:
        return 0
    width = len(work[0])
    if any(len(row) != width for row in work):
        raise ValueError("matrix rows must have equal width")
    pivot_row = 0
    for column in range(width):
        pivot = next(
            (row for row in range(pivot_row, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][column]
        work[pivot_row] = [value / scale for value in work[pivot_row]]
        for row in range(len(work)):
            if row == pivot_row or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [
                value - factor * pivot_value
                for value, pivot_value in zip(work[row], work[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == len(work):
            break
    return pivot_row


def validate_sizes(sizes: Sequence[int]) -> tuple[int, ...]:
    normalized = tuple(sizes)
    if len(normalized) < 2 or any(isinstance(size, bool) or not isinstance(size, int) for size in normalized):
        raise ValueError("sizes must contain at least two exact integers")
    if any(size <= 0 for size in normalized):
        raise ValueError("sizes must be positive")
    if any(right != 2 * left for left, right in zip(normalized, normalized[1:])):
        raise ValueError("sizes must form an ordered doubling lineage")
    return normalized


def parse_batches(
    records: Sequence[Mapping[str, Any]], sizes: Sequence[int]
) -> tuple[tuple[str, dict[int, tuple[Fraction, Fraction]]], ...]:
    lineage = validate_sizes(sizes)
    if len(records) < 2:
        raise ValueError("delete-one jackknife requires at least two aligned batches")
    parsed = []
    seen = set()
    for record in records:
        if set(record) != {"batch", "values"}:
            raise ValueError("batch fields must be exactly batch,values")
        batch = record["batch"]
        if not isinstance(batch, str) or not batch.strip() or batch in seen:
            raise ValueError("batch identifiers must be unique nonempty strings")
        seen.add(batch)
        values = record["values"]
        if not isinstance(values, Mapping) or set(values) != {str(size) for size in lineage}:
            raise ValueError("each batch must contain exactly the declared generation sizes")
        parsed_values = {}
        for size in lineage:
            pair = values[str(size)]
            if not isinstance(pair, Mapping) or set(pair) != {"A_M", "A_K"}:
                raise ValueError("generation fields must be exactly A_M,A_K")
            parsed_values[size] = (
                exact_fraction(pair["A_M"], field=f"{batch}.{size}.A_M"),
                exact_fraction(pair["A_K"], field=f"{batch}.{size}.A_K"),
            )
        parsed.append((batch, parsed_values))
    return tuple(parsed)


def mean_coordinates(
    batches: Sequence[tuple[str, Mapping[int, tuple[Fraction, Fraction]]]],
    sizes: Sequence[int],
) -> dict[int, tuple[Fraction, Fraction]]:
    count = len(batches)
    if not count:
        raise ValueError("cannot average zero batches")
    return {
        size: (
            sum(values[size][0] for _, values in batches) / count,
            sum(values[size][1] for _, values in batches) / count,
        )
        for size in sizes
    }


def adjacent_wedges(
    coordinates: Mapping[int, tuple[Fraction, Fraction]], sizes: Sequence[int]
) -> tuple[Fraction, ...]:
    wedges = []
    for left, right in zip(sizes, sizes[1:]):
        m_left, k_left = coordinates[left]
        m_right, k_right = coordinates[right]
        wedges.append(m_right * k_left - k_right * m_left)
    return tuple(wedges)


def principal_minors(matrix: Sequence[Sequence[Fraction]]) -> dict[str, Fraction]:
    result = {}
    for order in range(1, len(matrix) + 1):
        for indices in itertools.combinations(range(len(matrix)), order):
            result[",".join(str(index) for index in indices)] = determinant(
                [[matrix[row][column] for column in indices] for row in indices]
            )
    return result


def analyze(records: Sequence[Mapping[str, Any]], sizes_source: Sequence[int]) -> dict[str, Any]:
    sizes = validate_sizes(sizes_source)
    batches = parse_batches(records, sizes)
    full_coordinates = mean_coordinates(batches, sizes)
    full_wedges = adjacent_wedges(full_coordinates, sizes)

    replicates = []
    for omitted in range(len(batches)):
        kept = batches[:omitted] + batches[omitted + 1 :]
        coordinates = mean_coordinates(kept, sizes)
        replicates.append((batches[omitted][0], adjacent_wedges(coordinates, sizes)))

    width = len(full_wedges)
    replicate_center = tuple(
        sum(vector[index] for _, vector in replicates) / len(replicates)
        for index in range(width)
    )
    factor = Fraction(len(batches) - 1, len(batches))
    covariance = [
        [
            factor
            * sum(
                (vector[row] - replicate_center[row])
                * (vector[column] - replicate_center[column])
                for _, vector in replicates
            )
            for column in range(width)
        ]
        for row in range(width)
    ]
    minors = principal_minors(covariance)
    if any(value < 0 for value in minors.values()):
        raise ArithmeticError("exact jackknife covariance is not positive semidefinite")

    return {
        "sizes": list(sizes),
        "batch_count": len(batches),
        "full_coordinates": {
            str(size): {"A_M": fraction_text(pair[0]), "A_K": fraction_text(pair[1])}
            for size, pair in full_coordinates.items()
        },
        "full_wedges": {
            str(size): fraction_text(value) for size, value in zip(sizes[:-1], full_wedges)
        },
        "delete_one_replicates": [
            {
                "omitted_batch": batch,
                "wedges": {
                    str(size): fraction_text(value)
                    for size, value in zip(sizes[:-1], vector)
                },
            }
            for batch, vector in replicates
        ],
        "jackknife_center": [fraction_text(value) for value in replicate_center],
        "jackknife_covariance": [
            [fraction_text(value) for value in row] for row in covariance
        ],
        "covariance_rank": matrix_rank(covariance),
        "covariance_principal_minors": {
            key: fraction_text(value) for key, value in minors.items()
        },
        "exact_checks": {
            "all_replicates_recomputed_from_aligned_batches": True,
            "covariance_symmetric": all(
                covariance[row][column] == covariance[column][row]
                for row in range(width)
                for column in range(width)
            ),
            "covariance_psd": True,
        },
    }


def synthetic_batches(loadings: Sequence[ExactInput]) -> list[dict[str, Any]]:
    sizes = (85, 170, 340, 680)
    loading_values = [exact_fraction(value, field="loading") for value in loadings]
    if len(loading_values) != len(sizes):
        raise ValueError("one loading is required per generation")
    rows = []
    for index, amplitude in enumerate((1, 2, 3, 4)):
        rows.append(
            {
                "batch": f"b{index}",
                "values": {
                    str(size): {
                        "A_M": fraction_text(Fraction(amplitude, 2**generation)),
                        "A_K": fraction_text(
                            loading_values[generation] * Fraction(amplitude, 2**generation)
                        ),
                    }
                    for generation, size in enumerate(sizes)
                },
            }
        )
    return rows


def build_artifact() -> dict[str, Any]:
    sizes = (85, 170, 340, 680)
    common = analyze(synthetic_batches((2, 2, 2, 2)), sizes)
    drift = analyze(synthetic_batches((2, "5/2", 3, "7/2")), sizes)
    if any(value != "0" for value in common["full_wedges"].values()):
        raise AssertionError("common-ray control produced a nonzero wedge")
    if common["covariance_rank"] != 0:
        raise AssertionError("common-ray control produced jackknife variance")
    if any(value == "0" for value in drift["full_wedges"].values()):
        raise AssertionError("loading-drift control did not separate from the common ray")
    if drift["covariance_rank"] != 1:
        raise AssertionError("loading-drift covariance rank changed")
    return {
        "schema": SCHEMA,
        "issue": 439,
        "status": "exact_aligned_batch_wedge_jackknife_control",
        "wedge_definition": "D_N=A_M(2N)*A_K(N)-A_K(2N)*A_M(N)",
        "common_ray_control": common,
        "loading_drift_control": drift,
        "claim_boundary": {
            "included": "exact wedge and delete-one aligned-batch jackknife algebra",
            "excluded": "raw archives, empirical covariance, model selection, rejection, or physics claims",
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    expected = build_artifact()
    if artifact != expected:
        raise ValueError("common-ray wedge certificate does not exactly reproduce")
    return {
        "schema": SCHEMA,
        "status": "valid_exact_aligned_batch_wedge_jackknife_control",
        "common_wedge_count": len(expected["common_ray_control"]["full_wedges"]),
        "drift_wedge_count": len(expected["loading_drift_control"]["full_wedges"]),
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        artifact = json.loads(args.validate.read_text(encoding="utf-8"))
        print(json.dumps(validate_artifact(artifact), indent=2, sort_keys=True))
        return 0
    rendered = json.dumps(build_artifact(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
