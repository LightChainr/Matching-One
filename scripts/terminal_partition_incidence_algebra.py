#!/usr/bin/env python3
"""Exact incidence algebra of the four-terminal partition lattice."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from math import factorial
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

try:
    from scripts.terminal_partition_canonical import RGS, enumerate_rgs, rgs_to_blocks, validate_rgs
except ModuleNotFoundError:
    from terminal_partition_canonical import RGS, enumerate_rgs, rgs_to_blocks, validate_rgs


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "analysis" / "terminal_partition_incidence_algebra_certificate.json"
SCHEMA = "matching-one/terminal-partition-incidence-algebra/v1"


def refines(fine: Sequence[int], coarse: Sequence[int]) -> bool:
    """Return whether every block of ``fine`` is contained in a block of ``coarse``."""

    fine = validate_rgs(fine)
    coarse = validate_rgs(coarse, len(fine))
    owner = {}
    for block in rgs_to_blocks(fine):
        labels = {coarse[index] for index in block}
        if len(labels) != 1:
            return False
        owner[block[0]] = next(iter(labels))
    return True


def mobius_closed_form(fine: Sequence[int], coarse: Sequence[int]) -> int:
    fine = validate_rgs(fine)
    coarse = validate_rgs(coarse, len(fine))
    if not refines(fine, coarse):
        return 0
    fine_blocks = rgs_to_blocks(fine)
    value = 1
    for coarse_block in rgs_to_blocks(coarse):
        contained = sum(set(block).issubset(coarse_block) for block in fine_blocks)
        value *= (-1) ** (contained - 1) * factorial(contained - 1)
    return value


def matrices(partitions: Sequence[RGS]) -> tuple[list[list[int]], list[list[int]]]:
    if len(partitions) != len(set(partitions)):
        raise ValueError("partition catalog contains duplicates")
    normalized = tuple(validate_rgs(value) for value in partitions)
    if len({len(value) for value in normalized}) != 1:
        raise ValueError("partition catalog mixes terminal counts")
    zeta = [[int(refines(left, right)) for right in normalized] for left in normalized]
    mobius = [[mobius_closed_form(left, right) for right in normalized] for left in normalized]
    return zeta, mobius


def matmul(left: Sequence[Sequence[int]], right: Sequence[Sequence[int]]) -> list[list[int]]:
    if not left or not right or any(len(row) != len(right) for row in left):
        raise ValueError("incompatible matrix dimensions")
    width = len(right[0])
    if any(len(row) != width for row in right):
        raise ValueError("ragged matrix")
    return [
        [sum(left[i][k] * right[k][j] for k in range(len(right))) for j in range(width)]
        for i in range(len(left))
    ]


def recurrence_mobius(zeta: Sequence[Sequence[int]]) -> list[list[int]]:
    """Independent interval recurrence, without using the partition closed form."""

    size = len(zeta)
    if size == 0 or any(len(row) != size for row in zeta):
        raise ValueError("zeta matrix must be square")
    ranks = [sum(row) for row in zeta]
    order = sorted(range(size), key=lambda index: ranks[index])
    result = [[0] * size for _ in range(size)]
    for left in reversed(order):
        result[left][left] = 1
        above = [right for right in order if right != left and zeta[left][right]]
        above.sort(key=lambda index: ranks[index], reverse=True)
        for right in above:
            result[left][right] = -sum(
                result[left][middle]
                for middle in range(size)
                if middle != right and zeta[left][middle] and zeta[middle][right]
            )
    return result


def transform(values: Sequence[int], matrix: Sequence[Sequence[int]]) -> tuple[int, ...]:
    size = len(matrix)
    if len(values) != size or any(len(row) != size for row in matrix):
        raise ValueError("vector and matrix dimensions must agree")
    if any(isinstance(value, bool) or not isinstance(value, int) for value in values):
        raise TypeError("incidence transforms require exact integer values")
    return tuple(sum(matrix[i][j] * values[j] for j in range(size)) for i in range(size))


def determinant_bareiss(matrix: Sequence[Sequence[int]]) -> int:
    work = [list(row) for row in matrix]
    size = len(work)
    if size == 0 or any(len(row) != size for row in work):
        raise ValueError("matrix must be nonempty and square")
    sign = 1
    previous = 1
    for pivot_index in range(size - 1):
        if work[pivot_index][pivot_index] == 0:
            swap = next((row for row in range(pivot_index + 1, size) if work[row][pivot_index]), None)
            if swap is None:
                return 0
            work[pivot_index], work[swap] = work[swap], work[pivot_index]
            sign *= -1
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = work[row][column] * pivot - work[row][pivot_index] * work[pivot_index][column]
                if numerator % previous:
                    raise ArithmeticError("Bareiss division was not exact")
                work[row][column] = numerator // previous
        previous = pivot
    return sign * work[-1][-1]


def build_artifact() -> dict[str, Any]:
    partitions = enumerate_rgs(4)
    zeta, mobius = matrices(partitions)
    recurrence = recurrence_mobius(zeta)
    identity = [[int(i == j) for j in range(len(partitions))] for i in range(len(partitions))]
    probe = tuple(range(1, len(partitions) + 1))
    upper_events = transform(probe, zeta)
    recovered = transform(upper_events, mobius)
    histogram = Counter(value for row in mobius for value in row)
    return {
        "schema": SCHEMA,
        "issue": 13,
        "status": "exact_four_terminal_partition_incidence_algebra",
        "orientation": "fine_refines_coarse",
        "partition_catalog": [list(value) for value in partitions],
        "zeta_matrix": zeta,
        "mobius_matrix": mobius,
        "counts": {
            "states": len(partitions),
            "ordered_pairs": len(partitions) ** 2,
            "comparable_pairs": sum(sum(row) for row in zeta),
            "mobius_value_histogram": {str(key): count for key, count in sorted(histogram.items())},
        },
        "unimodularity": {
            "zeta_determinant": determinant_bareiss(zeta),
            "mobius_determinant": determinant_bareiss(mobius),
        },
        "integer_probe": {
            "partition_weights": list(probe),
            "principal_upset_sums": list(upper_events),
            "recovered_weights": list(recovered),
        },
        "exact_checks": {
            "closed_form_matches_interval_recurrence": mobius == recurrence,
            "zeta_times_mobius_is_identity": matmul(zeta, mobius) == identity,
            "mobius_times_zeta_is_identity": matmul(mobius, zeta) == identity,
            "integer_probe_round_trips": recovered == probe,
            "all_diagonal_entries_are_one": all(zeta[i][i] == mobius[i][i] == 1 for i in range(len(partitions))),
        },
        "claim_boundary": {
            "included": "four-terminal partition refinement, zeta/Mobius matrices, and exact event inversion",
            "excluded": "noncrossing restriction, planar duality, self-duality, reliability data, thresholds, or bounds",
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    expected = build_artifact()
    if artifact != expected:
        raise ValueError("incidence-algebra artifact does not exactly reproduce")
    if set(expected["exact_checks"].values()) != {True}:
        raise ValueError("all exact checks must pass")
    return {
        "schema": SCHEMA,
        "status": "valid_exact_four_terminal_partition_incidence_algebra",
        "states": expected["counts"]["states"],
        "comparable_pairs": expected["counts"]["comparable_pairs"],
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        print(json.dumps(validate_artifact(json.loads(args.validate.read_text())), indent=2, sort_keys=True))
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
