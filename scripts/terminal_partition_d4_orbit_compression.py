#!/usr/bin/env python3
"""Exact no-go certificate for deterministic D4-orbit serial composition."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

try:
    from scripts.terminal_partition_canonical import RGS, enumerate_rgs
    from scripts.terminal_partition_d4_equivariance import d4_group, d4_orbits
    from scripts.terminal_partition_serial_category import serial_compose
except ModuleNotFoundError:
    from terminal_partition_canonical import RGS, enumerate_rgs
    from terminal_partition_d4_equivariance import d4_group, d4_orbits
    from terminal_partition_serial_category import serial_compose


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "analysis" / "terminal_partition_d4_orbit_compression_certificate.json"
SCHEMA = "matching-one/terminal-partition-d4-orbit-compression/v1"


def orbit_output_counts(
    partitions: Sequence[RGS], orbits: Sequence[Sequence[int]]
) -> list[list[Counter[int]]]:
    orbit_index = {member: orbit for orbit, members in enumerate(orbits) for member in members}
    if len(orbit_index) != len(partitions):
        raise ValueError("orbits do not partition the state catalog")
    index = {partition: position for position, partition in enumerate(partitions)}
    return [
        [
            Counter(
                orbit_index[index[serial_compose(partitions[left], partitions[right])]]
                for left in left_members
                for right in right_members
            )
            for right_members in orbits
        ]
        for left_members in orbits
    ]


def normalized_kernel(counts: Sequence[Sequence[Counter[int]]], orbit_count: int) -> list[list[tuple[Fraction, ...]]]:
    result = []
    for row in counts:
        output_row = []
        for counter in row:
            denominator = sum(counter.values())
            if denominator <= 0:
                raise ValueError("orbit-pair output count is empty")
            output_row.append(tuple(Fraction(counter[index], denominator) for index in range(orbit_count)))
        result.append(output_row)
    return result


def compose_distributions(
    left: Sequence[Fraction], right: Sequence[Fraction], kernel: Sequence[Sequence[Sequence[Fraction]]]
) -> tuple[Fraction, ...]:
    size = len(kernel)
    if len(left) != size or len(right) != size or any(len(row) != size for row in kernel):
        raise ValueError("distribution and kernel dimensions must agree")
    if any(len(cell) != size for row in kernel for cell in row):
        raise ValueError("kernel outputs have the wrong dimension")
    result = [Fraction(0) for _ in range(size)]
    for left_index, left_weight in enumerate(left):
        for right_index, right_weight in enumerate(right):
            for output_index, probability in enumerate(kernel[left_index][right_index]):
                result[output_index] += left_weight * right_weight * probability
    return tuple(result)


def fraction_vector(values: Sequence[Fraction]) -> list[str]:
    return [str(value) for value in values]


def build_artifact() -> dict[str, Any]:
    partitions = enumerate_rgs(4)
    orbits = d4_orbits(partitions, d4_group())
    counts = orbit_output_counts(partitions, orbits)
    kernel = normalized_kernel(counts, len(orbits))
    ambiguous = []
    for left in range(len(orbits)):
        for right in range(len(orbits)):
            if len(counts[left][right]) > 1:
                ambiguous.append((len(orbits[left]) * len(orbits[right]), left, right, counts[left][right]))
    ambiguous.sort(key=lambda row: (row[0], row[1], row[2]))
    _, first_left, first_right, first_counts = ambiguous[0]

    labelled_witnesses = []
    output_orbit = {}
    orbit_index = {member: orbit for orbit, members in enumerate(orbits) for member in members}
    partition_index = {partition: index for index, partition in enumerate(partitions)}
    for left in orbits[first_left]:
        for right in orbits[first_right]:
            output = serial_compose(partitions[left], partitions[right])
            result_orbit = orbit_index[partition_index[output]]
            if result_orbit not in output_orbit:
                output_orbit[result_orbit] = True
                labelled_witnesses.append(
                    {
                        "left": list(partitions[left]),
                        "right": list(partitions[right]),
                        "output": list(output),
                        "output_orbit": result_orbit,
                    }
                )

    associativity_failures = 0
    first_associativity = None
    basis = [tuple(Fraction(int(index == selected)) for index in range(len(orbits))) for selected in range(len(orbits))]
    for left in range(len(orbits)):
        for middle in range(len(orbits)):
            for right in range(len(orbits)):
                left_grouped = compose_distributions(
                    compose_distributions(basis[left], basis[middle], kernel), basis[right], kernel
                )
                right_grouped = compose_distributions(
                    basis[left], compose_distributions(basis[middle], basis[right], kernel), kernel
                )
                if left_grouped != right_grouped:
                    associativity_failures += 1
                    if first_associativity is None:
                        first_associativity = {
                            "input_orbits": [left, middle, right],
                            "left_grouped": fraction_vector(left_grouped),
                            "right_grouped": fraction_vector(right_grouped),
                        }

    serialized_counts = [
        [
            {str(output): count for output, count in sorted(counter.items())}
            for counter in row
        ]
        for row in counts
    ]
    serialized_kernel = [
        [[str(value) for value in cell] for cell in row]
        for row in kernel
    ]
    return {
        "schema": SCHEMA,
        "issue": 13,
        "status": "exact_d4_orbit_compression_no_go",
        "partition_catalog": [list(value) for value in partitions],
        "d4_orbits": [[list(partitions[index]) for index in orbit] for orbit in orbits],
        "orbit_pair_output_counts": serialized_counts,
        "uniform_orbit_pair_kernel": serialized_kernel,
        "deterministic_quotient": {
            "orbit_pairs": len(orbits) ** 2,
            "unambiguous_pairs": len(orbits) ** 2 - len(ambiguous),
            "ambiguous_pairs": len(ambiguous),
            "output_support_size_histogram_on_ambiguous_pairs": {
                str(size): count
                for size, count in sorted(Counter(len(row[3]) for row in ambiguous).items())
            },
            "smallest_counterexample": {
                "left_orbit": first_left,
                "right_orbit": first_right,
                "labelled_pair_count": sum(first_counts.values()),
                "output_counts": {str(key): value for key, value in sorted(first_counts.items())},
                "labelled_witnesses": labelled_witnesses,
            },
        },
        "averaging_boundary": {
            "basis_triples_checked": len(orbits) ** 3,
            "associativity_failures": associativity_failures,
            "first_failure": first_associativity,
            "interpretation": "uniform orbit averaging is normalized but is not an associative exact quotient",
        },
        "exact_checks": {
            "seven_d4_orbits_cover_15_states": len(orbits) == 7 and sum(map(len, orbits)) == 15,
            "deterministic_orbit_quotient_fails": len(ambiguous) > 0,
            "smallest_counterexample_has_two_labelled_pairs": sum(first_counts.values()) == 2,
            "every_uniform_kernel_row_has_unit_mass": all(
                sum(cell) == 1 for row in kernel for cell in row
            ),
            "uniform_orbit_averaging_is_nonassociative": associativity_failures > 0,
        },
        "claim_boundary": {
            "included": "exact information loss under D4 orbit compression for typed serial composition",
            "excluded": "preferred gauge fixing, noncrossing reduction, planar duality, word search, thresholds, or bounds",
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    expected = build_artifact()
    if artifact != expected:
        raise ValueError("D4-orbit compression artifact does not exactly reproduce")
    if set(expected["exact_checks"].values()) != {True}:
        raise ValueError("all exact checks must pass")
    return {
        "schema": SCHEMA,
        "status": "valid_exact_d4_orbit_compression_no_go",
        "orbits": len(expected["d4_orbits"]),
        "ambiguous_pairs": expected["deterministic_quotient"]["ambiguous_pairs"],
        "averaged_associativity_failures": expected["averaging_boundary"]["associativity_failures"],
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
