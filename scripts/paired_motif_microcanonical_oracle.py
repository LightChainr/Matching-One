#!/usr/bin/env python3
"""Exact fixed-rank certificate for paired same-N motif controls."""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from math import comb
from pathlib import Path
from typing import Any, Iterable

from integer_period_torus import (
    IntegerTorusGeometry,
    axis_integer_torus,
    diamond_integer_torus,
    gaussian_integer_torus,
    integer_torus_geometry,
    matrix_product,
)


Embedding = tuple[int, ...]
MotifTable = dict[str, tuple[Embedding, ...]]


MOTIF_SIZES = {"nn_edge": 2, "diagonal_pair": 2, "face": 4, "right_angle": 3}
UNIMODULAR_CHANGE = ((1, 1), (0, 1))


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _embedding(vertices: Iterable[int], expected_size: int) -> Embedding:
    result = tuple(sorted(set(vertices)))
    if len(result) != expected_size:
        raise ValueError(f"motif collapsed to {len(result)} vertices; expected {expected_size}")
    return result


def motif_embeddings(geometry: IntegerTorusGeometry) -> MotifTable:
    edges = tuple(_embedding((edge.i, edge.j), 2) for edge in geometry.primal_edges)
    diagonals = []
    faces = []
    corners = []
    for vertex, (x, y) in enumerate(geometry.coordinates):
        diagonals.extend([
            _embedding((vertex, geometry.vertex((x + 1, y + 1))), 2),
            _embedding((vertex, geometry.vertex((x - 1, y + 1))), 2),
        ])
        faces.append(_embedding((
            vertex,
            geometry.vertex((x + 1, y)),
            geometry.vertex((x, y + 1)),
            geometry.vertex((x + 1, y + 1)),
        ), 4))
        corners.append(_embedding((
            vertex,
            geometry.vertex((x + 1, y)),
            geometry.vertex((x, y + 1)),
        ), 3))
    return {
        "nn_edge": edges,
        "diagonal_pair": tuple(diagonals),
        "face": tuple(faces),
        "right_angle": tuple(corners),
    }


def direct_count(mask: int, embeddings: tuple[Embedding, ...]) -> int:
    return sum(all(mask & (1 << vertex) for vertex in item) for item in embeddings)


def popcount(mask: int) -> int:
    """Return the number of set bits on every supported Python version."""
    return bin(mask).count("1")


def incremental_count(mask: int, n: int, embeddings: tuple[Embedding, ...]) -> int:
    incident: dict[int, list[Embedding]] = defaultdict(list)
    for item in embeddings:
        for vertex in item:
            incident[vertex].append(item)
    active = [False] * n
    count = 0
    for vertex in range(n):
        if not mask & (1 << vertex):
            continue
        active[vertex] = True
        count += sum(all(active[member] for member in item) for item in incident[vertex])
    return count


def expected_count(n: int, k: int, embeddings: int, motif_size: int) -> Fraction:
    if not 0 <= k <= n:
        raise ValueError("K must lie in [0,N]")
    if not 0 <= motif_size <= n:
        raise ValueError("motif size must lie in [0,N]")
    return Fraction(embeddings * comb(k, motif_size), comb(n, motif_size)) if k >= motif_size else Fraction(0)


def exact_summary(geometry: IntegerTorusGeometry) -> dict[str, Any]:
    motifs = motif_embeddings(geometry)
    sums = {name: [0] * (geometry.n + 1) for name in motifs}
    incremental_failures = 0
    for mask in range(1 << geometry.n):
        k = popcount(mask)
        for name, embeddings in motifs.items():
            direct = direct_count(mask, embeddings)
            sums[name][k] += direct
            incremental_failures += incremental_count(mask, geometry.n, embeddings) != direct

    formula_failures = 0
    for name, values in sums.items():
        size = MOTIF_SIZES[name]
        for k, total in enumerate(values):
            empirical = Fraction(total, comb(geometry.n, k))
            formula_failures += empirical != expected_count(geometry.n, k, len(motifs[name]), size)

    return {
        "name": geometry.name,
        "N": geometry.n,
        "embedding_multiplicities": {name: len(items) for name, items in motifs.items()},
        "checked_masks": 1 << geometry.n,
        "checked_K_values": geometry.n + 1,
        "formula_failures": formula_failures,
        "incremental_failures": incremental_failures,
        "microcanonical_sums": sums,
    }


def paired_summary(first: IntegerTorusGeometry, second: IntegerTorusGeometry) -> dict[str, Any]:
    if first.n != second.n:
        raise ValueError("paired geometries must have the same N")
    left = exact_summary(first)
    right = exact_summary(second)
    multiplicities_equal = left["embedding_multiplicities"] == right["embedding_multiplicities"]
    differences = {
        name: [a - b for a, b in zip(left["microcanonical_sums"][name], right["microcanonical_sums"][name])]
        for name in MOTIF_SIZES
    }
    configurationwise_nonzero_masks = 0
    left_motifs = motif_embeddings(first)
    right_motifs = motif_embeddings(second)
    for mask in range(1 << first.n):
        if any(
            direct_count(mask, left_motifs[name]) != direct_count(mask, right_motifs[name])
            for name in MOTIF_SIZES
        ):
            configurationwise_nonzero_masks += 1
    return {
        "pair": [first.name, second.name],
        "N": first.n,
        "embedding_multiplicities": left["embedding_multiplicities"],
        "multiplicities_equal": multiplicities_equal,
        "checked_shared_masks": 1 << first.n,
        "checked_K_values": first.n + 1,
        "configurationwise_nonzero_masks": configurationwise_nonzero_masks,
        "max_absolute_fixed_K_difference_sum": {
            name: max(map(abs, values), default=0) for name, values in differences.items()
        },
        "all_fixed_K_difference_sums_zero": all(not any(values) for values in differences.values()),
        "formula_failures": left["formula_failures"] + right["formula_failures"],
        "incremental_failures": left["incremental_failures"] + right["incremental_failures"],
    }


def declared_pair_gate(
    first_parameters: tuple[int, int],
    second_parameters: tuple[int, int],
    witness_vertices: tuple[int, ...],
) -> dict[str, Any]:
    first = gaussian_integer_torus(*first_parameters)
    second = gaussian_integer_torus(*second_parameters)
    if first.n != second.n:
        raise ValueError("declared pair must have the same norm")
    left = motif_embeddings(first)
    right = motif_embeddings(second)
    left_multiplicities = {name: len(items) for name, items in left.items()}
    right_multiplicities = {name: len(items) for name, items in right.items()}
    witness_mask = sum(1 << vertex for vertex in witness_vertices)
    witness_differences = {
        name: direct_count(witness_mask, left[name]) - direct_count(witness_mask, right[name])
        for name in MOTIF_SIZES
    }
    incremental_differences = {
        name: incremental_count(witness_mask, first.n, left[name])
        - incremental_count(witness_mask, second.n, right[name])
        for name in MOTIF_SIZES
    }
    return {
        "pair": [first.name, second.name],
        "N": first.n,
        "embedding_multiplicities": left_multiplicities,
        "multiplicities_equal": left_multiplicities == right_multiplicities,
        "fixed_K_mean_difference": {name: "0 for every K=0,...,N" for name in MOTIF_SIZES},
        "nontrivial_witness_vertices": witness_vertices,
        "witness_K": len(witness_vertices),
        "witness_differences": witness_differences,
        "incremental_witness_agrees": incremental_differences == witness_differences,
        "all_witness_differences_nonzero": all(witness_differences.values()),
    }


def joint_histogram(geometry: IntegerTorusGeometry) -> dict[int, Counter[tuple[int, ...]]]:
    motifs = motif_embeddings(geometry)
    order = tuple(MOTIF_SIZES)
    result: dict[int, Counter[tuple[int, ...]]] = defaultdict(Counter)
    for mask in range(1 << geometry.n):
        result[popcount(mask)][tuple(direct_count(mask, motifs[name]) for name in order)] += 1
    return result


def basis_invariance_summary(geometry: IntegerTorusGeometry) -> dict[str, Any]:
    changed_matrix = matrix_product(geometry.periods.matrix, UNIMODULAR_CHANGE)
    changed = integer_torus_geometry(changed_matrix, name=f"{geometry.name}-basis-changed")
    return {
        "name": geometry.name,
        "N": geometry.n,
        "change": UNIMODULAR_CHANGE,
        "determinant": 1,
        "joint_histograms_equal": joint_histogram(geometry) == joint_histogram(changed),
    }


def compact_control(summary: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in summary.items() if key != "microcanonical_sums"}


def serialize(value: Any) -> Any:
    if isinstance(value, Fraction):
        return fraction_text(value)
    if isinstance(value, dict):
        return {str(key): serialize(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize(item) for item in value]
    return value


@lru_cache(maxsize=1)
def build_artifact() -> dict[str, Any]:
    gaussian5 = (gaussian_integer_torus(2, 1), gaussian_integer_torus(1, 2))
    gaussian13 = (gaussian_integer_torus(3, 2), gaussian_integer_torus(2, 3))
    controls = [axis_integer_torus(3), diamond_integer_torus(2)]
    pairs = [paired_summary(*gaussian5), paired_summary(*gaussian13)]
    declared_pairs = [
        declared_pair_gate((8, 1), (7, 4), (0, 1, 8, 9)),
        declared_pair_gate((9, 2), (7, 6), (0, 2, 9, 11)),
        declared_pair_gate((11, 3), (9, 7), (0, 2, 9, 123)),
        declared_pair_gate((12, 1), (9, 8), (0, 1, 9, 137)),
        declared_pair_gate((13, 1), (11, 7), (0, 1, 13, 14)),
    ]
    control_summaries = [compact_control(exact_summary(item)) for item in controls]
    basis_checks = [basis_invariance_summary(item) for item in (*gaussian5, *gaussian13, *controls)]
    assert all(item["all_fixed_K_difference_sums_zero"] for item in pairs)
    assert all(item["multiplicities_equal"] and item["all_witness_differences_nonzero"] for item in declared_pairs)
    assert all(item["formula_failures"] == item["incremental_failures"] == 0 for item in pairs + control_summaries)
    assert all(item["joint_histograms_equal"] for item in basis_checks)
    return serialize({
        "schema": "matching-one/paired-motif-microcanonical-certificate/v1",
        "issue": 40,
        "data_class": "complete exact enumeration only",
        "motifs": MOTIF_SIZES,
        "paired_gaussian_certificates": pairs,
        "declared_production_pair_gates": declared_pairs,
        "independent_geometry_controls": control_summaries,
        "unimodular_basis_checks": basis_checks,
        "theorem_checked": "E[T_m|K=k]=c_m*binom(k,r_m)/binom(N,r_m); equal c_m implies E[Z_m|K=k]=0",
        "boundary": (
            "No covariance, variance reduction, fitted control coefficient, production sample, matching/Euler "
            "identity, or >=2x promotion claim is included."
        ),
    })


def render_markdown(artifact: dict[str, Any]) -> str:
    lines = [
        "# Paired motif microcanonical certificate", "",
        "All entries come from complete exact enumeration; no Monte Carlo samples are used.", "",
        "| Gaussian pair | N | masks | K values | equal multiplicities | maximum fixed-K difference sum |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for item in artifact["paired_gaussian_certificates"]:
        maximum = max(item["max_absolute_fixed_K_difference_sum"].values())
        lines.append(
            f"| `{' / '.join(item['pair'])}` | {item['N']} | {item['checked_shared_masks']} | "
            f"{item['checked_K_values']} | `{item['multiplicities_equal']}` | `{maximum}` |"
        )
    lines.extend([
        "", "The N=5 and N=13 conjugate/swapped controls are configurationwise degenerate under the",
        "deterministic quotient labelling (zero nontrivial masks); they certify the exhaustive counter and",
        "fixed-K algebra, not useful control covariance.", "", "Declared nontrivial same-N gates:", "",
    ])
    for item in artifact["declared_production_pair_gates"]:
        lines.append(
            f"- `{' / '.join(item['pair'])}` (N={item['N']}): equal multiplicities; K=4 witness "
            f"differences `{item['witness_differences']}`."
        )
    lines.extend(["", "Independent controls:", ""])
    for item in artifact["independent_geometry_controls"]:
        lines.append(
            f"- `{item['name']}` (N={item['N']}): {item['checked_masks']} masks, "
            f"formula failures `{item['formula_failures']}`, incremental failures `{item['incremental_failures']}`."
        )
    lines.extend([
        "", "All six unimodular-basis joint-histogram checks pass exactly.", "",
        "## Interpretation boundary", "", artifact["boundary"], "",
    ])
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    artifact = build_artifact()
    rendered = json.dumps(artifact, indent=2, sort_keys=True) + "\n" if args.format == "json" else render_markdown(artifact)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
