#!/usr/bin/env python3
"""Norm-2 H4/H8 production scorer for Issue #156.

The four square-bond geometries comprise two parent/child lineages related by
the active Gaussian multiplier 1+i. Every design receives a distinct random
seed block; parent/child cross-covariances are therefore zero by the declared
independent-stream design, not estimated as common-random-number covariance.
"""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
import json
from math import sqrt
import os
from pathlib import Path
from typing import Iterable, Optional, Sequence

from integer_period_torus import Matrix, integer_torus_geometry
from square_bond_primitive_pilot import (
    BatchResult,
    CATEGORIES,
    _run_batch,
    _splitmix64,
    analyze_design,
    exact_oracle,
)


ROTATE_ONE_PLUS_I: Matrix = ((1, -1), (1, 1))
LINEAGES = (
    {
        "lineage": "pell_Dminus2",
        "parent": ("pell_Dminus2_N30", ((6, 3), (0, 5))),
        "child": ("pell_Dminus2_N60_norm2", ((6, -2), (6, 8))),
    },
    {
        "lineage": "pell_Dplus1",
        "parent": ("pell_Dplus1_N56", ((8, 4), (0, 7))),
        "child": ("pell_Dplus1_N112_norm2", ((8, -3), (8, 11))),
    },
)
PRODUCTION_DESIGNS = tuple(
    (lineage[role][0], lineage[role][1], lineage["lineage"], role)
    for lineage in LINEAGES
    for role in ("parent", "child")
)


def matrix_multiply(left: Matrix, right: Matrix) -> Matrix:
    return (
        (
            left[0][0] * right[0][0] + left[0][1] * right[1][0],
            left[0][0] * right[0][1] + left[0][1] * right[1][1],
        ),
        (
            left[1][0] * right[0][0] + left[1][1] * right[1][0],
            left[1][0] * right[0][1] + left[1][1] * right[1][1],
        ),
    )


def validate_designs() -> None:
    for lineage in LINEAGES:
        parent_matrix = lineage["parent"][1]
        child_matrix = lineage["child"][1]
        if matrix_multiply(ROTATE_ONE_PLUS_I, parent_matrix) != child_matrix:
            raise ValueError(f"invalid norm-2 child in {lineage['lineage']}")
        parent_n = integer_torus_geometry(parent_matrix).n
        child_n = integer_torus_geometry(child_matrix).n
        if child_n != 2 * parent_n:
            raise ValueError(f"child order is not twice parent in {lineage['lineage']}")


def seed_blocks(seed: int) -> dict[str, int]:
    """Return four distinct deterministic seed blocks, one per design."""

    blocks = {
        identifier: _splitmix64(seed ^ _splitmix64(index + 0x156))
        for index, (identifier, _, _, _) in enumerate(PRODUCTION_DESIGNS)
    }
    if len(set(blocks.values())) != len(blocks):
        raise RuntimeError("seed block collision")
    return blocks


def run_production_batches(
    *, samples_per_design: int, batches: int, seed: int, workers: int
) -> tuple[list[BatchResult], dict[str, int]]:
    validate_designs()
    if samples_per_design <= 0 or batches <= 1:
        raise ValueError("samples_per_design must be positive and batches > 1")
    if samples_per_design % batches:
        raise ValueError("samples_per_design must be divisible by batches")
    if workers <= 0:
        raise ValueError("workers must be positive")
    per_batch = samples_per_design // batches
    blocks = seed_blocks(seed)
    tasks = [
        (identifier, matrix, batch, per_batch, blocks[identifier])
        for identifier, matrix, _, _ in PRODUCTION_DESIGNS
        for batch in range(batches)
    ]
    if workers == 1:
        return [_run_batch(task) for task in tasks], blocks
    with ProcessPoolExecutor(max_workers=workers) as executor:
        return list(executor.map(_run_batch, tasks)), blocks


def _block_diagonal(
    parent: Sequence[Sequence[float]], child: Sequence[Sequence[float]]
) -> list[list[float]]:
    return [
        [
            parent[i][j]
            if i < 3 and j < 3
            else child[i - 3][j - 3]
            if i >= 3 and j >= 3
            else 0.0
            for j in range(6)
        ]
        for i in range(6)
    ]


def analyze_lineage(parent: dict, child: dict) -> dict:
    parent_c = parent["contrasts"]["C_nontrivial_real"]["value"]
    child_c = child["contrasts"]["C_nontrivial_real"]["value"]
    parent_variance = parent["contrast_covariance_of_mean"][0][0]
    child_variance = child["contrast_covariance_of_mean"][0][0]
    common_variance = parent_variance + 4 * child_variance
    covariance_h4_h8 = -parent_variance + 4 * child_variance
    values = {
        "H4_null_2Cchild_plus_Cparent": 2 * child_c + parent_c,
        "H8_null_2Cchild_minus_Cparent": 2 * child_c - parent_c,
    }
    standard_error = sqrt(max(0.0, common_variance))
    scores = {
        name: {
            "value": value,
            "standard_error": standard_error,
            "z": value / standard_error if standard_error else None,
        }
        for name, value in values.items()
    }
    return {
        "lineage": parent["design"].rsplit("_N", 1)[0],
        "parent_design": parent["design"],
        "child_design": child["design"],
        "parent_child_sampling": "independent_seed_blocks_no_CRN",
        "C_parent": parent_c,
        "C_child": child_c,
        "C_child_over_C_parent": child_c / parent_c if parent_c else None,
        "C_parent_times_C_child": parent_c * child_c,
        "six_coordinate_order": [
            "C_parent",
            "Q_parent",
            "S_parent",
            "C_child",
            "Q_child",
            "S_child",
        ],
        "six_coordinate_covariance_of_mean": _block_diagonal(
            parent["contrast_covariance_of_mean"],
            child["contrast_covariance_of_mean"],
        ),
        "null_order": [
            "H4_null_2Cchild_plus_Cparent",
            "H8_null_2Cchild_minus_Cparent",
        ],
        "null_scores": scores,
        "null_covariance_of_mean": [
            [common_variance, covariance_h4_h8],
            [covariance_h4_h8, common_variance],
        ],
    }


def build_result(
    rows: Sequence[BatchResult],
    blocks: dict[str, int],
    *,
    samples_per_design: int,
    batches: int,
    seed: int,
    dps: int,
) -> dict:
    analyses = {
        identifier: analyze_design(identifier, matrix, rows, dps=dps)
        for identifier, matrix, _, _ in PRODUCTION_DESIGNS
    }
    designs = []
    for identifier, matrix, lineage, role in PRODUCTION_DESIGNS:
        payload = analyses[identifier]
        payload["lineage"] = lineage
        payload["lineage_role"] = role
        payload["seed_block"] = blocks[identifier]
        designs.append(payload)
    lineage_results = [
        analyze_lineage(
            analyses[lineage["parent"][0]], analyses[lineage["child"][0]]
        )
        for lineage in LINEAGES
    ]
    return {
        "schema": "p156-square-bond-primitive-norm2-v1",
        "issue": 156,
        "model": "square_bond_percolation",
        "p": 0.5,
        "samples_per_design": samples_per_design,
        "batches_per_design": batches,
        "master_seed": seed,
        "workers_affect_results": False,
        "sampling_contract": {
            "seed_blocks": blocks,
            "parent_child_common_random_numbers": False,
            "cross_design_covariance": "zero by independent RNG stream design",
            "batch_labels_are_not_paired_across_designs": True,
        },
        "geometry_contract": {
            "active_multiplier": "1+i",
            "multiplier_matrix": [list(row) for row in ROTATE_ONE_PLUS_I],
            "norm": 2,
            "child_order": "N_child=2*N_parent",
        },
        "exact_oracle": exact_oracle(),
        "designs": designs,
        "lineages": lineage_results,
        "primary_scores": {
            "H4": "2*C_child+C_parent=0",
            "H8": "2*C_child-C_parent=0",
            "radial_assumption": "same leading C~N^-1 transfer for the H4/H8 comparison",
        },
        "interpretation_boundary": [
            "H4 and H8 nulls are separate fixed hypotheses, not fitted signs.",
            "Parent and child use independent seed blocks; no CRN covariance is claimed.",
            "C/Q/S covariance is retained for every design and as a six-coordinate lineage block.",
            "Q remains a reflection/convention null and S remains a scalar diagnostic.",
            "The norm-2 score assumes the live same-radial C~N^-1 comparison.",
        ],
    }


def write_batches(
    path: Path, rows: Sequence[BatchResult], blocks: dict[str, int]
) -> None:
    metadata = {
        identifier: (lineage, role)
        for identifier, _, lineage, role in PRODUCTION_DESIGNS
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ("lineage", "role", "design", "seed_block", "batch", "samples")
            + CATEGORIES
        )
        for row in sorted(rows, key=lambda item: (item.design, item.batch)):
            lineage, role = metadata[row.design]
            writer.writerow(
                (lineage, role, row.design, blocks[row.design], row.batch, row.samples)
                + tuple(row.counts[category] for category in CATEGORIES)
            )


def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples-per-design", type=int, default=2_000_000)
    parser.add_argument("--batches", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20260829)
    parser.add_argument("--workers", type=int, default=max(1, os.cpu_count() or 1))
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args(list(argv) if argv is not None else None)
    if not exact_oracle()["passed"]:
        raise SystemExit("N=4 exact oracle failed")
    rows, blocks = run_production_batches(
        samples_per_design=args.samples_per_design,
        batches=args.batches,
        seed=args.seed,
        workers=args.workers,
    )
    payload = build_result(
        rows,
        blocks,
        samples_per_design=args.samples_per_design,
        batches=args.batches,
        seed=args.seed,
        dps=args.dps,
    )
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    Path(str(args.output_prefix) + ".json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_batches(Path(str(args.output_prefix) + ".batches.csv"), rows, blocks)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
