#!/usr/bin/env python3
"""Second norm-2 generation for the Issue #156 primitive C3 character."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
import json
from math import exp, sqrt
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
from square_bond_primitive_norm2 import ROTATE_ONE_PLUS_I, matrix_multiply


LINEAGES = (
    {
        "lineage": "pell_Dminus2_generation2",
        "parent": "pell_Dminus2_N60_norm2",
        "parent_matrix": ((6, -2), (6, 8)),
        "child": "pell_Dminus2_N120_norm2_generation2",
        "child_matrix": ((0, -10), (12, 6)),
    },
    {
        "lineage": "pell_Dplus1_generation2",
        "parent": "pell_Dplus1_N112_norm2",
        "parent_matrix": ((8, -3), (8, 11)),
        "child": "pell_Dplus1_N224_norm2_generation2",
        "child_matrix": ((0, -14), (16, 8)),
    },
)
MODELS = (
    ("rank4_H4", -0.5),
    ("even_nonlocal_character", 0.5),
    ("quadratic_H4", 0.25),
    ("local_H8_bound_saturated", 0.125),
)


def validate_designs() -> None:
    for lineage in LINEAGES:
        expected = matrix_multiply(ROTATE_ONE_PLUS_I, lineage["parent_matrix"])
        if expected != lineage["child_matrix"]:
            raise ValueError(f"bad child matrix for {lineage['lineage']}")
        parent_n = integer_torus_geometry(lineage["parent_matrix"]).n
        child_n = integer_torus_geometry(lineage["child_matrix"]).n
        if child_n != 2 * parent_n:
            raise ValueError(f"bad norm for {lineage['lineage']}")


def run_batches(
    *, samples_per_design: int, batches: int, seed: int, workers: int
) -> tuple[list[BatchResult], dict[str, int]]:
    validate_designs()
    if samples_per_design <= 0 or batches <= 1 or samples_per_design % batches:
        raise ValueError("samples_per_design must be positive and divisible by batches>1")
    if workers <= 0:
        raise ValueError("workers must be positive")
    per_batch = samples_per_design // batches
    blocks = {
        lineage["child"]: _splitmix64(seed ^ _splitmix64(index + 0x2156))
        for index, lineage in enumerate(LINEAGES)
    }
    tasks = [
        (lineage["child"], lineage["child_matrix"], batch, per_batch,
         blocks[lineage["child"]])
        for lineage in LINEAGES
        for batch in range(batches)
    ]
    if workers == 1:
        return [_run_batch(task) for task in tasks], blocks
    with ProcessPoolExecutor(max_workers=workers) as executor:
        return list(executor.map(_run_batch, tasks)), blocks


def _design_by_name(payload: dict, name: str) -> dict:
    matches = [row for row in payload["designs"] if row["design"] == name]
    if len(matches) != 1:
        raise ValueError(f"expected one parent design {name}, found {len(matches)}")
    return matches[0]


def score_model(parent: dict, child: dict, ratio: float) -> dict:
    parent_c = parent["contrasts"]["C_nontrivial_real"]["value"]
    child_c = child["contrasts"]["C_nontrivial_real"]["value"]
    parent_var = parent["contrast_covariance_of_mean"][0][0]
    child_var = child["contrast_covariance_of_mean"][0][0]
    residual = child_c - ratio * parent_c
    variance = child_var + ratio * ratio * parent_var
    standard_error = sqrt(max(0.0, variance))
    return {
        "fixed_child_over_parent_ratio": ratio,
        "residual_Cchild_minus_ratio_Cparent": residual,
        "standard_error": standard_error,
        "z": residual / standard_error if standard_error else None,
    }


def build_result(
    rows: Sequence[BatchResult],
    blocks: dict[str, int],
    parent_payload: dict,
    *,
    samples_per_design: int,
    batches: int,
    seed: int,
    dps: int,
) -> dict:
    children = {
        lineage["child"]: analyze_design(
            lineage["child"], lineage["child_matrix"], rows, dps=dps
        )
        for lineage in LINEAGES
    }
    lineage_results = []
    joint = {name: {"chi_square": 0.0, "df": 2} for name, _ in MODELS}
    for lineage in LINEAGES:
        parent = _design_by_name(parent_payload, lineage["parent"])
        child = children[lineage["child"]]
        scores = {}
        for name, ratio in MODELS:
            score = score_model(parent, child, ratio)
            scores[name] = score
            joint[name]["chi_square"] += score["z"] ** 2
        lineage_results.append(
            {
                "lineage": lineage["lineage"],
                "parent_design": lineage["parent"],
                "child_design": lineage["child"],
                "parent_source": "P156 norm2 generation-1 result.json",
                "parent_child_sampling": "independent_generation_seed_blocks",
                "C_parent": parent["contrasts"]["C_nontrivial_real"]["value"],
                "C_child": child["contrasts"]["C_nontrivial_real"]["value"],
                "C_child_over_C_parent": (
                    child["contrasts"]["C_nontrivial_real"]["value"]
                    / parent["contrasts"]["C_nontrivial_real"]["value"]
                ),
                "model_scores": scores,
            }
        )
    for model in joint.values():
        model["p_value_df2"] = exp(-model["chi_square"] / 2)
    return {
        "schema": "p156-square-bond-primitive-norm2-generation2-v1",
        "issue": 156,
        "model": "square_bond_percolation",
        "p": 0.5,
        "samples_per_child_design": samples_per_design,
        "batches_per_child_design": batches,
        "master_seed": seed,
        "sampling_contract": {
            "child_seed_blocks": blocks,
            "generation1_generation2_common_random_numbers": False,
            "cross_generation_covariance": "zero by independent RNG domains",
        },
        "geometry_contract": {
            "active_multiplier": "1+i",
            "multiplier_matrix": [list(row) for row in ROTATE_ONE_PLUS_I],
            "norm": 2,
        },
        "exact_oracle": exact_oracle(),
        "children": [children[lineage["child"]] for lineage in LINEAGES],
        "lineages": lineage_results,
        "joint_fixed_model_scores": joint,
        "interpretation_boundary": [
            "The fixed ratios were frozen before generation-2 targets were produced.",
            "The generation-1 parents and generation-2 children use independent RNG domains.",
            "The primary question is convergence toward the rank-4 H4 ratio -1/2.",
            "No exponent is fitted before the fixed-model scores are reported.",
        ],
    }


def write_batches(path: Path, rows: Sequence[BatchResult], blocks: dict[str, int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("design", "seed_block", "batch", "samples") + CATEGORIES)
        for row in sorted(rows, key=lambda item: (item.design, item.batch)):
            writer.writerow(
                (row.design, blocks[row.design], row.batch, row.samples)
                + tuple(row.counts[category] for category in CATEGORIES)
            )


def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent-result", type=Path, required=True)
    parser.add_argument("--samples-per-design", type=int, default=5_000_000)
    parser.add_argument("--batches", type=int, default=250)
    parser.add_argument("--seed", type=int, default=202608292)
    parser.add_argument("--workers", type=int, default=max(1, os.cpu_count() or 1))
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args(list(argv) if argv is not None else None)
    parent_payload = json.loads(args.parent_result.read_text(encoding="utf-8"))
    rows, blocks = run_batches(
        samples_per_design=args.samples_per_design,
        batches=args.batches,
        seed=args.seed,
        workers=args.workers,
    )
    payload = build_result(
        rows,
        blocks,
        parent_payload,
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
