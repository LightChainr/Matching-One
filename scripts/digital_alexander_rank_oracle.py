#!/usr/bin/env python3
"""Exact finite-torus oracle for the matching-pair ambient-rank identity."""

from __future__ import annotations

import argparse
import json
import random
from collections import Counter
from itertools import product
from pathlib import Path
from typing import Iterable, Sequence

from c4_self_matching_exact import c4_self_matching_torus
from integer_period_torus import (
    IntegerTorusGeometry,
    axis_integer_torus,
    classify_configuration,
    diamond_integer_torus,
    gaussian_integer_torus,
    integer_torus_geometry,
)


def active_from_mask(mask: int, n: int) -> tuple[bool, ...]:
    return tuple(bool(mask & (1 << vertex)) for vertex in range(n))


def rank_identity_record(
    geometry: IntegerTorusGeometry, active: Sequence[bool]
) -> dict[str, int]:
    black, _ = classify_configuration(geometry, active)
    white, _ = classify_configuration(
        geometry, tuple(not value for value in active), matching=True
    )
    rank_black = black.max_rank
    rank_white = white.max_rank
    q_either = int(black.either) - int(white.either)
    q_cross = int(black.cross) - int(white.cross)
    return {
        "rank_black": rank_black,
        "rank_white": rank_white,
        "q_either": q_either,
        "q_cross": q_cross,
        "common_channel_residual": q_either - q_cross,
        "weak_rank_residual": 2 * q_either - (rank_black - rank_white),
        "strong_rank_sum_residual": rank_black + rank_white - 2,
    }


def elementary_rank_lemma() -> list[dict[str, int | bool]]:
    """Exhaust the nine abstract rank pairs and expose the premise exactly."""

    rows = []
    for rank_black, rank_white in product(range(3), repeat=2):
        q_either = int(rank_black > 0) - int(rank_white > 0)
        q_cross = int(rank_black == 2) - int(rank_white == 2)
        premise = q_either == q_cross
        weak = 2 * q_either == rank_black - rank_white
        rows.append({
            "rank_black": rank_black,
            "rank_white": rank_white,
            "q_either": q_either,
            "q_cross": q_cross,
            "common_channel_premise": premise,
            "weak_identity": weak,
        })
    if any(row["common_channel_premise"] and not row["weak_identity"] for row in rows):
        raise AssertionError("the elementary rank implication failed")
    return rows


def geometry_from_spec(spec: dict) -> IntegerTorusGeometry:
    kind = spec["kind"]
    if kind == "axis":
        return axis_integer_torus(int(spec["L"]))
    if kind == "diamond":
        return diamond_integer_torus(int(spec["L"]))
    if kind == "gaussian":
        return gaussian_integer_torus(int(spec["a"]), int(spec["b"]))
    if kind == "c4_self_matching":
        return c4_self_matching_torus(int(spec["a"]), int(spec["b"]))
    raise ValueError(f"unsupported geometry kind {kind!r}")


def _joint_key(record: dict[str, int]) -> str:
    return (
        f"r{record['rank_black']}_r{record['rank_white']}"
        f"_q{record['q_either']:+d}"
    )


def summarize_records(
    records: Iterable[tuple[int, dict[str, int]]],
    *,
    keep_counterexamples: int = 8,
) -> dict[str, object]:
    joint: Counter[str] = Counter()
    q_counts: Counter[str] = Counter()
    common_failures = []
    weak_failures = []
    strong_failures = []
    raw_rank_equality_count = 0
    q_zero_count = 0
    raw_rank_equality_outside_q_zero_count = 0
    total = 0
    for mask, record in records:
        total += 1
        joint[_joint_key(record)] += 1
        q_counts[f"{record['q_either']:+d}"] += 1
        payload = {"mask": mask, **record}
        raw_equal = (
            record["q_either"]
            == record["rank_black"] - record["rank_white"]
        )
        if raw_equal:
            raw_rank_equality_count += 1
        if record["q_either"] == 0:
            q_zero_count += 1
        if raw_equal and record["q_either"] != 0:
            raw_rank_equality_outside_q_zero_count += 1
        if record["common_channel_residual"] and len(common_failures) < keep_counterexamples:
            common_failures.append(payload)
        if record["weak_rank_residual"] and len(weak_failures) < keep_counterexamples:
            weak_failures.append(payload)
        if record["strong_rank_sum_residual"] and len(strong_failures) < keep_counterexamples:
            strong_failures.append(payload)
    return {
        "configurations": total,
        "joint_rank_q_counts": dict(sorted(joint.items())),
        "q_counts": dict(sorted(q_counts.items())),
        "raw_rank_equality_count": raw_rank_equality_count,
        "q_zero_count": q_zero_count,
        "raw_rank_equality_outside_q_zero_count": (
            raw_rank_equality_outside_q_zero_count
        ),
        "common_channel_failure_count_capped": len(common_failures),
        "weak_rank_failure_count_capped": len(weak_failures),
        "strong_rank_sum_failure_count_capped": len(strong_failures),
        "common_channel_counterexamples": common_failures,
        "weak_rank_counterexamples": weak_failures,
        "strong_rank_sum_counterexamples": strong_failures,
    }


def enumerate_exact(spec: dict) -> dict[str, object]:
    geometry = geometry_from_spec(spec)
    if geometry.n > int(spec.get("max_exact_n", 20)):
        raise ValueError(f"N={geometry.n} exceeds exact enumeration cap")
    records = (
        (mask, rank_identity_record(geometry, active_from_mask(mask, geometry.n)))
        for mask in range(1 << geometry.n)
    )
    summary = summarize_records(records)
    summary.update({
        "id": spec["id"],
        "kind": spec["kind"],
        "N": geometry.n,
        "period_matrix": [list(row) for row in geometry.periods.matrix],
        "matching_edges_equal_primal_edges": (
            geometry.matching_edges == geometry.primal_edges
        ),
    })
    return summary


def _random_matrices(config: dict, rng: random.Random) -> list[tuple[tuple[int, int], tuple[int, int]]]:
    count = int(config["matrix_count"])
    entry_bound = int(config["entry_bound"])
    min_order = int(config["min_order"])
    max_order = int(config["max_order"])
    matrices = set()
    attempts = 0
    while len(matrices) < count:
        attempts += 1
        if attempts > count * 10000:
            raise ValueError("could not generate enough distinct period matrices")
        matrix = (
            (rng.randint(-entry_bound, entry_bound), rng.randint(-entry_bound, entry_bound)),
            (rng.randint(-entry_bound, entry_bound), rng.randint(-entry_bound, entry_bound)),
        )
        determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
        if min_order <= abs(determinant) <= max_order:
            matrices.add(matrix)
    return sorted(matrices, key=lambda m: (abs(m[0][0] * m[1][1] - m[0][1] * m[1][0]), m))


def deterministic_counterexample_search(config: dict) -> dict[str, object]:
    rng = random.Random(int(config["seed"]))
    matrices = _random_matrices(config, rng)
    samples_per_matrix = int(config["samples_per_matrix"])
    weak_counterexamples = []
    common_counterexamples = []
    strong_counterexamples = []
    evaluated = 0
    rank_pair_counts: Counter[str] = Counter()

    for matrix in matrices:
        geometry = integer_torus_geometry(matrix, name="counterexample-search")
        full_mask = (1 << geometry.n) - 1
        masks = {0, full_mask}
        while len(masks) < min(samples_per_matrix, 1 << geometry.n):
            masks.add(rng.getrandbits(geometry.n))
        for mask in sorted(masks):
            record = rank_identity_record(geometry, active_from_mask(mask, geometry.n))
            evaluated += 1
            rank_pair_counts[
                f"r{record['rank_black']}_r{record['rank_white']}"
            ] += 1
            payload = {
                "period_matrix": [list(row) for row in matrix],
                "N": geometry.n,
                "mask": mask,
                **record,
            }
            if record["common_channel_residual"] and len(common_counterexamples) < 8:
                common_counterexamples.append(payload)
            if record["weak_rank_residual"] and len(weak_counterexamples) < 8:
                weak_counterexamples.append(payload)
            if record["strong_rank_sum_residual"] and len(strong_counterexamples) < 8:
                strong_counterexamples.append(payload)

    return {
        "seed": int(config["seed"]),
        "matrix_count": len(matrices),
        "samples_per_matrix": samples_per_matrix,
        "configurations_evaluated": evaluated,
        "order_range": [int(config["min_order"]), int(config["max_order"])],
        "rank_pair_counts": dict(sorted(rank_pair_counts.items())),
        "common_channel_counterexamples": common_counterexamples,
        "weak_rank_counterexamples": weak_counterexamples,
        "strong_rank_sum_counterexamples": strong_counterexamples,
    }


def analyze(config: dict) -> dict[str, object]:
    lemma = elementary_rank_lemma()
    exhaustive = [enumerate_exact(spec) for spec in config["exhaustive_geometries"]]
    search = deterministic_counterexample_search(config["counterexample_search"])
    return {
        "schema_version": 1,
        "issue": 269,
        "elementary_rank_lemma": {
            "rows": lemma,
            "premise_implies_weak_identity": all(
                not row["common_channel_premise"] or row["weak_identity"]
                for row in lemma
            ),
            "premise_rank_pairs": [
                [row["rank_black"], row["rank_white"]]
                for row in lemma
                if row["common_channel_premise"]
            ],
        },
        "exhaustive_geometries": exhaustive,
        "counterexample_search": search,
        "conclusion": {
            "common_channel_holds_on_all_exhaustive": all(
                not row["common_channel_counterexamples"] for row in exhaustive
            ),
            "weak_identity_holds_on_all_exhaustive": all(
                not row["weak_rank_counterexamples"] for row in exhaustive
            ),
            "strong_rank_sum_holds_on_all_exhaustive": all(
                not row["strong_rank_sum_counterexamples"] for row in exhaustive
            ),
            "random_search_found_weak_counterexample": bool(
                search["weak_rank_counterexamples"]
            ),
            "claim_level": "finite_exact_plus_deterministic_counterexample_search_not_proof",
        },
        "scientific_boundary": config["scientific_boundary"],
    }


def render_markdown(result: dict[str, object]) -> str:
    conclusion = result["conclusion"]
    lines = [
        "# Digital Alexander rank oracle",
        "",
        "Finite exact enumeration plus deterministic counterexample search; not a general-topology proof.",
        "",
        "| geometry | N | configurations | joint `(r_b,r_w,q)` counts | weak failures | strong failures |",
        "|---|---:|---:|---|---:|---:|",
    ]
    for row in result["exhaustive_geometries"]:
        joint = ", ".join(f"{key}:{value}" for key, value in row["joint_rank_q_counts"].items())
        lines.append(
            f"| {row['id']} | {row['N']} | {row['configurations']} | `{joint}` | "
            f"{len(row['weak_rank_counterexamples'])} | {len(row['strong_rank_sum_counterexamples'])} |"
        )
    search = result["counterexample_search"]
    lines += [
        "",
        "## Verdict",
        "",
        f"- common either/cross channel premise on every exhaustive configuration: `{conclusion['common_channel_holds_on_all_exhaustive']}`",
        f"- `2q = r_black-r_white` on every exhaustive configuration: `{conclusion['weak_identity_holds_on_all_exhaustive']}`",
        f"- stronger `r_black+r_white=2` on every exhaustive configuration: `{conclusion['strong_rank_sum_holds_on_all_exhaustive']}`",
        "- on every declared exhaustive geometry, the archived equality `q = r_black-r_white` occurs exactly when `q=0`, hence at rank pair `(1,1)`",
        f"- deterministic search: {search['configurations_evaluated']} configurations on {search['matrix_count']} integer-period tori; weak counterexamples: {len(search['weak_rank_counterexamples'])}",
        "",
        "The nine-case rank lemma proves that equality of the `either` and `cross` differences implies the weak rank identity. The finite oracle verifies the premise; it does not replace a digital Alexander-duality proof for arbitrary tori.",
        "",
        "## Boundary",
        "",
    ]
    lines.extend(f"- {item}" for item in result["scientific_boundary"])
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--json-output", type=Path)
    parser.add_argument("--markdown-output", type=Path)
    args = parser.parse_args()
    result = analyze(json.loads(args.manifest.read_text(encoding="utf-8")))
    if args.json_output:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    else:
        print(json.dumps(result, indent=2, sort_keys=True))
    if args.markdown_output:
        args.markdown_output.parent.mkdir(parents=True, exist_ok=True)
        args.markdown_output.write_text(render_markdown(result), encoding="utf-8")


if __name__ == "__main__":
    main()
