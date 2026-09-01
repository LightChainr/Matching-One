#!/usr/bin/env python3
"""Exact modular-scalar classification of torus homology channels.

The repository represents every wrapping component by an integer winding
subgroup of H_1(T^2, Z).  A modular change of homology generators acts by an
SL(2,Z) matrix on every winding vector.  This oracle separates channel labels
that depend only on the rational rank of that subgroup from labels that depend
on the chosen generators.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable, Sequence


Vector = tuple[int, int]
Matrix = tuple[tuple[int, int], tuple[int, int]]


def determinant(matrix: Matrix) -> int:
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def parse_matrix(payload: Sequence[Sequence[int]]) -> Matrix:
    if len(payload) != 2 or any(len(row) != 2 for row in payload):
        raise ValueError("an SL(2,Z) matrix must be 2 by 2")
    matrix = (
        (int(payload[0][0]), int(payload[0][1])),
        (int(payload[1][0]), int(payload[1][1])),
    )
    if determinant(matrix) != 1:
        raise ValueError(f"matrix is not in SL(2,Z): determinant={determinant(matrix)}")
    return matrix


def apply_matrix(matrix: Matrix, vector: Vector) -> Vector:
    return (
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1],
    )


def subgroup_rank(basis: Sequence[Vector]) -> int:
    nonzero = [vector for vector in basis if vector != (0, 0)]
    if not nonzero:
        return 0
    first = nonzero[0]
    if any(first[0] * other[1] != first[1] * other[0] for other in nonzero[1:]):
        return 2
    return 1


def channel_flags(basis: Sequence[Vector]) -> dict[str, object]:
    """Mirror the exact semantics in scripts/torus_homology.py."""

    rank = subgroup_rank(basis)
    direction_0 = any(x != 0 for x, _ in basis)
    direction_1 = any(y != 0 for _, y in basis)
    return {
        "rank": rank,
        "direction_0": direction_0,
        "direction_1": direction_1,
        "either": rank > 0,
        "both": direction_0 and direction_1,
        "cross": rank == 2,
    }


def transform_basis(matrix: Matrix, basis: Sequence[Vector]) -> tuple[Vector, ...]:
    return tuple(apply_matrix(matrix, vector) for vector in basis)


def changed_channels(before: dict[str, object], after: dict[str, object]) -> list[str]:
    return sorted(key for key in before if before[key] != after[key])


def counterexample(matrix: Matrix, basis: Sequence[Vector]) -> dict[str, object]:
    transformed = transform_basis(matrix, basis)
    before = channel_flags(basis)
    after = channel_flags(transformed)
    return {
        "matrix": [list(row) for row in matrix],
        "determinant": determinant(matrix),
        "basis_before": [list(vector) for vector in basis],
        "basis_after": [list(vector) for vector in transformed],
        "flags_before": before,
        "flags_after": after,
        "changed_channels": changed_channels(before, after),
    }


def allowed_by_stabilizer(spin: int, order: int) -> bool:
    """A homogeneous scalar-channel spin response survives iff exp(2pi*i*s/n)=1."""

    if order <= 0:
        raise ValueError("stabilizer order must be positive")
    return spin % order == 0


def _matrix_records(config: dict) -> list[dict[str, object]]:
    records = []
    for row in config["sl2_transformations"]:
        matrix = parse_matrix(row["matrix"])
        records.append({"id": row["id"], "matrix": matrix})
    return records


def _subgroup_records(config: dict) -> list[dict[str, object]]:
    records = []
    for row in config["representative_subgroups"]:
        basis = tuple((int(v[0]), int(v[1])) for v in row["basis"])
        records.append({"id": row["id"], "basis": basis})
    return records


def _finite_regression(
    transformations: Iterable[dict[str, object]],
    subgroups: Iterable[dict[str, object]],
) -> list[dict[str, object]]:
    rows = []
    for subgroup in subgroups:
        basis = subgroup["basis"]
        before = channel_flags(basis)
        for transformation in transformations:
            matrix = transformation["matrix"]
            after_basis = transform_basis(matrix, basis)
            after = channel_flags(after_basis)
            rows.append({
                "subgroup": subgroup["id"],
                "transformation": transformation["id"],
                "rank_preserved": before["rank"] == after["rank"],
                "either_preserved": before["either"] == after["either"],
                "cross_preserved": before["cross"] == after["cross"],
                "changed_channels": changed_channels(before, after),
            })
    return rows


def analyze(config: dict) -> dict:
    transformations = _matrix_records(config)
    subgroups = _subgroup_records(config)
    regression = _finite_regression(transformations, subgroups)
    if not all(
        row["rank_preserved"] and row["either_preserved"] and row["cross_preserved"]
        for row in regression
    ):
        raise AssertionError("an invertible homology change altered a rank channel")

    basis = ((1, 1),)
    direction_0_example = counterexample(parse_matrix(((1, -1), (0, 1))), basis)
    direction_1_example = counterexample(parse_matrix(((1, 0), (-1, 1))), basis)

    candidates = [int(spin) for spin in config["candidate_spins"]]
    square_order = int(config["elliptic_stabilizers"]["square_tau_i"])
    hexagonal_order = int(config["elliptic_stabilizers"]["hexagonal_rho"])
    spin_filter = [
        {
            "spin": spin,
            "square_tau_i_allowed": allowed_by_stabilizer(spin, square_order),
            "hexagonal_rho_allowed": allowed_by_stabilizer(spin, hexagonal_order),
            "square_lattice_and_hexagonal_allowed": (
                allowed_by_stabilizer(spin, square_order)
                and allowed_by_stabilizer(spin, hexagonal_order)
            ),
        }
        for spin in candidates
    ]

    scalar_channels = {
        "rank": "rank(MW)=rank(W) because M in SL(2,Z) is invertible over Q",
        "either": "either is exactly rank>0",
        "cross": "cross is exactly rank=2",
    }
    non_scalar_channels = {
        "direction_0": "depends on the selected homology generator",
        "direction_1": "depends on the selected homology generator",
        "both": "a rank-1 spiral can use both generators and shear to one generator",
    }
    lifted = {
        base: {
            combination: "modular_scalar"
            for combination in ("primal", "matching", "even", "odd")
        }
        for base in ("either", "cross")
    }

    return {
        "schema_version": 1,
        "issue": 114,
        "theorem": {
            "action": "W -> M W for M in SL(2,Z)",
            "scalar_channels": scalar_channels,
            "basis_dependent_channels": non_scalar_channels,
            "matching_combination_lift": lifted,
            "lift_reason": (
                "occupation complement commutes with geometric relabelling; sums and "
                "differences of two scalar rank channels remain scalar"
            ),
        },
        "counterexamples": {
            "direction_0_and_both": direction_0_example,
            "direction_1_and_both": direction_1_example,
        },
        "finite_regression": {
            "case_count": len(regression),
            "all_rank_channels_preserved": True,
            "cases": regression,
        },
        "elliptic_spin_filter": {
            "premise": "homogeneous first-order spin response in a modular-scalar channel",
            "square_stabilizer_order": square_order,
            "hexagonal_stabilizer_order": hexagonal_order,
            "intersection_period": 12,
            "candidates": spin_filter,
        },
        "scientific_boundary": config["scientific_boundary"],
    }


def render_markdown(result: dict) -> str:
    lines = [
        "# Modular homology-channel oracle",
        "",
        "Exact integer-homology classification; no Monte Carlo or fitted continuum field is used.",
        "",
        "| channel | classification | reason |",
        "|---|---|---|",
    ]
    theorem = result["theorem"]
    for channel, reason in theorem["scalar_channels"].items():
        lines.append(f"| `{channel}` | modular scalar | {reason} |")
    for channel, reason in theorem["basis_dependent_channels"].items():
        lines.append(f"| `{channel}` | basis-dependent | {reason} |")

    lines += [
        "",
        "The scalar classification lifts to primal, matching, even, and odd combinations of `cross` or `either`.",
        "Complement commutes with geometric relabelling, and linear combinations of scalar channels remain scalar.",
        "",
        "## Exact counterexample",
        "",
        "The rank-1 spiral basis `(1,1)` has `both=true`.  The determinant-one shear",
        "`[[1,-1],[0,1]]` maps it to `(0,1)`, so `direction_0` and `both` become false while",
        "`rank/either/cross` stay unchanged.  The transpose shear similarly changes `direction_1`.",
        "",
        "## Scalar-channel elliptic filter",
        "",
        "| spin | tau=i | hexagonal rho | square and hexagonal |",
        "|---:|:---:|:---:|:---:|",
    ]
    for row in result["elliptic_spin_filter"]["candidates"]:
        mark = lambda value: "yes" if value else "no"
        lines.append(
            f"| H{row['spin']} | {mark(row['square_tau_i_allowed'])} | "
            f"{mark(row['hexagonal_rho_allowed'])} | "
            f"{mark(row['square_lattice_and_hexagonal_allowed'])} |"
        )
    lines += [
        "",
        "Thus a homogeneous first-order response in scalar `cross/either` channels kills H4/H8 at the hexagonal point and permits H12.",
        "This conclusion does not apply to a primitive winding character, `both`, or another vector-valued channel.",
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
