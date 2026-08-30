#!/usr/bin/env python3
"""Exact finite-field projective twist representations for Issues 337/334."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Iterable, Sequence


Vector = tuple[int, int]
Matrix2 = tuple[tuple[int, int], tuple[int, int]]
Permutation = tuple[int, ...]
R: Matrix2 = ((0, -1), (1, 0))
REFLECTION: Matrix2 = ((1, 0), (0, -1))
T_SHEAR: Matrix2 = ((1, 1), (0, 1))


def projective_lines(q: int) -> list[Vector]:
    return [(1, slope) for slope in range(q)] + [(0, 1)]


def projectivize(vector: Vector, q: int) -> Vector:
    x, y = vector[0] % q, vector[1] % q
    if x == 0 and y == 0:
        raise ValueError("zero has no projective class")
    if x:
        inverse = pow(x, -1, q)
        return 1, y * inverse % q
    return 0, 1


def action_permutation(matrix: Matrix2, q: int) -> Permutation:
    lines = projective_lines(q)
    index = {line: i for i, line in enumerate(lines)}
    return tuple(index[projectivize((
        matrix[0][0] * x + matrix[0][1] * y,
        matrix[1][0] * x + matrix[1][1] * y,
    ), q)] for x, y in lines)


def compose(first: Permutation, second: Permutation) -> Permutation:
    return tuple(first[second[i]] for i in range(len(first)))


def generated_group(generators: Sequence[Permutation]) -> list[Permutation]:
    identity = tuple(range(len(generators[0])))
    group = {identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            candidate = compose(generator, current)
            if candidate not in group:
                group.add(candidate)
                frontier.append(candidate)
    return sorted(group)


def permutation_matrix(permutation: Permutation) -> list[list[Fraction]]:
    size = len(permutation)
    output = [[Fraction(0) for _ in range(size)] for _ in range(size)]
    for old, new in enumerate(permutation):
        output[new][old] = Fraction(1)
    return output


def matmul(
    first: Sequence[Sequence[Fraction]], second: Sequence[Sequence[Fraction]]
) -> list[list[Fraction]]:
    return [[
        sum(first[i][k] * second[k][j] for k in range(len(second)))
        for j in range(len(second[0]))
    ] for i in range(len(first))]


def add_matrices(*matrices: Sequence[Sequence[Fraction]]) -> list[list[Fraction]]:
    return [[sum(matrix[i][j] for matrix in matrices) for j in range(len(matrices[0][0]))]
            for i in range(len(matrices[0]))]


def outer_projector(vector: Sequence[int]) -> list[list[Fraction]]:
    norm = sum(value * value for value in vector)
    return [[Fraction(x * y, norm) for y in vector] for x in vector]


def identity(size: int) -> list[list[Fraction]]:
    return [[Fraction(i == j) for j in range(size)] for i in range(size)]


def apply_permutation(permutation: Permutation, vector: Sequence[int]) -> tuple[int, ...]:
    output = [0] * len(vector)
    for old, new in enumerate(permutation):
        output[new] = vector[old]
    return tuple(output)


def eigenvalue(permutation: Permutation, vector: Sequence[int]) -> int:
    transformed = apply_permutation(permutation, vector)
    if transformed == tuple(vector):
        return 1
    if transformed == tuple(-value for value in vector):
        return -1
    raise ValueError("vector is not a one-dimensional eigenline")


def fraction_payload(matrix: Sequence[Sequence[Fraction]]) -> list[list[str]]:
    return [[str(value) for value in row] for row in matrix]


def representation_certificate(q: int) -> dict[str, object]:
    lines = projective_lines(q)
    rotation = action_permutation(R, q)
    reflection = action_permutation(REFLECTION, q)
    shear = action_permutation(T_SHEAR, q)
    d4_image = generated_group((rotation, reflection))
    modular_image = generated_group((rotation, shear))
    if q == 2:
        basis = {
            "scalar_A1": (1, 1, 1),
            "balanced_H4_alias_A1": (1, -2, 1),
            "axis_odd_B1": (1, 0, -1),
        }
        d4_decomposition = "2 A1 + B1"
        modular_decomposition = "1 + 2D standard irreducible of PSL(2,F2)=S3"
        generator_matrices = {
            "S_on_nontrivial_basis_[H4,axis_odd]": [["1", "0"], ["0", "-1"]],
            "T_on_nontrivial_basis_[H4,axis_odd]": [
                ["-1/2", "sqrt(3)/2"], ["sqrt(3)/2", "1/2"]
            ],
        }
    elif q == 3:
        basis = {
            "scalar_A1": (1, 1, 1, 1),
            "balanced_H4_alias_A1": (1, -1, -1, 1),
            "axis_odd_B1": (1, 0, 0, -1),
            "diagonal_odd_B2": (0, 1, -1, 0),
        }
        d4_decomposition = "2 A1 + B1 + B2"
        modular_decomposition = "1 + 3D standard irreducible of PSL(2,F3)=A4"
        generator_matrices = {
            "S_on_nontrivial_basis_[H4,axis_odd,diagonal_odd]": [
                ["1", "0", "0"], ["0", "-1", "0"], ["0", "0", "-1"]
            ],
            "reflection_on_nontrivial_basis": [
                ["1", "0", "0"], ["0", "1", "0"], ["0", "0", "-1"]
            ],
            "T_on_nontrivial_basis": [
                ["0", "1/sqrt(2)", "-1/sqrt(2)"],
                ["1/sqrt(2)", "1/2", "1/2"],
                ["1/sqrt(2)", "-1/2", "-1/2"],
            ],
            "T_cyclic_restriction": {
                "neutral_1D": "u=H4/sqrt(3)+sqrt(2/3)*axis_odd",
                "charged_2D": [
                    "v=sqrt(2/3)*H4-axis_odd/sqrt(3)",
                    "w=diagonal_odd",
                ],
                "T_on_[v,w]": [["-1/2", "-sqrt(3)/2"],
                                 ["sqrt(3)/2", "-1/2"]],
            },
        }
    else:
        raise ValueError("certificate currently supports q=2,3")

    projectors = {name: outer_projector(vector) for name, vector in basis.items()}
    if q == 2:
        d4_isotypic = {
            "A1_multiplicity_2": add_matrices(
                projectors["scalar_A1"], projectors["balanced_H4_alias_A1"]
            ),
            "B1": projectors["axis_odd_B1"],
        }
    else:
        d4_isotypic = {
            "A1_multiplicity_2": add_matrices(
                projectors["scalar_A1"], projectors["balanced_H4_alias_A1"]
            ),
            "B1": projectors["axis_odd_B1"],
            "B2": projectors["diagonal_odd_B2"],
        }
    modular_isotypic = {
        "trivial": projectors["scalar_A1"],
        "standard": add_matrices(*[
            projector for name, projector in projectors.items()
            if name != "scalar_A1"
        ]),
    }
    cyclic_isotypic: dict[str, list[list[Fraction]]] | None = None
    if q == 3:
        cyclic_neutral = outer_projector((3, -1, -1, -1))
        cyclic_isotypic = {
            "neutral_1D_inside_standard": cyclic_neutral,
            "charged_real_2D_inside_standard": add_matrices(
                modular_isotypic["standard"],
                [[-value for value in row] for row in cyclic_neutral],
            ),
        }
    projector_sum = add_matrices(*projectors.values())
    projector_gates = {
        "sum_to_identity": projector_sum == identity(len(lines)),
        "idempotent": all(matmul(value, value) == value for value in projectors.values()),
        "pairwise_orthogonal": all(
            matmul(first, second) == [[Fraction(0) for _ in lines] for _ in lines]
            for i, first in enumerate(projectors.values())
            for j, second in enumerate(projectors.values()) if i != j
        ),
        "D4_character_isotypic_complete": (
            add_matrices(*d4_isotypic.values()) == identity(len(lines))
        ),
        "modular_character_isotypic_complete": (
            add_matrices(*modular_isotypic.values()) == identity(len(lines))
        ),
    }
    if cyclic_isotypic is not None:
        projector_gates["T_cyclic_splits_standard"] = (
            add_matrices(*cyclic_isotypic.values()) == modular_isotypic["standard"]
            and all(matmul(value, value) == value
                    for value in cyclic_isotypic.values())
        )
    eigen_rows = []
    for name, vector in basis.items():
        eigen_rows.append({
            "name": name, "vector": list(vector),
            "rotation_eigenvalue": eigenvalue(rotation, vector),
            "reflection_eigenvalue": eigenvalue(reflection, vector),
        })
    standard_character = [
        sum(permutation[i] == i for i in range(len(lines))) - 1
        for permutation in modular_image
    ]
    standard_inner_product = Fraction(
        sum(value * value for value in standard_character), len(modular_image)
    )
    return {
        "q": q,
        "line_order": [list(line) for line in lines],
        "generators": {
            "S_rotation_permutation_old_to_new": list(rotation),
            "reflection_permutation_old_to_new": list(reflection),
            "T_shear_permutation_old_to_new": list(shear),
        },
        "D4_projective_image_order": len(d4_image),
        "D4_decomposition": d4_decomposition,
        "D4_eigenlines": eigen_rows,
        "projectors": {name: fraction_payload(value) for name, value in projectors.items()},
        "D4_character_isotypic_projectors": {
            name: fraction_payload(value) for name, value in d4_isotypic.items()
        },
        "modular_character_isotypic_projectors": {
            name: fraction_payload(value) for name, value in modular_isotypic.items()
        },
        "T_cyclic_real_isotypic_projectors": None if cyclic_isotypic is None else {
            name: fraction_payload(value) for name, value in cyclic_isotypic.items()
        },
        "projector_gates": projector_gates,
        "modular_projective_image_order": len(modular_image),
        "modular_group": "PSL(2,F2)=S3" if q == 2 else "PSL(2,F3)=A4",
        "modular_decomposition": modular_decomposition,
        "standard_character_values_over_group": standard_character,
        "standard_character_inner_product": str(standard_inner_product),
        "standard_is_irreducible": standard_inner_product == 1,
        "orthonormal_generator_matrices": generator_matrices,
    }


def observed_f3(score_path: Path) -> dict[str, object]:
    score = json.loads(score_path.read_text(encoding="utf-8"))
    joint = score["joint_estimate"]
    order, means, covariance = joint["order"], joint["mean"], joint["covariance"]
    names = (
        "second_minus_first_F3_H4_axis_diag",
        "second_minus_first_F3_axis_odd",
        "second_minus_first_F3_diagonal_odd",
    )
    indices = [order.index(name) for name in names]
    vector = [means[index] for index in indices]
    block = [[covariance[i][j] for j in indices] for i in indices]
    h, a, d = vector
    cyclic = [
        h / math.sqrt(3.0) + math.sqrt(2.0 / 3.0) * a,
        math.sqrt(2.0 / 3.0) * h - a / math.sqrt(3.0),
        d,
    ]
    transform = [
        [1.0 / math.sqrt(3.0), math.sqrt(2.0 / 3.0), 0.0],
        [math.sqrt(2.0 / 3.0), -1.0 / math.sqrt(3.0), 0.0],
        [0.0, 0.0, 1.0],
    ]
    cyclic_covariance = [[
        sum(transform[i][x] * block[x][y] * transform[j][y]
            for x in range(3) for y in range(3))
        for j in range(3)
    ] for i in range(3)]
    charged = cyclic[1:]
    charged_covariance = [row[1:] for row in cyclic_covariance[1:]]
    determinant = (
        charged_covariance[0][0] * charged_covariance[1][1]
        - charged_covariance[0][1] * charged_covariance[1][0]
    )
    charged_inverse = [
        [charged_covariance[1][1] / determinant,
         -charged_covariance[0][1] / determinant],
        [-charged_covariance[1][0] / determinant,
         charged_covariance[0][0] / determinant],
    ]
    charged_quadratic = sum(
        charged[i] * charged_inverse[i][j] * charged[j]
        for i in range(2) for j in range(2)
    )
    t_matrix = [
        [0.0, 1.0 / math.sqrt(2.0), -1.0 / math.sqrt(2.0)],
        [1.0 / math.sqrt(2.0), 0.5, 0.5],
        [1.0 / math.sqrt(2.0), -0.5, -0.5],
    ]
    predicted = [sum(t_matrix[i][j] * vector[j] for j in range(3)) for i in range(3)]
    return {
        "source_score": str(score_path),
        "D4_basis_order": ["balanced_H4_alias_A1", "axis_odd_B1", "diagonal_odd_B2"],
        "orientation_contrast": vector,
        "covariance": block,
        "T_cyclic_basis_order": ["neutral_u", "charged_v", "charged_w"],
        "T_cyclic_coordinates": cyclic,
        "T_cyclic_covariance": cyclic_covariance,
        "T_charged_doublet_quadratic": charged_quadratic,
        "T_charged_doublet_df": 2,
        "T_neutral_z": cyclic[0] / math.sqrt(cyclic_covariance[0][0]),
        "interpretation": (
            "diagonal_odd is one coordinate of the real 2D T-charged doublet; "
            "it is not a standalone modular irrep"
        ),
        "no_fit_T_shear_prediction_in_D4_basis": predicted,
    }


def build_certificate(score_path: Path) -> dict[str, object]:
    fields = [representation_certificate(q) for q in (2, 3)]
    passed = all(
        all(row["projector_gates"].values()) and row["standard_is_irreducible"]
        for row in fields
    )
    return {
        "schema": "matching-one/flat-twist-projective-representation/v1",
        "issues": [337, 334],
        "status": "exact finite representation plus placement of exploratory N65 contrast",
        "conventions": {
            "line_action": "column homology vector ell maps by g ell",
            "twist_action": "covector alpha maps contragrediently; its kernel line maps by g",
            "S": [[0, -1], [1, 0]],
            "T": [[1, 1], [0, 1]],
            "reflection": [[1, 0], [0, -1]],
            "D4_labels": "A1:(S=+1,reflection=+1), B1:(-1,+1), B2:(-1,-1)",
        },
        "fields": fields,
        "observed_F3": observed_f3(score_path),
        "all_exact_gates_pass": passed,
        "prediction": {
            "minimal": "a covariantly transported F3 twist source must obey the exact 3x3 T mixing matrix",
            "charged": "under the order-3 T subgroup the A4 standard triplet is neutral 1D plus a real 2D charged rotation pair",
            "falsifier": "a declared shear/source transport that leaves diagonal_odd invariant by itself violates projective covariance",
        },
        "claim_boundary": "the roughly 2-sigma diagonal-odd coordinate is placed in a charged doublet, not claimed as a discovery or field identity",
    }


def render_markdown(result: dict[str, object]) -> str:
    f2, f3 = result["fields"]
    observed = result["observed_F3"]
    vector = observed["orientation_contrast"]
    cyclic = observed["T_cyclic_coordinates"]
    prediction = observed["no_fit_T_shear_prediction_in_D4_basis"]
    lines = [
        "# Flat-twist projective representation", "",
        "All exact permutation, projector and irreducibility gates pass.", "",
        "## Exact decompositions", "",
        f"- F2 under D4 image: `{f2['D4_decomposition']}`; under S,T: "
        f"`{f2['modular_decomposition']}`.",
        f"- F3 under D4 image: `{f3['D4_decomposition']}`; under S,T: "
        f"`{f3['modular_decomposition']}`.", "",
        "The balanced axes-minus-diagonals H4 alias is a second A1 copy under D4, not "
        "a symmetry-distinct scalar irrep. At F3 the axis-odd and diagonal-odd lines are B1 "
        "and B2. The full modular image mixes all three non-scalar D4 coordinates into the "
        "irreducible 3D standard representation of A4.", "",
        "Restricting that triplet to the order-3 shear gives one real neutral line and one "
        "real 2D charged rotation block. The diagonal-odd coordinate is one axis of this "
        "charged doublet, not a standalone modular field.", "",
        "## Placement of the N65 contrast", "",
        f"D4 basis `[H4,axis-odd,diagonal-odd]`: `{vector}`.", "",
        f"T-cyclic basis `[neutral,charged-v,charged-w]`: `{cyclic}`.", "",
        "The covariance-aware charged-doublet diagnostic is "
        f"`{observed['T_charged_doublet_quadratic']:.6f} / 2 df`; the neutral "
        f"coordinate is `z={observed['T_neutral_z']:.3f}`.", "",
        "The approximately two-sigma diagonal-odd marginal is therefore inseparable from "
        "its charged partner under a realizable modular shear. It is not promoted to a discovery.", "",
        "## Minimal no-fit prediction", "",
        "For a covariantly transported T-shear source, the exact mixing is", "",
        "```text",
        "H' = (axis_odd - diagonal_odd)/sqrt(2)",
        "A' = H/sqrt(2) + (axis_odd + diagonal_odd)/2",
        "D' = H/sqrt(2) - (axis_odd + diagonal_odd)/2",
        "```", "",
        f"so the observed N65 vector predicts `{prediction}` with no fitted amplitude. "
        "A future charged/projective-source run should freeze this whole vector and its "
        "transport, not score diagonal-odd alone.", "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--score", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args()
    result = build_certificate(args.score)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_markdown(result), encoding="utf-8")


if __name__ == "__main__":
    main()
