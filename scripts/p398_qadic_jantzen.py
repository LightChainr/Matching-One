#!/usr/bin/env python3
"""Exact width-three/four Q-adic Gram filtration and mark projections."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import comb
from pathlib import Path
from typing import Any, Sequence

from noncrossing_connectivity_codec import noncrossing_states
from p321_homology_trace_certificate import action_matrix, rotate_state
from p333_charge2_landing import landing_reference_state as charge2_reference_state
from p333_gram_source_intertwiner import (
    encode_fraction,
    join_block_count,
    matrix_residual_rank,
    multiply,
    rref_solve,
    subtract,
    transpose,
)
from p333_minimal_multimark_jet import falling_factorial
from p333_source_landing_doublet import (
    landing_pair_state,
    landing_rotation as c3_rotation,
)
from p333_source_landing_doublet_width4 import (
    landing_reference_state as c4_reference_state,
    landing_rotation as c4_rotation,
)


SCHEMA = "matching-one/p398-qadic-jantzen/v1"
ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "analysis/p398_qadic_jantzen_protocol.json"
DEFAULT_JSON = ROOT / "results/p398-qadic-jantzen/latest.json"


def identity(size: int) -> list[list[Fraction]]:
    return [
        [Fraction(int(row == column)) for column in range(size)]
        for row in range(size)
    ]


def encode_vector(vector: Sequence[Fraction]) -> list[int | str]:
    return [encode_fraction(value) for value in vector]


def encode_matrix(matrix: Sequence[Sequence[Fraction]]) -> list[list[int | str]]:
    return [encode_vector(row) for row in matrix]


def exact_determinant(matrix: Sequence[Sequence[int | Fraction]]) -> Fraction:
    work = [[Fraction(value) for value in row] for row in matrix]
    size = len(work)
    determinant = Fraction(1)
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column]), None
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            determinant = -determinant
        value = work[column][column]
        determinant *= value
        work[column] = [entry / value for entry in work[column]]
        for row in range(column + 1, size):
            scale = work[row][column]
            if not scale:
                continue
            work[row] = [
                entry - scale * pivot_entry
                for entry, pivot_entry in zip(work[row], work[column])
            ]
    return determinant


def radical_basis(size: int) -> list[list[Fraction]]:
    """Columns e_i-e_last, a canonical basis of ker(all-ones)."""

    return [
        [Fraction(int(row == column) - int(row == size - 1)) for column in range(size - 1)]
        for row in range(size)
    ]


def gram_coefficient(states: Sequence[Sequence[int]], order: int) -> list[list[Fraction]]:
    """Coefficient of t^order in G(1+t)."""

    return [
        [
            Fraction(comb(join_block_count(left, right), order))
            if join_block_count(left, right) >= order
            else Fraction(0)
            for right in states
        ]
        for left in states
    ]


def solve_unique(matrix: Sequence[Sequence[Fraction]], target: Sequence[Fraction]) -> list[Fraction]:
    solution = rref_solve(matrix, target, len(matrix))
    if not solution["consistent"] or solution["dimension"] != 0:
        raise AssertionError("leading radical form must give a unique H-dual")
    return solution["particular"]


def response_columns(width: int) -> dict[str, dict[str, Any]]:
    states = noncrossing_states(width)
    scalar_orders = (1, 2) if width == 3 else (1, 2, 3)
    families: dict[str, dict[str, Any]] = {
        "scalar_falling_factorial": {
            "sector": "C3 trivial" if width == 3 else "C4 trivial",
            "provenance": "081a5ed/e7e6c80",
            "columns": [
                [Fraction(falling_factorial(len(set(state)), order)) for order in scalar_orders]
                for state in states
            ],
            "mark_action": identity(len(scalar_orders)),
            "orders": list(scalar_orders),
        }
    }
    if width == 3:
        references = [landing_pair_state(site) for site in range(width)]
        raw = [
            [join_block_count(state, reference) for reference in references]
            for state in states
        ]
        families["C3_charge1_landing"] = {
            "sector": "C3 charge-one rational doublet",
            "provenance": "b82e8cc",
            "columns": [
                [Fraction(row[0] - row[2]), Fraction(row[1] - row[2])]
                for row in raw
            ],
            "mark_action": [[Fraction(value) for value in row] for row in c3_rotation()],
        }
    else:
        references = [c4_reference_state(site) for site in range(width)]
        charge2_references = [charge2_reference_state(site) for site in range(width)]
        if references != charge2_references:
            raise AssertionError("charge-one and charge-two landing registries differ")
        raw = [
            [join_block_count(state, reference) for reference in references]
            for state in states
        ]
        families["C4_charge1_landing"] = {
            "sector": "C4 charge-one rational doublet",
            "provenance": "7b40ec7",
            "columns": [
                [Fraction(row[0] - row[2]), Fraction(row[1] - row[3])]
                for row in raw
            ],
            "mark_action": [[Fraction(value) for value in row] for row in c4_rotation()],
        }
        families["C4_charge2_landing"] = {
            "sector": "C4 charge-two sign character",
            "provenance": "ab3eed8",
            "columns": [
                [Fraction(row[0] - row[1] + row[2] - row[3])]
                for row in raw
            ],
            "mark_action": [[Fraction(-1)]],
        }
    return families


def translation_sector_dimensions(width: int, traces: Sequence[int]) -> dict[str, int]:
    if width == 3:
        trivial_full = sum(traces) // 3
        charge1_copies = (traces[0] - traces[1]) // 3
        return {
            "trivial": trivial_full - 1,
            "charge1_rational": 2 * charge1_copies,
        }
    trivial_full = sum(traces) // 4
    charge2 = (traces[0] - traces[1] + traces[2] - traces[3]) // 4
    charge1_copies = (traces[0] - traces[2]) // 4
    return {
        "trivial": trivial_full - 1,
        "charge1_rational": 2 * charge1_copies,
        "charge2": charge2,
    }


def family_projection(
    family: dict[str, Any],
    basis: Sequence[Sequence[Fraction]],
    leading_form: Sequence[Sequence[Fraction]],
    translation: Sequence[Sequence[Fraction]],
) -> dict[str, Any]:
    columns = family["columns"]
    mark_action = family["mark_action"]
    functional = multiply(transpose(basis), columns)
    dual_columns = [solve_unique(leading_form, column) for column in zip(*functional)]
    dual_coordinates = [list(row) for row in zip(*dual_columns)]
    full_vectors = multiply(basis, dual_coordinates)
    direct_covariance = subtract(
        multiply(transpose(translation), multiply(columns, mark_action)), columns
    )
    projected_covariance = subtract(
        multiply(translation, full_vectors), multiply(full_vectors, mark_action)
    )
    output = {
        "sector": family["sector"],
        "provenance": family["provenance"],
        "mark_dimension": len(mark_action),
        "grade_one_functional_rank": matrix_residual_rank(functional),
        "grade_one_functional_B_coordinates": encode_matrix(functional),
        "grade_one_H_dual_B_coordinates": encode_matrix(dual_coordinates),
        "translation_covariance": {
            "raw_response_residual_rank": matrix_residual_rank(direct_covariance),
            "H_dual_state_residual_rank": matrix_residual_rank(projected_covariance),
        },
    }
    if "orders" in family:
        output["derivative_orders"] = family["orders"]
    return output


def width_result(width: int) -> dict[str, Any]:
    states = noncrossing_states(width)
    size = len(states)
    g0 = gram_coefficient(states, 0)
    g1 = gram_coefficient(states, 1)
    basis = radical_basis(size)
    leading_form = multiply(transpose(basis), multiply(g1, basis))
    leading_rank = matrix_residual_rank(leading_form)
    leading_det = exact_determinant(leading_form)
    if matrix_residual_rank(g0) != 1:
        raise AssertionError("Q=1 Gram is not the expected all-ones rank-one form")
    if leading_rank != size - 1:
        raise AssertionError("first-order stop failed; higher t-order is required")

    translation = [
        [Fraction(value) for value in row]
        for row in action_matrix(width, lambda state: rotate_state(state, 1))
    ]
    traces = [
        sum(rotate_state(state, power) == state for state in states)
        for power in range(width)
    ]
    sectors = translation_sector_dimensions(width, traces)
    projections = {
        name: family_projection(family, basis, leading_form, translation)
        for name, family in response_columns(width).items()
    }
    functional_columns: list[list[Fraction]] = []
    for family in response_columns(width).values():
        restricted = multiply(transpose(basis), family["columns"])
        functional_columns.extend([list(column) for column in zip(*restricted)])
    combined_rank = matrix_residual_rank(functional_columns)
    if width == 3:
        sector_coverage = {
            "trivial": {
                "layer_dimension": sectors["trivial"],
                "tested_projection_rank": projections["scalar_falling_factorial"]["grade_one_functional_rank"],
            },
            "charge1_rational": {
                "layer_dimension": sectors["charge1_rational"],
                "tested_projection_rank": projections["C3_charge1_landing"]["grade_one_functional_rank"],
            },
        }
    else:
        sector_coverage = {
            "trivial": {
                "layer_dimension": sectors["trivial"],
                "tested_projection_rank": projections["scalar_falling_factorial"]["grade_one_functional_rank"],
            },
            "charge1_rational": {
                "layer_dimension": sectors["charge1_rational"],
                "tested_projection_rank": projections["C4_charge1_landing"]["grade_one_functional_rank"],
            },
            "charge2": {
                "layer_dimension": sectors["charge2"],
                "tested_projection_rank": projections["C4_charge2_landing"]["grade_one_functional_rank"],
            },
        }
    for entry in sector_coverage.values():
        entry["uncovered_dimension"] = entry["layer_dimension"] - entry["tested_projection_rank"]

    return {
        "width": width,
        "basis_dimension": size,
        "exact_polynomial_degree_available": width,
        "G0_rank": 1,
        "radical_dimension": size - 1,
        "leading_radical_form": {
            "rank": leading_rank,
            "determinant": encode_fraction(leading_det),
            "unimodular": abs(leading_det) == 1,
        },
        "local_invariant_factor_valuations": [0] + [1] * (size - 1),
        "jantzen_filtration_dimensions": {"J0": size, "J1": size - 1, "J2": 0},
        "associated_graded_dimensions": {"grade_0": 1, "grade_1": size - 1},
        "highest_nonzero_layer": 1,
        "higher_coefficients_needed": False,
        "proof": (
            "rank G(1)=1 gives exactly dim(V)-1 positive local Smith valuations; "
            "the nonsingular leading radical form forces their sum to dim(V)-1, "
            "so every positive valuation is exactly one and J2=0."
        ),
        "translation_traces_on_V": traces,
        "grade_one_sector_dimensions": sectors,
        "mark_projections": projections,
        "sector_coverage": sector_coverage,
        "combined_tested_projection_rank": combined_rank,
        "combined_uncovered_dimension": size - 1 - combined_rank,
    }


def build_result() -> dict[str, Any]:
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    if protocol["status"] != "frozen_before_exact_filtration":
        raise AssertionError("protocol is not frozen")
    widths = [width_result(width) for width in (3, 4)]
    return {
        "schema": SCHEMA,
        "status": "exact_q_adic_filtration_certificate",
        "issues": [398, 333],
        "protocol": {
            "path": str(PROTOCOL.relative_to(ROOT)),
            "sha256": hashlib.sha256(PROTOCOL.read_bytes()).hexdigest(),
            "commit": "916ea59",
        },
        "widths": widths,
        "global_decision": (
            "The unmarked width-three and width-four Gram degenerations have Loewy/Jantzen "
            "length two only: all positive local invariant-factor valuations equal one and "
            "J2 vanishes. Width three's tested scalar plus C3-charge-one responses span the "
            "entire grade-one layer. At width four the tested scalar, C4-charge-one and "
            "C4-charge-two responses cover only 6 of 13 grade-one dimensions, leaving exact "
            "deficits 2, 2 and 3 in the trivial, charge-one and charge-two sectors. The missing "
            "datum is therefore multiplicity/rooted-connectivity within existing sectors, not "
            "a deeper Q-adic order or another terminal C4 irrep."
        ),
        "nilpotent_separation": {
            "base_parameter_nilpotent": {
                "construction": "multiplication by t on V tensor Q[t]/(t^2)",
                "square_zero": True,
                "rank_before_specialization": {str(row["width"]): row["basis_dimension"] for row in widths},
                "rank_after_Q_equals_1_specialization": 0,
                "automatic_for_every_differentiable_family": True,
            },
            "fixed_Q_extension": {
                "construction": "a separate marked/rooted module after t acts trivially",
                "not_implied_by_nonzero_Jantzen_projection": True,
                "existing_countercontrols": [
                    "width-three scalar responses project with full trivial-sector rank but their frozen scalar endpoint extensions fail",
                    "all three width-four tested sector families project nontrivially but their individual fixed-Q closure tests fail",
                    "the width-three C3 charge-one doublet succeeds only after its separate affine/endpoint/Gram/source test"
                ],
            },
        },
        "claim_boundary": [
            "Exact local-Q Gram algebra for the repository's declared unmarked connectivity modules at widths three and four only.",
            "The associated-graded projections classify already-tested response covectors; they do not construct a fixed-Q extension.",
            "J2=0 rules out a deeper base-Q filtration in these modules, not nonsemisimple rooted, direct-sum or physical transfer modules.",
            "No continuum LCFT field, periodic-TL cell-module dictionary or physical Jordan block is identified."
        ],
    }


def render_markdown(result: dict[str, Any]) -> str:
    lines = [
        "# P398 canonical Q-adic/Jantzen control at widths three and four",
        "",
        result["global_decision"],
        "",
        "| width | dim V | valuations | dim gr0 | dim gr1 | dim J2 | det leading radical form | tested rank | uncovered |",
        "|---:|---:|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in result["widths"]:
        counts = {value: row["local_invariant_factor_valuations"].count(value) for value in (0, 1)}
        lines.append(
            f"| {row['width']} | {row['basis_dimension']} | `0^{counts[0]},1^{counts[1]}` | "
            f"{row['associated_graded_dimensions']['grade_0']} | {row['associated_graded_dimensions']['grade_1']} | "
            f"{row['jantzen_filtration_dimensions']['J2']} | {row['leading_radical_form']['determinant']} | "
            f"{row['combined_tested_projection_rank']} | {row['combined_uncovered_dimension']} |"
        )
    lines += ["", "## Associated-grade coverage", ""]
    for row in result["widths"]:
        lines.append(f"### Width {row['width']}")
        lines.append("")
        lines.append("| sector | gr1 dimension | tested projection rank | uncovered |")
        lines.append("|---|---:|---:|---:|")
        for sector, entry in row["sector_coverage"].items():
            lines.append(
                f"| {sector} | {entry['layer_dimension']} | {entry['tested_projection_rank']} | {entry['uncovered_dimension']} |"
            )
        lines.append("")
    lines += [
        "## Interpretation",
        "",
        "- The full polynomial Gram family is available through degree `width`, but the exact unimodular leading radical form proves that first order already identifies every nonzero layer.",
        "- At width three the tested responses cover the complete grade-one representation. At width four every C4 sector still has uncovered multiplicity, so changing only the terminal character cannot close the gap.",
        "- The automatic dual-number `t` action is square-zero before specialization and zero after `Q=1`; it is a base-parameter nilpotent, not a fixed-Q marked or physical extension.",
        "- Nonzero associated-grade projection is necessary descriptive information but not a closure certificate: the existing width-four fixed-Q tests fail despite nonzero projection in every tested sector.",
        "",
        "## Boundary",
        "",
    ]
    lines += [f"- {item}" for item in result["claim_boundary"]]
    lines.append("")
    return "\n".join(lines)


def render_card(result: dict[str, Any]) -> str:
    width3, width4 = result["widths"]
    return "\n".join(
        [
            "# Scientific card: canonical Q-adic/Jantzen control",
            "",
            "- **Mechanism space changed:** the base connectivity Gram degeneration has no hidden higher Q-adic layer at widths three or four: the valuations are `0^1,1^4` and `0^1,1^13`, with `J2=0` exactly.",
            f"- **Width-three projection:** the frozen scalar and C3 charge-one responses span all `{width3['combined_tested_projection_rank']}/4` grade-one dimensions.",
            f"- **Width-four projection:** the frozen scalar/charge-one/charge-two responses span only `{width4['combined_tested_projection_rank']}/13`; uncovered sector dimensions are trivial `2`, charge one `2`, charge two `3`.",
            "- **New discriminator:** the missing width-four information is multiplicity/rooted-connectivity inside already present C4 sectors, not another terminal irrep or a deeper power of `Q-1`.",
            "- **Nilpotent type:** the square-zero dual-number action is base-parameter bookkeeping and vanishes on fixed-Q specialization; a physical or fixed-Q extension still needs its own module and closure certificate.",
            "- **Not proved:** no LCFT field, physical transfer Jordan block or periodic-TL cell-module identification.",
            "",
        ]
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--markdown", type=Path, required=True)
    parser.add_argument("--card", type=Path, required=True)
    args = parser.parse_args(argv)
    result = build_result()
    for path in (args.json, args.markdown, args.card):
        path.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.markdown.write_text(render_markdown(result), encoding="utf-8")
    args.card.write_text(render_card(result), encoding="utf-8")
    print(
        json.dumps(
            {
                f"width_{row['width']}": {
                    "valuations": row["local_invariant_factor_valuations"],
                    "tested_rank": row["combined_tested_projection_rank"],
                }
                for row in result["widths"]
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
