#!/usr/bin/env python3
"""Exact minimal context-Hankel oracle for Gaussian cover sectors."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

from gaussian_semigroup_design import Gaussian


Q = Fraction
Matrix = list[list[Q]]


def rank(matrix: Matrix) -> int:
    values = [row[:] for row in matrix]
    if not values:
        return 0
    row = 0
    for column in range(len(values[0])):
        pivot = next(
            (candidate for candidate in range(row, len(values)) if values[candidate][column]),
            None,
        )
        if pivot is None:
            continue
        values[row], values[pivot] = values[pivot], values[row]
        scale = values[row][column]
        values[row] = [value / scale for value in values[row]]
        for other in range(len(values)):
            if other == row:
                continue
            scale = values[other][column]
            values[other] = [
                value - scale * pivot_value
                for value, pivot_value in zip(values[other], values[row])
            ]
        row += 1
        if row == len(values):
            break
    return row


def determinant(matrix: Matrix) -> Q:
    if any(len(row) != len(matrix) for row in matrix):
        raise ValueError("determinant requires a square matrix")
    values = [row[:] for row in matrix]
    result = Q(1)
    for column in range(len(values)):
        pivot = next(
            (candidate for candidate in range(column, len(values)) if values[candidate][column]),
            None,
        )
        if pivot is None:
            return Q(0)
        if pivot != column:
            values[column], values[pivot] = values[pivot], values[column]
            result = -result
        scale = values[column][column]
        result *= scale
        values[column] = [value / scale for value in values[column]]
        for other in range(column + 1, len(values)):
            scale = values[other][column]
            values[other] = [
                value - scale * pivot_value
                for value, pivot_value in zip(values[other], values[column])
            ]
    return result


def diagonal(values: tuple[int, ...]) -> Matrix:
    return [
        [Q(value if row == column else 0) for column, value in enumerate(values)]
        for row in range(len(values))
    ]


def multiply(left: Matrix, right: Matrix) -> Matrix:
    return [
        [
            sum(left[row][inner] * right[inner][column] for inner in range(len(right)))
            for column in range(len(right[0]))
        ]
        for row in range(len(left))
    ]


def subtract(left: Matrix, right: Matrix) -> Matrix:
    return [
        [left[row][column] - right[row][column] for column in range(len(left[0]))]
        for row in range(len(left))
    ]


def text_matrix(matrix: Matrix) -> list[list[str]]:
    return [[str(value) for value in row] for row in matrix]


def build_artifact() -> dict[str, object]:
    identity = Gaussian(1, 0)
    a = Gaussian(1, 1)
    aa = a.multiply(a)
    direct = Gaussian(0, 2)
    if aa != direct:
        raise AssertionError("(1+i)^2 must equal 2i")

    words = ("epsilon", "a", "aa", "direct_2i")
    endpoints = (identity, a, aa, direct)
    endpoint_hankel: Matrix = [
        [Q(1) for _ in words],
        [Q(value.norm) for value in endpoints],
    ]
    endpoint_minor = [[endpoint_hankel[row][column] for column in (0, 1)] for row in range(2)]

    # Sector-resolved minimal witness.  The endpoint block is the exact
    # constant/norm character table above; the third entry is one allowed
    # opposite-character marked matrix element.  Selection rules make the
    # off-diagonal blocks exactly zero.
    sector_block: Matrix = [
        [endpoint_minor[0][0], endpoint_minor[0][1], Q(0)],
        [endpoint_minor[1][0], endpoint_minor[1][1], Q(0)],
        [Q(0), Q(0), Q(1)],
    ]

    # A three-state exact witness showing why endpoint composition does not
    # constrain the charged coordinate.  On the first two coordinates the
    # direct degree-four map equals two degree-two steps.  On the marked
    # coordinate it may differ, and endpoint covectors annihilate that defect.
    A_a = diagonal((1, 2, -1))
    A_direct = diagonal((1, 4, 0))
    composition_defect = subtract(A_direct, multiply(A_a, A_a))

    return {
        "schema": "matching-one/context-hankel-exact-oracle/v1",
        "issue": 249,
        "status": "exact_witness_plus_frozen_acquisition_prediction",
        "exact_gaussian_words": {
            "a": "1+i",
            "aa": "(1+i)^2=2i",
            "direct": "2i",
            "same_endpoint": aa == direct,
            "endpoint_norms": dict(zip(words, (value.norm for value in endpoints))),
            "endpoint_smith": {
                word: list(value.smith_invariants()) for word, value in zip(words, endpoints)
            },
        },
        "endpoint_product_sector": {
            "declared_exact_outputs": ["constant_character", "norm_character"],
            "columns": list(words),
            "hankel_block": text_matrix(endpoint_hankel),
            "rank": rank(endpoint_hankel),
            "aa_equals_direct_column": [
                endpoint_hankel[row][2] == endpoint_hankel[row][3] for row in range(2)
            ],
            "rank_two_minor_columns": ["epsilon", "a"],
            "rank_two_minor_determinant": str(determinant(endpoint_minor)),
            "meaning": (
                "Any endpoint-only series factors through the Gaussian product and therefore "
                "identifies aa with direct_2i; this exact block is a rank-two witness, not a fit "
                "to the physical thermal amplitudes."
            ),
        },
        "sector_resolved_enrichment": {
            "rows_and_columns": ["endpoint_constant", "endpoint_norm", "charged_marked"],
            "block_hankel_minor": text_matrix(sector_block),
            "rank": rank(sector_block),
            "determinant": str(determinant(sector_block)),
            "off_sector_zeros": "exact deck-character selection rule",
            "interpretation": "rank 3 = rank 2 trivial endpoint sector direct-sum rank 1 charged sector",
        },
        "morphism_sensitive_exact_witness": {
            "A_a": text_matrix(A_a),
            "A_direct_2i": text_matrix(A_direct),
            "A_direct_minus_A_a_squared": text_matrix(composition_defect),
            "defect_rank": rank(composition_defect),
            "endpoint_projection_of_defect": [["0", "0"], ["0", "0"]],
            "claim": (
                "Product composition is exact on the endpoint subspace while one marked/charged "
                "coordinate can retain a rank-one context defect invisible to endpoint covectors."
            ),
            "boundary": "This is an exact identifiable model class, not a measured lattice defect.",
        },
        "frozen_physical_prediction": {
            "endpoint_rank": 2,
            "first_allowed_charged_sector_rank": 1,
            "combined_unstructured_rank": 3,
            "sector_statement": (
                "The third state is predicted to live in a nontrivial deck/marked sector, not to be "
                "a third scalar bulk correction field."
            ),
            "falsifiers": [
                "endpoint held-out contexts reject every rank-two realization",
                "the allowed charged block is zero after exact labeling/null controls",
                "the charged block itself has predictive rank greater than one",
            ],
        },
        "archive_gate": {
            "rank_scored_from_current_archive": False,
            "reason": (
                "Current unmarked N65/N130/N260 endpoints and norm-four pilots do not contain a "
                "common-replica intermediate Z2 charged context for the aa path and its direct-2i comparator."
            ),
            "do_not_do": "Do not fill the missing charged entry from endpoint amplitudes or a fitted latent model.",
        },
        "minimal_acquisition": {
            "source": "one N65 parent with exact common parent/fiber labels",
            "words": ["epsilon", "a=1+i", "aa=(1+i)(1+i)", "direct d=2i"],
            "endpoint_rows": [
                "frozen normalized P4[D] response",
                "frozen normalized P4[Dprime] or the already-declared second endpoint mode",
            ],
            "charged_row": (
                "R_chi=<O_chibar S_chi>/[p(1-p)] with the intermediate Z2 flag for aa and the "
                "corresponding flagged character of the Z2xZ2 direct cover"
            ),
            "required_raw_structure": (
                "same replicas for every available row; delete-one blocks sufficient to rebuild the "
                "whole partial Hankel and its covariance"
            ),
            "first_scores": [
                "endpoint rank-2 held-out context prediction",
                "charged-sector nonzero gate",
                "sector-resolved rank (2,1) versus enriched rank 3",
                "A_direct-A_a^2 residual only inside the charged block",
            ],
            "minimum_new_measurement": (
                "one flagged charged aa/direct-2i pair; unmarked aa and direct-2i are the same endpoint "
                "and must not be counted as independent data"
            ),
        },
        "claim_boundary": {
            "exact": "endpoint product blindness, character-sector block zeros, and ranks of the displayed oracle",
            "prediction": "physical endpoint rank 2, charged rank 1, combined rank 3",
            "not_claimed": "that the current archive has already measured any of these predictive ranks",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = build_artifact()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
