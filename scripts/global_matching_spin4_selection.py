#!/usr/bin/env python3
"""Exact singlet/[2] selection rule and a one-seam spin-4 discriminator."""

from __future__ import annotations

import argparse
import importlib.util
import json
from fractions import Fraction
from itertools import combinations
from pathlib import Path


HERE = Path(__file__).resolve().parent
PROJECTOR_PATH = HERE / "p262_confluent_potts_projectors.py"
SPEC = importlib.util.spec_from_file_location("p262_projectors", PROJECTOR_PATH)
P262 = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(P262)

Matrix = list[list[Fraction]]
Vector = list[Fraction]


def fraction_record(value: Fraction) -> dict[str, int | str]:
    return {"numerator": value.numerator, "denominator": value.denominator, "text": str(value)}


def row_times_matrix(row: Vector, matrix: Matrix) -> Vector:
    return [sum(row[i] * matrix[i][j] for i in range(len(row))) for j in range(len(matrix[0]))]


def matrix_times_vector(matrix: Matrix, vector: Vector) -> Vector:
    return [sum(matrix[i][j] * vector[j] for j in range(len(vector))) for i in range(len(matrix))]


def dot(left: Vector, right: Vector) -> Fraction:
    return sum(a * b for a, b in zip(left, right))


def pair_permutation_matrix(q: int, permutation: tuple[int, ...]) -> Matrix:
    """Permutation action on unordered colour pairs; columns are input states."""
    if sorted(permutation) != list(range(q)):
        raise ValueError("permutation must contain each colour exactly once")
    pairs = list(combinations(range(q), 2))
    index = {pair: i for i, pair in enumerate(pairs)}
    matrix = [[Fraction(0) for _ in pairs] for _ in pairs]
    for column, pair in enumerate(pairs):
        image = tuple(sorted((permutation[pair[0]], permutation[pair[1]])))
        matrix[index[image]][column] = Fraction(1)
    return matrix


def projector_character(projector: Matrix, group_action: Matrix) -> Fraction:
    return P262.trace(P262.multiply(projector, group_action))


def integer_selection_check(q: int) -> dict:
    projectors = P262.unordered_pair_projectors(q)
    size = q * (q - 1) // 2
    invariant = [Fraction(1) for _ in range(size)]
    seed = [Fraction(i == 0) for i in range(size)]
    charged = matrix_times_vector(projectors["two_row_2"], seed)
    return {
        "Q": q,
        "invariant_covector_matrix_elements": {
            name: [fraction_record(value) for value in row_times_matrix(invariant, projector)]
            for name, projector in projectors.items()
        },
        "exact_linear_null": all(value == 0 for value in row_times_matrix(invariant, projectors["two_row_2"])),
        "charged_positive_control": {
            "definition": "v=P_[2] e_{01}",
            "invariant_one_point": fraction_record(dot(invariant, charged)),
            "two_point_norm": fraction_record(dot(charged, charged)),
        },
    }


def s4_transposition_oracle() -> dict:
    q = 4
    projectors = P262.unordered_pair_projectors(q)
    identity_action = pair_permutation_matrix(q, tuple(range(q)))
    transposition_action = pair_permutation_matrix(q, (1, 0, 2, 3))
    rows = {}
    for name, projector in projectors.items():
        dimension = projector_character(projector, identity_action)
        twisted = projector_character(projector, transposition_action)
        rows[name] = {
            "S4_irrep": {"singlet": "[4]", "standard": "[3,1]", "two_row_2": "[2,2]"}[name],
            "dimension": fraction_record(dimension),
            "transposition_character": fraction_record(twisted),
            "twist_to_identity_trace_ratio": fraction_record(twisted / dimension),
        }
    return rows


def analyze() -> dict:
    return {
        "schema_version": 1,
        "issue": 257,
        "scope": "exact regular one-insertion selection rule plus frozen discriminator; no MC run",
        "selection_rule": {
            "exact_statement": (
                "For integer Q>=4, an unlabeled S_Q-invariant covector on the unordered-pair carrier "
                "annihilates the [1] and [2] projectors exactly and preserves only the singlet projector. "
                "Any regular analytic continuation of this one-insertion matrix element is therefore zero at Q=1."
            ),
            "observable_hypothesis": (
                "The ordinary fixed-Q global matching observable is an unlabeled zero-leg scalar before its C4 "
                "orientation projection: it introduces neither an FK cluster-colour character nor defect endpoints."
            ),
            "channel_table": [
                {
                    "candidate": "V_(2,+/-2)",
                    "dimension": "17/4",
                    "spin": "+/-4",
                    "S_Q_sector": "[2]",
                    "legs": 4,
                    "regular_global_matching_linear_coupling": "exact_zero",
                    "failed_gates": ["S_Q singlet", "zero-leg/topological endpoint"],
                },
                {
                    "candidate": "thermal Q4 epsilon",
                    "dimension": "21/4",
                    "spin": "+/-4",
                    "S_Q_sector": "singlet",
                    "legs": 0,
                    "regular_global_matching_linear_coupling": "allowed_not_proved_nonzero",
                    "passed_gates": ["S_Q singlet", "zero-leg", "C4 spin"],
                },
            ],
            "integer_oracles": [integer_selection_check(q) for q in (4, 5, 6)],
        },
        "claim_boundary": {
            "minimal_counterexample": (
                "An invariant torus trace need not exclude [2]: tr(P_[2])=Q(Q-3)/2, whose analytic value at Q=1 is -1. "
                "Thus symmetry proves the ordinary linear matrix-element null, not absence from every intermediate-state trace."
            ),
            "categorical_trace_at_Q1": fraction_record(Fraction(-1)),
            "live_loopholes": [
                "an explicit cluster-colour/character insertion",
                "a twisted or defect-sector trace",
                "a singular Q-derivative/confluent residue such as the #262 J direction",
                "failure of the assumed regular generic-Q continuation of the global matching endpoint",
            ],
            "positive_control": "P_[2] e_{01} has zero invariant one-point but a nonzero invariant quadratic norm.",
        },
        "one_shot_discriminator": {
            "status": "frozen_not_run",
            "model": "critical Q=4 square-lattice Potts/FK torus",
            "comparison": "identity seam versus one colour-transposition seam",
            "quantity": (
                "unnormalized orientation-resolved transfer/numerator amplitude; restore the partition-function "
                "factor if simulation records normalized expectations"
            ),
            "target_formula": "A_g/A_e=chi_lambda(g)/dim(lambda)",
            "exact_character_oracle": s4_transposition_oracle(),
            "primary_targets": {
                "thermal_Q4_singlet": fraction_record(Fraction(1)),
                "V22_[2,2]": fraction_record(Fraction(0)),
            },
            "decision": "a single handed H4 amplitude ratio separates 1 from 0; no exponent fit is needed",
        },
        "evidence_layers": {
            "exact": [
                "pair-space projector matrix elements",
                "charged positive control",
                "S4 transposition characters",
                "categorical-trace counterexample to the overbroad claim",
            ],
            "derived_under_explicit_hypothesis": [
                "global matching is a regular unlabeled zero-leg generic-Q endpoint",
                "therefore its ordinary linear overlap with V_(2,+/-2) vanishes",
            ],
            "conjectural": [
                "the observed global matching spin-4 amplitude is dominated by thermal Q4 epsilon",
            ],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = json.dumps(analyze(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
