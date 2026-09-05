#!/usr/bin/env python3
"""Exact witnesses for the finite-noise Jordan/ordinary-mode boundary.

These are synthetic finite-word realizations, not fits to Matching One data.
The general radial-semigroup proof is in the companion note. No numerical
matrix diagonalization or floating-point arithmetic enters the verifier.
"""
from __future__ import annotations

import argparse
import itertools
import json
from fractions import Fraction as F
from pathlib import Path
from typing import Sequence

Matrix = tuple[tuple[F, F], tuple[F, F]]
IDENTITY: Matrix = ((F(1), F(0)), (F(0), F(1)))


def matrix(rows: Sequence[Sequence[object]]) -> Matrix:
    if len(rows) != 2 or any(len(row) != 2 for row in rows):
        raise ValueError("matrix must have exactly two rows and two columns")
    if any(isinstance(x, (float, bool)) for row in rows for x in row):
        raise TypeError("certificate entries must be exact rational strings or integers")
    return tuple(tuple(F(x) for x in row) for row in rows)  # type: ignore[return-value]


def multiply(a: Matrix, b: Matrix) -> Matrix:
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))  # type: ignore[return-value]


def difference(a: Matrix, b: Matrix) -> Matrix:
    return tuple(tuple(a[i][j] - b[i][j] for j in range(2))
                 for i in range(2))  # type: ignore[return-value]


def norm_inf(a: Matrix) -> F:
    return max(sum(abs(x) for x in row) for row in a)


def serialize(a: Matrix) -> list[list[str]]:
    return [[str(x) for x in row] for row in a]


def transfer(cocycle: F, epsilon: F) -> Matrix:
    return ((F(1), cocycle), (cocycle * epsilon ** 2, F(1)))


def words(alphabet: Sequence[str], depth: int) -> list[tuple[str, ...]]:
    if not isinstance(depth, int) or isinstance(depth, bool) or depth < 0:
        raise ValueError("word depth must be a nonnegative integer")
    return [word for size in range(depth + 1)
            for word in itertools.product(alphabet, repeat=size)]


def word_matrix(word: Sequence[str], generators: dict[str, Matrix]) -> Matrix:
    result = IDENTITY
    for letter in word:
        result = multiply(result, generators[letter])
    return result


def word_error_bound(word: Sequence[str], cocycles: dict[str, F]) -> F:
    """C_w with ||A_w(eps)-A_w(0)||_inf <= eps^2 C_w for |eps|<=1.

    Telescoping the product, each changed factor has norm eps^2 |lambda|;
    every other factor has norm at most 1+|lambda|. This bound is rational.
    """
    result = F(0)
    for j, letter in enumerate(word):
        term = abs(cocycles[letter])
        for k, other in enumerate(word):
            if k != j:
                term *= 1 + abs(cocycles[other])
        result += term
    return result


def build_certificate() -> dict:
    cocycles = {"a": F(1), "b": F(3)}
    depth = 5
    contexts = words(tuple(cocycles), depth)
    jordan = {letter: transfer(value, F(0)) for letter, value in cocycles.items()}
    witnesses = []
    for exponent in range(1, 7):
        epsilon = F(1, 10 ** exponent)
        ordinary = {letter: transfer(value, epsilon)
                    for letter, value in cocycles.items()}
        diagonalizer = ((F(1), F(1)), (epsilon, -epsilon))
        inverse = ((F(1, 2), 1 / (2 * epsilon)),
                   (F(1, 2), -1 / (2 * epsilon)))
        max_error = max(norm_inf(difference(word_matrix(word, ordinary),
                                           word_matrix(word, jordan)))
                        for word in contexts)
        max_bound = max(epsilon ** 2 * word_error_bound(word, cocycles)
                        for word in contexts)
        witnesses.append({
            "epsilon": str(epsilon),
            "ordinary_generators": {letter: serialize(value)
                                    for letter, value in ordinary.items()},
            "common_diagonalizer": serialize(diagonalizer),
            "diagonalizer_inverse": serialize(inverse),
            "eigenvalues": {letter: [str(1 + value * epsilon),
                                      str(1 - value * epsilon)]
                            for letter, value in cocycles.items()},
            "max_context_error_inf": str(max_error),
            "max_telescoping_error_bound": str(max_bound),
            "diagonalizer_condition_inf": str(norm_inf(diagonalizer) * norm_inf(inverse)),
        })
    return {
        "schema": "matching-one/jordan-nonseparation/v1",
        "baseline_commit": "4d70c1787ff97dbb98cb5e96022f947bb8fad97e",
        "related_issues": [218, 249, 370],
        "kind": "exact_synthetic_counterfamily_not_a_data_fit",
        "cocycles": {key: str(value) for key, value in cocycles.items()},
        "cocycle_scope": "rational finite-word control; 1 and 3 are NOT log(2) and log(5)",
        "word_order": "A_(uv)=A_u A_v",
        "max_word_length": depth,
        "context_count_per_witness": len(contexts),
        "jordan_generators": {letter: serialize(value) for letter, value in jordan.items()},
        "witnesses": witnesses,
        "claim_boundary": {
            "verified_here": "common real diagonalizability, exact commutation, and finite-word approximation bounds",
            "proved_in_note": "closure obstruction also for a common ordinary-power radial semigroup",
            "not_claimed": "new Monte Carlo evidence, a physical Jordan identification, or exclusion of a spectrally separated ordinary model",
        },
    }


def verify_certificate(payload: dict) -> dict:
    """Check stored witnesses by exact identities and re-evaluated contexts.

    Does not call the certificate builder or a numerical optimizer. A missing,
    changed, or false witness raises ValueError rather than returning a score.
    """
    def require(condition: bool, message: str) -> None:
        if not condition:
            raise ValueError(message)

    require(payload["schema"] == "matching-one/jordan-nonseparation/v1", "schema")
    require(payload["cocycles"] == {"a": "1", "b": "3"}, "declared rational control")
    cocycles = {key: F(value) for key, value in payload["cocycles"].items()}
    require(payload["max_word_length"] == 5, "declared word depth")
    contexts = words(tuple(cocycles), payload["max_word_length"])
    require(len(contexts) == payload["context_count_per_witness"], "context count")
    jordan = {key: matrix(value) for key, value in payload["jordan_generators"].items()}
    require(set(jordan) == set(cocycles), "Jordan alphabet")
    for key, value in cocycles.items():
        require(jordan[key] == ((F(1), value), (F(0), F(1))), "Jordan control")
        nilpotent = difference(jordan[key], IDENTITY)
        require(norm_inf(nilpotent) > 0 and norm_inf(multiply(nilpotent, nilpotent)) == 0,
                "nonzero rank-one nilpotent")
    require(len(payload["witnesses"]) == 6, "witness count")
    comparisons = 0
    previous_error = None
    for index, witness in enumerate(payload["witnesses"], start=1):
        epsilon = F(witness["epsilon"])
        require(epsilon == F(1, 10 ** index), "declared epsilon sequence")
        ordinary = {key: matrix(value) for key, value in witness["ordinary_generators"].items()}
        require(set(ordinary) == set(cocycles), "ordinary alphabet")
        s, inverse = matrix(witness["common_diagonalizer"]), matrix(witness["diagonalizer_inverse"])
        require(multiply(s, inverse) == IDENTITY and multiply(inverse, s) == IDENTITY, "inverse")
        require(F(witness["diagonalizer_condition_inf"]) == norm_inf(s) * norm_inf(inverse),
                "condition number")
        for key, value in cocycles.items():
            require(ordinary[key] == ((F(1), value), (value * epsilon ** 2, F(1))),
                    "declared approximating family")
            eigenvalues = tuple(F(x) for x in witness["eigenvalues"][key])
            require(eigenvalues == (1 + value * epsilon, 1 - value * epsilon), "eigenvalues")
            require(eigenvalues[0] != eigenvalues[1] and min(eigenvalues) > 0,
                    "distinct positive real eigenvalues")
            diagonal = ((eigenvalues[0], F(0)), (F(0), eigenvalues[1]))
            require(multiply(multiply(inverse, ordinary[key]), s) == diagonal, "diagonalization")
        require(multiply(ordinary["a"], ordinary["b"]) == multiply(ordinary["b"], ordinary["a"]),
                "common-generator commutation")
        errors, bounds = [], []
        for word in contexts:
            error = norm_inf(difference(word_matrix(word, ordinary), word_matrix(word, jordan)))
            bound = epsilon ** 2 * word_error_bound(word, cocycles)
            require(error <= bound, "finite-word error bound")
            errors.append(error)
            bounds.append(bound)
            comparisons += 1
        require(max(errors) == F(witness["max_context_error_inf"]), "stored context error")
        require(max(bounds) == F(witness["max_telescoping_error_bound"]), "stored bound")
        require(previous_error is None or max(errors) < previous_error, "convergence sequence")
        previous_error = max(errors)
    return {"status": "verified_exactly", "ordinary_witnesses": 6,
            "contexts_per_witness": len(contexts), "context_error_inequalities": comparisons,
            "floating_point_operations": 0}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", type=Path, help="verify a stored JSON certificate")
    args = parser.parse_args()
    output = (verify_certificate(json.loads(args.verify.read_text(encoding="utf-8")))
              if args.verify else build_certificate())
    print(json.dumps(output, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
