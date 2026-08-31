#!/usr/bin/env python3
"""Compile a bounded typed matrix model into canonical sparse rational equations."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT / "results/model-certificates/synthetic/m2d-vs-m2j/latest.json"
DEFAULT_OUTPUT = ROOT / "results/model-certificates/framework/algebraic-problem/latest.json"
SCHEMA = "matching-one/canonical-algebraic-model-problem/v1"
Polynomial = dict[tuple[int, ...], Fraction]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _constant(value: Any, count: int) -> Polynomial:
    coefficient = Fraction(value)
    return {} if coefficient == 0 else {(0,) * count: coefficient}


def _variable(index: int, count: int) -> Polynomial:
    exponents = [0] * count
    exponents[index] = 1
    return {tuple(exponents): Fraction(1)}


def _add(left: Polynomial, right: Polynomial) -> Polynomial:
    result = dict(left)
    for monomial, coefficient in right.items():
        result[monomial] = result.get(monomial, Fraction()) + coefficient
        if result[monomial] == 0:
            del result[monomial]
    return result


def _scale(value: Any, polynomial: Polynomial) -> Polynomial:
    factor = Fraction(value)
    return {monomial: factor * coefficient for monomial, coefficient in polynomial.items() if factor * coefficient}


def _multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for left_monomial, left_coefficient in left.items():
        for right_monomial, right_coefficient in right.items():
            monomial = tuple(a + b for a, b in zip(left_monomial, right_monomial))
            result[monomial] = result.get(monomial, Fraction()) + left_coefficient * right_coefficient
            if result[monomial] == 0:
                del result[monomial]
    return result


def _poly_matrix_multiply(left: Sequence[Sequence[Polynomial]], right: Sequence[Sequence[Polynomial]]) -> list[list[Polynomial]]:
    _require(left and right and len(left[0]) == len(right), "polynomial matrix product shape mismatch")
    result: list[list[Polynomial]] = []
    for row in range(len(left)):
        output_row = []
        for column in range(len(right[0])):
            value: Polynomial = {}
            for inner in range(len(right)):
                value = _add(value, _multiply(left[row][inner], right[inner][column]))
            output_row.append(value)
        result.append(output_row)
    return result


def _identity(dimension: int, count: int) -> list[list[Polynomial]]:
    return [[_constant(row == column, count) for column in range(dimension)] for row in range(dimension)]


def _serialize(polynomial: Polynomial) -> list[dict[str, Any]]:
    return [
        {"coefficient": str(polynomial[monomial]), "exponents": list(monomial)}
        for monomial in sorted(polynomial)
    ]


def evaluate_equation(equation: Mapping[str, Any], variable_order: Sequence[str], assignment: Mapping[str, Any]) -> Fraction:
    _require(set(assignment) == set(variable_order), "assignment variables drift")
    values = [Fraction(assignment[name]) for name in variable_order]
    result = Fraction()
    for term in equation["terms"]:
        value = Fraction(term["coefficient"])
        for base, exponent in zip(values, term["exponents"]):
            value *= base ** exponent
        result += value
    return result


def compile_problem(descriptor: Mapping[str, Any]) -> Mapping[str, Any]:
    _require(set(descriptor) == {"dimension", "matrix_symbols", "relations"}, "descriptor fields drift")
    dimension = descriptor["dimension"]
    symbols = descriptor["matrix_symbols"]
    relations = descriptor["relations"]
    _require(isinstance(dimension, int) and dimension >= 1, "invalid dimension")
    _require(isinstance(symbols, list) and symbols and len(symbols) == len(set(symbols)), "matrix symbols must be unique")
    _require(all(isinstance(name, str) and name.isidentifier() for name in symbols), "invalid matrix symbol")
    variable_order = [f"{name}_{row}_{column}" for name in symbols for row in range(dimension) for column in range(dimension)]
    variable_index = {name: index for index, name in enumerate(variable_order)}
    count = len(variable_order)
    matrices = {
        name: [[_variable(variable_index[f"{name}_{row}_{column}"], count) for column in range(dimension)] for row in range(dimension)]
        for name in symbols
    }
    equations: list[dict[str, Any]] = []

    def append_matrix(label: str, matrix: Sequence[Sequence[Polynomial]]) -> None:
        for row in range(len(matrix)):
            for column in range(len(matrix[row])):
                if matrix[row][column]:
                    equations.append({"label": f"{label}[{row},{column}]", "terms": _serialize(matrix[row][column])})

    for relation_index, relation in enumerate(relations):
        relation_type = relation.get("type")
        label = relation.get("id", f"relation-{relation_index}")
        if relation_type == "matrix_affine_identity":
            _require(set(relation) == {"id", "type", "left", "identity_coefficient", "right", "right_coefficient"}, "affine relation fields drift")
            _require(relation["left"] in matrices and relation["right"] in matrices, "unknown matrix symbol")
            residual = []
            for row in range(dimension):
                residual_row = []
                for column in range(dimension):
                    value = matrices[relation["left"]][row][column]
                    value = _add(value, _scale(-Fraction(relation["identity_coefficient"]), _constant(row == column, count)))
                    value = _add(value, _scale(-Fraction(relation["right_coefficient"]), matrices[relation["right"]][row][column]))
                    residual_row.append(value)
                residual.append(residual_row)
            append_matrix(label, residual)
        elif relation_type == "matrix_product_zero":
            _require(set(relation) == {"id", "type", "factors"}, "product relation fields drift")
            factors = relation["factors"]
            _require(isinstance(factors, list) and factors and all(name in matrices for name in factors), "unknown product factor")
            product = _identity(dimension, count)
            for name in factors:
                product = _poly_matrix_multiply(product, matrices[name])
            append_matrix(label, product)
        elif relation_type == "typed_observation":
            _require(set(relation) == {"id", "type", "readout", "word", "source", "expected"}, "observation relation fields drift")
            readout = [Fraction(value) for value in relation["readout"]]
            source = [Fraction(value) for value in relation["source"]]
            word = relation["word"]
            _require(len(readout) == dimension and len(source) == dimension, "observation dimension mismatch")
            _require(isinstance(word, list) and all(name in matrices for name in word), "unknown observation word")
            product = _identity(dimension, count)
            for name in word:
                product = _poly_matrix_multiply(product, matrices[name])
            value: Polynomial = _constant(-Fraction(relation["expected"]), count)
            for row in range(dimension):
                for column in range(dimension):
                    value = _add(value, _scale(readout[row] * source[column], product[row][column]))
            if value:
                equations.append({"label": label, "terms": _serialize(value)})
        else:
            raise ValueError(f"unsupported relation type: {relation_type}")
    max_degree = max((sum(term["exponents"]) for equation in equations for term in equation["terms"]), default=0)
    return {
        "coefficient_field": "Q",
        "variable_order": variable_order,
        "equations": equations,
        "equation_count": len(equations),
        "maximum_total_degree": max_degree,
        "canonical_order": "relation order, matrix row-major, monomial exponent lexicographic",
        "status": "canonical_sparse_rational_problem_compiled",
    }


def descriptor_from_source(source: Mapping[str, Any]) -> dict[str, Any]:
    realization = source["m2j_extracted_realization"]
    moments = source["synthetic_input"]["moments"][:3]
    relations: list[dict[str, Any]] = [
        {"id": "shared_nilpotent_form", "type": "matrix_affine_identity", "left": "A", "identity_coefficient": "1", "right": "K", "right_coefficient": "1"},
        {"id": "nilpotent_square", "type": "matrix_product_zero", "factors": ["K", "K"]},
    ]
    for exponent, expected in enumerate(moments):
        relations.append({
            "id": f"moment_{exponent}",
            "type": "typed_observation",
            "readout": realization["readout"],
            "word": ["A"] * exponent,
            "source": realization["source"],
            "expected": expected,
        })
    return {"dimension": 2, "matrix_symbols": ["A", "K"], "relations": relations}


def _solution_assignment(source: Mapping[str, Any]) -> dict[str, str]:
    realization = source["m2j_extracted_realization"]
    return {
        f"{name}_{row}_{column}": matrix[row][column]
        for name, matrix in (("A", realization["generator"]), ("K", realization["nilpotent"]))
        for row in range(2)
        for column in range(2)
    }


def build_result(source_path: Path = DEFAULT_SOURCE) -> dict[str, Any]:
    source = json.loads(source_path.read_text(encoding="utf-8"))
    descriptor = descriptor_from_source(source)
    problem = compile_problem(descriptor)
    assignment = _solution_assignment(source)
    residuals = [evaluate_equation(equation, problem["variable_order"], assignment) for equation in problem["equations"]]
    _require(all(value == 0 for value in residuals), "frozen supplied realization does not satisfy compiled problem")
    return {
        "schema": SCHEMA,
        "issue": 370,
        "claim_level": "exact_compilation",
        "source": {
            "path": str(source_path.relative_to(ROOT)),
            "sha256": _sha256_file(source_path),
            "dependency_group": source["synthetic_input"]["dependency_group"],
        },
        "descriptor": descriptor,
        "problem": problem,
        "supplied_solution_check": {"residuals": [str(value) for value in residuals], "all_zero": True},
        "claim_boundary": {
            "included": "deterministic sparse-rational compilation of affine matrix, product-zero, and typed observation equalities for one bounded model descriptor",
            "excluded": "inequalities, parameter bounds, gauge-atlas construction, SOS relaxation generation, witness search, SDP solving, noisy data, or physical model validation",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any], source_path: Path = DEFAULT_SOURCE) -> Mapping[str, Any]:
    expected = build_result(source_path)
    _require(result == expected, "algebraic problem artifact does not exactly reproduce")
    return {
        "schema": result["schema"],
        "status": "valid_canonical_algebraic_model_problem",
        "variable_count": len(result["problem"]["variable_order"]),
        "equation_count": result["problem"]["equation_count"],
        "maximum_total_degree": result["problem"]["maximum_total_degree"],
        "source_sha256": result["source"]["sha256"],
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        value = json.loads(args.validate.read_text(encoding="utf-8"))
        print(json.dumps(validate_result(value, args.source), indent=2, sort_keys=True))
        return 0
    rendered = json.dumps(build_result(args.source), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
