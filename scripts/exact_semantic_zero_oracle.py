#!/usr/bin/env python3
"""Exact synthetic infeasibility witness for a contradictory semantic-zero row."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results/model-certificates/synthetic/semantic-zero-contradiction/latest.json"
SCHEMA = "matching-one.synthetic-model-certificate.semantic-zero.v1"


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")


def _sha256(value: Any) -> str:
    return hashlib.sha256(_canonical_bytes(value)).hexdigest()


def _add_polynomials(*polynomials: Sequence[Fraction]) -> list[Fraction]:
    size = max(map(len, polynomials))
    return [
        sum((poly[index] if index < len(poly) else Fraction() for poly in polynomials), Fraction())
        for index in range(size)
    ]


def _scale_polynomial(scale: Fraction, polynomial: Sequence[Fraction]) -> list[Fraction]:
    return [scale * coefficient for coefficient in polynomial]


def build_result() -> dict[str, Any]:
    semantic_zero = [Fraction(0), Fraction(1)]  # z
    observed_row = [Fraction(-1), Fraction(7)]  # 7z-1
    multipliers = [Fraction(7), Fraction(-1)]
    bezout_identity = _add_polynomials(
        _scale_polynomial(multipliers[0], semantic_zero),
        _scale_polynomial(multipliers[1], observed_row),
    )
    observed_value = Fraction(1, 7)
    if bezout_identity != [Fraction(1), Fraction(0)]:
        raise AssertionError("primitive exact-zero Bezout witness did not equal one")
    if observed_value == 0:
        raise AssertionError("negative-control observation must violate the semantic zero")

    synthetic_input = {
        "row_descriptor": {
            "observable": "synthetic_unmarked_global_response",
            "character": "nontrivial_deck_character",
            "perturbation_order": 1,
            "semantic_constraint": "exact_zero",
        },
        "observed_value": str(observed_value),
        "dependency_group": "synthetic-exact-semantic-zero",
    }
    return {
        "schema": SCHEMA,
        "issue": 370,
        "claim_level": "exact_synthetic_control",
        "synthetic_input": synthetic_input,
        "synthetic_input_sha256": _sha256(synthetic_input),
        "polynomial_system": {
            "variable": "z",
            "semantic_zero_coefficients_constant_first": ["0", "1"],
            "observed_row_coefficients_constant_first": ["-1", "7"],
            "equations": ["z=0", "7*z-1=0"],
        },
        "primitive_bezout_infeasibility_witness": {
            "multipliers": [str(value) for value in multipliers],
            "identity": "7*(z)-1*(7*z-1)=1",
            "result_coefficients_constant_first": [str(value) for value in bezout_identity],
            "gcd_of_integer_multipliers": 1,
            "status": "exactly_infeasible",
        },
        "gate_result": {
            "semantic_zero_enforced_as_hard_equality": True,
            "contradictory_row_caught_before_optimization": True,
            "solver_invoked": False,
        },
        "claim_boundary": {
            "included": "one exact linear ideal-membership contradiction for a typed synthetic zero row",
            "excluded": "derivation of a physical selection rule, noisy near-zero handling, general polynomial infeasibility, SOS completeness, or model validation",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any]) -> Mapping[str, Any]:
    expected = build_result()
    if result != expected:
        raise ValueError("semantic-zero certificate does not exactly reproduce")
    return {
        "schema": result["schema"],
        "status": "valid_exact_synthetic_control",
        "bezout_identity": result["primitive_bezout_infeasibility_witness"]["identity"],
        "contradiction_status": result["primitive_bezout_infeasibility_witness"]["status"],
        "solver_invoked": result["gate_result"]["solver_invoked"],
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        print(json.dumps(validate_result(json.loads(args.validate.read_text())), indent=2, sort_keys=True))
        return 0
    rendered = json.dumps(build_result(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
