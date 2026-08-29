#!/usr/bin/env python3
"""Exact spin-4 Q velocities and a gated two-size covariance scorer.

At Q=1, d/d(log Q)=d/dQ.  The scorer deliberately accepts only a total
generic-Q field derivative.  A fixed-observable FK measure-score response is
an input term, not a substitute for the explicit Q derivative of the field
projector/normalization.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Dict, Iterable, Mapping, Sequence


REQUIRED_DERIVATIVE_SEMANTICS = (
    "total_dQ_field=measure_score_response+explicit_field_definition_derivative"
)


def exact_targets() -> Dict[str, object]:
    sqrt3_over_pi = math.sqrt(3.0) / math.pi
    fields = {
        "four_leg_V_2_2": {
            "x": Fraction(17, 4),
            "spin": -4,
            "absolute_spin": 4,
            "leg_count": 4,
            "dx_du": Fraction(-15, 8),
            "velocity_coefficient_sqrt3_over_pi": Fraction(-5, 16),
            "generic_Q_family": "V_(2,2) dense/FK loop primary",
            "potts_multiplicity": "unresolved",
        },
        "thermal_Q4_epsilon": {
            "x": Fraction(21, 4),
            "spin": 4,
            "absolute_spin": 4,
            "leg_count": 0,
            "dx_du": Fraction(-27, 8),
            "velocity_coefficient_sqrt3_over_pi": Fraction(-9, 16),
            "generic_Q_family": "level-4 chiral descendant of Potts energy",
            "potts_multiplicity": "thermal family known; lattice Q4 overlap unresolved",
        },
    }
    rendered = {}
    for identifier, row in fields.items():
        coefficient = row["velocity_coefficient_sqrt3_over_pi"]
        rendered[identifier] = {
            key: ({"text": str(value), "numerator": value.numerator,
                   "denominator": value.denominator, "decimal": float(value)}
                  if isinstance(value, Fraction) else value)
            for key, value in row.items()
        }
        rendered[identifier]["dx_dQ_at_Q1"] = {
            "exact": f"{coefficient}*sqrt(3)/pi",
            "coefficient_sqrt3_over_pi": str(coefficient),
            "decimal": float(coefficient) * sqrt3_over_pi,
        }
    gap = Fraction(1, 4)
    return {
        "potts_coordinate": {
            "u": "beta^2",
            "sqrt_Q": "-2*cos(pi*u)",
            "Q1_u": "2/3",
            "dQ_du_at_Q1": "2*pi*sqrt(3)",
            "d_dt_equals_d_dQ_at_Q1": True,
            "t": "log(Q)",
        },
        "fields": rendered,
        "velocity_gap_four_leg_minus_thermal_Q4": {
            "exact": "sqrt(3)/(4*pi)",
            "coefficient_sqrt3_over_pi": str(gap),
            "decimal": float(gap) * sqrt3_over_pi,
        },
        "dimension_gap_thermal_minus_four_leg": "1",
    }


def solve(matrix: Sequence[Sequence[float]], rhs: Sequence[float]) -> list[float]:
    size = len(rhs)
    a = [list(map(float, matrix[row])) + [float(rhs[row])] for row in range(size)]
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(a[row][column]))
        if abs(a[pivot][column]) < 1e-22:
            raise ArithmeticError("singular covariance")
        a[column], a[pivot] = a[pivot], a[column]
        scale = a[column][column]
        a[column] = [value / scale for value in a[column]]
        for row in range(size):
            if row == column:
                continue
            factor = a[row][column]
            a[row] = [a[row][j] - factor * a[column][j] for j in range(size + 1)]
    return [a[row][-1] for row in range(size)]


def validate_covariance(covariance: Sequence[Sequence[float]]) -> None:
    if len(covariance) != 4 or any(len(row) != 4 for row in covariance):
        raise ValueError("velocity scorer requires a full 4x4 covariance")
    for i in range(4):
        if covariance[i][i] <= 0.0:
            raise ValueError("covariance diagonal must be positive")
        for j in range(4):
            if abs(covariance[i][j] - covariance[j][i]) > 1e-12:
                raise ValueError("covariance must be symmetric")
    # A solve is a compact positive-definite/nonsingular gate for this 4x4 use.
    solve(covariance, [1.0, 0.0, 0.0, 0.0])


def velocity_estimate(sizes: Sequence[float], point: Sequence[float],
                      covariance: Sequence[Sequence[float]]) -> Dict[str, object]:
    """Score [O1,dQO1,O2,dQO2] with the full covariance."""
    if len(sizes) != 2 or sizes[0] <= 0 or sizes[1] <= sizes[0]:
        raise ValueError("sizes must be two increasing positive lengths")
    if len(point) != 4:
        raise ValueError("point order must be [O1,dQO1,O2,dQO2]")
    validate_covariance(covariance)
    observable1, derivative1, observable2, derivative2 = map(float, point)
    if observable1 == 0.0 or observable2 == 0.0:
        raise ValueError("log-response velocity requires nonzero field means")
    log_ratio = math.log(float(sizes[1]) / float(sizes[0]))
    response1 = derivative1 / observable1
    response2 = derivative2 / observable2
    velocity = -(response2 - response1) / log_ratio
    gradient = [
        -derivative1 / (observable1 * observable1 * log_ratio),
        1.0 / (observable1 * log_ratio),
        derivative2 / (observable2 * observable2 * log_ratio),
        -1.0 / (observable2 * log_ratio),
    ]
    variance = sum(gradient[i] * covariance[i][j] * gradient[j]
                   for i in range(4) for j in range(4))
    if variance <= 0.0:
        raise ValueError("propagated velocity variance must be positive")
    return {
        "size_order": list(sizes),
        "point_order": ["O_L1", "dQ_O_L1", "O_L2", "dQ_O_L2"],
        "log_response": [response1, response2],
        "velocity": velocity,
        "gradient": gradient,
        "variance": variance,
        "standard_error": math.sqrt(variance),
    }


def score_targets(estimate: Mapping[str, object], targets: Mapping[str, object]) -> Dict[str, object]:
    result = {}
    velocity = float(estimate["velocity"])
    standard_error = float(estimate["standard_error"])
    for identifier, field in targets["fields"].items():
        target = float(field["dx_dQ_at_Q1"]["decimal"])
        result[identifier] = {
            "target": target,
            "residual": velocity - target,
            "signed_z": (velocity - target) / standard_error,
        }
    return result


def score_input(payload: Mapping[str, object], targets: Mapping[str, object]) -> Dict[str, object]:
    if payload.get("schema") != "matching-one.q-velocity-two-size-input.v1":
        raise ValueError("unexpected Q-velocity input schema")
    if payload.get("derivative_semantics") != REQUIRED_DERIVATIVE_SEMANTICS:
        return {
            "status": "NOT_SCOREABLE",
            "reason": "input is not the total derivative of a declared generic-Q field",
            "required_derivative_semantics": REQUIRED_DERIVATIVE_SEMANTICS,
        }
    if not payload.get("explicit_field_definition_derivative_included"):
        return {
            "status": "NOT_SCOREABLE",
            "reason": "measure score is present but explicit field-definition derivative is absent",
            "required_derivative_semantics": REQUIRED_DERIVATIVE_SEMANTICS,
        }
    estimate = velocity_estimate(payload["sizes"], payload["point"], payload["covariance"])
    return {
        "status": "SCORED_TOTAL_FIELD_DERIVATIVE",
        "field_family": payload.get("field_family", "unspecified"),
        "estimate": estimate,
        "fixed_target_scores": score_targets(estimate, targets),
        "evidence_rule": "the two sizes and both response components are one covariance block",
    }


def synthetic_oracle(targets: Mapping[str, object]) -> Dict[str, object]:
    sizes = [8.0, 16.0]
    y = Fraction(9, 4)
    amplitude = 2.75
    normalization_velocity = -0.43
    covariance = [
        [2.0e-8, 0.3e-8, 0.2e-8, -0.1e-8],
        [0.3e-8, 4.0e-8, 0.1e-8, 0.2e-8],
        [0.2e-8, 0.1e-8, 3.0e-8, 0.4e-8],
        [-0.1e-8, 0.2e-8, 0.4e-8, 5.0e-8],
    ]
    records = {}
    for identifier, field in targets["fields"].items():
        velocity = float(field["dx_dQ_at_Q1"]["decimal"])
        point = []
        for length in sizes:
            observable = amplitude * length ** (-float(y))
            derivative = observable * (normalization_velocity - velocity * math.log(length))
            point.extend([observable, derivative])
        estimate = velocity_estimate(sizes, point, covariance)
        records[identifier] = {
            "synthetic_model": "A(Q)*L^[-y(Q)] with arbitrary A(1), dQ log A",
            "point": point,
            "recovered_velocity": estimate["velocity"],
            "target_velocity": velocity,
            "absolute_error": abs(estimate["velocity"] - velocity),
            "normalization_velocity_cancelled": normalization_velocity,
        }
    return records


def archive_audit(root: Path) -> Dict[str, object]:
    fk = root / "results/fk-q-score/latest.json"
    local_site = root / "results/local-20260829/P225-multiradius-pivotal/raw/n130_n170_20k.batches.csv"
    return {
        "scoreable_now": False,
        "fk_q_score_archive": {
            "path": str(fk.relative_to(root)),
            "present": fk.exists(),
            "sizes": [2] if fk.exists() else [],
            "observable_semantics": "fixed topology functions" if fk.exists() else None,
            "contains_measure_score": fk.exists(),
            "contains_explicit_field_derivative": False,
        },
        "site_percolation_score_archive": {
            "path": str(local_site.relative_to(root)),
            "present": local_site.exists(),
            "coordinate": "Bernoulli site t/lambda; not Potts critical-manifold Q",
            "usable_for_Q_velocity": False,
        },
        "missing": [
            "two sizes of one declared generic-Q spin-4 field family",
            "explicit derivative of its Q-dependent projector/normalization",
            "full covariance of [O1,dQO1,O2,dQO2]",
        ],
    }


def render(root: Path, input_path: Path = None) -> Dict[str, object]:
    targets = exact_targets()
    result = {
        "schema": "matching-one.q1-spin4-velocity-oracle.v1",
        "issue": 261,
        "exact_targets": targets,
        "archive_audit": archive_audit(root),
        "synthetic_two_size_oracle": synthetic_oracle(targets),
        "scorer_contract": {
            "input_schema": "matching-one.q-velocity-two-size-input.v1",
            "point_order": ["O_L1", "dQ_O_L1", "O_L2", "dQ_O_L2"],
            "required_derivative_semantics": REQUIRED_DERIVATIVE_SEMANTICS,
            "formula": "-[(dQO2/O2)-(dQO1/O1)]/log(L2/L1)",
            "covariance": "full 4x4; propagated with the exact gradient",
        },
        "claim_boundary": {
            "exact": "continuum dimensions and Q velocities",
            "synthetic": "amplitude cancellation and covariance propagation",
            "unresolved": "physical Potts multiplicity, lattice overlap, and explicit field derivative",
            "forbidden": "substitute the FK measure score for the total field derivative",
        },
    }
    if input_path is None:
        result["lattice_score"] = {
            "status": "NOT_SCOREABLE",
            "reason": "no archive contains two sizes plus an explicit generic-Q field derivative",
        }
    else:
        payload = json.loads(input_path.read_text(encoding="utf-8"))
        result["lattice_score"] = score_input(payload, targets)
        result["lattice_score"]["input"] = str(input_path)
    return result


def main(argv: Iterable[str] = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    root = Path(__file__).resolve().parents[1]
    result = render(root, args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
