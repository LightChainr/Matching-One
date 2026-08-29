#!/usr/bin/env python3
"""One-amplitude retrospective score of the Issue #231 KdV sector vector."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import mpmath as mp


DESIGNS = ("pell_Dminus2_N30", "pell_Dplus1_N56")
CONTRAST_ORDER = ("C_nontrivial_real", "Q_reflection_null", "S_scalar")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def block_covariance(pilot_by_id) -> mp.matrix:
    covariance = mp.zeros(6)
    for block, identifier in enumerate(DESIGNS):
        local = pilot_by_id[identifier]["contrast_covariance_of_mean"]
        for row in range(3):
            for column in range(3):
                covariance[3 * block + row, 3 * block + column] = local[row][column]
    return covariance


def one_amplitude_gls(
    observed: mp.matrix,
    theory: mp.matrix,
    covariance: mp.matrix,
) -> dict[str, object]:
    inverse = covariance**-1
    information = (theory.T * inverse * theory)[0]
    amplitude = (theory.T * inverse * observed)[0] / information
    residual = observed - amplitude * theory
    chi_square = (residual.T * inverse * residual)[0]
    degrees = len(observed) - 1
    survival = mp.gammainc(
        mp.mpf(degrees) / 2, chi_square / 2, mp.inf
    ) / mp.gamma(mp.mpf(degrees) / 2)
    return {
        "amplitude": float(amplitude),
        "amplitude_standard_error": float(mp.sqrt(1 / information)),
        "residual": [float(value) for value in residual],
        "chi_square": float(chi_square),
        "degrees_of_freedom": degrees,
        "chi_square_survival_p": float(survival),
    }


def build_score(pilot_path: Path, oracle_path: Path) -> dict[str, object]:
    pilot = json.loads(pilot_path.read_text(encoding="utf-8"))
    oracle = json.loads(oracle_path.read_text(encoding="utf-8"))
    if pilot.get("schema") != "p156-square-bond-primitive-pilot-v1":
        raise ValueError("unexpected PR #222 pilot schema")
    if oracle.get("schema") != "matching-one/p231-pinson-arguin-vacuum-kdv-oracle/v1":
        raise ValueError("unexpected Issue #231 oracle schema")
    pilot_by_id = {record["design"]: record for record in pilot["pilot"]}
    oracle_by_id = {record["id"]: record for record in oracle["records"]}
    if set(DESIGNS) - pilot_by_id.keys() or set(DESIGNS) - oracle_by_id.keys():
        raise ValueError("N30/N56 design missing from pilot or oracle")

    observed_values = []
    theory_values = []
    for identifier in DESIGNS:
        source = pilot_by_id[identifier]
        deterministic = oracle_by_id[identifier]
        if tuple(source["contrast_order"]) != CONTRAST_ORDER:
            raise ValueError(f"{identifier} pilot contrast order changed")
        if tuple(deterministic["contrast_order"]) != CONTRAST_ORDER:
            raise ValueError(f"{identifier} oracle contrast order changed")
        if source["N_vertices"] != deterministic["N"]:
            raise ValueError(f"{identifier} size mismatch")
        (a, _), (c, _) = source["period_matrix_rows"]
        if a * a + c * c != deterministic["omega1_length_squared"]:
            raise ValueError(f"{identifier} omega1 scale mismatch")
        observed_values.extend(
            source["contrasts"][name]["value"] for name in CONTRAST_ORDER
        )
        theory_values.extend(
            float(value) for value in deterministic["finite_size_design_vector"]
        )

    observed = mp.matrix(observed_values)
    theory = mp.matrix(theory_values)
    covariance = block_covariance(pilot_by_id)
    full = one_amplitude_gls(observed, theory, covariance)

    # The non-scalar C projection is reported separately because S visibly
    # contains a scalar correction independent of a spin-4 C3 response.
    selected = (0, 3)
    c_observed = mp.matrix([observed[index] for index in selected])
    c_theory = mp.matrix([theory[index] for index in selected])
    c_covariance = mp.matrix(
        [[covariance[first, second] for second in selected] for first in selected]
    )
    non_scalar = one_amplitude_gls(c_observed, c_theory, c_covariance)

    return {
        "schema": "matching-one/p231-vacuum-kdv-sector-score/v1",
        "issue": 231,
        "analysis_class": "retrospective_mechanism_diagnostic_on_PR222",
        "joint_order": [
            f"{identifier}_{name}" for identifier in DESIGNS for name in CONTRAST_ORDER
        ],
        "observed": observed_values,
        "theory_vector_per_unit_g4": theory_values,
        "covariance": [
            [float(covariance[row, column]) for column in range(6)]
            for row in range(6)
        ],
        "full_CQS_one_amplitude_score": full,
        "non_scalar_C_only_diagnostic": non_scalar,
        "structural_predictions": {
            "Q_reflection_null": [theory_values[1], theory_values[4]],
            "C_same_sign_across_N30_N56": theory_values[0] * theory_values[3] > 0,
            "C_theory_ratio_N30_over_N56": theory_values[0] / theory_values[3],
            "C_observed_ratio_N30_over_N56": observed_values[0] / observed_values[3],
            "S_response_changes_sign": theory_values[2] * theory_values[5] < 0,
        },
        "decision": {
            "full_sector_vector": "FAIL" if full["chi_square_survival_p"] < 0.05 else "PASS",
            "non_scalar_C_direction": (
                "COMPATIBLE"
                if non_scalar["chi_square_survival_p"] >= 0.05
                else "INCOMPATIBLE"
            ),
            "summary": (
                "The normalized KdV sector response predicts the observed same-sign C "
                "transport and its two-size ratio, but one KdV amplitude cannot explain "
                "the large positive scalar S residual. The PR222 non-scalar direction "
                "remains KdV-compatible; the complete C/Q/S vector is not KdV-pure."
            ),
        },
        "inputs": {
            "pilot_result": str(pilot_path),
            "pilot_result_sha256": sha256(pilot_path),
            "kdv_oracle": str(oracle_path),
            "kdv_oracle_sha256": sha256(oracle_path),
        },
        "governance": {
            "new_monte_carlo": False,
            "new_independent_evidence": False,
            "overall_amplitude_fitted_after_reveal": True,
            "do_not_count_separately_from_PR222": True,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot-result", type=Path, required=True)
    parser.add_argument("--oracle", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = build_score(args.pilot_result, args.oracle)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
