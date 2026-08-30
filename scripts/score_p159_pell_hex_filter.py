#!/usr/bin/env python3
"""Score the Issue #159 Pell/hex phase gate from PR #222 batch counts.

This is a post-reveal operational bridge audit.  It reuses the exact oracle,
continuum baselines, and 200k primitive-sector sufficient statistics already
produced by PR #222; it does not create a second simulation or a new
preregistered evidence block.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp

from hexagonal_pell_spin_filter import eisenstein_e4


DESIGNS = {
    "pell_Dminus2_N30": {
        "N": 30, "D": -2, "period_matrix_rows": [[6, 3], [0, 5]],
        "task_phase_coordinate_approx": -0.150698,
    },
    "pell_Dplus1_N56": {
        "N": 56, "D": 1, "period_matrix_rows": [[8, 4], [0, 7]],
        "task_phase_coordinate_approx": 0.040410,
    },
}
CATEGORIES = (
    "rank0", "l0", "l1", "l2", "rank1_other", "rank2",
    "invariant_failure",
)
TARGET_LINES = [[1, 0], [0, 1], [1, -1]]
CONTRAST_ORDER = ("C_nontrivial_real", "Q_reflection_null", "S_scalar")
TRANSFORM = (
    (1.0, -0.5, -0.5),
    (0.0, -math.sqrt(3.0) / 2.0, math.sqrt(3.0) / 2.0),
    (1.0, 1.0, 1.0),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def covariance_of_mean(rows):
    count = len(rows)
    means = [math.fsum(row[index] for row in rows) / count
             for index in range(len(rows[0]))]
    return [[
        math.fsum((row[i] - means[i]) * (row[j] - means[j]) for row in rows)
        / (count * (count - 1))
        for j in range(len(means))
    ] for i in range(len(means))]


def transform_covariance(transform, covariance):
    return [[
        math.fsum(transform[i][k] * covariance[k][ell] * transform[j][ell]
                  for k in range(3) for ell in range(3))
        for j in range(3)
    ] for i in range(3)]


def read_batches(path: Path):
    grouped = {identifier: [] for identifier in DESIGNS}
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            identifier = raw["design"]
            if identifier not in grouped:
                raise ValueError(f"unexpected pilot design {identifier}")
            row = {
                "batch": int(raw["batch"]), "samples": int(raw["samples"]),
                "counts": {name: int(raw[name]) for name in CATEGORIES},
            }
            if sum(row["counts"].values()) != row["samples"]:
                raise ValueError(f"{identifier} batch {row['batch']} counts do not partition samples")
            if row["counts"]["invariant_failure"] != 0:
                raise ValueError(f"{identifier} has a primitive-sector invariant failure")
            grouped[identifier].append(row)
    for identifier, rows in grouped.items():
        rows.sort(key=lambda row: row["batch"])
        if [row["batch"] for row in rows] != list(range(100)):
            raise ValueError(f"{identifier} must contain 100 contiguous batches")
        if any(row["samples"] != 2000 for row in rows):
            raise ValueError(f"{identifier} must contain 2000 samples per batch")
    return grouped


def score_design(identifier, rows, source_design):
    spec = DESIGNS[identifier]
    if source_design["N_vertices"] != spec["N"]:
        raise ValueError(f"{identifier} N changed")
    if source_design["period_matrix_rows"] != spec["period_matrix_rows"]:
        raise ValueError(f"{identifier} period matrix changed")
    if source_design["engine_windings"] != TARGET_LINES:
        raise ValueError(f"{identifier} target-line basis changed")
    baselines = [float(value) for value in source_design["continuum_baselines"]]
    batch_probabilities = [
        [row["counts"][name] / row["samples"] for name in ("l0", "l1", "l2")]
        for row in rows
    ]
    probabilities = [
        math.fsum(row[index] for row in batch_probabilities) / len(rows)
        for index in range(3)
    ]
    residuals = [value - baseline for value, baseline in zip(probabilities, baselines)]
    contrasts = [
        math.fsum(TRANSFORM[i][j] * residuals[j] for j in range(3))
        for i in range(3)
    ]
    probability_covariance = covariance_of_mean(batch_probabilities)
    contrast_covariance = transform_covariance(TRANSFORM, probability_covariance)
    output_contrasts = {}
    for index, name in enumerate(CONTRAST_ORDER):
        standard_error = math.sqrt(max(0.0, contrast_covariance[index][index]))
        output_contrasts[name] = {
            "value": contrasts[index],
            "standard_error": standard_error,
            "z": contrasts[index] / standard_error,
        }
        registered = source_design["contrasts"][name]
        if not math.isclose(contrasts[index], registered["value"], abs_tol=2e-15):
            raise ValueError(f"{identifier}/{name} disagrees with PR #222 result")
    if any(
        not math.isclose(contrast_covariance[i][j],
                         source_design["contrast_covariance_of_mean"][i][j],
                         abs_tol=2e-18)
        for i in range(3) for j in range(3)
    ):
        raise ValueError(f"{identifier} covariance disagrees with PR #222 result")
    return {
        "N": spec["N"], "D": spec["D"],
        "period_matrix_rows": spec["period_matrix_rows"],
        "probabilities": probabilities, "continuum_baselines": baselines,
        "residuals": residuals, "contrasts": output_contrasts,
        "probability_covariance_of_mean": probability_covariance,
        "contrast_order": list(CONTRAST_ORDER),
        "contrast_covariance_of_mean": contrast_covariance,
    }


def block_diagonal_covariance(scored):
    output = [[0.0] * 6 for _ in range(6)]
    for block, identifier in enumerate(DESIGNS):
        covariance = scored[identifier]["contrast_covariance_of_mean"]
        for i in range(3):
            for j in range(3):
                output[3 * block + i][3 * block + j] = covariance[i][j]
    return output


def joint_zero_score(scored, contrast_name):
    chi_square = 0.0
    values = []
    for identifier in DESIGNS:
        row = scored[identifier]["contrasts"][contrast_name]
        values.append(row["value"])
        chi_square += row["z"] ** 2
    # The two design seeds are independently derived in PR #222.  For two
    # degrees of freedom the chi-square survival is exactly exp(-x/2).
    return {
        "values": values, "chi_square": chi_square, "dof": 2,
        "chi_square_survival_p": math.exp(-chi_square / 2.0),
    }


def e4_phase_coordinates():
    with mp.workdps(80):
        reference = eisenstein_e4(1j, 150)
        output = {}
        for identifier, spec in DESIGNS.items():
            matrix = spec["period_matrix_rows"]
            tau = mp.mpc(matrix[0][1], matrix[1][1]) / matrix[0][0]
            ratio = mp.re(eisenstein_e4(tau, 150) / reference)
            output[identifier] = {
                "tau_real": float(mp.re(tau)), "tau_imag": float(mp.im(tau)),
                "task_phase_coordinate_approx": spec["task_phase_coordinate_approx"],
                "repository_E4_over_E4_i": float(ratio),
                "signs_agree": (ratio > 0) == (spec["task_phase_coordinate_approx"] > 0),
            }
        return output


def validate_basis_transport(source):
    pilot = {row["design"]: row for row in source["pilot"]}
    for identifier, spec in DESIGNS.items():
        row = pilot[identifier]
        matrix = spec["period_matrix_rows"]
        if matrix[1][0] != 0 or 2 * matrix[0][1] != matrix[0][0]:
            raise ValueError(f"{identifier} is not in the positive-rho vertical basis")
        if row["engine_windings"] != TARGET_LINES:
            raise ValueError(f"{identifier} line registry changed")
    # A=[[0,-1],[1,1]] acts on column homology coordinates.  Up to sign it
    # cycles l0 -> l1 -> l2 -> l0 in the common transported registry.
    action = ((0, -1), (1, 1))
    def canonical(vector):
        x, y = vector
        divisor = math.gcd(abs(x), abs(y))
        x, y = x // divisor, y // divisor
        return [-x, -y] if x < 0 or (x == 0 and y < 0) else [x, y]
    cycled = [
        canonical((action[0][0] * line[0] + action[0][1] * line[1],
                   action[1][0] * line[0] + action[1][1] * line[1]))
        for line in TARGET_LINES
    ]
    if cycled != [TARGET_LINES[1], TARGET_LINES[2], TARGET_LINES[0]]:
        raise AssertionError("positive-rho C3 action does not cycle the frozen lines")
    return {
        "passed": True,
        "normalized_period_basis": "omega1=1, omega2=1/2+i*y with y>0",
        "transport_Dminus2_to_Dplus1": [[1, 0], [0, 1]],
        "positive_rho_C3_action": [list(row) for row in action],
        "ordered_unoriented_line_cycle": ["l0->l1", "l1->l2", "l2->l0"],
        "C_row_transport": [1.0, -0.5, -0.5],
    }


def build_score(batch_path: Path, source_path: Path):
    source = json.loads(source_path.read_text(encoding="utf-8"))
    if source.get("schema") != "p156-square-bond-primitive-pilot-v1":
        raise ValueError("unexpected PR #222 source schema")
    if source.get("samples_per_design") != 200000 or source.get("batches_per_design") != 100:
        raise ValueError("source pilot is not the frozen 200k/100-batch run")
    oracle = source.get("exact_oracle", {})
    expected_counts = {
        "rank0": 75, "l0": 57, "l1": 24, "l2": 24,
        "rank1_other": 1, "rank2": 75, "invariant_failure": 0,
    }
    if not oracle.get("passed") or oracle.get("counts") != expected_counts:
        raise ValueError("N=4 exact primitive-sector oracle failed")
    basis_transport = validate_basis_transport(source)
    grouped = read_batches(batch_path)
    source_designs = {row["design"]: row for row in source["pilot"]}
    scored = {
        identifier: score_design(identifier, grouped[identifier], source_designs[identifier])
        for identifier in DESIGNS
    }
    phase = e4_phase_coordinates()

    c_minus = scored["pell_Dminus2_N30"]["contrasts"]["C_nontrivial_real"]
    c_plus = scored["pell_Dplus1_N56"]["contrasts"]["C_nontrivial_real"]
    observed_same_sign = c_minus["value"] * c_plus["value"] > 0.0
    phase_opposite_sign = (
        phase["pell_Dminus2_N30"]["repository_E4_over_E4_i"] *
        phase["pell_Dplus1_N56"]["repository_E4_over_E4_i"] < 0.0
    )
    h4_sign_pass = observed_same_sign != phase_opposite_sign

    # Historically frozen conditional simple-zero rule (commit 46f3a6f):
    # for the square-bond identity-like H4 sector, A_D=N^2*C_D and
    # A_minus/A_plus -> -2.  Applying that old rule to C is post-reveal.
    amplitude_minus = 30.0**2 * c_minus["value"]
    amplitude_plus = 56.0**2 * c_plus["value"]
    amplitude_minus_se = 30.0**2 * c_minus["standard_error"]
    amplitude_plus_se = 56.0**2 * c_plus["standard_error"]
    minus2_residual = amplitude_minus + 2.0 * amplitude_plus
    minus2_residual_se = math.hypot(amplitude_minus_se, 2.0 * amplitude_plus_se)

    detection = joint_zero_score(scored, "C_nontrivial_real")
    reflection = joint_zero_score(scored, "Q_reflection_null")
    scalar = joint_zero_score(scored, "S_scalar")
    detection_pass = all(
        scored[identifier]["contrasts"]["C_nontrivial_real"]["z"] >= 3.0
        for identifier in DESIGNS
    )
    reflection_pass = all(
        abs(scored[identifier]["contrasts"]["Q_reflection_null"]["z"]) <= 2.0
        for identifier in DESIGNS
    ) and reflection["chi_square_survival_p"] >= 0.05

    return {
        "schema": "matching-one/p159-pell-hex-filter-score/v1",
        "issue": 159,
        "analysis_class": "post_reveal_operational_bridge_audit",
        "observable_descriptor": {
            "model": "critical square-bond percolation",
            "topology_channel": "full-configuration primitive rank-1 homology subgroup",
            "occupation_convention": "open bonds at fixed p=1/2",
            "combination": "Pinson-Arguin continuum-subtracted C3 character coordinates",
            "orientation_order": "positive-rho basis l0=(1,0),l1=(0,1),l2=(1,-1)",
            "probability_coordinate": "fixed_p",
            "normalization": "raw sector probabilities; C,Q,S linear transform",
            "quantity": "C primary, Q reflection null, S scalar control",
            "period_basis_convention": "period-matrix columns; engine (u,v)=u*omega1+v*omega2",
        },
        "exact_oracle": {"passed": True, "counts": expected_counts},
        "basis_transport": basis_transport,
        "phase_coordinates": phase,
        "design_scores": scored,
        "joint_contrast_order": [
            "N30_C", "N30_Q", "N30_S", "N56_C", "N56_Q", "N56_S",
        ],
        "joint_contrast_covariance_of_mean": block_diagonal_covariance(scored),
        "zero_scores": {
            "C_nontrivial_primary": detection,
            "Q_reflection_null": reflection,
            "S_scalar_control": scalar,
        },
        "gates": {
            "exact_oracle": {"passed": True},
            "basis_parallel_transport": {"passed": basis_transport["passed"]},
            "nontrivial_character_detection": {
                "rule": "C z>=3 with the transported positive convention on both sides",
                "passed": detection_pass,
            },
            "reflection_null": {
                "rule": "abs(Q z)<=2 on each side and joint zero p>=0.05",
                "passed": reflection_pass,
            },
            "ordinary_H4_simple_zero_phase": {
                "rule": "transported C must reverse sign when the E4 phase coordinate reverses",
                "phase_coordinates_have_opposite_sign": phase_opposite_sign,
                "observed_C_has_same_sign": observed_same_sign,
                "passed": h4_sign_pass,
            },
            "H4_specific_filter": {
                "passed": False,
                "reason": (
                    "the simple-zero H4 sign gate fails; the observed same-sign C pattern is "
                    "compatible with an even-in-E4/H8-like phase and the three-line real "
                    "character cannot otherwise distinguish H4 from H8"
                ),
            },
        },
        "exploratory_historical_minus2_score": {
            "status": "post_reveal_application_not_preregistered_for_C",
            "historical_rule_provenance": {
                "commit": "46f3a6fc175a620c9be2763781ed5820581e2982",
                "path": "predictions/elliptic_point_pell_spin_projector_20260828.yaml",
                "conditional_target": "A_D=N^2*C_D; A_Dminus2/A_Dplus1=-2",
            },
            "A_Dminus2": amplitude_minus,
            "A_Dminus2_standard_error": amplitude_minus_se,
            "A_Dplus1": amplitude_plus,
            "A_Dplus1_standard_error": amplitude_plus_se,
            "observed_ratio_Dminus2_over_Dplus1": amplitude_minus / amplitude_plus,
            "target_residual_Aminus_plus_2_Aplus": minus2_residual,
            "residual_standard_error": minus2_residual_se,
            "residual_z": minus2_residual / minus2_residual_se,
            "passes_abs_z_le_2": abs(minus2_residual / minus2_residual_se) <= 2.0,
        },
        "decision": {
            "primitive_character_bridge": (
                "PASS" if detection_pass and reflection_pass else "FAIL"
            ),
            "ordinary_H4_simple_zero_bridge": "PASS" if h4_sign_pass else "FAIL",
            "square_site_H4_promotion": "BLOCKED",
            "summary": (
                "The primitive C3 observable is real and well controlled, but it fails the "
                "transported ordinary-H4 phase reversal. It is not an H4-specific Pell filter."
            ),
        },
        "inputs": {
            "batches": str(batch_path), "batches_sha256": sha256(batch_path),
            "source_result": str(source_path), "source_result_sha256": sha256(source_path),
        },
        "governance": {
            "new_simulation": False,
            "new_primary_evidence": False,
            "do_not_count_separately_from_PR222": True,
            "minus2_score_is_preregistered_for_this_C_observable": False,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batches", type=Path, required=True)
    parser.add_argument("--source-result", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = build_score(args.batches, args.source_result)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
