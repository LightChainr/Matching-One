#!/usr/bin/env python3
"""Geometry-aware crosswalk of the N65/N85/N145 natural A current."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Mapping, Optional, Sequence


AREA_POWER = -13.0 / 8.0
GEOMETRIES = {
    65: ((8, 1), (7, 4)),
    85: ((9, 2), (7, 6)),
    145: ((12, 1), (9, 8)),
}
SOURCE_FILES = {
    65: "analysis/p337_natural_current_scale_preregistration.json",
    85: "results/server-20260830/P337-natural-current-scale-N85/score.json",
    145: "results/server-20260830/P337-natural-current-third-scale-N145/score.json",
}


def solve(matrix: Sequence[Sequence[float]], vector: Sequence[float]) -> list[float]:
    """Solve a small dense system with partial pivoting."""

    size = len(vector)
    rows = [list(map(float, matrix[i])) + [float(vector[i])] for i in range(size)]
    if any(len(row) != size + 1 for row in rows):
        raise ValueError("matrix must be square")
    tolerance = max(abs(value) for row in rows for value in row[:-1]) * 1e-14
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(rows[row][column]))
        if abs(rows[pivot][column]) <= max(tolerance, 1e-300):
            raise ValueError("matrix is singular")
        rows[column], rows[pivot] = rows[pivot], rows[column]
        scale = rows[column][column]
        rows[column] = [value / scale for value in rows[column]]
        for row in range(size):
            if row == column:
                continue
            scale = rows[row][column]
            rows[row] = [left - scale * right
                         for left, right in zip(rows[row], rows[column])]
    return [rows[i][-1] for i in range(size)]


def inverse(matrix: Sequence[Sequence[float]]) -> list[list[float]]:
    size = len(matrix)
    columns = [solve(matrix, [float(i == j) for i in range(size)])
               for j in range(size)]
    return [[columns[j][i] for j in range(size)] for i in range(size)]


def matvec(matrix: Sequence[Sequence[float]], vector: Sequence[float]) -> list[float]:
    return [math.fsum(a * b for a, b in zip(row, vector)) for row in matrix]


def matmul(
    left: Sequence[Sequence[float]], right: Sequence[Sequence[float]]
) -> list[list[float]]:
    columns = list(zip(*right))
    return [[math.fsum(a * b for a, b in zip(row, column)) for column in columns]
            for row in left]


def transpose(matrix: Sequence[Sequence[float]]) -> list[list[float]]:
    return [list(column) for column in zip(*matrix)]


def add(left: Sequence[Sequence[float]], right: Sequence[Sequence[float]]) -> list[list[float]]:
    return [[left[i][j] + right[i][j] for j in range(len(left))]
            for i in range(len(left))]


def quadratic(vector: Sequence[float], covariance: Sequence[Sequence[float]]) -> float:
    return math.fsum(a * b for a, b in zip(vector, solve(covariance, vector)))


def block_diagonal(blocks: Sequence[Sequence[Sequence[float]]]) -> list[list[float]]:
    width = sum(len(block) for block in blocks)
    output = [[0.0] * width for _ in range(width)]
    offset = 0
    for block in blocks:
        for i, row in enumerate(block):
            for j, value in enumerate(row):
                output[offset + i][offset + j] = float(value)
        offset += len(block)
    return output


def gaussian_chi4(a: int, b: int) -> tuple[Fraction, Fraction]:
    radius2 = a * a + b * b
    denominator = radius2 * radius2
    return (
        Fraction(a**4 - 6 * a * a * b * b + b**4, denominator),
        Fraction(4 * a * b * (a * a - b * b), denominator),
    )


def exact_descriptor(n: int, representation: tuple[int, int]) -> dict[str, object]:
    a, b = representation
    if a * a + b * b != n:
        raise ValueError(f"{representation} does not have norm {n}")
    real, imaginary = gaussian_chi4(a, b)
    return {
        "representation": [a, b],
        "tau": "i",
        "E4_modulus_class": "E4(i), common to every Gaussian square quotient",
        "z_axis_chi4": {
            "real": str(real), "imaginary": str(imaginary),
            "unit_norm_exact": real * real + imaginary * imaginary == 1,
        },
        "H4_reflection_even_covector": float(real),
        "H4_reflection_odd_covector_audit_only": float(imaginary),
        "charged_projective_A_scalar": 0.5,
        "charged_identity": (
            "q_A^2=(1,1,0,0)=(u+H_F3)/2 on P1(F3); the scalar u/2 "
            "is the one additional charged/projective descriptor"
        ),
    }


def load_scale(root: Path, n: int) -> dict[str, object]:
    path = root / SOURCE_FILES[n]
    payload = json.loads(path.read_text(encoding="utf-8"))
    natural = payload["natural_coordinate"]
    values = [float(value) for value in natural["value"][:2]]
    covariance = [[float(value) for value in row[:2]]
                  for row in natural["covariance"][:2]]
    return {
        "N": n,
        "path": SOURCE_FILES[n],
        "values": values,
        "covariance": covariance,
        "representations": [list(row) for row in GEOMETRIES[n]],
        "descriptors": [exact_descriptor(n, row) for row in GEOMETRIES[n]],
    }


def design_row(scale: Mapping[str, object], columns: Sequence[str]) -> list[list[float]]:
    radial = float(scale["N"]) ** AREA_POWER
    rows = []
    for descriptor in scale["descriptors"]:
        values = {
            "pure_N": 1.0,
            "H4_geometry": descriptor["H4_reflection_even_covector"],
            "A_projective_scalar": descriptor["charged_projective_A_scalar"],
        }
        rows.append([radial * values[name] for name in columns])
    return rows


def fit_model(
    name: str,
    columns: Sequence[str],
    training: Sequence[Mapping[str, object]],
    heldout: Mapping[str, object],
) -> dict[str, object]:
    y = [value for scale in training for value in scale["values"]]
    covariance = block_diagonal([scale["covariance"] for scale in training])
    x = [row for scale in training for row in design_row(scale, columns)]
    weight = inverse(covariance)
    xt = transpose(x)
    normal = matmul(matmul(xt, weight), x)
    parameter_covariance = inverse(normal)
    rhs = matvec(matmul(xt, weight), y)
    parameters = matvec(parameter_covariance, rhs)
    fitted = matvec(x, parameters)
    training_residual = [actual - expected for actual, expected in zip(y, fitted)]

    x_test = design_row(heldout, columns)
    predicted = matvec(x_test, parameters)
    fit_covariance = matmul(matmul(x_test, parameter_covariance), transpose(x_test))
    predictive_covariance = add(heldout["covariance"], fit_covariance)
    heldout_residual = [actual - expected
                        for actual, expected in zip(heldout["values"], predicted)]
    pair = [-1.0, 1.0]
    pair_observed = math.fsum(a * b for a, b in zip(pair, heldout["values"]))
    pair_predicted = math.fsum(a * b for a, b in zip(pair, predicted))
    pair_variance = math.fsum(pair[i] * predictive_covariance[i][j] * pair[j]
                              for i in range(2) for j in range(2))
    return {
        "name": name,
        "columns": list(columns),
        "area_power": AREA_POWER,
        "fit_scales": [scale["N"] for scale in training],
        "parameters": parameters,
        "parameter_standard_errors": [
            math.sqrt(max(0.0, parameter_covariance[i][i]))
            for i in range(len(parameters))
        ],
        "parameter_covariance": parameter_covariance,
        "training_prediction": fitted,
        "training_residual": training_residual,
        "training_quadratic": quadratic(training_residual, covariance),
        "training_df": len(y) - len(parameters),
        "heldout_N": heldout["N"],
        "heldout_observed": heldout["values"],
        "heldout_prediction": predicted,
        "heldout_residual": heldout_residual,
        "heldout_measurement_covariance": heldout["covariance"],
        "heldout_fit_covariance": fit_covariance,
        "heldout_predictive_covariance": predictive_covariance,
        "heldout_predictive_quadratic": quadratic(heldout_residual, predictive_covariance),
        "heldout_df": 2,
        "heldout_pair_contrast": {
            "observed": pair_observed,
            "predicted": pair_predicted,
            "residual": pair_observed - pair_predicted,
            "predictive_standard_error": math.sqrt(pair_variance),
            "z": (pair_observed - pair_predicted) / math.sqrt(pair_variance),
        },
    }


def crosswalk(root: Path) -> dict[str, object]:
    scales = {n: load_scale(root, n) for n in GEOMETRIES}
    training = [scales[65], scales[85]]
    heldout = scales[145]
    models = [
        fit_model("pure_N_law", ["pure_N"], training, heldout),
        fit_model("one_H4_geometry_covector", ["H4_geometry"], training, heldout),
        fit_model(
            "H4_geometry_plus_A_projective_scalar",
            ["H4_geometry", "A_projective_scalar"],
            training,
            heldout,
        ),
    ]
    by_name = {model["name"]: model for model in models}
    h4 = by_name["one_H4_geometry_covector"]
    extended = by_name["H4_geometry_plus_A_projective_scalar"]
    pure = by_name["pure_N_law"]
    pair = [-1.0, 1.0]

    def pair_value(scale: Mapping[str, object]) -> float:
        return math.fsum(a * b for a, b in zip(pair, scale["values"]))

    def pair_variance(scale: Mapping[str, object]) -> float:
        return math.fsum(
            pair[i] * scale["covariance"][i][j] * pair[j]
            for i in range(2) for j in range(2)
        )

    def delta_h4(scale: Mapping[str, object]) -> float:
        values = [row["H4_reflection_even_covector"] for row in scale["descriptors"]]
        return values[1] - values[0]

    radial_85_to_145 = (145.0 / 85.0) ** AREA_POWER
    angular_85_to_145 = delta_h4(scales[145]) / delta_h4(scales[85])
    radial_only_factor = radial_85_to_145
    geometry_factor = radial_85_to_145 * angular_85_to_145
    observed_pair_145 = pair_value(scales[145])
    radial_only_target = pair_value(scales[85]) * radial_only_factor
    geometry_target = pair_value(scales[85]) * geometry_factor
    total_rebound = observed_pair_145 - radial_only_target
    geometry_increment = geometry_target - radial_only_target
    curvature_remainder = observed_pair_145 - geometry_target

    def anchored_score(target: float, transfer: float) -> dict[str, float]:
        variance = pair_variance(scales[145]) + transfer * transfer * pair_variance(scales[85])
        residual = observed_pair_145 - target
        return {
            "target": target,
            "residual": residual,
            "predictive_standard_error": math.sqrt(variance),
            "z": residual / math.sqrt(variance),
            "quadratic": residual * residual / variance,
        }

    rebound_decomposition = {
        "N85_to_N145_radial_factor": radial_85_to_145,
        "N85_to_N145_H4_angular_ratio": angular_85_to_145,
        "N85_to_N145_geometry_aware_factor": geometry_factor,
        "radial_only_N85_anchored": anchored_score(radial_only_target, radial_only_factor),
        "geometry_aware_N85_anchored": anchored_score(geometry_target, geometry_factor),
        "central_geometry_increment": geometry_increment,
        "central_scale_curvature_remainder": curvature_remainder,
        "central_geometry_fraction_of_apparent_rebound": geometry_increment / total_rebound,
        "central_curvature_fraction_of_apparent_rebound": curvature_remainder / total_rebound,
        "interpretation_boundary": (
            "the fractions are a central-value accounting, not a covariance-weighted model probability"
        ),
    }
    next_geometry = {
        "selection": "N170 exact angle-flip child of the N85 lineage",
        "parent": {"N": 85, "pair": [[9, 2], [7, 6]]},
        "common_Gaussian_multiplier": [1, 1],
        "child": {"N": 170, "pair": [[11, 7], [13, 1]]},
        "exact_transport": {
            "H4_geometry_covector": "child=-parent for both orientations",
            "A_projective_scalar": "child=parent=1/2",
            "area_ratio": 2,
            "project_H4_radial_ratio": 2 ** AREA_POWER,
            "F3_reduction": "nondegenerate; norm multiplier 2 mod 3",
        },
        "why": (
            "it follows the better-measured N85 source, lies near N145 in area, "
            "and flips H4 while preserving the charged scalar, separating the two "
            "descriptors without selecting another arbitrary Gaussian pair"
        ),
        "archive_status": "no N170 projective-birth archive with tau1,ell,tau2 exists",
        "minimum_missing_fields": [
            "orientation", "batch", "samples", "tau1", "tau2",
            "kind", "ell_x", "ell_y", "count",
        ],
        "production_status": "selected only; no new simulation authorized or run",
    }
    covariance_6d = block_diagonal([scales[n]["covariance"] for n in (65, 85, 145)])
    values_6d = [value for n in (65, 85, 145) for value in scales[n]["values"]]
    return {
        "schema": "matching-one/P337-natural-current-geometry-crosswalk/v1",
        "status": "existing-data geometry-aware crosswalk; N145 held out of direction selection",
        "observable": "K_A=p(1-p)Jminus_A/W_A=d_eta log W_A",
        "scale_contract": {
            "area_power": AREA_POWER,
            "name": "project H4 N^-13/8",
            "no_exponent_fit": True,
        },
        "descriptor_contract": {
            "tau": "all six quotients have exact square modulus tau=i; E4(i) is common and cannot distinguish them",
            "H4": "reflection-even Re[((a+ib)/sqrt(N))^4], stored exactly",
            "charged_projective": "q_A^2=(u+H_F3)/2, hence one scalar u/2 descriptor equal to 1/2",
            "reflection_odd_boundary": "Im chi4 is stored as an audit but excluded from the reflection-even K_A models",
        },
        "data": {
            "order": [
                "N65_first", "N65_second", "N85_first", "N85_second",
                "N145_first", "N145_second",
            ],
            "values": values_6d,
            "covariance": covariance_6d,
            "scales": [scales[n] for n in (65, 85, 145)],
            "independence": "block diagonal across N; full paired orientation covariance retained within N",
        },
        "models": models,
        "diagnosis": {
            "pure_N_training_quadratic": pure["training_quadratic"],
            "H4_training_quadratic": h4["training_quadratic"],
            "H4_plus_projective_training_quadratic": extended["training_quadratic"],
            "H4_heldout_predictive_quadratic": h4["heldout_predictive_quadratic"],
            "H4_plus_projective_heldout_predictive_quadratic": extended["heldout_predictive_quadratic"],
            "projective_increment_training_delta_quadratic": (
                h4["training_quadratic"] - extended["training_quadratic"]
            ),
            "central_rebound_decomposition": rebound_decomposition,
            "reading": (
                "the pure area law without geometry fails already on N65/N85. A single "
                "exact H4 covector fits their direction and predicts both held-out N145 "
                "orientations; adding the charged scalar changes little. Exact geometry "
                "rotation explains 23 percent of the central rebound from the N85-anchored "
                "radial target, leaving 77 percent as a central scale-curvature remainder. "
                "That remainder is not resolved: the joint H4 geometry model scores 1.52/2 "
                "on N145. Geometry is necessary, while extra scale curvature is suggested "
                "by central values but not required at current precision."
            ),
            "boundary": (
                "geometry sufficiency is not asymptotic H4 proof; N65 is noisy, each N uses "
                "a different pair, and the one extra charged descriptor is only weakly identified"
            ),
        },
        "next_same_lineage_geometry": next_geometry,
    }


def render_markdown(payload: Mapping[str, object]) -> str:
    lines = [
        "# Geometry-aware crosswalk of the natural A current",
        "",
        "The model direction is selected from N65/N85 only. N145 is a held-out two-component diagnostic.",
        "",
        "All Gaussian quotients have exact `tau=i`, so `E4(i)` is common. The varying exact descriptor is `Re z_axis=cos(4 theta)`. The charged source obeys `q_A^2=(u+H_F3)/2`, giving one additional projective scalar `u/2=1/2`.",
        "",
        "| model | training chi2 / df | N145 predictive chi2 / 2 | predicted N145 pair | pair residual / SE |",
        "|---|---:|---:|---:|---:|",
    ]
    for model in payload["models"]:
        pair = model["heldout_pair_contrast"]
        lines.append(
            f"| {model['name']} | {model['training_quadratic']:.4f} / {model['training_df']} | "
            f"{model['heldout_predictive_quadratic']:.4f} | {pair['predicted']:.8f} | "
            f"{pair['z']:.3f} |"
        )
    h4 = payload["models"][1]
    extended = payload["models"][2]
    lines += [
        "",
        f"The exact H4 model predicts N145 components `{h4['heldout_prediction'][0]:+.8f}, {h4['heldout_prediction'][1]:+.8f}` against observed `{h4['heldout_observed'][0]:+.8f}, {h4['heldout_observed'][1]:+.8f}`. Its full predictive score is `{h4['heldout_predictive_quadratic']:.3f}/2`; the pair residual is only `{h4['heldout_pair_contrast']['z']:.3f}` predictive SE.",
        "",
        f"Adding the charged scalar improves training chi-square by only `{payload['diagnosis']['projective_increment_training_delta_quadratic']:.3f}` and gives held-out `{extended['heldout_predictive_quadratic']:.3f}/2`. It is not selected by these data.",
        "",
        f"In central-value accounting, exact angle rotation contributes `{payload['diagnosis']['central_rebound_decomposition']['central_geometry_fraction_of_apparent_rebound']:.1%}` of the apparent rebound relative to the N85-anchored radial target; the remaining `{payload['diagnosis']['central_rebound_decomposition']['central_curvature_fraction_of_apparent_rebound']:.1%}` is a scale-curvature remainder. But that remainder is not resolved: the joint H4 geometry model already gives an acceptable held-out score. Geometry is necessary; extra curvature is suggested centrally but not required statistically.",
        "",
        "The next clean geometry is the N170 angle-flip child `(11+7i,13+i)` of N85 under common multiplier `1+i`: H4 flips exactly, the charged scalar stays fixed, and the area is close to N145. No matching projective-birth archive exists, so it is selected but not run.",
        "",
    ]
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args(argv)
    payload = crosswalk(args.root)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_markdown(payload), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
