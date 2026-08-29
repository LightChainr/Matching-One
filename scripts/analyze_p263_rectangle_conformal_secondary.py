#!/usr/bin/env python3
"""Post-reveal fixed-rectangle conformal secondary for Issue #263.

The frozen primary scorer treats the normalized bottom-boundary coordinates as
upper-half-plane coordinates.  The actual lattice domains instead converge,
under mesh refinement, to the fixed rectangle [-2,4] x [0,4].  This script
maps that rectangle to the upper half-plane, recomputes the four cross ratios,
and applies the boundary-primary Jacobian and conformal-prefactor correction.
It never changes or replaces the frozen primary score.
"""

from __future__ import annotations

import argparse
import csv
from fractions import Fraction
import json
import math
from pathlib import Path

import mpmath as mp

from p263_boundary_tangent_ode import frobenius_high_branch
from score_p263_boundary_qscore_pilot import GEOMETRY_ORDER, _gls


LAMBDA_ORDER = (
    Fraction(1, 4),
    Fraction(1, 3),
    Fraction(2, 3),
    Fraction(3, 4),
)
ANCHOR_INDEX = 1
ACTIVE_INDICES = (0, 2, 3)
CHANNEL_FIELDS = {
    "1234": "count_1234",
    "12|34": "count_12_34",
    "14|23": "count_14_23",
}


def _mp_fraction(value: Fraction) -> mp.mpf:
    return mp.mpf(value.numerator) / value.denominator


def rectangle_modulus() -> dict[str, mp.mpf]:
    """Elliptic parameter for width 6, height 4.

    For z=sn(u|m), the rectangle u in [-K,K]+i[0,K'] maps to the
    upper half-plane.  Conformality of u=(K/3)(w-1) requires K'/K=4/3.
    """

    nome = mp.exp(-4 * mp.pi / 3)
    theta2 = mp.jtheta(2, 0, nome)
    theta3 = mp.jtheta(3, 0, nome)
    parameter = (theta2 / theta3) ** 4
    complete = mp.ellipk(parameter)
    complementary = mp.ellipk(1 - parameter)
    return {
        "nome": nome,
        "parameter_m": parameter,
        "modulus_k": mp.sqrt(parameter),
        "K": complete,
        "K_prime": complementary,
    }


def cai_prefactor(points: list[mp.mpf]) -> mp.mpf:
    x1, x2, x3, x4 = points
    return (x3 - x1) * (x4 - x2) / (
        (x2 - x1) * (x3 - x2) * (x4 - x3) * (x4 - x1)
    )


def cross_ratio(points: list[mp.mpf]) -> mp.mpf:
    x1, x2, x3, x4 = points
    return (x2 - x1) * (x4 - x3) / ((x3 - x1) * (x4 - x2))


def mapped_geometry(lam: Fraction, modulus: dict[str, mp.mpf]) -> dict:
    complete = modulus["K"]
    parameter = modulus["parameter_m"]
    second = 2 * _mp_fraction(lam) / (1 + _mp_fraction(lam))
    original = [mp.mpf(0), second, mp.mpf(1), mp.mpf(2)]
    images: list[mp.mpf] = []
    derivatives: list[mp.mpf] = []
    for point in original:
        u = complete * (point - 1) / 3
        images.append(mp.ellipfun("sn", u, parameter))
        derivatives.append(
            (complete / 3)
            * mp.ellipfun("cn", u, parameter)
            * mp.ellipfun("dn", u, parameter)
        )
    original_prefactor = cai_prefactor(original)
    image_prefactor = cai_prefactor(images)
    invariant_log_factor = (
        sum(mp.log(abs(value)) for value in derivatives)
        + 2 * mp.log(image_prefactor)
        - 2 * mp.log(original_prefactor)
    )
    h_prime = mp.sqrt(3) / (3 * mp.pi)
    return {
        "euclidean_lambda": _mp_fraction(lam),
        "bottom_points": original,
        "uhp_images": images,
        "map_derivatives": derivatives,
        "effective_lambda": cross_ratio(images),
        "euclidean_K": original_prefactor,
        "mapped_K": image_prefactor,
        "invariant_log_factor": invariant_log_factor,
        "z_correction": h_prime * invariant_log_factor,
    }


def high_branch_log_tangent(lam: mp.mpf, order: int = 140) -> mp.mpf:
    """d_Q log V at Q=1 modulo its lambda-independent amplitude."""

    _, coefficients = frobenius_high_branch(order)
    ordinary = mp.fsum(
        _mp_fraction(coefficient.value) * lam**index
        for index, coefficient in enumerate(coefficients)
    )
    regular_unit = mp.fsum(
        _mp_fraction(-Fraction(3, 2) * coefficient.derivative) * lam**index
        for index, coefficient in enumerate(coefficients)
    )
    return (mp.sqrt(3) / mp.pi) * (mp.log(lam) + regular_unit / ordinary)


def _anchored(values: list[mp.mpf]) -> list[mp.mpf]:
    return [value - values[ANCHOR_INDEX] for value in values]


def _score_residual(
    residual: list[mp.mpf], covariance: list[list[float]]
) -> dict:
    gls = _gls([float(value) for value in residual], covariance)
    chi_square = mp.mpf(str(gls["chi_square"]))
    degrees = int(gls["degrees_of_freedom"])
    survival = mp.gammainc(
        mp.mpf(degrees) / 2, chi_square / 2, mp.inf
    ) / mp.gamma(mp.mpf(degrees) / 2)
    return {
        "residual": [float(value) for value in residual],
        "joint_gls": gls,
        "chi_square_survival": float(survival),
    }


def conformal_rescore(score_payload: dict, geometries: list[dict]) -> dict:
    target_raw = [
        high_branch_log_tangent(geometry["effective_lambda"])
        for geometry in geometries
    ]
    target = _anchored(target_raw)
    old_anchored = [
        mp.mpf(str(value)) for value in score_payload["anchored_lattice_tangent"]
    ]
    old_z = [
        mp.mpf(str(row["z_before_amplitude_projection"]))
        for row in score_payload["estimates"]
    ]
    corrections = [geometry["z_correction"] for geometry in geometries]
    corrected_z = [value - correction for value, correction in zip(old_z, corrections)]
    corrected_anchored = _anchored(corrected_z)
    covariance = score_payload["residual_covariance"]

    lambda_only_residual = [
        old_anchored[index] - target[index] for index in ACTIVE_INDICES
    ]
    full_residual = [
        corrected_anchored[index] - target[index] for index in ACTIVE_INDICES
    ]
    return {
        "original_primary": {
            "residual": score_payload["residual"],
            "joint_gls": score_payload["joint_gls"],
        },
        "effective_target": [float(value) for value in target],
        "lambda_only_secondary": {
            "anchored_lattice_tangent": [float(value) for value in old_anchored],
            **_score_residual(lambda_only_residual, covariance),
        },
        "full_rectangle_conformal_secondary": {
            "anchored_lattice_tangent": [
                float(value) for value in corrected_anchored
            ],
            "anchored_deterministic_z_correction": [
                float(value - corrections[ANCHOR_INDEX]) for value in corrections
            ],
            **_score_residual(full_residual, covariance),
        },
    }


def event_totals(path: Path) -> dict[str, dict[str, int]]:
    totals = {
        geometry: {"samples": 0, **{channel: 0 for channel in CHANNEL_FIELDS}}
        for geometry in GEOMETRY_ORDER
    }
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            geometry = row["geometry_id"]
            totals[geometry]["samples"] += int(row["samples"])
            for channel, field in CHANNEL_FIELDS.items():
                totals[geometry][channel] += int(row[field])
    return totals


def ordinary_scaling(level1_path: Path, level2_path: Path) -> dict:
    level1 = event_totals(level1_path)
    level2 = event_totals(level2_path)
    expected_probability_ratio = 2 ** (-4 / 3)
    rows = []
    for geometry in GEOMETRY_ORDER:
        for channel in CHANNEL_FIELDS:
            count1 = level1[geometry][channel]
            count2 = level2[geometry][channel]
            probability1 = count1 / level1[geometry]["samples"]
            probability2 = count2 / level2[geometry]["samples"]
            ratio = probability2 / probability1
            rows.append(
                {
                    "geometry_id": geometry,
                    "channel": channel,
                    "level1_count": count1,
                    "level2_count": count2,
                    "level1_probability": probability1,
                    "level2_probability": probability2,
                    "probability_ratio_level2_over_level1": ratio,
                    "ratio_over_two_to_minus_four_thirds": (
                        ratio / expected_probability_ratio
                    ),
                    "effective_decay_exponent": -math.log(ratio, 2),
                }
            )
    pooled = []
    for channel in CHANNEL_FIELDS:
        count1 = sum(level1[geometry][channel] for geometry in GEOMETRY_ORDER)
        count2 = sum(level2[geometry][channel] for geometry in GEOMETRY_ORDER)
        samples1 = sum(level1[geometry]["samples"] for geometry in GEOMETRY_ORDER)
        samples2 = sum(level2[geometry]["samples"] for geometry in GEOMETRY_ORDER)
        ratio = (count2 / samples2) / (count1 / samples1)
        pooled.append(
            {
                "channel": channel,
                "level1_count": count1,
                "level2_count": count2,
                "probability_ratio_level2_over_level1": ratio,
                "ratio_over_two_to_minus_four_thirds": (
                    ratio / expected_probability_ratio
                ),
            }
        )
    return {
        "h_at_Q1": "1/3",
        "predicted_probability_ratio_for_doubled_span": expected_probability_ratio,
        "level2_over_level1_sample_ratio": 2.5,
        "predicted_level2_over_level1_event_count_ratio": (
            2.5 * expected_probability_ratio
        ),
        "by_geometry_and_channel": rows,
        "pooled_by_channel": pooled,
    }


def render(
    level1_batches: Path,
    level1_score: Path,
    level2_batches: Path,
    level2_score: Path,
) -> dict:
    mp.mp.dps = 80
    modulus = rectangle_modulus()
    geometries = [mapped_geometry(lam, modulus) for lam in LAMBDA_ORDER]
    scores = {}
    for name, path in (("level1_200k", level1_score), ("level2_500k", level2_score)):
        payload = json.loads(path.read_text(encoding="utf-8"))
        scores[name] = conformal_rescore(payload, geometries)
    return {
        "schema": "matching-one.p263-rectangle-conformal-secondary.v1",
        "issue": 263,
        "status": "post_reveal_secondary_not_primary",
        "frozen_primary_unchanged": True,
        "rectangle": {
            "normalized_domain": "[-2,4] x [0,4]",
            "map": "z=sn((K/3)*(w-1)|m)",
            "aspect_equation": "K'(m)/K(m)=4/3",
            "nome": mp.nstr(modulus["nome"], 60),
            "parameter_m": mp.nstr(modulus["parameter_m"], 60),
            "modulus_k": mp.nstr(modulus["modulus_k"], 60),
            "K": mp.nstr(modulus["K"], 60),
            "K_prime": mp.nstr(modulus["K_prime"], 60),
            "verified_K_prime_over_K": mp.nstr(
                modulus["K_prime"] / modulus["K"], 60
            ),
        },
        "geometries": [
            {
                key: [mp.nstr(item, 50) for item in value]
                if isinstance(value, list)
                else mp.nstr(value, 50)
                for key, value in geometry.items()
            }
            for geometry in geometries
        ],
        "scores": scores,
        "ordinary_probability_scaling": ordinary_scaling(
            level1_batches, level2_batches
        ),
        "interpretation_boundary": [
            "This is a post-reveal secondary; it does not replace or tune the frozen primary scorer.",
            "Span doubling refines a fixed rectangle and does not send the top or side boundaries to infinity.",
            "The full secondary includes both effective cross ratios and the primary-field Jacobian/prefactor correction.",
            "The two acquisition levels are scored separately and are not combined as independent continuum evidence.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--level1-batches", type=Path, required=True)
    parser.add_argument("--level1-score", type=Path, required=True)
    parser.add_argument("--level2-batches", type=Path, required=True)
    parser.add_argument("--level2-score", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = render(
        args.level1_batches,
        args.level1_score,
        args.level2_batches,
        args.level2_score,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
