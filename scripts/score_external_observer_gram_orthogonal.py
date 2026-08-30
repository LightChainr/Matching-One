#!/usr/bin/env python3
"""Frozen Gram-orthogonal continuation for P267 Target 1."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

import mpmath as mp

import score_marked_birth_path as base
import score_external_observer_transfer as transfer


METRICS = (
    "P4_perp_re", "P4_perp_im",
    "P4_far_D_re", "P4_far_D_im",
    "P4_far_S_re", "P4_far_S_im",
    "beta_first", "beta_second",
    "Gram_imag_first", "Gram_imag_second",
    "Gram_residual_re_first", "Gram_residual_im_first",
    "Gram_residual_re_second", "Gram_residual_im_second",
)


def orthogonalize(observables: dict[str, mp.mpf]) -> dict[str, mp.mpf | mp.mpc]:
    denominator = observables["Gram_abs_J_S2"]
    if denominator <= 0:
        raise ValueError("non-positive JS Gram denominator")
    gram = mp.mpc(
        observables["Gram_J_D_conj_J_S_re"],
        observables["Gram_J_D_conj_J_S_im"],
    )
    # The lattice sources have a real exact Gram.  Its imaginary entry is an
    # audit, not a fit direction.
    beta = mp.re(gram) / denominator
    far_d = mp.mpc(
        observables["connected_O_far_J_D_re"],
        observables["connected_O_far_J_D_im"],
    )
    far_s = mp.mpc(
        observables["connected_O_far_J_S_re"],
        observables["connected_O_far_J_S_im"],
    )
    return {
        "beta": beta,
        "gram_imag": mp.im(gram),
        "gram_residual": gram - beta * denominator,
        "far_D": far_d,
        "far_S": far_s,
        "perp": far_d - beta * far_s,
    }


def projected(first: Sequence[base.PathRow], second: Sequence[base.PathRow]) -> tuple[mp.mpf, dict[str, mp.mpf]]:
    center = base.intrinsic_center(first, second)
    left = orthogonalize(base.orientation_observables(first, center))
    right = orthogonalize(base.orientation_observables(second, center))
    leverage = base.cos4(first[0].a, first[0].b) - base.cos4(second[0].a, second[0].b)
    if leverage == 0:
        raise ValueError("zero H4 leverage")

    def p4(name: str) -> mp.mpc:
        return (left[name] - right[name]) / leverage

    perp, far_d, far_s = p4("perp"), p4("far_D"), p4("far_S")
    left_residual = left["gram_residual"]
    right_residual = right["gram_residual"]
    return center, {
        "P4_perp_re": mp.re(perp), "P4_perp_im": mp.im(perp),
        "P4_far_D_re": mp.re(far_d), "P4_far_D_im": mp.im(far_d),
        "P4_far_S_re": mp.re(far_s), "P4_far_S_im": mp.im(far_s),
        "beta_first": left["beta"], "beta_second": right["beta"],
        "Gram_imag_first": left["gram_imag"],
        "Gram_imag_second": right["gram_imag"],
        "Gram_residual_re_first": mp.re(left_residual),
        "Gram_residual_im_first": mp.im(left_residual),
        "Gram_residual_re_second": mp.re(right_residual),
        "Gram_residual_im_second": mp.im(right_residual),
    }


def score_prefix(prefix: Path) -> dict:
    groups = base.read_path(Path(str(prefix) + ".path.csv"))
    sizes = {key[0] for key in groups}
    if len(sizes) != 1:
        raise ValueError("score input contains multiple sizes")
    n = sizes.pop()
    batches = sorted(
        set(key[2] for key in groups if key[1] == "first")
        & set(key[2] for key in groups if key[1] == "second")
    )
    if len(batches) < 2:
        raise ValueError("delete-one score needs aligned batches")

    def combined(kept: Sequence[int], side: str) -> list[base.PathRow]:
        return base.combine([groups[(n, side, batch)] for batch in kept])

    center, point = projected(combined(batches, "first"), combined(batches, "second"))
    delete_centers = []
    delete_rows = []
    for omitted in batches:
        kept = [batch for batch in batches if batch != omitted]
        leave_center, leave_point = projected(
            combined(kept, "first"), combined(kept, "second")
        )
        delete_centers.append(leave_center)
        delete_rows.append(leave_point)
    covariance = base.covariance(delete_rows)
    standard_errors = {
        name: mp.sqrt(max(mp.mpf(0), covariance[index][index]))
        for index, name in enumerate(METRICS)
    }
    center_mean = mp.fsum(delete_centers) / len(delete_centers)
    center_se = mp.sqrt(
        mp.mpf(len(delete_centers) - 1) / len(delete_centers)
        * mp.fsum((value - center_mean) ** 2 for value in delete_centers)
    )
    return {
        "schema": "matching-one/external-observer-gram-orthogonal-single/v1",
        "N": n,
        "prefix": str(prefix),
        "batches": len(batches),
        "intrinsic_center": base._text(center),
        "intrinsic_center_delete_one_se": base._text(center_se),
        "metric_order": list(METRICS),
        "point": {name: base._text(point[name]) for name in METRICS},
        "standard_error": {
            name: base._text(standard_errors[name]) for name in METRICS
        },
        "delete_one_covariance": [
            [base._text(value) for value in row] for row in covariance
        ],
    }


def complex_point(report: dict) -> complex:
    return complex(float(report["point"]["P4_perp_re"]),
                   float(report["point"]["P4_perp_im"]))


def complex_covariance(report: dict) -> list[list[float]]:
    order = report["metric_order"]
    cov = report["delete_one_covariance"]
    indices = [order.index("P4_perp_re"), order.index("P4_perp_im")]
    return [[float(cov[i][j]) for j in indices] for i in indices]


def build(first_prefix: Path, second_prefix: Path) -> dict:
    first = score_prefix(first_prefix)
    second = score_prefix(second_prefix)
    if first["N"] >= second["N"]:
        raise ValueError("first prefix must be the smaller size")
    z1, z2 = complex_point(first), complex_point(second)
    c1, c2 = complex_covariance(first), complex_covariance(second)
    return {
        "schema": "matching-one/external-observer-gram-orthogonal-transfer/v1",
        "frozen_rule": "per-size/orientation beta=Re<JD,JS>/<|JS|^2>; delete-one recomputes root and beta; O and transfer never train beta",
        "sizes": [first["N"], second["N"]],
        "single_size": {str(first["N"]): first, str(second["N"]): second},
        "primary": {
            str(first["N"]): {
                "complex": [z1.real, z1.imag],
                "covariance_re_im": c1,
                "nonzero_mahalanobis_chi2_2d": transfer.mahalanobis(z1, c1),
            },
            str(second["N"]): {
                "complex": [z2.real, z2.imag],
                "covariance_re_im": c2,
                "nonzero_mahalanobis_chi2_2d": transfer.mahalanobis(z2, c2),
            },
            "transfer_second_over_first": transfer.ratio_summary(
                z1, z2, transfer.block_diagonal(c1, c2)
            ),
        },
        "sufficiency": {
            "coupling_and_transfer": "complete from existing per-batch path aggregates",
            "missing": "sum_abs_J_D2",
            "impact": "cannot score residual source norm or L2 energy fraction; coupling is unaffected",
        },
        "claim_boundary": "orthogonal in the recorded same-next-site Gram metric, not a CFT field inner product",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first-prefix", type=Path, required=True)
    parser.add_argument("--second-prefix", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--dps", type=int, default=50)
    args = parser.parse_args()
    mp.mp.dps = args.dps
    result = build(args.first_prefix, args.second_prefix)
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()

