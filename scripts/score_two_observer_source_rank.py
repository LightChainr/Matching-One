#!/usr/bin/env python3
"""Score the frozen 2-observer x 2-source projective-rank lane.

Rows are the bulk Euler observer O_far and the separated two-orbit arm
projection O_sep4.  Columns are the same-birth-path sources JD_perp and JS.
The primary null is rank one over C, equivalently det(C)=0.  No individual
matrix entry is interpreted as a continuum-field identification.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

import mpmath as mp

import score_marked_birth_path as base


METRICS = (
    "C_far_perp_re", "C_far_perp_im",
    "C_far_S_re", "C_far_S_im",
    "C_sep4_perp_re", "C_sep4_perp_im",
    "C_sep4_S_re", "C_sep4_S_im",
    "det_re", "det_im", "normalized_wedge",
    "beta_first", "beta_second",
)


def orientation_lane(observables: dict[str, mp.mpf]) -> dict[str, mp.mpf | mp.mpc]:
    denominator = observables["Gram_abs_J_S2"]
    if denominator <= 0:
        raise ValueError("non-positive JS Gram denominator")
    beta = observables["Gram_J_D_conj_J_S_re"] / denominator

    def coupling(observer: str, source: str) -> mp.mpc:
        return mp.mpc(
            observables[f"connected_{observer}_{source}_re"],
            observables[f"connected_{observer}_{source}_im"],
        )

    far_d = coupling("O_far", "J_D")
    far_s = coupling("O_far", "J_S")
    sep_d = coupling("O_sep4", "J_D")
    sep_s = coupling("O_sep4", "J_S")
    return {
        "beta": beta,
        "far_perp": far_d - beta * far_s,
        "far_S": far_s,
        "sep4_perp": sep_d - beta * sep_s,
        "sep4_S": sep_s,
    }


def projected(first: Sequence[base.PathRow], second: Sequence[base.PathRow]) -> tuple[mp.mpf, dict[str, mp.mpf]]:
    center = base.intrinsic_center(first, second)
    left = orientation_lane(base.orientation_observables(first, center))
    right = orientation_lane(base.orientation_observables(second, center))
    leverage = base.cos4(first[0].a, first[0].b) - base.cos4(second[0].a, second[0].b)
    if leverage == 0:
        raise ValueError("zero H4 leverage")

    def p4(name: str) -> mp.mpc:
        return (left[name] - right[name]) / leverage

    far_perp = p4("far_perp")
    far_s = p4("far_S")
    sep_perp = p4("sep4_perp")
    sep_s = p4("sep4_S")
    determinant = far_perp * sep_s - far_s * sep_perp
    far_norm = mp.sqrt(abs(far_perp) ** 2 + abs(far_s) ** 2)
    sep_norm = mp.sqrt(abs(sep_perp) ** 2 + abs(sep_s) ** 2)
    wedge = abs(determinant) / (far_norm * sep_norm) if far_norm and sep_norm else mp.nan
    return center, {
        "C_far_perp_re": mp.re(far_perp),
        "C_far_perp_im": mp.im(far_perp),
        "C_far_S_re": mp.re(far_s),
        "C_far_S_im": mp.im(far_s),
        "C_sep4_perp_re": mp.re(sep_perp),
        "C_sep4_perp_im": mp.im(sep_perp),
        "C_sep4_S_re": mp.re(sep_s),
        "C_sep4_S_im": mp.im(sep_s),
        "det_re": mp.re(determinant),
        "det_im": mp.im(determinant),
        "normalized_wedge": wedge,
        "beta_first": left["beta"],
        "beta_second": right["beta"],
    }


def covariance(rows: Sequence[dict[str, mp.mpf]]) -> list[list[mp.mpf]]:
    count = len(rows)
    means = {name: mp.fsum(row[name] for row in rows) / count for name in METRICS}
    factor = mp.mpf(count - 1) / count
    return [[
        factor * mp.fsum(
            (row[left] - means[left]) * (row[right] - means[right]) for row in rows
        )
        for right in METRICS
    ] for left in METRICS]


def determinant_chi2(point: dict[str, mp.mpf], cov: list[list[mp.mpf]]) -> mp.mpf:
    i = METRICS.index("det_re")
    j = METRICS.index("det_im")
    a, b, c = cov[i][i], cov[i][j], cov[j][j]
    det = a * c - b * b
    if det <= 0:
        return mp.nan
    x, y = point["det_re"], point["det_im"]
    return (c * x * x - 2 * b * x * y + a * y * y) / det


def score_prefix(prefix: Path) -> dict:
    groups = base.read_path(Path(str(prefix) + ".path.csv"))
    sizes = {key[0] for key in groups}
    if len(sizes) != 1:
        raise ValueError("score input contains multiple sizes")
    n = sizes.pop()
    batches = sorted(
        {key[2] for key in groups if key[1] == "first"}
        & {key[2] for key in groups if key[1] == "second"}
    )
    if len(batches) < 2:
        raise ValueError("rank score needs at least two aligned batches")

    def combined(kept: Sequence[int], side: str) -> list[base.PathRow]:
        return base.combine([groups[(n, side, batch)] for batch in kept])

    center, point = projected(combined(batches, "first"), combined(batches, "second"))
    delete_rows = []
    delete_centers = []
    for omitted in batches:
        kept = [batch for batch in batches if batch != omitted]
        one_center, one_point = projected(
            combined(kept, "first"), combined(kept, "second")
        )
        delete_centers.append(one_center)
        delete_rows.append(one_point)
    cov = covariance(delete_rows)
    standard_error = {
        name: mp.sqrt(max(mp.mpf(0), cov[index][index]))
        for index, name in enumerate(METRICS)
    }
    center_mean = mp.fsum(delete_centers) / len(delete_centers)
    center_se = mp.sqrt(
        mp.mpf(len(delete_centers) - 1) / len(delete_centers)
        * mp.fsum((value - center_mean) ** 2 for value in delete_centers)
    )
    return {
        "schema": "matching-one/two-observer-two-source-rank-score/v1",
        "N": n,
        "prefix": str(prefix),
        "batches": len(batches),
        "intrinsic_center": base._text(center),
        "intrinsic_center_delete_one_se": base._text(center_se),
        "matrix_rows": ["O_far", "O_sep4_axis_minus_diagonal"],
        "matrix_columns": ["JD_perp", "JS"],
        "orthogonalization_rule": "per orientation beta=Re<JD,JS>/<|JS|^2>, recomputed in every delete-one replicate",
        "primary_null": "complex coupling matrix has rank one, equivalently det(C)=0",
        "metric_order": list(METRICS),
        "point": {name: base._text(point[name]) for name in METRICS},
        "standard_error": {name: base._text(standard_error[name]) for name in METRICS},
        "determinant_nonzero_mahalanobis_chi2_2d": base._text(
            determinant_chi2(point, cov)
        ),
        "delete_one_covariance": [
            [base._text(value) for value in row] for row in cov
        ],
        "claim_boundary": "individual nonzero entries are not field labels; only projective rank is primary",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--dps", type=int, default=50)
    args = parser.parse_args()
    mp.mp.dps = args.dps
    result = score_prefix(args.prefix)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
