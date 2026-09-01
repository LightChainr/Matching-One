#!/usr/bin/env python3
"""Profile the production (A_top, E_top) state plane along radial lineages.

For a parent observation ``y0`` and child observation ``y1`` the rank-one
null is

    E[y0] = mu,       E[y1] = lambda * mu,

with both ``mu in R^2`` and ``lambda in R`` nuisance parameters.  At fixed
``lambda`` the two-dimensional truth profiles out analytically.  The remaining
Mahalanobis discrepancy is

    (y1-lambda*y0)' (S1+lambda^2*S0)^-1 (y1-lambda*y0).

Its derivative is a degree-at-most-six polynomial divided by a strictly
positive denominator.  Enumerating all real polynomial roots and the compact
boundary ``|lambda|=infinity`` therefore gives a global, not local, minimum.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Iterable, Sequence

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
CROSSWALK = ROOT / "results/rank-plane-crosswalk/latest.json"

LINEAGES = (
    ("P49-N130-to-N170", "P49-N130", "P49-N170", "P49"),
    ("P43-N185-to-N265", "P43-N185", "P43-N265", "P43"),
    ("P50-N145-to-N290", "P50-N145", "P50-N290", "P50"),
    ("P57-N325-to-N425", "P57-N325", "P57-N425", "P57"),
)

PRIMARY_COVARIANCE = "covariance_intrinsic_center_first_order_influence"
SENSITIVITY_COVARIANCE = "covariance_fixed_center_exact_batch_estimator"
STATE_METRICS = ("P4_A_top", "P4_E_top")


def _number(value: mp.mpf | int | float, digits: int = 17) -> float:
    return float(mp.nstr(value, digits))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _poly_add(*polynomials: Sequence[mp.mpf]) -> list[mp.mpf]:
    size = max(len(poly) for poly in polynomials)
    return [
        mp.fsum(poly[index] for poly in polynomials if index < len(poly))
        for index in range(size)
    ]


def _poly_scale(poly: Sequence[mp.mpf], scale: mp.mpf) -> list[mp.mpf]:
    return [scale * value for value in poly]


def _poly_multiply(
    left: Sequence[mp.mpf], right: Sequence[mp.mpf]
) -> list[mp.mpf]:
    result = [mp.mpf(0)] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            result[i + j] += left_value * right_value
    return result


def _poly_derivative(poly: Sequence[mp.mpf]) -> list[mp.mpf]:
    return [index * poly[index] for index in range(1, len(poly))]


def _poly_value(poly: Sequence[mp.mpf], value: mp.mpf) -> mp.mpf:
    result = mp.mpf(0)
    for coefficient in reversed(poly):
        result = result * value + coefficient
    return result


def _inverse_2x2(matrix: Sequence[Sequence[mp.mpf]]) -> list[list[mp.mpf]]:
    determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    if determinant <= 0 or matrix[0][0] <= 0 or matrix[1][1] <= 0:
        raise ValueError("state covariance must be positive definite")
    return [
        [matrix[1][1] / determinant, -matrix[0][1] / determinant],
        [-matrix[1][0] / determinant, matrix[0][0] / determinant],
    ]


def _matvec(
    matrix: Sequence[Sequence[mp.mpf]], vector: Sequence[mp.mpf]
) -> list[mp.mpf]:
    return [mp.fsum(a * b for a, b in zip(row, vector)) for row in matrix]


def _quadratic_form_2x2(
    vector: Sequence[mp.mpf], covariance: Sequence[Sequence[mp.mpf]]
) -> mp.mpf:
    inverse = _inverse_2x2(covariance)
    transformed = _matvec(inverse, vector)
    return mp.fsum(a * b for a, b in zip(vector, transformed))


def _chi2_survival(chi2: mp.mpf, degrees_of_freedom: int) -> mp.mpf:
    if chi2 < 0 or degrees_of_freedom <= 0:
        raise ValueError("invalid chi-square request")
    half_df = mp.mpf(degrees_of_freedom) / 2
    return mp.gammainc(half_df, chi2 / 2, mp.inf) / mp.gamma(half_df)


def _state_row(dataset: dict[str, Any], covariance_key: str) -> dict[str, Any]:
    order = dataset["covariance_metric_order"]
    indices = [order.index(metric) for metric in STATE_METRICS]
    covariance = dataset[covariance_key]
    return {
        "id": dataset["id"],
        "N": dataset["N"],
        "orientation_pair": dataset["orientations"],
        "batches": dataset["batches"],
        "samples_per_orientation": dataset["samples_per_orientation"],
        "p0": dataset["p0"],
        "estimate": [mp.mpf(str(dataset["point"][metric])) for metric in STATE_METRICS],
        "covariance": [
            [mp.mpf(str(covariance[i][j])) for j in indices] for i in indices
        ],
    }


def _profile_polynomials(
    parent: dict[str, Any], child: dict[str, Any]
) -> tuple[list[mp.mpf], list[mp.mpf], list[mp.mpf]]:
    y0 = parent["estimate"]
    y1 = child["estimate"]
    s0 = parent["covariance"]
    s1 = child["covariance"]

    v00 = [s1[0][0], mp.mpf(0), s0[0][0]]
    v01 = [s1[0][1], mp.mpf(0), s0[0][1]]
    v11 = [s1[1][1], mp.mpf(0), s0[1][1]]
    d0 = [y1[0], -y0[0]]
    d1 = [y1[1], -y0[1]]

    numerator = _poly_add(
        _poly_multiply(v11, _poly_multiply(d0, d0)),
        _poly_scale(
            _poly_multiply(v01, _poly_multiply(d0, d1)), mp.mpf(-2)
        ),
        _poly_multiply(v00, _poly_multiply(d1, d1)),
    )
    denominator = _poly_add(
        _poly_multiply(v00, v11),
        _poly_scale(_poly_multiply(v01, v01), mp.mpf(-1)),
    )
    derivative_numerator = _poly_add(
        _poly_multiply(_poly_derivative(numerator), denominator),
        _poly_scale(
            _poly_multiply(numerator, _poly_derivative(denominator)),
            mp.mpf(-1),
        ),
    )
    scale = max(abs(value) for value in derivative_numerator)
    if scale == 0:
        return numerator, denominator, [mp.mpf(0)]
    # Structural cancellation removes the nominal degree-seven term.  Trim
    # only values far below the precision of the decimal input certificate.
    while (
        len(derivative_numerator) > 1
        and abs(derivative_numerator[-1]) < scale * mp.mpf("1e-50")
    ):
        derivative_numerator.pop()
    return numerator, denominator, derivative_numerator


def _fit_parent_truth(
    parent: dict[str, Any], child: dict[str, Any], lam: mp.mpf
) -> list[mp.mpf]:
    inverse_parent = _inverse_2x2(parent["covariance"])
    inverse_child = _inverse_2x2(child["covariance"])
    precision = [
        [
            inverse_parent[i][j] + lam * lam * inverse_child[i][j]
            for j in range(2)
        ]
        for i in range(2)
    ]
    right = [
        value_parent + lam * value_child
        for value_parent, value_child in zip(
            _matvec(inverse_parent, parent["estimate"]),
            _matvec(inverse_child, child["estimate"]),
        )
    ]
    return _matvec(_inverse_2x2(precision), right)


def profile_rank_one(parent: dict[str, Any], child: dict[str, Any]) -> dict[str, Any]:
    """Return the globally profiled one-degree-of-freedom rank-one test."""

    # Validate covariance before constructing the rational profile.
    _inverse_2x2(parent["covariance"])
    _inverse_2x2(child["covariance"])
    numerator, denominator, derivative = _profile_polynomials(parent, child)
    derivative_scale = max(abs(value) for value in derivative)
    if derivative_scale == 0:
        raise ArithmeticError("rank-one profile is constant")
    normalized = [value / derivative_scale for value in derivative]
    roots = mp.polyroots(
        normalized, maxsteps=2000, extraprec=80, error=False, asc=True
    )
    imaginary_tolerance = mp.power(10, -max(25, mp.mp.dps // 2))
    real_roots = sorted(
        mp.re(root) for root in roots if abs(mp.im(root)) < imaginary_tolerance
    )

    def discrepancy(lam: mp.mpf) -> mp.mpf:
        return _poly_value(numerator, lam) / _poly_value(denominator, lam)

    candidates = [(discrepancy(root), root) for root in real_roots]
    infinity_discrepancy = _quadratic_form_2x2(
        parent["estimate"], parent["covariance"]
    )
    candidates.append((infinity_discrepancy, mp.inf))
    chi2, lam = min(candidates, key=lambda row: row[0])
    if not mp.isfinite(lam):
        fitted_parent = [mp.mpf(0), mp.mpf(0)]
        fitted_child = child["estimate"]
    else:
        fitted_parent = _fit_parent_truth(parent, child, lam)
        fitted_child = [lam * value for value in fitted_parent]

    residuals = [abs(_poly_value(normalized, root)) for root in real_roots]
    return {
        "lambda": _number(lam) if mp.isfinite(lam) else "infinity",
        "fitted_parent_state": [_number(value) for value in fitted_parent],
        "fitted_child_state": [_number(value) for value in fitted_child],
        "min_chi2": _number(chi2),
        "degrees_of_freedom": 1,
        "p_value": _number(_chi2_survival(chi2, 1)),
        "optimizer_certificate": {
            "method": "all real roots of analytic rational-profile derivative plus |lambda|=infinity",
            "derivative_polynomial_degree": len(derivative) - 1,
            "number_of_complex_roots": len(roots) - len(real_roots),
            "real_stationary_points": [_number(value) for value in real_roots],
            "stationary_chi2": [_number(discrepancy(value)) for value in real_roots],
            "infinity_chi2": _number(infinity_discrepancy),
            "max_normalized_root_residual": _number(max(residuals, default=mp.mpf(0))),
        },
    }


def _determinant_diagnostic(parent: dict[str, Any], child: dict[str, Any]) -> dict[str, Any]:
    a0, e0 = parent["estimate"]
    a1, e1 = child["estimate"]
    determinant = a0 * e1 - e0 * a1
    gradient_parent = [e1, -a1]
    gradient_child = [-e0, a0]

    def covariance_quadratic(
        gradient: Sequence[mp.mpf], covariance: Sequence[Sequence[mp.mpf]]
    ) -> mp.mpf:
        return mp.fsum(
            gradient[i] * covariance[i][j] * gradient[j]
            for i in range(2)
            for j in range(2)
        )

    variance = covariance_quadratic(
        gradient_parent, parent["covariance"]
    ) + covariance_quadratic(gradient_child, child["covariance"])
    standard_error = mp.sqrt(variance)
    z_score = determinant / standard_error
    return {
        "determinant": _number(determinant),
        "delta_method_standard_error": _number(standard_error),
        "delta_method_z": _number(z_score),
        "boundary": "descriptive first-order diagnostic; the profiled Mahalanobis statistic is the primary rank-one test",
    }


def _zero_even_baseline(rows: Sequence[dict[str, Any]]) -> dict[str, Any]:
    contributions = [
        row["estimate"][1] ** 2 / row["covariance"][1][1] for row in rows
    ]
    chi2 = mp.fsum(contributions)
    return {
        "null": "P4_E_top=0 independently at every listed size",
        "chi2": _number(chi2),
        "degrees_of_freedom": len(rows),
        "p_value": _number(_chi2_survival(chi2, len(rows))),
        "per_size_chi2": {
            row["id"]: _number(value) for row, value in zip(rows, contributions)
        },
    }


def analyze_crosswalk(
    crosswalk: dict[str, Any], covariance_key: str, alpha: mp.mpf
) -> dict[str, Any]:
    datasets = {dataset["id"]: dataset for dataset in crosswalk["datasets"]}
    lineage_results = []
    all_rows: list[dict[str, Any]] = []
    for lineage_id, parent_id, child_id, dependency_group in LINEAGES:
        parent = _state_row(datasets[parent_id], covariance_key)
        child = _state_row(datasets[child_id], covariance_key)
        all_rows.extend([parent, child])
        rank_one = profile_rank_one(parent, child)
        zero_even = _zero_even_baseline([parent, child])
        lineage_results.append(
            {
                "lineage": lineage_id,
                "dependency_group": dependency_group,
                "parent": {
                    "id": parent["id"],
                    "N": parent["N"],
                    "orientation_pair": parent["orientation_pair"],
                    "batches": parent["batches"],
                    "samples_per_orientation": parent["samples_per_orientation"],
                    "p0": parent["p0"],
                    "estimate": [_number(value) for value in parent["estimate"]],
                    "covariance": [
                        [_number(value) for value in row]
                        for row in parent["covariance"]
                    ],
                },
                "child": {
                    "id": child["id"],
                    "N": child["N"],
                    "orientation_pair": child["orientation_pair"],
                    "batches": child["batches"],
                    "samples_per_orientation": child["samples_per_orientation"],
                    "p0": child["p0"],
                    "estimate": [_number(value) for value in child["estimate"]],
                    "covariance": [
                        [_number(value) for value in row]
                        for row in child["covariance"]
                    ],
                },
                "determinant_diagnostic": _determinant_diagnostic(parent, child),
                "rank_one_common_ray": {
                    **rank_one,
                    "decision_at_alpha": "eliminated"
                    if mp.mpf(str(rank_one["p_value"])) < alpha
                    else "survives",
                },
                "zero_even_baseline": {
                    **zero_even,
                    "decision_at_alpha": "eliminated"
                    if mp.mpf(str(zero_even["p_value"])) < alpha
                    else "survives",
                },
            }
        )

    joint_rank_one_chi2 = mp.fsum(
        mp.mpf(str(row["rank_one_common_ray"]["min_chi2"]))
        for row in lineage_results
    )
    joint_rank_one_p = _chi2_survival(joint_rank_one_chi2, len(lineage_results))
    zero_even = _zero_even_baseline(all_rows)
    return {
        "covariance_key": covariance_key,
        "metric_order": list(STATE_METRICS),
        "lineages": lineage_results,
        "joint_rank_one": {
            "combination": "sum across four block-diagonal production dependency groups",
            "dependency_groups": [row[3] for row in LINEAGES],
            "chi2": _number(joint_rank_one_chi2),
            "degrees_of_freedom": len(lineage_results),
            "p_value": _number(joint_rank_one_p),
            "decision_at_alpha": "eliminated" if joint_rank_one_p < alpha else "survives",
        },
        "joint_zero_even_baseline": {
            **zero_even,
            "decision_at_alpha": "eliminated"
            if mp.mpf(str(zero_even["p_value"])) < alpha
            else "survives",
        },
    }


def _source_hashes(crosswalk: dict[str, Any]) -> list[dict[str, Any]]:
    target_ids = {parent for _, parent, _, _ in LINEAGES} | {
        child for _, _, child, _ in LINEAGES
    }
    rows = []
    for dataset in crosswalk["datasets"]:
        if dataset["id"] not in target_ids:
            continue
        row = {"id": dataset["id"]}
        for kind in ("histogram", "moments"):
            path = ROOT / dataset[kind]
            row[kind] = str(path.relative_to(ROOT))
            row[f"{kind}_sha256"] = _sha256(path)
        rows.append(row)
    return sorted(rows, key=lambda row: row["id"])


def build_report(
    source: str = "json", alpha: str = "0.01", dps: int = 80
) -> dict[str, Any]:
    mp.mp.dps = max(50, dps)
    checked = json.loads(CROSSWALK.read_text(encoding="utf-8"))
    raw_matches_checked: bool | None = None
    if source == "raw":
        # Import lazily so synthetic/unit use does not pay the archive rebuild.
        sys.path.insert(0, str(ROOT / "scripts"))
        import rank_plane_crosswalk  # pylint: disable=import-outside-toplevel

        # The v1 crosswalk has a frozen 35-dps serialization contract.  Rebuild
        # under that contract, then restore the higher precision used by the
        # profile-root certificate below.
        requested_dps = mp.mp.dps
        mp.mp.dps = 35
        crosswalk = rank_plane_crosswalk.build_report()
        mp.mp.dps = requested_dps
        raw_matches_checked = crosswalk == checked
        if not raw_matches_checked:
            raise ValueError("raw same-batch reconstruction differs from checked crosswalk")
    elif source == "json":
        crosswalk = checked
    else:
        raise ValueError(f"unknown source {source}")

    alpha_value = mp.mpf(alpha)
    primary = analyze_crosswalk(crosswalk, PRIMARY_COVARIANCE, alpha_value)
    sensitivity = analyze_crosswalk(crosswalk, SENSITIVITY_COVARIANCE, alpha_value)
    return {
        "schema": "matching-one/etop-production-rank1-elimination/v1",
        "issues": [337, 370],
        "status": "production_archive_model_elimination",
        "alpha": _number(alpha_value),
        "source": {
            "mode": source,
            "checked_crosswalk": str(CROSSWALK.relative_to(ROOT)),
            "checked_crosswalk_sha256": _sha256(CROSSWALK),
            "raw_reconstruction_matches_checked_crosswalk": raw_matches_checked,
            "same_batch_archive_hashes": _source_hashes(crosswalk),
        },
        "estimand": {
            "state": "(P4_A_top, P4_E_top) at each archive's intrinsic matching center p0",
            "rank_one_null": "child=lambda*parent with unrestricted real lambda and unrestricted two-component parent truth",
            "cross_size_covariance": "block diagonal, as specified for the independent production size blocks",
            "primary_covariance": PRIMARY_COVARIANCE,
            "sensitivity_covariance": SENSITIVITY_COVARIANCE,
            "boundary": "intrinsic-center covariance is first-order influence; fixed-center covariance is exact for batch estimators conditional on the plug-in p0",
        },
        "primary": primary,
        "fixed_center_sensitivity": sensitivity,
        "scientific_card": [
            "MECHANISM SPACE: test whether radial evolution of the two-component A_top/E_top state needs a rotating second direction, rather than whether E_top exists.",
            "RESULT: all four common-ray rank-one lineages survive separately and jointly, while the global E_top=0 baseline is eliminated.",
            "NOT PROVED: survival does not identify an asymptotic field or make the center influence exact; no stored cross-size CRN covariance is claimed.",
            "OBSERVER-SECTOR-SOURCE-GEOMETRY: P4(A_top,E_top) | Alexander odd/even rank plane | threshold-rank source | four explicit Gaussian parent-child lineages.",
            "DEPENDENCY GROUPS: P49, P43, P50 and P57 are kept as four named production blocks; within each size A_top/E_top covariance is full.",
            "UPWEIGHT OBSERVATION: a future child with a frozen, large ray-rotation prediction or a stored cross-size CRN covariance can turn this surviving nuisance-ray test into a sharper discriminator.",
        ],
    }


def render_markdown(report: dict[str, Any]) -> str:
    primary = report["primary"]
    fixed = report["fixed_center_sensitivity"]
    lines = [
        "# Production E_top rank-one elimination",
        "",
        "## Answer",
        "",
        "The production state plane does **not** eliminate a common radial ray.",
        f"The four-lineage joint profile gives `chi2={primary['joint_rank_one']['chi2']:.6g}`",
        f"on `{primary['joint_rank_one']['degrees_of_freedom']}` df",
        f"(`p={primary['joint_rank_one']['p_value']:.6g}`). In contrast, the",
        f"zero-even baseline gives `chi2={primary['joint_zero_even_baseline']['chi2']:.6g}`",
        f"on `{primary['joint_zero_even_baseline']['degrees_of_freedom']}` df",
        f"(`p={primary['joint_zero_even_baseline']['p_value']:.6g}`) and is eliminated.",
        "Thus E_top is resolved, but its observed parent-child evolution is compatible",
        "with the same two-component ray rescaled by a lineage-specific nuisance lambda.",
        "",
        "## Lineage profiles",
        "",
        "| lineage | determinant | lambda | min chi2 / 1 df | p | rank-one | E_top=0 p / 2 df |",
        "|---|---:|---:|---:|---:|---|---:|",
    ]
    for row in primary["lineages"]:
        rank = row["rank_one_common_ray"]
        zero = row["zero_even_baseline"]
        lines.append(
            f"| {row['lineage']} | {row['determinant_diagnostic']['determinant']:.6g} "
            f"| {rank['lambda']:.6g} | {rank['min_chi2']:.6g} "
            f"| {rank['p_value']:.6g} | {rank['decision_at_alpha']} "
            f"| {zero['p_value']:.6g} |"
        )
    lines.extend(
        [
            "",
            "Each minimum is global: the scorer enumerates every real root of the",
            "degree-at-most-six analytic profile derivative and the compact boundary",
            "`|lambda|=infinity`. Parent truth and lambda are nuisance parameters, leaving",
            "one degree of freedom per lineage.",
            "",
            "## Covariance and sensitivity",
            "",
            f"Primary covariance: `{report['estimand']['primary_covariance']}`.",
            "It is the full same-batch A_top/E_top covariance with the displayed",
            "first-order influence correction for each fitted matching center.",
            f"The fixed-center exact-batch sensitivity gives joint rank-one `chi2={fixed['joint_rank_one']['chi2']:.6g}`",
            f"(`p={fixed['joint_rank_one']['p_value']:.6g}`), so the decision is unchanged.",
            "Cross-size blocks are block diagonal by the production independence contract;",
            "the four dependency-group labels remain explicit rather than being treated as",
            "eight unrelated evidence rows.",
            "",
            "## Scientific card",
            "",
        ]
    )
    lines.extend(f"- {item}" for item in report["scientific_card"])
    lines.extend(
        [
            "",
            "## Reproduction",
            "",
            "```bash",
            "python3 scripts/etop_rank1_elimination.py --source raw --format json --output results/etop-rank1-elimination/latest.json",
            "python3 scripts/etop_rank1_elimination.py --format markdown --output results/etop-rank1-elimination/REPORT.md",
            "python3 -m unittest discover -s tests -p 'test_etop_rank1_elimination.py'",
            "```",
            "",
            "The first command rebuilds the rank-plane values and covariance from the",
            "same-batch raw histograms/moments, checks equality with the committed",
            "crosswalk, and records SHA-256 hashes for all 16 input files.",
            "",
            "## Claim boundary",
            "",
            report["estimand"]["boundary"] + ". The rank-one result is a model",
            "survival statement at these four production edges, not an asymptotic field",
            "identification. The E_top=0 rejection establishes a resolved companion",
            "direction, not that the companion evolves independently of A_top.",
            "",
        ]
    )
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", choices=("json", "raw"), default="json")
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--alpha", default="0.01")
    parser.add_argument("--dps", type=int, default=80)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    report = build_report(source=args.source, alpha=args.alpha, dps=args.dps)
    payload = (
        json.dumps(report, indent=2, sort_keys=True) + "\n"
        if args.format == "json"
        else render_markdown(report)
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
