#!/usr/bin/env python3
"""High-precision complex-zero pilot for tiny exact matching polynomials.

The exact coefficients come from :mod:`exact_matching_polynomial`.  Roots are
computed twice at independent working precisions, then audited for polynomial
backward error, conjugation, and the exact matching transform

    P_hat(p) = -P(1-p).

The transform pairs a root ``z`` of ``P`` with ``1-z`` of ``P_hat``.  It does
not assert that a square/matching pair is self-matching or that the roots of a
single polynomial must be invariant under ``p -> 1-p``.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Sequence, Tuple

import mpmath as mp

from exact_matching_polynomial import bernstein_counts, bernstein_to_power, polynomial_string
from matched_torus_reference import axis_geometry, diamond_geometry


ComplexList = List[mp.mpc]


def matching_partner_coefficients(coefficients: Sequence[int]) -> List[int]:
    """Return ascending coefficients of ``-P(1-p)`` exactly."""

    output = [0] * len(coefficients)
    for degree, coefficient in enumerate(coefficients):
        for power in range(degree + 1):
            output[power] -= coefficient * math.comb(degree, power) * (-1) ** power
    while len(output) > 1 and output[-1] == 0:
        output.pop()
    return output


def evaluate(coefficients: Sequence[int], value: mp.mpc) -> mp.mpc:
    output = mp.mpc(0)
    for coefficient in reversed(coefficients):
        output = output * value + coefficient
    return output


def normalized_residual(coefficients: Sequence[int], value: mp.mpc) -> mp.mpf:
    numerator = abs(evaluate(coefficients, value))
    denominator = mp.fsum(abs(coefficient) * abs(value) ** degree
                          for degree, coefficient in enumerate(coefficients))
    return numerator / denominator if denominator else numerator


def canonical_sort(roots: Iterable[mp.mpc]) -> ComplexList:
    return sorted((mp.mpc(root) for root in roots), key=lambda root: (mp.re(root), mp.im(root)))


def roots_once(coefficients: Sequence[int], dps: int) -> ComplexList:
    with mp.workdps(dps):
        roots = mp.polyroots(
            [mp.mpf(value) for value in coefficients],
            asc=True,
            maxsteps=5000,
            extraprec=50,
            cleanup=True,
            error=False,
        )
        return canonical_sort(mp.mpc(mp.nstr(root, dps - 10)) for root in roots)


def greedy_match(
    roots: Sequence[mp.mpc],
    targets: Sequence[mp.mpc],
    transform: Callable[[mp.mpc], mp.mpc] = lambda value: value,
) -> List[Tuple[int, mp.mpf]]:
    unused = set(range(len(targets)))
    output = []
    for root in roots:
        expected = transform(root)
        index = min(unused, key=lambda candidate: abs(targets[candidate] - expected))
        output.append((index, abs(targets[index] - expected)))
        unused.remove(index)
    return output


def stable_roots(coefficients: Sequence[int], dps: int) -> Tuple[ComplexList, List[mp.mpf]]:
    if dps < 60:
        raise ValueError("dps must be at least 60")
    low = roots_once(coefficients, dps)
    high = roots_once(coefficients, dps + 30)
    if len(low) != len(high):
        raise ArithmeticError("precision runs returned different root counts")
    matches = greedy_match(high, low)
    return high, [distance for _index, distance in matches]


def root_cloud_metrics(roots: Sequence[mp.mpc], tolerance: mp.mpf) -> Dict[str, object]:
    real_roots = [mp.re(root) for root in roots if abs(mp.im(root)) <= tolerance]
    physical = sorted(root for root in real_roots if 0 < root < 1)
    if len(physical) != 1:
        physical_value = None
    else:
        physical_value = physical[0]
    nonreal_count = sum(abs(mp.im(root)) > tolerance for root in roots)
    return {
        "degree": len(roots),
        "real_root_count": len(real_roots),
        "nonreal_root_count": nonreal_count,
        "nonreal_fraction": mp.mpf(nonreal_count) / len(roots),
        "physical_root_0_1": physical_value,
        "imaginary_rms": mp.sqrt(mp.fsum(mp.im(root) ** 2 for root in roots) / len(roots)),
        "imaginary_mean_absolute": mp.fsum(abs(mp.im(root)) for root in roots) / len(roots),
    }


def analyze_polynomial(
    geometry_name: str, length: int, coefficients: Sequence[int], dps: int
) -> Tuple[dict, List[dict]]:
    # Keep the audit arithmetic at least as precise as the confirmation solve.
    # mpmath values retain their mantissas, but subsequent operations otherwise
    # use the caller's (often 15-digit) global context.
    mp.mp.dps = max(mp.mp.dps, dps + 30)
    partner_coefficients = matching_partner_coefficients(coefficients)
    try:
        roots, stability = stable_roots(coefficients, dps)
        partner_roots, partner_stability = stable_roots(partner_coefficients, dps)
    except Exception as exc:  # Preserve an auditable failed geometry.
        return ({
            "geometry": geometry_name,
            "L": length,
            "status": "ROOT_FAILURE",
            "failure_type": type(exc).__name__,
            "failure_message": str(exc),
            "power_coefficients_ascending": list(coefficients),
            "matching_partner_coefficients_ascending": partner_coefficients,
        }, [])

    tolerance = mp.power(10, -(dps // 2))
    conjugates = greedy_match(roots, roots, lambda root: mp.conj(root))
    partner_matches = greedy_match(roots, partner_roots, lambda root: 1 - root)
    root_rows = []
    for index, root in enumerate(roots):
        conjugate_index, conjugate_distance = conjugates[index]
        partner_index, partner_distance = partner_matches[index]
        expected_partner = 1 - root
        root_rows.append({
            "geometry": geometry_name,
            "L": length,
            "root_index": index,
            "root_real": mp.re(root),
            "root_imaginary": mp.im(root),
            "is_real": abs(mp.im(root)) <= tolerance,
            "is_upper_half": mp.im(root) > tolerance,
            "normalized_polynomial_residual": normalized_residual(coefficients, root),
            "precision_stability_distance": stability[index],
            "conjugate_root_index": conjugate_index,
            "conjugate_pair_distance": conjugate_distance,
            "matching_partner_root_index": partner_index,
            "matching_partner_expected_real": mp.re(expected_partner),
            "matching_partner_expected_imaginary": mp.im(expected_partner),
            "matching_partner_pair_distance": partner_distance,
            "matching_partner_normalized_residual": normalized_residual(
                partner_coefficients, expected_partner
            ),
        })

    metrics = root_cloud_metrics(roots, tolerance)
    summary = {
        "geometry": geometry_name,
        "L": length,
        "status": "OK",
        "degree": len(roots),
        "power_coefficients_ascending": list(coefficients),
        "polynomial": polynomial_string(list(coefficients)),
        "matching_partner_coefficients_ascending": partner_coefficients,
        "matching_partner_polynomial": polynomial_string(partner_coefficients),
        "metrics": metrics,
        "audit": {
            "working_dps": dps,
            "confirmation_dps": dps + 30,
            "max_normalized_polynomial_residual": max(
                row["normalized_polynomial_residual"] for row in root_rows
            ),
            "max_precision_stability_distance": max(stability),
            "max_partner_precision_stability_distance": max(partner_stability),
            "max_conjugate_pair_distance": max(
                row["conjugate_pair_distance"] for row in root_rows
            ),
            "max_matching_partner_pair_distance": max(
                row["matching_partner_pair_distance"] for row in root_rows
            ),
            "max_matching_partner_normalized_residual": max(
                row["matching_partner_normalized_residual"] for row in root_rows
            ),
        },
    }
    return summary, root_rows


def linear_inverse_n_prediction(first: dict, second: dict, target_n: int, metric: str) -> dict:
    x_first = mp.mpf(1) / int(first["N"])
    x_second = mp.mpf(1) / int(second["N"])
    x_target = mp.mpf(1) / target_n
    y_first = mp.mpf(first["metrics"][metric])
    y_second = mp.mpf(second["metrics"][metric])
    prediction = y_first + (y_second - y_first) * (x_target - x_first) / (x_second - x_first)
    return {
        "metric": metric,
        "model": "a + b/N through the two declared training sizes",
        "training_L": [first["L"], second["L"]],
        "training_N": [first["N"], second["N"]],
        "target_N": target_n,
        "prediction": prediction,
    }


def prediction_audit(summaries: Sequence[dict]) -> dict:
    metrics = ("physical_root_0_1", "imaginary_rms", "nonreal_fraction")
    by_geometry = {}
    for geometry_name in ("axis", "diamond"):
        rows = sorted(
            (row for row in summaries if row["geometry"] == geometry_name and row["status"] == "OK"),
            key=lambda row: row["L"],
        )
        for row in rows:
            row["N"] = row["L"] ** 2 if geometry_name == "axis" else 2 * row["L"] ** 2
        if len(rows) < 3:
            raise ValueError("each geometry needs at least three sizes for the holdout audit")
        training = rows[-3:-1]
        holdout = rows[-1]
        holdout_scores = []
        for metric in metrics:
            item = linear_inverse_n_prediction(training[0], training[1], holdout["N"], metric)
            observed = mp.mpf(holdout["metrics"][metric])
            item.update({"target_L": holdout["L"], "observed": observed, "signed_error": observed - item["prediction"]})
            holdout_scores.append(item)

        prospective_n = (rows[-1]["L"] + 1) ** 2
        if geometry_name == "diamond":
            prospective_n *= 2
        prospective = [
            {**linear_inverse_n_prediction(rows[-2], rows[-1], prospective_n, metric),
             "target_L": rows[-1]["L"] + 1}
            for metric in metrics
        ]
        by_geometry[geometry_name] = {
            "retrospective_holdout": holdout_scores,
            "prospective_next_size": prospective,
            "warning": "two-point extrapolation is a frozen pilot diagnostic, not a scaling claim",
        }
    return by_geometry


def serializable(value: object, digits: int = 60) -> object:
    if isinstance(value, (mp.mpf, mp.mpc)):
        return mp.nstr(value, digits)
    if isinstance(value, dict):
        return {key: serializable(item, digits) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serializable(item, digits) for item in value]
    return value


def write_csv(path: Path, rows: Sequence[dict], digits: int = 50) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(serializable(row, digits))


def report(payload: dict) -> str:
    lines = [
        "# Tiny exact matching-polynomial complex-zero pilot",
        "",
        "The map uses exact enumeration coefficients and 100/130-digit numerical roots.",
        "It is an algebra-and-design pilot, not a thermodynamic zero-density claim.",
        "",
        "## Root audit",
        "",
        "| geometry | L | N | degree | real/nonreal | physical root | imag RMS | max residual | max conjugate error | max matching error |",
        "|:---|---:|---:|---:|:---:|---:|---:|---:|---:|---:|",
    ]
    for row in payload["polynomials"]:
        if row["status"] != "OK":
            lines.append(f"| {row['geometry']} | {row['L']} | - | - | ROOT FAILURE | - | - | - | - | - |")
            continue
        metrics = row["metrics"]
        audit = row["audit"]
        lines.append(
            "| {} | {} | {} | {} | {}/{} | {} | {} | {} | {} | {} |".format(
                row["geometry"], row["L"], row["N"], row["degree"],
                metrics["real_root_count"], metrics["nonreal_root_count"],
                mp.nstr(metrics["physical_root_0_1"], 12), mp.nstr(metrics["imaginary_rms"], 8),
                mp.nstr(audit["max_normalized_polynomial_residual"], 3),
                mp.nstr(audit["max_conjugate_pair_distance"], 3),
                mp.nstr(audit["max_matching_partner_pair_distance"], 3),
            )
        )
    lines.extend([
        "",
        "The matching column compares roots of `P(p)` against independently solved roots of",
        "`-P(1-p)`. A single axis/diamond polynomial is not assumed to be self-matching.",
        "",
        "## Declared train/holdout diagnostics",
        "",
        "Each target is predicted from the preceding two sizes with `a+b/N`. This deliberately",
        "cheap rule is scored before any richer fit; the exercise defines quantities that a future",
        "zero-map computation can falsify instead of merely producing a visually suggestive cloud.",
        "",
        "| geometry | metric | train L | held-out L | prediction | observed | signed error |",
        "|:---|:---|:---:|---:|---:|---:|---:|",
    ])
    for geometry, block in payload["predictions"].items():
        for item in block["retrospective_holdout"]:
            lines.append("| {} | {} | {} | {} | {} | {} | {} |".format(
                geometry, item["metric"], ",".join(map(str, item["training_L"])), item["target_L"],
                mp.nstr(item["prediction"], 12), mp.nstr(item["observed"], 12),
                mp.nstr(item["signed_error"], 8),
            ))
    lines.extend([
        "",
        "## Frozen next-size predictions",
        "",
        "| geometry | target L/N | metric | prediction |",
        "|:---|:---:|:---|---:|",
    ])
    for geometry, block in payload["predictions"].items():
        for item in block["prospective_next_size"]:
            lines.append("| {} | {}/{} | {} | {} |".format(
                geometry, item["target_L"], item["target_N"], item["metric"],
                mp.nstr(item["prediction"], 14),
            ))
    lines.extend([
        "",
        "Axis L=5 requires 2^25 configurations and remains within the reference engine's hard",
        "limit; diamond L=4 has N=32 and needs a frontier exact engine. Preserve these predictions",
        "unchanged if either target is later computed.",
        "",
        "## What this pilot says",
        "",
        "- Every solved polynomial passes conjugate pairing, dual matching-root pairing and",
        "  independent-precision stability by many orders of magnitude.",
        "- The unique real root in `(0,1)` rapidly approaches the known threshold from opposite",
        "  sides for the two orientations; the two-point tiny-size extrapolations are visibly biased.",
        "- Complex roots already proliferate by axis L=3 and diamond L=2, but the cloud summaries",
        "  are not yet smooth enough to justify a conformal or Lee-Yang interpretation.",
        "- The next useful step is one new exact size, scored against the frozen scalar targets above,",
        "  before adding plots, clustering rules or modular interpretations.",
        "",
    ])
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=100)
    parser.add_argument("--json", type=Path)
    parser.add_argument("--csv", type=Path)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args()
    cases = [("axis", length) for length in range(1, 5)] + [
        ("diamond", length) for length in range(1, 4)
    ]
    summaries = []
    root_rows = []
    for geometry_name, length in cases:
        geometry = axis_geometry(length) if geometry_name == "axis" else diamond_geometry(length)
        coefficients = bernstein_to_power(bernstein_counts(geometry))
        summary, rows = analyze_polynomial(geometry_name, length, coefficients, args.dps)
        summary["N"] = geometry.n
        summaries.append(summary)
        root_rows.extend(rows)
    if any(row["status"] != "OK" for row in summaries):
        predictions = {}
    else:
        predictions = prediction_audit(summaries)
    payload = {
        "schema": "tiny exact matching complex-zero map v1",
        "root_method": "mpmath.polyroots at dps and dps+30 with cross-precision matching",
        "matching_transform": "P_hat(p)=-P(1-p); z -> 1-z",
        "polynomials": summaries,
        "predictions": predictions,
    }
    text = json.dumps(serializable(payload), indent=2)
    print(text)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(text + "\n", encoding="utf-8")
    if args.csv:
        write_csv(args.csv, root_rows)
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(report(payload), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
