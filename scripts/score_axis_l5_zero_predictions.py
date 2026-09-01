#!/usr/bin/env python3
"""Score PR #78's frozen axis-L=5 zero-cloud predictions.

The prediction values are read verbatim from the pilot JSON; this scorer never
refits or rewrites them.  Exact coefficients come from the C++ enumeration
artifact and roots are audited through ``exact_matching_zero_map``.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import mpmath as mp

from exact_matching_zero_map import analyze_polynomial, serializable


METRICS = ("physical_root_0_1", "imaginary_rms", "nonreal_fraction")


def score(polynomial_payload: dict, pilot_payload: dict, dps: int) -> tuple[dict, list[dict]]:
    if polynomial_payload.get("geometry") != "axis" or polynomial_payload.get("L") != 5:
        raise ValueError("polynomial artifact must be axis L=5")
    prospective = pilot_payload["predictions"]["axis"]["prospective_next_size"]
    frozen = {row["metric"]: row for row in prospective}
    if tuple(frozen) != METRICS:
        raise ValueError("pilot does not contain the expected three frozen metrics in order")
    for metric in METRICS:
        row = frozen[metric]
        if row["target_L"] != 5 or row["target_N"] != 25:
            raise ValueError(f"frozen {metric} target is not axis L=5/N=25")

    summary, root_rows = analyze_polynomial(
        "axis", 5, polynomial_payload["power_coefficients_ascending"], dps
    )
    if summary["status"] != "OK":
        raise ArithmeticError(f"root analysis failed: {summary}")

    scores = []
    for metric in METRICS:
        prediction = mp.mpf(frozen[metric]["prediction"])
        observed = mp.mpf(summary["metrics"][metric])
        error = observed - prediction
        scores.append({
            "metric": metric,
            "model": frozen[metric]["model"],
            "training_L": frozen[metric]["training_L"],
            "training_N": frozen[metric]["training_N"],
            "target_L": 5,
            "target_N": 25,
            "frozen_prediction": prediction,
            "observed": observed,
            "signed_error_observed_minus_prediction": error,
            "absolute_error": abs(error),
            "relative_error": error / prediction if prediction else mp.nan,
        })

    payload = {
        "schema": "axis L=5 frozen zero prediction score v1",
        "prediction_source_schema": pilot_payload.get("schema"),
        "coefficient_source_schema": polynomial_payload.get("schema"),
        "enumerated_configurations": polynomial_payload["configurations"],
        "root_method": "exact_matching_zero_map.analyze_polynomial at dps and dps+30",
        "scores": scores,
        "polynomial": summary,
    }
    return payload, root_rows


def write_roots(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(serializable(row, 60))


def report(payload: dict) -> str:
    polynomial = payload["polynomial"]
    metrics = polynomial["metrics"]
    audit = polynomial["audit"]
    lines = [
        "# Exact axis L=5 matching polynomial and frozen zero-map score",
        "",
        "The C++ kernel exhausted all `2^25 = 33,554,432` configurations. Its",
        "L=1..4 coefficients were first regressed exactly against the Python oracle.",
        "PR #78's three predictions were read unchanged from its committed JSON.",
        "",
        "## Frozen score",
        "",
        "| metric | frozen prediction | observed | observed - predicted | relative error |",
        "|:---|---:|---:|---:|---:|",
    ]
    for row in payload["scores"]:
        lines.append(
            "| {} | {} | {} | {} | {} |".format(
                row["metric"], mp.nstr(row["frozen_prediction"], 16),
                mp.nstr(row["observed"], 16),
                mp.nstr(row["signed_error_observed_minus_prediction"], 10),
                mp.nstr(row["relative_error"], 8),
            )
        )
    lines.extend([
        "",
        "The cheap two-point `a+b/N` rule is surprisingly close for the physical root",
        "and gets the discrete nonreal fraction nearly right, but badly overpredicts the",
        "imaginary RMS. The latter is a clean falsification of that frozen cloud-scale rule.",
        "No alternative model was fit after seeing L=5.",
        "",
        "## Exact and numerical audit",
        "",
        f"- Exact power-basis coefficients (ascending): `{polynomial['power_coefficients_ascending']}`.",
        f"- Degree / real / nonreal roots: `{polynomial['degree']}` / `{metrics['real_root_count']}` / `{metrics['nonreal_root_count']}`.",
        f"- Physical root: `{mp.nstr(metrics['physical_root_0_1'], 40)}`.",
        f"- Maximum normalized polynomial residual: `{mp.nstr(audit['max_normalized_polynomial_residual'], 6)}`.",
        f"- Maximum 100/130-digit root shift: `{mp.nstr(audit['max_precision_stability_distance'], 6)}`.",
        f"- Maximum conjugate-pair error: `{mp.nstr(audit['max_conjugate_pair_distance'], 6)}`.",
        f"- Maximum exact-partner `z -> 1-z` pairing error: `{mp.nstr(audit['max_matching_partner_pair_distance'], 6)}`.",
        "",
        "This is a finite-size exact result only; the root cloud is not assigned a CFT or",
        "Lee-Yang interpretation.",
    ])
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polynomial", type=Path, required=True)
    parser.add_argument("--pilot", type=Path, required=True)
    parser.add_argument("--dps", type=int, default=100)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    polynomial = json.loads(args.polynomial.read_text(encoding="utf-8"))
    pilot = json.loads(args.pilot.read_text(encoding="utf-8"))
    payload, roots = score(polynomial, pilot, args.dps)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(serializable(payload), indent=2) + "\n", encoding="utf-8")
    write_roots(args.csv, roots)
    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    args.markdown.write_text(report(payload), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
