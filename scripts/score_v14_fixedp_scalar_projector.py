#!/usr/bin/env python3
"""Project the orientation-independent matching-odd scalar from fixed-p pairs.

Consumes the standard `gaussian_orientation_mc` analysis CSV. For a same-N
orientation pair with c_i=cos(4 theta_i) and matching-function estimates M_i,

    M0 = (c1*M2 - c2*M1)/(c1-c2)

cancels a pure H4 contribution and retains the orientation-independent sector.
The covariance between M1 and M2 is reconstructed exactly from the reported
standard errors of M1, M2, and their difference:

    cov12 = (var1 + var2 - var(M1-M2))/2.

The V_<1,4> candidate predicts M0(pc) ~ N^-25/8.  A fixed p_ref sufficiently
close to pc adds a common thermal displacement and must therefore be reported
explicitly; this script is a discovery/power diagnostic, not a proof of the
operator assignment.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def project_row(row: dict) -> dict:
    c1 = float(row["cos4_first"])
    c2 = float(row["cos4_second"])
    dc = c1 - c2
    if abs(dc) < 1e-15:
        raise ValueError("vanishing DeltaCos4")
    m1 = float(row["first_estimate"])
    m2 = float(row["second_estimate"])
    s1 = float(row["first_batch_se"])
    s2 = float(row["second_batch_se"])
    sd = float(row["difference_batch_se"])
    cov = (s1 * s1 + s2 * s2 - sd * sd) / 2.0
    w1 = -c2 / dc
    w2 = c1 / dc
    value = w1 * m1 + w2 * m2
    variance = w1 * w1 * s1 * s1 + w2 * w2 * s2 * s2 + 2.0 * w1 * w2 * cov
    if variance < -1e-20:
        raise ValueError("reconstructed scalar variance is negative")
    se = math.sqrt(max(variance, 0.0))
    n = int(row["N"])
    scale = n ** (25.0 / 8.0)
    return {
        "N": n,
        "first": [int(row["a1"]), int(row["b1"])],
        "second": [int(row["a2"]), int(row["b2"])],
        "cos4_first": c1,
        "cos4_second": c2,
        "delta_cos4": dc,
        "M_first": m1,
        "M_second": m2,
        "M_scalar_H4_null": value,
        "M_scalar_se": se,
        "M_scalar_z": value / se if se else None,
        "within_pair_correlation": cov / (s1 * s2),
        "N25_8_scaled_scalar": value * scale,
        "N25_8_scaled_se": se * scale,
    }


def load_rows(path: Path, channel: str) -> list[dict]:
    selected = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["channel"] == channel and row["sector"] == "matching_function":
                selected.append(project_row(row))
    if not selected:
        raise ValueError("no matching_function rows found")
    return sorted(selected, key=lambda row: row["N"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("analysis_csv", type=Path)
    parser.add_argument("--channel", default="direction_1")
    parser.add_argument("--p-ref", type=float, default=0.592746050790)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    rows = load_rows(args.analysis_csv, args.channel)
    payload = {
        "format_version": 1,
        "classification": "retrospective discovery/power diagnostic",
        "p_ref": args.p_ref,
        "channel": args.channel,
        "hypothesis": "V_<1,4> scalar: M0(pc) proportional to N^(-25/8)",
        "rows": rows,
        "rules": [
            "Do not infer absence from a non-significant scalar projector.",
            "Do not fit a radial exponent when all per-size scalar z scores are underpowered.",
            "Use smaller N/high statistics or the leading-H4 annihilator for production tests.",
        ],
    }
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
