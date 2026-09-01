#!/usr/bin/env python3
"""Metric-cancelling amplitude ratios for the two-sector derivative ladder.

The construction uses scaling dimensions already frozen in the P48
four-channel program:

    P4[S]    ~ A_S    N^{-1}
    P4[D]    ~ A_D    N^{-13/8}
    P4[D']   ~ A_Dp   N^{-5/8}
    P4[S']   ~ A_Sp   N^{-5/4}
    Mbar'    ~ B      N^{3/8}

The thermal metric b that converts probability to the scaling field z sits
in every thermal derivative and in Mbar'. It therefore cancels in

    R_I = P4[D'] / (P4[S]  * Mbar') = A_Dp / (A_S  * B)
    R_T = P4[S'] / (P4[D]  * Mbar') = A_Sp / (A_D  * B)

Raw A_T / A_I does not cancel two independent lattice couplings and is not a
universality candidate. Descriptive P48 reconstructions are development-only;
this module does not retroactively define an Issue #57 target.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Dict, Mapping, Tuple

ROOT = Path(__file__).resolve().parents[1]
P48_SCORE = ROOT / "results" / "server-20260828" / "P48-new-geometry-score" / "score.json"
P35_CLOSURE = ROOT / "results" / "server-20260828" / "P35-amplitude-closure" / "closure.csv"

EXPONENTS: Dict[str, Fraction] = {
    "P4_S": Fraction(-1),
    "P4_D": Fraction(-13, 8),
    "P4_D_prime": Fraction(-5, 8),
    "P4_S_prime": Fraction(-5, 4),
    "Mbar_prime": Fraction(3, 8),
}


def ratio_exponent(numerator: str, denominator: Tuple[str, ...]) -> Fraction:
    power = EXPONENTS[numerator]
    for name in denominator:
        power -= EXPONENTS[name]
    return power


def metric_free_exponents() -> Dict[str, Fraction]:
    return {
        "R_I": ratio_exponent("P4_D_prime", ("P4_S", "Mbar_prime")),
        "R_T": ratio_exponent("P4_S_prime", ("P4_D", "Mbar_prime")),
        "raw_A_T_over_A_I": ratio_exponent("P4_D", ("P4_S",)),
    }


def scaled_ratio(numerator: float, left: float, right: float) -> float:
    return float(numerator) / (float(left) * float(right))


def p35_last_B(path: Path = P35_CLOSURE) -> Tuple[int, float]:
    rows = path.read_text(encoding="utf-8").strip().splitlines()
    header = rows[0].split(",")
    index_n = header.index("N")
    index_b = header.index("B")
    last = rows[-1].split(",")
    return int(last[index_n]), float(last[index_b])


def descriptive_p48_ratios(
    score_path: Path = P48_SCORE, b_path: Path = P35_CLOSURE
) -> dict[str, object]:
    payload = json.loads(score_path.read_text(encoding="utf-8"))
    source_n, source_b = p35_last_B(b_path)
    target_sizes = payload["target_sizes"]
    channels: Mapping[str, Mapping[str, object]] = payload["channels"]
    rows = []
    for index, size in enumerate(target_sizes):
        a_s = float(channels["P4_S"]["observed_scaled"][index])  # type: ignore[index]
        a_d = float(channels["P4_D"]["observed_scaled"][index])  # type: ignore[index]
        a_dp = float(channels["P4_D_prime"]["observed_scaled"][index])  # type: ignore[index]
        a_sp = float(channels["P4_S_prime"]["observed_scaled"][index])  # type: ignore[index]
        rows.append(
            {
                "N": size,
                "R_I": scaled_ratio(a_dp, a_s, source_b),
                "R_T": scaled_ratio(a_sp, a_d, source_b),
                "raw_A_D_over_A_S": a_d / a_s,
                "A_S_prime": a_sp,
            }
        )
    return {
        "role": "development_descriptive_only",
        "not_a_numerical_target_for_issue_57": True,
        "B_approximation": {
            "source": "P35 last row",
            "source_N": source_n,
            "B": source_b,
            "caveat": (
                "common-metric approximation transported from P35; "
                "not a delete-one B at the P48 target sizes"
            ),
        },
        "p48_score": str(score_path.relative_to(ROOT)),
        "rows": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    payload = {
        "exponents": {name: str(power) for name, power in metric_free_exponents().items()},
        "descriptive": descriptive_p48_ratios(),
    }
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print("net N-powers (must be 0 for metric-free ratios):")
        for name, power in metric_free_exponents().items():
            print(f"  {name:18} {power}")
        descriptive = payload["descriptive"]
        print(
            "descriptive P48 reconstruction "
            f"(B from P35 N={descriptive['B_approximation']['source_N']}):"
        )
        for row in descriptive["rows"]:
            print(
                f"  N={row['N']}  R_I={row['R_I']:.4f}  "
                f"R_T={row['R_T']:.4f}  A_D/A_S={row['raw_A_D_over_A_S']:.4f}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
