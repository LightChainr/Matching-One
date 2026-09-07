#!/usr/bin/env python3
"""Decide, in one command, what this repository already says about a claimed p_c.

Built for the case where someone -- an external model, a preprint, a referee --
hands us a closed form or an exact value for the square-site percolation
threshold.  The repository holds three exhaustive certified exclusions and four
mutually disjoint published intervals.  Between them, most claims can be refuted
or placed within seconds, and the ones that cannot are exactly the interesting
ones.

The tool refutes.  It never confirms: passing every check here means only that
the claim is not already dead, and the report says what would still have to be
done.

Usage
-----
    # a claimed minimal polynomial, coefficients ascending
    threshold_claim_intake.py --polynomial -1,0,3,0,-1,0,1

    # a claimed numeric value, as many digits as are being claimed
    threshold_claim_intake.py --decimal 0.5927460508

    # a claimed closed form, evaluated with mpmath at high precision
    threshold_claim_intake.py --expression "(sqrt(5)-1)/2 + 0.0"
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis" / "pslq_search_contract.json"
RESULTS = ROOT / "results"

SCHEMA = "matching-one.threshold-claim-intake.v1"

# The censused classes.  Whether a class is actually *excluded* varies by
# interval -- the wider Mertens intervals have surviving quartics -- so the
# verdict reads each committed artifact rather than assuming.
EXCLUSIONS = (
    {
        "id": "degree4-height100",
        "degree_min": 1,
        "degree_max": 4,
        "height_max": 100,
        "artifact": "results/pslq-degree4-{interval}/latest.json",
        "note": "primitive sign-normalized integer quartics, exhaustive",
    },
    {
        "id": "degree6-height3",
        "degree_min": 1,
        "degree_max": 6,
        "height_max": 3,
        "artifact": "results/pslq-degree6-low-height-{interval}/latest.json",
        "note": "the historical complexity range as first read, at height 3; "
                "superseded by degree6-height4, which is the class that "
                "actually covers the exactly-known thresholds",
    },
    {
        "id": "degree6-height4",
        "degree_min": 1,
        "degree_max": 6,
        "height_max": 4,
        "artifact": "results/pslq-degree6-height4-{interval}/latest.json",
        "note": "the corrected historical complexity range: every exactly-known "
                "planar percolation threshold has degree <= 6 and height <= 4, "
                "the height being set by the Ziff 2006 A-lattice quintic",
    },
)


def normalize(coefficients: Sequence[int]) -> list[int]:
    """Primitive, leading coefficient positive -- the census normalization.

    A claim written as 2x^4-4x^2+2 or as -(x^4-2x^2+1) is the same polynomial as
    far as its roots go, and the census enumerated one representative of each
    class.  Comparing without normalizing would silently miss a match.
    """
    from math import gcd

    values = list(coefficients)
    divisor = 0
    for value in values:
        divisor = gcd(divisor, abs(value))
    if divisor > 1:
        values = [value // divisor for value in values]
    if values[-1] < 0:
        values = [-value for value in values]
    return values


def exclusion_status(exclusion: dict, interval_id: str) -> dict:
    """Read the committed census artifact for this class and interval."""
    path = ROOT / exclusion["artifact"].format(interval=interval_id)
    if not path.is_file():
        return {"artifact": None, "excluded": None, "witnesses": []}
    payload = json.loads(path.read_text(encoding="utf-8"))
    result = payload.get("interval_result", {})
    witnesses = [
        normalize(row["coefficients_ascending"])
        for row in result.get("root_witnesses", [])
    ]
    if not witnesses:
        for row in result.get("by_degree", []):
            witnesses.extend(
                normalize(item["coefficients_ascending"])
                for item in row.get("root_witnesses", [])
            )
    return {
        "artifact": exclusion["artifact"].format(interval=interval_id),
        "excluded": result.get("excluded"),
        "witnesses": witnesses,
    }


def load_intervals() -> list[dict]:
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    return contract["intervals"]


def parse_polynomial(text: str) -> list[int]:
    values = [int(part.strip()) for part in text.split(",")]
    while values and values[-1] == 0:
        values.pop()
    if len(values) < 2:
        raise ValueError("a minimal polynomial needs degree at least 1")
    return values


def height(coefficients: Sequence[int]) -> int:
    return max(abs(value) for value in coefficients)


def evaluate(coefficients: Sequence[int], point: Fraction) -> Fraction:
    total = Fraction(0)
    power = Fraction(1)
    for value in coefficients:
        total += value * power
        power *= point
    return total


def roots_in(coefficients: Sequence[int], lower: Fraction, upper: Fraction) -> list:
    import sys

    sys.path.insert(0, str(ROOT / "scripts"))
    from exact_polynomial_root_certificate import isolate_roots  # noqa: E402

    return isolate_roots([Fraction(v) for v in coefficients], lower, upper, bits=120)


def decimal_window(text: str) -> tuple[Fraction, Fraction]:
    """A decimal written to d places asserts a half-ulp window, not a point.

    "0.5927460508" is a claim about ten digits, so it is the window
    0.5927460508 +/- 5e-11.  Treating it as an exact rational would make a
    correct-to-ten-digits claim look like it contradicts an estimate quoted to
    fourteen, which is the opposite of the truth.
    """
    cleaned = text.strip()
    if "." not in cleaned:
        return Fraction(cleaned), Fraction(cleaned)
    places = len(cleaned.split(".", 1)[1])
    value = Fraction(cleaned)
    half_ulp = Fraction(1, 2 * 10**places)
    return value - half_ulp, value + half_ulp


def value_from_arguments(arguments) -> tuple[tuple[Fraction, Fraction] | None, str, str]:
    """Return (window, how it was obtained, the text as given)."""
    if arguments.decimal is not None:
        return decimal_window(arguments.decimal), "decimal_half_ulp_window", arguments.decimal
    if arguments.expression is not None:
        import mpmath as mp

        with mp.workdps(80):
            value = mp.mpf(eval(arguments.expression, {"__builtins__": {}}, vars(mp)))
            text = mp.nstr(value, 50)
        return decimal_window(text), "mpmath_80_digits_reported_to_50", arguments.expression
    return None, "none", ""


def containing_intervals(window: tuple[Fraction, Fraction], intervals: Iterable[dict]) -> list[str]:
    """Published intervals the claim's own window overlaps."""
    low, high = window
    overlapping = []
    for interval in intervals:
        if Fraction(interval["lower"]) <= high and low <= Fraction(interval["upper"]):
            overlapping.append(interval["id"])
    return overlapping


def polynomial_report(coefficients: Sequence[int], intervals) -> dict:
    canonical = normalize(coefficients)
    degree = len(canonical) - 1
    h = height(canonical)
    per_interval = []
    for interval in intervals:
        lower = Fraction(interval["lower"])
        upper = Fraction(interval["upper"])
        found = roots_in(canonical, lower, upper)
        classes = []
        for exclusion in EXCLUSIONS:
            if not (exclusion["degree_min"] <= degree <= exclusion["degree_max"]):
                continue
            if h > exclusion["height_max"]:
                continue
            status = exclusion_status(exclusion, interval["id"])
            classes.append(
                {
                    "class": exclusion["id"],
                    "artifact": status["artifact"],
                    "class_is_excluded_on_this_interval": status["excluded"],
                    "claim_is_a_recorded_survivor": canonical in status["witnesses"],
                    "recorded_survivors_in_this_class": len(status["witnesses"]),
                }
            )
        per_interval.append(
            {
                "interval": interval["id"],
                "has_root_in_interval": bool(found),
                "isolating_windows": [[str(a), str(b)] for a, b in found],
                "censused_classes_covering_this_claim": classes,
            }
        )
    return {
        "degree": degree,
        "height": h,
        "coefficients_ascending_as_given": list(coefficients),
        "coefficients_ascending_normalized": canonical,
        "per_interval": per_interval,
    }


def verdict(report: dict) -> dict:
    contradictions = []
    known_survivor = []
    uncensused = []

    for row in report.get("polynomial", {}).get("per_interval", []):
        if not row["has_root_in_interval"]:
            continue
        classes = row["censused_classes_covering_this_claim"]
        if not classes:
            uncensused.append(row["interval"])
            continue
        for entry in classes:
            if entry["class_is_excluded_on_this_interval"]:
                contradictions.append(
                    f"{row['interval']}: the claimed polynomial has a root there, "
                    f"but {entry['class']} is an exhaustive census of that class "
                    f"on that interval and certified that none does "
                    f"({entry['artifact']})"
                )
            elif entry["claim_is_a_recorded_survivor"]:
                known_survivor.append(
                    f"{row['interval']}: already one of the "
                    f"{entry['recorded_survivors_in_this_class']} recorded "
                    f"{entry['class']} survivors there, catalogued as a width "
                    f"artifact rather than a candidate formula"
                )

    inside = report.get("value", {}).get("intervals_containing_it")

    polynomial = report.get("polynomial")
    if polynomial and not any(
        row["has_root_in_interval"] for row in polynomial["per_interval"]
    ):
        return {
            "outcome": "the_polynomial_has_no_root_in_any_published_interval",
            "because": [
                "certified by exact Sturm isolation at 120 bits on all four "
                "intervals: this polynomial has no real root in any of them"
            ],
            "what_this_means": (
                "it cannot be a minimal polynomial for the square-site threshold "
                "unless every published estimate is wrong by more than its own "
                "quoted uncertainty"
            ),
        }

    if contradictions:
        return {
            "outcome": "refuted_by_a_committed_certificate",
            "because": contradictions,
            "what_this_means": (
                "the claim and one of our exhaustive censuses cannot both be "
                "right; the census is exact integer arithmetic with a second "
                "independent implementation agreeing cell by cell"
            ),
        }

    if known_survivor:
        return {
            "outcome": "already_catalogued_as_a_width_artifact",
            "because": known_survivor,
            "what_this_means": (
                "having a root inside one of the wider published intervals is "
                "not evidence: the census found these by enumeration, each "
                "survives exactly one interval, and the four intervals are "
                "pairwise disjoint so at least three do not contain p_c"
            ),
        }

    if inside is not None and not inside:
        return {
            "outcome": "contradicts_every_published_interval",
            "because": [
                "the claimed value's own precision window overlaps none of the "
                "four published estimates"
            ],
            "what_this_means": (
                "not automatically wrong -- the four intervals are pairwise "
                "disjoint, so at least three already fail to contain p_c -- but "
                "it must explain why every published estimate is off"
            ),
        }

    return {
        "outcome": "survives_our_checks",
        "because": (
            [f"{name}: root lies in a class never censused here" for name in uncensused]
            or ["nothing committed here refutes it"]
        ),
        "what_this_means": "this is the interesting case, and it is not a confirmation",
        "what_would_still_have_to_be_done": [
            "reproduce the claimed value to more digits than any published "
            "estimate, by a method independent of the one that produced it",
            "state the minimal polynomial or closed form exactly, so that degree "
            "and height are checkable rather than asserted",
            "if algebraic: census the claim's own (degree, height) class on the "
            "interval and confirm it is the unique root there",
            "reconcile it with the four disjoint published intervals, at least "
            "three of which must be wrong",
            "supply the mechanism -- every exactly-known planar threshold comes "
            "from a duality, star-triangle or transfer-matrix identity, not from "
            "a number that happens to fit",
        ],
    }


def render(arguments) -> dict:
    intervals = load_intervals()
    report: dict = {
        "schema": SCHEMA,
        "claim_level": "not_applicable_this_is_a_filter",
        "input": {
            "polynomial": arguments.polynomial,
            "decimal": arguments.decimal,
            "expression": arguments.expression,
        },
        "published_intervals": [
            {
                "id": row["id"],
                "lower": row["lower"],
                "upper": row["upper"],
            }
            for row in intervals
        ],
        "censused_classes_held_here": [
            {
                "id": row["id"],
                "degrees": [row["degree_min"], row["degree_max"]],
                "height_max": row["height_max"],
                "note": row["note"],
                "excluded_on": [
                    interval["id"]
                    for interval in intervals
                    if exclusion_status(row, interval["id"])["excluded"] is True
                ],
                "has_recorded_survivors_on": [
                    interval["id"]
                    for interval in intervals
                    if exclusion_status(row, interval["id"])["excluded"] is False
                ],
            }
            for row in EXCLUSIONS
        ],
    }

    window, how, text = value_from_arguments(arguments)
    if window is not None:
        report["value"] = {
            "as_given": text,
            "obtained_by": how,
            "window": [str(window[0]), str(window[1])],
            "intervals_containing_it": containing_intervals(window, intervals),
        }

    if arguments.polynomial:
        coefficients = parse_polynomial(arguments.polynomial)
        report["polynomial"] = polynomial_report(coefficients, intervals)

    report["verdict"] = verdict(report)
    report["this_tool_cannot"] = [
        "confirm a claim; it can only refute one or fail to",
        "check a proof, a mechanism, or an argument -- only a number or a "
        "polynomial",
        "rule out an algebraic value of degree 5 or 6 with height above 3, or of "
        "degree above 6 at any height: those classes have never been censused here",
    ]
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polynomial", help="integer coefficients, ascending, comma separated")
    parser.add_argument("--decimal", help="claimed value as a decimal string")
    parser.add_argument("--expression", help="closed form, evaluated by mpmath")
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    if not any((arguments.polynomial, arguments.decimal, arguments.expression)):
        parser.error("give at least one of --polynomial, --decimal, --expression")
    text = json.dumps(render(arguments), indent=2, sort_keys=True) + "\n"
    if arguments.output is None:
        print(text, end="")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
