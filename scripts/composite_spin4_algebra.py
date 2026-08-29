#!/usr/bin/env python3
"""Exact parity and harmonic algebra for composite spin-4 corrections.

The declared spin-4 generators all carry the angular factor cos(4 theta).
Products are expanded by Laurent convolution in z=exp(4 i theta), so every
Fourier coefficient and radial exponent remains rational.  This module checks
the selection rules discussed in Issue 58 without fitting any data.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from itertools import product
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Tuple


GENERATOR_ORDER = ("T4", "I4", "V4", "S0")
GENERATORS = {
    "T4": {"parity": -1, "spin": 4, "omega": Fraction(13, 4), "status": "observed_candidate"},
    "I4": {"parity": 1, "spin": 4, "omega": Fraction(2), "status": "observed_candidate"},
    "V4": {"parity": 1, "spin": 4, "omega": Fraction(8, 3), "status": "conditional_parity"},
    "S0": {"parity": 1, "spin": 0, "omega": Fraction(2), "status": "optional_scalar"},
}
LEADING_T4_OMEGA = Fraction(13, 4)


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return "%d/%d" % (value.numerator, value.denominator)


def laurent_cos_power(power: int) -> Dict[int, Fraction]:
    """Expand ``cos(4 theta)**power`` in powers of ``z=exp(4 i theta)``."""

    if power < 0:
        raise ValueError("power must be nonnegative")
    coefficients = {0: Fraction(1)}
    kernel = {-1: Fraction(1, 2), 1: Fraction(1, 2)}
    for _ in range(power):
        next_coefficients: Dict[int, Fraction] = {}
        for left_power, left_value in coefficients.items():
            for right_power, right_value in kernel.items():
                exponent = left_power + right_power
                next_coefficients[exponent] = (
                    next_coefficients.get(exponent, Fraction(0))
                    + left_value * right_value
                )
        coefficients = next_coefficients
    return coefficients


def cosine_harmonics(power: int) -> Dict[int, Fraction]:
    """Return real-cosine coefficients keyed by physical harmonic spin."""

    laurent = laurent_cos_power(power)
    harmonics: Dict[int, Fraction] = {}
    if 0 in laurent:
        harmonics[0] = laurent[0]
    for exponent in sorted(k for k in laurent if k > 0):
        if laurent.get(-exponent) != laurent[exponent]:
            raise AssertionError("reflection symmetry was lost")
        harmonics[4 * exponent] = 2 * laurent[exponent]
    return harmonics


def normalize_counts(counts: Mapping[str, int]) -> Dict[str, int]:
    normalized = {name: int(counts.get(name, 0)) for name in GENERATOR_ORDER}
    if any(value < 0 for value in normalized.values()):
        raise ValueError("generator multiplicities must be nonnegative")
    if not any(normalized.values()):
        raise ValueError("the empty monomial is not a correction")
    return normalized


def monomial_name(counts: Mapping[str, int]) -> str:
    factors = []
    for name in GENERATOR_ORDER:
        count = counts.get(name, 0)
        if count == 1:
            factors.append(name)
        elif count > 1:
            factors.append("%s^%d" % (name, count))
    return "*".join(factors)


def classify_monomial(counts: Mapping[str, int]) -> dict[str, Any]:
    """Classify one product relative to the leading T4 correction."""

    normalized = normalize_counts(counts)
    parity = 1
    total_omega = Fraction(0)
    spin4_factor_count = 0
    conditional = []
    for name, count in normalized.items():
        generator = GENERATORS[name]
        parity *= int(generator["parity"]) ** count
        total_omega += generator["omega"] * count
        if generator["spin"] == 4:
            spin4_factor_count += count
        if count and generator["status"] in ("conditional_parity", "optional_scalar"):
            conditional.append(name)

    relative_q = total_omega - LEADING_T4_OMEGA
    harmonics = cosine_harmonics(spin4_factor_count)
    return {
        "name": monomial_name(normalized),
        "counts": normalized,
        "degree": sum(normalized.values()),
        "matching_parity": parity,
        "spin4_factor_count": spin4_factor_count,
        "total_omega": fraction_text(total_omega),
        "relative_q": fraction_text(relative_q),
        "accelerated_w": fraction_text(Fraction(4) + relative_q),
        "harmonics": {"H%d" % spin: fraction_text(value) for spin, value in harmonics.items()},
        "has_H4": 4 in harmonics,
        "has_H12": 12 in harmonics,
        "H12_over_H4": (
            fraction_text(harmonics[12] / harmonics[4])
            if 4 in harmonics and 12 in harmonics
            else None
        ),
        "conditional_generators": conditional,
    }


def exponent_tuples(max_degree: int) -> Iterable[Tuple[int, int, int, int]]:
    for counts in product(range(max_degree + 1), repeat=len(GENERATOR_ORDER)):
        if 0 < sum(counts) <= max_degree:
            yield counts


def enumerate_matching_odd_h4(max_degree: int) -> list[dict[str, Any]]:
    rows = []
    for values in exponent_tuples(max_degree):
        counts = dict(zip(GENERATOR_ORDER, values))
        if counts["T4"] == 0:
            continue
        row = classify_monomial(counts)
        if row["matching_parity"] == -1 and row["has_H4"]:
            rows.append(row)
    rows.sort(
        key=lambda row: (
            Fraction(row["relative_q"]),
            row["degree"],
            row["name"],
        )
    )
    return rows


def q3_diophantine_certificate() -> dict[str, Any]:
    """Prove q=3 is absent from the full declared generator semiring.

    A matching-odd monomial has an odd positive number t of T4 factors.  Its
    exponent relative to one leading T4 is

      q = 13(t-1)/4 + 2(i+s) + 8v/3.

    Multiplication by 12 gives

      12q = 39(t-1) + 24(i+s) + 32v.

    For q=3 the right side must equal 36.  If t>=3, its first term is at
    least 78.  If t=1, division by four gives 6(i+s)+8v=9, impossible since
    the left side is even.
    """

    return {
        "target": "q=3",
        "equation": "39*(t-1)+24*(i+s)+32*v=36",
        "matching_odd_constraint": "t is a positive odd integer",
        "case_t_at_least_3": "39*(t-1)>=78>36",
        "case_t_equals_1": "6*(i+s)+8*v=9 after division by 4; even=odd is impossible",
        "conclusion": "no monomial at any degree in the declared generator semiring has q=3",
    }


def build_artifact(max_degree: int = 5) -> dict[str, Any]:
    named_counts = {
        "leading_T4": {"T4": 1},
        "conditional_scalar": {"T4": 1, "S0": 1},
        "one_I4_insertion": {"T4": 1, "I4": 1},
        "two_I4_insertions": {"T4": 1, "I4": 2},
        "mixed_I4_V4": {"T4": 1, "I4": 1, "V4": 1},
        "two_V4_insertions": {"T4": 1, "V4": 2},
    }
    named = {name: classify_monomial(counts) for name, counts in named_counts.items()}
    rows = enumerate_matching_odd_h4(max_degree)
    q3_rows = [row for row in rows if row["relative_q"] == "3"]
    q6_rows = [row["name"] for row in rows if row["relative_q"] == "6"]

    assert named["conditional_scalar"]["relative_q"] == "2"
    assert not named["one_I4_insertion"]["has_H4"]
    assert named["two_I4_insertions"]["harmonics"] == {"H4": "3/4", "H12": "1/4"}
    assert named["mixed_I4_V4"]["accelerated_w"] == "26/3"
    assert named["two_V4_insertions"]["accelerated_w"] == "28/3"
    assert not q3_rows

    return {
        "schema": "matching-one/composite-spin4-algebra/v1",
        "issue": 58,
        "status": "exact_selection_rule_oracle",
        "generator_order": list(GENERATOR_ORDER),
        "generators": {
            name: {
                "matching_parity": data["parity"],
                "spin": data["spin"],
                "omega_length": fraction_text(data["omega"]),
                "status": data["status"],
            }
            for name, data in GENERATORS.items()
        },
        "angular_rule": {
            "method": "exact Laurent convolution of ((z+z^-1)/2)^n",
            "z": "exp(4 i theta)",
            "H4_condition": "the number of spin-4 factors is odd",
            "matching_odd_condition": "the number of T4 factors is odd",
        },
        "named_cases": named,
        "q3_no_go": q3_diophantine_certificate(),
        "enumeration": {
            "max_total_degree": max_degree,
            "matching_odd_H4_row_count": len(rows),
            "q3_rows": q3_rows,
            "q6_collision_rows": q6_rows,
            "rows": rows,
        },
        "external_same_family_row": {
            "name": "next ordinary thermal spin-4 quasiprimary",
            "relative_q": "6",
            "accelerated_w": "10",
            "warning": "q=6 is not unique once optional analytic composites are admitted",
        },
        "claim_boundary": {
            "proved": (
                "parity, exact harmonic support and coefficients, rational exponent sums, "
                "and the all-degree q=3 no-go inside the declared generator semiring"
            ),
            "conditional": "the parity assignment of V4 and the existence/nonzero coupling of S0",
            "not_proved": (
                "operator existence, nonzero amplitudes, the absence of other generators, "
                "or the literal H12/H4 amplitude ratio after continuum mixing"
            ),
        },
    }


def render_markdown(artifact: dict[str, Any]) -> str:
    named = artifact["named_cases"]
    case_order = [
        "leading_T4",
        "conditional_scalar",
        "one_I4_insertion",
        "two_I4_insertions",
        "mixed_I4_V4",
        "two_V4_insertions",
    ]
    lines = [
        "# Exact composite spin-4 algebra",
        "",
        "## Named cases",
        "",
        "| case | parity | q | accelerated w | harmonics | H12/H4 |",
        "|---|---:|---:|---:|---|---:|",
    ]
    for key in case_order:
        row = named[key]
        harmonics = ", ".join("%s:%s" % item for item in row["harmonics"].items())
        lines.append(
            "| `%s` | `%+d` | `%s` | `%s` | `%s` | `%s` |"
            % (
                row["name"],
                row["matching_parity"],
                row["relative_q"],
                row["accelerated_w"],
                harmonics,
                row["H12_over_H4"] or "--",
            )
        )
    lines.extend(
        [
            "",
            "`T4*I4` is matching odd but has only H0/H8, so one even spin-4 insertion cannot",
            "correct the H4 channel. Every listed cubic spin-4 product has exact support",
            "`(3/4) H4 + (1/4) H12`, giving the elementary ratio `H12/H4=1/3`.",
            "",
            "## q=3 exclusion",
            "",
            "For generator counts `(t,i,v,s)`, the relative exponent obeys",
            "",
            "```text",
            "12q = 39(t-1) + 24(i+s) + 32v.",
            "```",
            "",
            "Matching oddness requires positive odd `t`. At `q=3`, `t>=3` already contributes",
            "at least 78 to a target of 36. For `t=1`, division by four requires",
            "`6(i+s)+8v=9`, an even/odd contradiction. Thus `q=3` is absent at every degree",
            "inside the declared generator semiring, not merely in the finite enumeration.",
            "",
            "## Additional exact warning",
            "",
            "The ordinary thermal-tower row `q=6, w=10` is not exponent-unique after optional",
            "analytic composites are admitted. The degree-limited artifact lists the colliding",
            "composite monomials explicitly; harmonic sidebands and independent amplitude controls",
            "are required to distinguish them.",
            "",
            "## Evidence boundary",
            "",
            "The algebra does not prove that any generator exists with nonzero lattice coupling.",
            "`V4` parity and `S0` existence remain conditional. Continuum response tensors or",
            "mixing may change the literal H12/H4 amplitude ratio, although the elementary harmonic",
            "support and rational exponent arithmetic are exact.",
            "",
            "## Reproduction",
            "",
            "```bash",
            "python scripts/composite_spin4_algebra.py --format json",
            "python scripts/composite_spin4_algebra.py --format markdown",
            "python -m unittest tests.test_composite_spin4_algebra",
            "```",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-degree", type=int, default=5)
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.max_degree < 1:
        parser.error("--max-degree must be positive")

    artifact = build_artifact(args.max_degree)
    if args.format == "json":
        rendered = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    else:
        rendered = render_markdown(artifact)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
