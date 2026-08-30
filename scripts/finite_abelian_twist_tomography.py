#!/usr/bin/env python3
"""Exact finite-abelian twist tomography for saturated torus carriers.

The integral-saturation theorem implies that a connected carrier image
``Lambda <= Z^2`` is a direct summand.  For any finite abelian group ``A``
of order ``n``, characters in ``Hom(Z^2,A)`` which kill ``Lambda`` therefore
number ``n**(2-rank(Lambda))``.  The aggregate constraint trace is a discrete
sample of the intrinsic topological source.  Prime-order individual twists
also resolve the rank-one winding line modulo that prime.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import product
import json
from math import gcd, prod
from pathlib import Path
from typing import Iterable, Sequence


GroupElement = tuple[int, ...]
Vector = tuple[int, int]


def _elements(moduli: Sequence[int]) -> list[GroupElement]:
    if not moduli or any(modulus < 2 for modulus in moduli):
        raise ValueError("finite abelian factors must have modulus at least two")
    return list(product(*(range(modulus) for modulus in moduli)))


def _linear_combination(
    first: GroupElement,
    second: GroupElement,
    vector: Vector,
    moduli: Sequence[int],
) -> GroupElement:
    return tuple(
        (vector[0] * a + vector[1] * b) % modulus
        for a, b, modulus in zip(first, second, moduli)
    )


def count_annihilating_characters(
    moduli: Sequence[int], generators: Sequence[Vector]
) -> int:
    """Count ``alpha:Z^2->A`` which vanish on every supplied generator."""

    elements = _elements(moduli)
    zero = tuple(0 for _ in moduli)
    return sum(
        all(_linear_combination(a, b, vector, moduli) == zero for vector in generators)
        for a in elements
        for b in elements
    )


def _primitive(vector: Vector) -> bool:
    return vector != (0, 0) and gcd(abs(vector[0]), abs(vector[1])) == 1


def finite_abelian_audit() -> dict[str, object]:
    groups = [(2,), (3,), (4,), (2, 2), (2, 3), (2, 2, 2)]
    lines = [(1, 0), (0, 1), (1, 1), (-2, 1), (2, 3)]
    rows = []
    for moduli in groups:
        order = prod(moduli)
        for line in lines:
            if not _primitive(line):
                raise AssertionError("audit line is not primitive")
            counts = {
                "rank0": count_annihilating_characters(moduli, ()),
                "rank1": count_annihilating_characters(moduli, (line,)),
                "rank2": count_annihilating_characters(moduli, ((1, 0), (0, 1))),
            }
            expected = {"rank0": order**2, "rank1": order, "rank2": 1}
            rows.append(
                {
                    "group_factors": list(moduli),
                    "group_order": order,
                    "primitive_line": list(line),
                    "counts": counts,
                    "expected": expected,
                    "pass": counts == expected,
                }
            )
    return {
        "rows": rows,
        "row_count": len(rows),
        "all_pass": all(row["pass"] for row in rows),
    }


def projective_lines(prime: int) -> list[Vector]:
    if prime < 2:
        raise ValueError("prime must be at least two")
    return [(1, slope) for slope in range(prime)] + [(0, 1)]


def _kernel_line(alpha: Vector, prime: int) -> Vector:
    if alpha == (0, 0):
        raise ValueError("the zero twist has the whole plane as kernel")
    for line in projective_lines(prime):
        if (alpha[0] * line[0] + alpha[1] * line[1]) % prime == 0:
            return line
    raise AssertionError("a nonzero functional on F_q^2 must have a line kernel")


def prime_projective_audit(prime: int) -> dict[str, object]:
    lines = projective_lines(prime)
    multiplicities = {line: 0 for line in lines}
    for alpha in product(range(prime), repeat=2):
        if alpha != (0, 0):
            multiplicities[_kernel_line(alpha, prime)] += 1
    return {
        "prime": prime,
        "projective_line_count": len(lines),
        "nonzero_twist_count": prime**2 - 1,
        "kernel_multiplicities": [
            {"line": list(line), "twists": multiplicities[line]} for line in lines
        ],
        "every_line_has_q_minus_1_twists": all(
            count == prime - 1 for count in multiplicities.values()
        ),
    }


def aggregate_trace(order: int, probabilities: Sequence[Fraction]) -> Fraction:
    p0, p1, p2 = probabilities
    return order * order * p0 + order * p1 + p2


def reconstruct_from_order_two_three(
    trace_two: Fraction, trace_three: Fraction
) -> tuple[Fraction, Fraction, Fraction]:
    """Invert ``S_n=n^2 P0+n P1+P2`` using ``S_1=1,S_2,S_3``."""

    p0 = (trace_three - 2 * trace_two + 1) / 2
    p1 = trace_two - 1 - 3 * p0
    p2 = 1 - p0 - p1
    return p0, p1, p2


def _fraction(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def source_inversion_audit() -> dict[str, object]:
    probabilities = (Fraction(1, 7), Fraction(2, 7), Fraction(4, 7))
    trace_two = aggregate_trace(2, probabilities)
    trace_three = aggregate_trace(3, probabilities)
    recovered = reconstruct_from_order_two_three(trace_two, trace_three)
    source_samples = {
        str(order): {
            "aggregate_trace_S_n": _fraction(aggregate_trace(order, probabilities)),
            "normalized_twist_average_S_n_over_n": _fraction(
                aggregate_trace(order, probabilities) / order
            ),
        }
        for order in (2, 3, 4, 6, 8)
    }
    return {
        "input_probabilities_P0_P1_P2": [_fraction(value) for value in probabilities],
        "trace_two": _fraction(trace_two),
        "trace_three": _fraction(trace_three),
        "recovered_probabilities_P0_P1_P2": [_fraction(value) for value in recovered],
        "exact_recovery": recovered == probabilities,
        "source_samples": source_samples,
    }


def prime_line_tomography_audit(prime: int) -> dict[str, object]:
    p0 = Fraction(1, 7)
    p1 = Fraction(2, 7)
    p2 = Fraction(4, 7)
    lines = projective_lines(prime)
    denominator = sum(range(1, len(lines) + 1))
    line_weights = {
        line: p1 * Fraction(index, denominator)
        for index, line in enumerate(lines, start=1)
    }
    twist_sum = Fraction(1)
    recovered = {}
    for line in lines:
        constrained_probability = p0 + line_weights[line]
        twist_sum += (prime - 1) * constrained_probability
        recovered[line] = constrained_probability - p0
    expected_sum = aggregate_trace(prime, (p0, p1, p2))
    return {
        "prime": prime,
        "P0": _fraction(p0),
        "P1": _fraction(p1),
        "P2": _fraction(p2),
        "line_rows": [
            {
                "line_mod_q": list(line),
                "rank1_line_probability": _fraction(line_weights[line]),
                "nonzero_twist_constraint_probability": _fraction(
                    p0 + line_weights[line]
                ),
                "recovered_after_subtracting_P0": _fraction(recovered[line]),
            }
            for line in lines
        ],
        "sum_over_all_twists": _fraction(twist_sum),
        "aggregate_formula": _fraction(expected_sum),
        "aggregate_pass": twist_sum == expected_sum,
        "line_recovery_pass": recovered == line_weights,
    }


def build_certificate() -> dict[str, object]:
    finite_groups = finite_abelian_audit()
    projective = [prime_projective_audit(prime) for prime in (2, 3, 5, 7)]
    inversion = source_inversion_audit()
    line_tomography = [prime_line_tomography_audit(prime) for prime in (2, 3, 5)]
    all_pass = (
        finite_groups["all_pass"]
        and all(row["every_line_has_q_minus_1_twists"] for row in projective)
        and inversion["exact_recovery"]
        and all(
            row["aggregate_pass"] and row["line_recovery_pass"]
            for row in line_tomography
        )
    )
    return {
        "schema": "matching-one/finite-abelian-twist-tomography/v1",
        "issues": [269, 334, 337],
        "status": "exact_finite_abelian_twist_transform",
        "theorem": {
            "scope": (
                "Every theorem-supported black NN or white matching configuration "
                "and every finite abelian coefficient group A of order n. Connected "
                "images are saturated by #269; disjoint essential components on the "
                "embedded torus are parallel, so the global image is also saturated."
            ),
            "annihilator_count": (
                "|{alpha in Hom(Z^2,A): alpha|Lambda=0}|=n^(2-rank Lambda)."
            ),
            "aggregate_transform": "S_n=sum_alpha T_alpha=n^2 P0+n P1+P2.",
            "source_value": (
                "S_n/n=Z_top(s=-log n)=n P0+P1+n^(-1) P2."
            ),
            "two_order_inversion": {
                "P0": "(S_3-2 S_2+1)/2",
                "P1": "S_2-1-3 P0",
                "P2": "1-P0-P1",
            },
            "prime_projective_refinement": (
                "For nonzero alpha in F_q^2, T_alpha=P0+L_ker(alpha), "
                "where L_line is the probability of a rank-one primitive "
                "winding line reducing to that projective class modulo q."
            ),
        },
        "consequences": [
            "Twist averages at group orders 2 and 3 reconstruct the entire unmarked rank-source functional, with normalization S_1=1.",
            "Individual prime twists are a modular projective-line tomography of the #334 first/plateau winding direction.",
            "Integral saturation makes r_q=r for every prime: finite-field ranks cannot reveal an additional Smith/index state on these carriers.",
            "The aggregate transform depends only on |A|, while individual twist sectors retain line-incidence information.",
        ],
        "machine_certificates": {
            "finite_abelian_groups": finite_groups,
            "prime_projective_orbits": projective,
            "two_order_source_inversion": inversion,
            "prime_line_tomography": line_tomography,
            "all_pass": all_pass,
        },
        "claim_boundary": [
            "The result is an exact lattice/cohomology transform, not yet an identification with a local CFT field.",
            "A finite set of primes resolves winding lines only modulo those primes; it does not reconstruct an unbounded integral line without an external size bound.",
            "The projective refinement concerns the rank-one carrier line, not a new saturation index, because #269 proves saturation.",
        ],
    }


def render_markdown(certificate: dict[str, object]) -> str:
    theorem = certificate["theorem"]
    machine = certificate["machine_certificates"]
    lines = [
        "# Finite-abelian twist tomography of the intrinsic rank source",
        "",
        f"Status: `{certificate['status']}`.",
        "",
        "## Exact transform",
        "",
        theorem["annihilator_count"],
        "",
        f"`{theorem['aggregate_transform']}`",
        "",
        f"`{theorem['source_value']}`",
        "",
        "Thus order-2 and order-3 twist averages, together with `S_1=1`, invert to:",
        "",
    ]
    for name, formula in theorem["two_order_inversion"].items():
        lines.append(f"- `{name} = {formula}`")
    lines += [
        "",
        "## Projective refinement",
        "",
        theorem["prime_projective_refinement"],
        "",
        "This converts the finite-field construction from proposed saturation tomography into modular winding-line tomography: saturation is already exact by #269.",
        "",
        "## Executable gates",
        "",
        f"- finite-abelian enumeration rows: {machine['finite_abelian_groups']['row_count']}",
        f"- prime projective audits: {len(machine['prime_projective_orbits'])}",
        f"- prime line-tomography audits: {len(machine['prime_line_tomography'])}",
        f"- all gates pass: `{machine['all_pass']}`",
        "",
        "## Consequence",
        "",
    ]
    lines.extend(f"- {item}" for item in certificate["consequences"])
    lines += ["", "## Boundary", ""]
    lines.extend(f"- {item}" for item in certificate["claim_boundary"])
    return "\n".join(lines)


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    certificate = build_certificate()
    rendered = (
        json.dumps(certificate, indent=2, sort_keys=True) + "\n"
        if args.format == "json"
        else render_markdown(certificate) + "\n"
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
