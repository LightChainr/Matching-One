#!/usr/bin/env python3
"""Exact tiny oracle for typed ambient-H1 rank-birth insertions.

For a site v and an environment on all other sites, the insertion records the
rank transition, the two elementary rank gates it crosses, the unique
rank-one plateau line when it exists, its integral saturation index, and an
optional local landing-sector H4 mark.  Exact polynomial arithmetic certifies
the Russo decomposition M'(p)=f_01(p)+f_12(p).
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
from math import comb
from pathlib import Path
from typing import Any, Optional, Sequence

from digital_alexander_filtration_oracle import rank_mark
from integer_period_torus import (
    IntegerTorusGeometry,
    Vector,
    axis_integer_torus,
    gaussian_integer_torus,
)
from marked_pivotal_h4_reference import landing_mark


Polynomial = list[Fraction]
RankState = tuple[int, Optional[Vector], Optional[int]]


def _fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _polynomial_text(values: Sequence[Fraction]) -> list[str]:
    return [_fraction_text(value) for value in values]


def _add_polynomial(left: Sequence[Fraction], right: Sequence[Fraction]) -> Polynomial:
    degree = max(len(left), len(right))
    return [
        (left[i] if i < len(left) else Fraction(0))
        + (right[i] if i < len(right) else Fraction(0))
        for i in range(degree)
    ]


def _scale_polynomial(value: int, polynomial: Sequence[Fraction]) -> Polynomial:
    return [value * coefficient for coefficient in polynomial]


def _bernstein_term(occupied: int, total: int) -> Polynomial:
    """Return p^occupied (1-p)^(total-occupied) in ascending power basis."""

    result = [Fraction(0)] * (total + 1)
    for extra in range(total - occupied + 1):
        result[occupied + extra] = Fraction(((-1) ** extra) * comb(total - occupied, extra))
    return result


def _differentiate(polynomial: Sequence[Fraction]) -> Polynomial:
    if len(polynomial) <= 1:
        return [Fraction(0)]
    return [degree * polynomial[degree] for degree in range(1, len(polynomial))]


def _evaluate(polynomial: Sequence[Fraction], value: Fraction) -> Fraction:
    result = Fraction(0)
    for coefficient in reversed(polynomial):
        result = result * value + coefficient
    return result


def _mask_active(mask: int, n: int) -> list[bool]:
    return [bool(mask & (1 << vertex)) for vertex in range(n)]


def _rank_state_cache(geometry: IntegerTorusGeometry) -> list[RankState]:
    return [
        rank_mark(geometry, _mask_active(mask, geometry.n), matching=False)
        for mask in range(1 << geometry.n)
    ]


def homology_h4_mark(geometry: IntegerTorusGeometry, line: Vector) -> dict[str, str]:
    """Exact spin-4 harmonic of a primitive line in physical lattice coordinates."""

    x, y = geometry.periods.period_vector(line)
    radius_squared = x * x + y * y
    denominator = radius_squared * radius_squared
    cos4 = Fraction(x**4 - 6 * x * x * y * y + y**4, denominator)
    sin4 = Fraction(4 * x * y * (x * x - y * y), denominator)
    return {
        "physical_vector": f"{x},{y}",
        "cos4": _fraction_text(cos4),
        "sin4": _fraction_text(sin4),
    }


def rank_birth_insertion(
    geometry: IntegerTorusGeometry,
    active_without: Sequence[bool],
    vertex: int,
    *,
    local_radius: Optional[int] = None,
) -> dict[str, Any]:
    """Return the complete typed insertion for one absent site.

    A direct 0->2 jump is represented as two simultaneous elementary births.
    It has no canonical intermediate projective line.  For 0->1 the line is
    the new rank-one image; for 1->2 it is the rank-one plateau line destroyed
    at the second birth.  This is the canonical line shared by the two
    endpoints of a nonempty essential-H1 interval.
    """

    if len(active_without) != geometry.n:
        raise ValueError("active mask length does not match geometry")
    if not 0 <= vertex < geometry.n:
        raise ValueError("vertex outside geometry")
    without = list(active_without)
    without[vertex] = False
    before_rank, before_line, before_index = rank_mark(
        geometry, without, matching=False
    )
    with_site = list(without)
    with_site[vertex] = True
    after_rank, after_line, after_index = rank_mark(
        geometry, with_site, matching=False
    )
    if not 0 <= before_rank <= after_rank <= 2:
        raise AssertionError("ambient rank is not monotone under site insertion")

    first_birth = int(before_rank == 0 and after_rank >= 1)
    second_birth = int(before_rank <= 1 and after_rank == 2)
    if after_rank - before_rank != first_birth + second_birth:
        raise AssertionError("rank jump did not decompose into the two elementary gates")

    line: Optional[Vector] = None
    index: Optional[int] = None
    line_role: Optional[str] = None
    if before_rank == 0 and after_rank == 1:
        line, index = after_line, after_index
        line_role = "new_rank_one_image_after_0_to_1"
    elif before_rank == 1 and after_rank == 2:
        line, index = before_line, before_index
        line_role = "rank_one_plateau_line_before_1_to_2"

    if line is not None and index is None:
        raise AssertionError("rank-one line is missing its saturation index")
    if line is None and index is not None:
        raise AssertionError("saturation index was recorded without a line")

    births = []
    for label, present in (("0_to_1", first_birth), ("1_to_2", second_birth)):
        if present:
            births.append(
                {
                    "type": label,
                    "simultaneous_0_to_2": bool(first_birth and second_birth),
                    "ell": list(line) if line is not None else None,
                    "ell_role": line_role,
                    "iota": index,
                    "homology_h4": homology_h4_mark(geometry, line) if line else None,
                }
            )

    local_mark: dict[str, Any]
    if local_radius is None:
        local_mark = {
            "status": "not_evaluated_on_short_period_control",
            "radius": None,
            "axis": None,
            "diagonal": None,
            "landed": None,
            "h4": None,
        }
    else:
        if vertex != geometry.vertex((0, 0)):
            raise ValueError("landing_mark currently uses the registered origin root")
        local_mark = {
            "status": "evaluated",
            "radius": local_radius,
            **landing_mark(
                geometry,
                without,
                local_radius,
                open_matching=False,
            ),
        }

    return {
        "vertex": vertex,
        "rank_before": before_rank,
        "rank_after": after_rank,
        "delta_rank": after_rank - before_rank,
        "gate_0_to_1": first_birth,
        "gate_1_to_2": second_birth,
        "births": births,
        "local_geometry_mark": local_mark,
    }


def _rank_polynomial(
    geometry: IntegerTorusGeometry,
    states: Sequence[RankState],
) -> Polynomial:
    result = [Fraction(0)] * (geometry.n + 1)
    for mask, state in enumerate(states):
        value = state[0] - 1
        term = _scale_polynomial(value, _bernstein_term(mask.bit_count(), geometry.n))
        result = _add_polynomial(result, term)
    return result


def _line_key(birth: dict[str, Any]) -> str:
    line = birth["ell"]
    if line is None:
        return "none(simultaneous-0-to-2)"
    return f"({line[0]},{line[1]});iota={birth['iota']}"


def _geometry_summary(
    name: str,
    geometry: IntegerTorusGeometry,
    *,
    roots: Sequence[int],
    root_multiplicity: int,
    local_radius: Optional[int],
) -> dict[str, Any]:
    states = _rank_state_cache(geometry)
    transition_counts: Counter[str] = Counter()
    birth_counts = {"0_to_1": [0] * geometry.n, "1_to_2": [0] * geometry.n}
    line_counts: dict[str, Counter[str]] = {
        "0_to_1": Counter(),
        "1_to_2": Counter(),
    }
    local_mark_sums = {
        "0_to_1": Counter(),
        "1_to_2": Counter(),
    }
    homology_cos4_sums = {"0_to_1": Fraction(0), "1_to_2": Fraction(0)}
    homology_sin4_sums = {"0_to_1": Fraction(0), "1_to_2": Fraction(0)}
    representatives: dict[str, Any] = {}

    for vertex in roots:
        bit = 1 << vertex
        for mask in range(1 << geometry.n):
            if mask & bit:
                continue
            before_rank = states[mask][0]
            after_rank = states[mask | bit][0]
            transition = f"{before_rank}->{after_rank}"
            transition_counts[transition] += root_multiplicity
            insertion = rank_birth_insertion(
                geometry,
                _mask_active(mask, geometry.n),
                vertex,
                local_radius=local_radius,
            )
            if insertion["delta_rank"] != after_rank - before_rank:
                raise AssertionError("cached and direct rank transitions disagree")
            representatives.setdefault(
                transition,
                {"environment_mask": mask, "insertion": insertion},
            )
            occupied = mask.bit_count()
            for birth in insertion["births"]:
                label = birth["type"]
                birth_counts[label][occupied] += root_multiplicity
                line_counts[label][_line_key(birth)] += root_multiplicity
                h4 = birth["homology_h4"]
                if h4 is not None:
                    homology_cos4_sums[label] += root_multiplicity * Fraction(h4["cos4"])
                    homology_sin4_sums[label] += root_multiplicity * Fraction(h4["sin4"])
                local = insertion["local_geometry_mark"]
                if local["status"] == "evaluated":
                    for key in ("axis", "diagonal", "both", "landed", "h4"):
                        local_mark_sums[label][key] += root_multiplicity * local[key]

    influence_polynomials: dict[str, Polynomial] = {}
    for label, counts in birth_counts.items():
        polynomial = [Fraction(0)] * geometry.n
        for occupied, count in enumerate(counts):
            polynomial = _add_polynomial(
                polynomial,
                _scale_polynomial(count, _bernstein_term(occupied, geometry.n - 1)),
            )
        influence_polynomials[label] = polynomial

    matching_polynomial = _rank_polynomial(geometry, states)
    matching_derivative = _differentiate(matching_polynomial)
    gate_sum = _add_polynomial(
        influence_polynomials["0_to_1"], influence_polynomials["1_to_2"]
    )
    if matching_derivative != gate_sum:
        raise AssertionError("exact Russo gate decomposition failed")

    half = Fraction(1, 2)
    return {
        "id": name,
        "N": geometry.n,
        "period_matrix": [list(row) for row in geometry.periods.matrix],
        "root_sampling": {
            "roots": list(roots),
            "multiplicity_per_environment": root_multiplicity,
            "exact_translation_reduction": root_multiplicity > 1,
        },
        "environment_insertions_weighted": sum(transition_counts.values()),
        "transition_counts": dict(sorted(transition_counts.items())),
        "birth_bernstein_counts_by_other_occupancy": birth_counts,
        "line_index_counts": {
            label: dict(sorted(counts.items())) for label, counts in line_counts.items()
        },
        "homology_line_h4_raw_sums": {
            label: {
                "cos4": _fraction_text(homology_cos4_sums[label]),
                "sin4": _fraction_text(homology_sin4_sums[label]),
            }
            for label in ("0_to_1", "1_to_2")
        },
        "local_landing_h4_raw_sums": {
            label: dict(local_mark_sums[label])
            for label in ("0_to_1", "1_to_2")
        },
        "power_basis_coefficients_ascending": {
            "M": _polynomial_text(matching_polynomial),
            "M_prime_direct": _polynomial_text(matching_derivative),
            "f_01": _polynomial_text(influence_polynomials["0_to_1"]),
            "f_12": _polynomial_text(influence_polynomials["1_to_2"]),
            "f_01_plus_f_12": _polynomial_text(gate_sum),
        },
        "p_equals_half": {
            "M_prime": _fraction_text(_evaluate(matching_derivative, half)),
            "f_01": _fraction_text(_evaluate(influence_polynomials["0_to_1"], half)),
            "f_12": _fraction_text(_evaluate(influence_polynomials["1_to_2"], half)),
        },
        "representatives": representatives,
    }


def geometry_specs() -> list[dict[str, Any]]:
    axis_l2 = axis_integer_torus(2)
    gaussian = gaussian_integer_torus(2, 1)
    axis_l4 = axis_integer_torus(4)
    return [
        {
            "name": "axis-L2-degenerate",
            "geometry": axis_l2,
            "roots": tuple(range(axis_l2.n)),
            "root_multiplicity": 1,
            "local_radius": None,
        },
        {
            "name": "gaussian-2-1",
            "geometry": gaussian,
            "roots": tuple(range(gaussian.n)),
            "root_multiplicity": 1,
            "local_radius": None,
        },
        {
            "name": "axis-L4-fixed-root",
            "geometry": axis_l4,
            "roots": (axis_l4.vertex((0, 0)),),
            "root_multiplicity": axis_l4.n,
            "local_radius": 1,
        },
    ]


def build_artifact() -> dict[str, Any]:
    geometries = [_geometry_summary(**spec) for spec in geometry_specs()]
    for geometry in geometries:
        coefficients = geometry["power_basis_coefficients_ascending"]
        if coefficients["M_prime_direct"] != coefficients["f_01_plus_f_12"]:
            raise AssertionError("serialized Russo identity failed")
    return {
        "schema": "matching-one/homology-rank-birth-insertion/v1",
        "issues": [215, 269, 276],
        "status": "tiny_exact_typed_insertion_oracle",
        "pointwise_theorem": {
            "definitions": {
                "I_01": "1[r(v=0)=0 and r(v=1)>=1]",
                "I_12": "1[r(v=0)<=1 and r(v=1)=2]",
            },
            "identity": "Delta_v r = I_01 + I_12",
            "russo": "M'(p)=d E_p[r]/dp=f_01(p)+f_12(p)",
            "direct_0_to_2": "counted once in each gate; no canonical intermediate ell",
        },
        "typed_line_semantics": {
            "0_to_1": "ell is the new primitive rank-one image; iota is its saturation index after insertion",
            "1_to_2": "ell is the canonical rank-one plateau line immediately before insertion; iota is its endpoint index",
            "shared_interval_mark": "on a nonempty rank-one plateau the two endpoint insertions carry the same ell",
        },
        "geometries": geometries,
        "claim_boundary": {
            "exact": (
                "pointwise gate decomposition for every monotone rank in {0,1,2}; exact exhaustive "
                "polynomial certificates on the declared geometries"
            ),
            "not_claimed": (
                "a canonical ell for a simultaneous 0->2 jump, production SNF streaming, "
                "or a continuum operator identification"
            ),
        },
    }


def render_markdown(artifact: dict[str, Any]) -> str:
    lines = [
        "# Canonical homology-rank-birth insertion",
        "",
        "The two elementary gates are typed before any CFT interpretation:",
        "",
        "```text",
        "I_01 = 1[r0=0 and r1>=1]",
        "I_12 = 1[r0<=1 and r1=2]",
        "Delta_v r = I_01 + I_12.",
        "```",
        "",
        "Thus Russo differentiation gives exactly `M'(p)=f_01(p)+f_12(p)`. A direct",
        "`0->2` jump contributes once to each gate and needs no artificial intermediate state.",
        "",
        "| geometry | weighted insertions | transitions | M'(1/2) | f01(1/2) | f12(1/2) |",
        "|---|---:|---|---:|---:|---:|",
    ]
    for row in artifact["geometries"]:
        transitions = ", ".join(f"{key}:{value}" for key, value in row["transition_counts"].items())
        half = row["p_equals_half"]
        lines.append(
            f"| {row['id']} | {row['environment_insertions_weighted']} | {transitions} | "
            f"{half['M_prime']} | {half['f_01']} | {half['f_12']} |"
        )
    lines.extend(
        [
            "",
            "For `0->1`, `ell` is the new primitive rank-one image and `iota` its integral",
            "saturation index. For `1->2`, the canonical mark is the rank-one plateau line",
            "immediately before the second birth. The same `ell` labels both endpoints of every",
            "nonempty essential-H1 interval. A simultaneous `0->2` jump has no canonical",
            "intermediate line and is recorded with `ell=null`.",
            "",
            "The axis-L4 control also attaches the existing radius-one landing-sector H4 mark",
            "and the exact spin-4 harmonic of the physical winding line. These remain typed",
            "coordinates of the same insertion, not extra evidence rows.",
            "",
        ]
    )
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    artifact = build_artifact()
    rendered = (
        json.dumps(artifact, indent=2, sort_keys=True) + "\n"
        if args.format == "json"
        else render_markdown(artifact)
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
