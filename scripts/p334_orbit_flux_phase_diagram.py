#!/usr/bin/env python3
"""Exact cross-p orbit-resolved source/sink phase diagram for Issue #334."""

from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import json
from math import comb
from pathlib import Path
from typing import Iterable, Optional, Sequence

from p334_n13_multiorbit_flux import P_REF, exact_census


Poly = tuple[Fraction, ...]  # ascending monomial coefficients
LABELS = ("axis_orbit", "diagonal_orbit")


def _trim(poly: Iterable[Fraction | int]) -> Poly:
    values = [Fraction(value) for value in poly]
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values or [Fraction(0)])


def _eval(poly: Poly, x: Fraction) -> Fraction:
    value = Fraction(0)
    for coefficient in reversed(poly):
        value = value * x + coefficient
    return value


def _derivative(poly: Poly) -> Poly:
    return _trim(index * value for index, value in enumerate(poly[1:], start=1))


def _scale(poly: Poly, scalar: Fraction | int) -> Poly:
    return _trim(Fraction(scalar) * value for value in poly)


def _sub(left: Poly, right: Poly) -> Poly:
    size = max(len(left), len(right))
    return _trim(
        (left[index] if index < len(left) else 0)
        - (right[index] if index < len(right) else 0)
        for index in range(size)
    )


def _divmod(numerator: Poly, denominator: Poly) -> tuple[Poly, Poly]:
    if denominator == (0,):
        raise ZeroDivisionError("polynomial division by zero")
    remainder = list(numerator)
    quotient = [Fraction(0)] * max(1, len(numerator) - len(denominator) + 1)
    while len(remainder) >= len(denominator) and any(remainder):
        degree = len(remainder) - len(denominator)
        factor = remainder[-1] / denominator[-1]
        quotient[degree] += factor
        for index, value in enumerate(denominator):
            remainder[degree + index] -= factor * value
        while len(remainder) > 1 and remainder[-1] == 0:
            remainder.pop()
    return _trim(quotient), _trim(remainder)


def _monic(poly: Poly) -> Poly:
    poly = _trim(poly)
    if poly == (0,):
        return poly
    return _scale(poly, 1 / poly[-1])


def _gcd(left: Poly, right: Poly) -> Poly:
    left, right = _trim(left), _trim(right)
    while right != (0,):
        _, remainder = _divmod(left, right)
        left, right = right, remainder
    return _monic(left)


def _square_free(poly: Poly) -> Poly:
    common = _gcd(poly, _derivative(poly))
    quotient, remainder = _divmod(poly, common)
    if remainder != (0,):
        raise AssertionError("polynomial gcd did not divide exactly")
    return _trim(quotient)


def _sturm(poly: Poly) -> list[Poly]:
    sequence = [_square_free(poly)]
    sequence.append(_derivative(sequence[0]))
    while sequence[-1] != (0,):
        _, remainder = _divmod(sequence[-2], sequence[-1])
        if remainder == (0,):
            break
        sequence.append(_scale(remainder, -1))
    return sequence


def _variations(sequence: Sequence[Poly], x: Fraction) -> int:
    signs = []
    for poly in sequence:
        value = _eval(poly, x)
        if value:
            signs.append(1 if value > 0 else -1)
    return sum(left != right for left, right in zip(signs, signs[1:]))


def _bernstein_edge_poly(
    census: dict[str, object], label: str, key: str
) -> Poly:
    degree = census["geometry"]["N"] - 1
    coefficients = [Fraction(0)] * (degree + 1)
    for row in census["coefficient_rows"]:
        k = row["lower_subset_size"]
        count = row[label][key]
        for offset in range(degree - k + 1):
            coefficients[k + offset] += count * (-1) ** offset * comb(
                degree - k, offset
            )
    return _trim(coefficients)


def _strip_endpoint_roots(poly: Poly) -> tuple[Poly, int, int]:
    zero_order = 0
    while len(poly) > 1 and poly[0] == 0:
        poly = _trim(poly[1:])
        zero_order += 1
    one_order = 0
    factor = (Fraction(-1), Fraction(1))
    while len(poly) > 1 and _eval(poly, Fraction(1)) == 0:
        quotient, remainder = _divmod(poly, factor)
        if remainder != (0,):
            raise AssertionError("p=1 factor did not divide")
        poly = quotient
        one_order += 1
    return poly, zero_order, one_order


class _ExactDyadicRoot(RuntimeError):
    def __init__(self, root: Fraction) -> None:
        self.root = root


def _isolate_no_exact_midpoint(poly: Poly, bits: int) -> list[tuple[Fraction, Fraction]]:
    sequence = _sturm(poly)

    def count(left: Fraction, right: Fraction) -> int:
        return _variations(sequence, left) - _variations(sequence, right)

    pending = [(Fraction(0), Fraction(1), count(Fraction(0), Fraction(1)))]
    isolated = []
    target = Fraction(1, 1 << bits)
    while pending:
        left, right, roots = pending.pop()
        if roots == 0:
            continue
        if roots == 1 and right - left <= target:
            isolated.append((left, right))
            continue
        midpoint = (left + right) / 2
        if _eval(poly, midpoint) == 0:
            raise _ExactDyadicRoot(midpoint)
        left_roots = count(left, midpoint)
        right_roots = roots - left_roots
        pending.append((midpoint, right, right_roots))
        pending.append((left, midpoint, left_roots))
    return sorted(isolated)


def _decimal(value: Fraction, digits: int = 30) -> str:
    with localcontext() as context:
        context.prec = digits + 10
        return format(
            Decimal(value.numerator) / Decimal(value.denominator), f".{digits}g"
        )


def isolate_open_unit_roots(poly: Poly, *, bits: int = 112) -> dict[str, object]:
    reduced, zero_order, one_order = _strip_endpoint_roots(poly)
    exact_roots: list[Fraction] = []
    while True:
        try:
            intervals = _isolate_no_exact_midpoint(reduced, bits)
            break
        except _ExactDyadicRoot as found:
            exact_roots.append(found.root)
            factor = (-found.root, Fraction(1))
            while _eval(reduced, found.root) == 0:
                reduced, remainder = _divmod(reduced, factor)
                if remainder != (0,):
                    raise AssertionError("exact dyadic root did not divide")
            reduced, _, _ = _strip_endpoint_roots(reduced)

    rows = [
        {
            "kind": "exact_dyadic",
            "root": str(root),
            "decimal": _decimal(root),
            "lower": str(root),
            "upper": str(root),
        }
        for root in exact_roots
        if 0 < root < 1
    ]
    rows.extend(
        {
            "kind": "rational_isolating_interval",
            "decimal": _decimal((left + right) / 2),
            "lower": str(left),
            "upper": str(right),
            "width": str(right - left),
        }
        for left, right in intervals
    )
    rows.sort(key=lambda row: Fraction(row["lower"]))
    return {
        "degree_after_trim": len(_trim(poly)) - 1,
        "endpoint_zero_order": zero_order,
        "endpoint_one_order": one_order,
        "interior_root_count": len(rows),
        "roots": rows,
    }


def _root_midpoint(row: dict[str, str]) -> Fraction:
    return (Fraction(row["lower"]) + Fraction(row["upper"])) / 2


def _sign(value: Fraction) -> str:
    return "positive" if value > 0 else "negative" if value < 0 else "zero"


def _phase_intervals(polys: dict[str, Poly], roots: dict[str, dict]) -> list[dict]:
    events = []
    for name in ("axis_net", "diagonal_net", "total_net"):
        for row in roots[name]["roots"]:
            events.append((_root_midpoint(row), name))
    events.sort()
    clusters: list[tuple[Fraction, list[str]]] = []
    for value, name in events:
        if clusters and abs(value - clusters[-1][0]) < Fraction(1, 1 << 90):
            clusters[-1][1].append(name)
        else:
            clusters.append((value, [name]))
    boundaries = [Fraction(0)] + [value for value, _ in clusters] + [Fraction(1)]
    rows = []
    for left, right in zip(boundaries, boundaries[1:]):
        probe = (left + right) / 2
        axis = _eval(polys["axis_net"], probe)
        diagonal = _eval(polys["diagonal_net"], probe)
        total = _eval(polys["total_net"], probe)
        rows.append(
            {
                "lower": _decimal(left),
                "upper": _decimal(right),
                "axis_net_sign": _sign(axis),
                "diagonal_net_sign": _sign(diagonal),
                "total_net_sign": _sign(total),
                "orbit_contributions": (
                    "reinforce" if axis * diagonal < 0 else "cancel"
                ),
                "axis_signed_share_at_midpoint": (
                    _decimal(axis / total) if total else None
                ),
            }
        )
    return rows


def _point_metrics(polys: dict[str, Poly], p: Fraction) -> dict[str, object]:
    values = {name: _eval(poly, p) for name, poly in polys.items()}
    slopes = {name: _eval(_derivative(poly), p) for name, poly in polys.items()}
    axis_share = values["axis_net"] / values["total_net"]
    axis_share_slope = (
        slopes["axis_net"] * values["total_net"]
        - values["axis_net"] * slopes["total_net"]
    ) / values["total_net"] ** 2
    return {
        "p": str(p),
        "values": {name: _decimal(value) for name, value in values.items()},
        "slopes_dp": {name: _decimal(value) for name, value in slopes.items()},
        "axis_signed_share": _decimal(axis_share),
        "axis_signed_share_slope_dp": _decimal(axis_share_slope),
        "activity_cancellation": {
            label: _decimal(
                abs(values[f"{label}_net"])
                / (values[f"{label}_birth"] + values[f"{label}_exit"])
            )
            for label in ("axis", "diagonal")
        },
        "orbit_composition": {
            "axis_fraction_of_birth": _decimal(
                values["axis_birth"]
                / (values["axis_birth"] + values["diagonal_birth"])
            ),
            "axis_fraction_of_exit": _decimal(
                values["axis_exit"]
                / (values["axis_exit"] + values["diagonal_exit"])
            ),
        },
    }


def _geometry_payload(a: int, b: int) -> dict[str, object]:
    census = exact_census(a, b, include_direct_rank2=True)
    axis_birth = _bernstein_edge_poly(census, "axis_orbit", "birth_edges")
    axis_exit = _bernstein_edge_poly(census, "axis_orbit", "exit_edges")
    diagonal_birth = _bernstein_edge_poly(census, "diagonal_orbit", "birth_edges")
    diagonal_exit = _bernstein_edge_poly(census, "diagonal_orbit", "exit_edges")
    polys = {
        "axis_birth": axis_birth,
        "axis_exit": axis_exit,
        "axis_net": _sub(axis_birth, axis_exit),
        "diagonal_birth": diagonal_birth,
        "diagonal_exit": diagonal_exit,
        "diagonal_net": _sub(diagonal_birth, diagonal_exit),
        "birth_character_total": _sub(axis_birth, diagonal_birth),
        "exit_character_total": _sub(axis_exit, diagonal_exit),
    }
    polys["total_net"] = _sub(polys["axis_net"], polys["diagonal_net"])
    roots = {name: isolate_open_unit_roots(poly) for name, poly in polys.items()}
    nearest = {}
    for name in ("axis_net", "diagonal_net", "total_net"):
        row = min(
            roots[name]["roots"],
            key=lambda candidate: abs(_root_midpoint(candidate) - P_REF),
        )
        root = _root_midpoint(row)
        nearest[name] = {
            "root": row["decimal"],
            "p_ref_minus_root": _decimal(P_REF - root),
            "slope_at_root": _decimal(_eval(_derivative(polys[name]), root)),
        }
    return {
        "id": census["geometry"]["id"],
        "N": census["geometry"]["N"],
        "period_matrix": census["geometry"]["period_matrix"],
        "chi4_axis": census["orbits"]["axis_orbit"]["chi4"],
        "polynomial_basis": "ascending monomial expansion of exact degree-(N-1) subset-boundary Bernstein sums",
        "roots": roots,
        "signed_share_singularities": {
            "axis_share_zero": roots["axis_net"]["roots"],
            "diagonal_share_zero": roots["diagonal_net"]["roots"],
            "common_pole": roots["total_net"]["roots"],
        },
        "phase_intervals": _phase_intervals(polys, roots),
        "p_ref_metrics": _point_metrics(polys, P_REF),
        "nearest_roots_to_p_ref": nearest,
    }


def build_certificate() -> dict[str, object]:
    n13 = _geometry_payload(3, 2)
    n17 = _geometry_payload(4, 1)
    return {
        "schema": "matching-one/p334-orbit-flux-phase-diagram/v1",
        "issue": 334,
        "parents": ["b8e286e", "e34140d"],
        "status": "exact_orbit_source_sink_phase_diagram",
        "geometries": {"N13": n13, "N17": n17},
        "mechanism_classification": {
            "class": "common_activity_with_orbit_composition_counterflow",
            "statement": (
                "The 76/24 split is not birth-dominated. At p_ref each orbit net "
                "is only a small residual of two large positive birth/exit currents. "
                "Birth and exit have similar axis/diagonal composition, and their "
                "small composition skew reverses between N13 and N17, making both "
                "orbit nets reverse together while their reinforcing signed share "
                "stays close. The signed-share slope itself reverses and roughly "
                "doubles, so 76/24 is not a geometry-independent constant."
            ),
            "falsifiable_prediction": (
                "On the next independently chosen two-orbit Gaussian quotient, the "
                "axis and diagonal net-flux zeros should remain a close ordered pair; "
                "between them the two chi4-weighted orbit contributions reinforce, "
                "while outside them they cancel. A quotient that lacks this paired-zero "
                "window falsifies the common-activity/counterflow classification."
            ),
        },
        "exactness_gates": {
            "birth_and_exit_have_no_interior_roots": all(
                geometry["roots"][f"{label}_{flow}"]["interior_root_count"] == 0
                for geometry in (n13, n17)
                for label in ("axis", "diagonal")
                for flow in ("birth", "exit")
            ),
            "character_source_and_sink_have_no_interior_roots": all(
                geometry["roots"][name]["interior_root_count"] == 0
                for geometry in (n13, n17)
                for name in ("birth_character_total", "exit_character_total")
            ),
            "all_root_intervals_width_below_2^-111": all(
                Fraction(row.get("width", "0")) <= Fraction(1, 1 << 111)
                for geometry in (n13, n17)
                for root_set in geometry["roots"].values()
                for row in root_set["roots"]
            ),
        },
        "claim_boundary": [
            "All root counts and isolating intervals come from exact rational polynomial arithmetic.",
            "Decimal root locations and slopes summarize exact rational brackets; they are not fitted values.",
            "The next-quotient statement is a mechanism prediction, not evidence already supplied by N13/N17.",
            "No new size enumeration, Monte Carlo sample, Huawei production, PR, or merge is used.",
        ],
    }


def render_markdown(payload: dict[str, object]) -> str:
    lines = [
        "# Orbit-resolved source/sink phase diagram",
        "",
        f"Status: `{payload['status']}`.",
        "",
        payload["mechanism_classification"]["statement"],
        "",
    ]
    for key in ("N13", "N17"):
        geometry = payload["geometries"][key]
        point = geometry["p_ref_metrics"]
        lines += [
            f"## {key}: {geometry['id']}",
            "",
            "| exact curve | interior roots in (0,1) | locations |",
            "|---|---:|---|",
        ]
        for name in (
            "axis_birth",
            "axis_exit",
            "diagonal_birth",
            "diagonal_exit",
            "axis_net",
            "diagonal_net",
            "birth_character_total",
            "exit_character_total",
            "total_net",
        ):
            root_set = geometry["roots"][name]
            locations = ", ".join(row["decimal"] for row in root_set["roots"])
            lines.append(
                f"| {name} | {root_set['interior_root_count']} | {locations or 'none'} |"
            )
        lines += [
            "",
            f"At `p_ref`, axis signed share = `{point['axis_signed_share']}` with slope "
            f"`{point['axis_signed_share_slope_dp']}`.",
            "",
            f"Activity residual fractions `|birth-exit|/(birth+exit)`: axis "
            f"`{point['activity_cancellation']['axis']}`, diagonal "
            f"`{point['activity_cancellation']['diagonal']}`.",
            "",
            f"Axis composition: birth `{point['orbit_composition']['axis_fraction_of_birth']}`, "
            f"exit `{point['orbit_composition']['axis_fraction_of_exit']}`.",
            "",
            "Thus the signed-share zero/pole map is: axis zero at "
            f"`{geometry['signed_share_singularities']['axis_share_zero'][0]['decimal']}`, "
            "diagonal zero at "
            f"`{geometry['signed_share_singularities']['diagonal_share_zero'][0]['decimal']}`, "
            "and common pole at "
            f"`{geometry['signed_share_singularities']['common_pole'][0]['decimal']}`.",
            "",
            "Nearest roots to `p_ref`:",
            "",
        ]
        for name, row in geometry["nearest_roots_to_p_ref"].items():
            lines.append(
                f"- {name}: `{row['root']}`; `p_ref-root={row['p_ref_minus_root']}`; "
                f"slope `{row['slope_at_root']}`."
            )
        lines += ["", "Phase intervals:", ""]
        for row in geometry["phase_intervals"]:
            lines.append(
                f"- `({row['lower']}, {row['upper']})`: axis {row['axis_net_sign']}, "
                f"diagonal {row['diagonal_net_sign']}, total {row['total_net_sign']}; "
                f"orbit contributions **{row['orbit_contributions']}**."
            )
        lines.append("")
    lines += [
        "## Next falsifiable prediction",
        "",
        payload["mechanism_classification"]["falsifiable_prediction"],
        "",
        "## Boundary",
        "",
    ]
    lines.extend(f"- {item}" for item in payload["claim_boundary"])
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    payload = build_certificate()
    rendered = (
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
        if args.format == "json"
        else render_markdown(payload) + "\n"
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
