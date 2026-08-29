#!/usr/bin/env python3
"""Exact Q-velocity fingerprints for Q=1 spin-sector competitors."""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path


def fraction_record(value: Fraction) -> dict:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "text": str(value),
        "decimal": float(value),
    }


def radical_velocity_record(coefficient: Fraction) -> dict:
    """Represent coefficient * sqrt(3) / pi exactly up to the named constants."""
    numerator = coefficient.numerator
    denominator = coefficient.denominator
    sign = "-" if numerator < 0 else ""
    magnitude = abs(numerator)
    factor = "" if magnitude == 1 else f"{magnitude}*"
    if denominator == 1:
        expression = f"{sign}{factor}sqrt(3)/pi"
    else:
        expression = f"{sign}{factor}sqrt(3)/({denominator}*pi)"
    return {
        "coefficient_of_sqrt3_over_pi": fraction_record(coefficient),
        "expression": expression,
        "decimal": float(coefficient) * math.sqrt(3) / math.pi,
    }


def loop_dimension(r: Fraction, s: Fraction, u: Fraction) -> Fraction:
    return 1 + ((r * r - 1) * u + (s * s - 1) / u) / 2


def loop_dx_du(r: Fraction, s: Fraction, u: Fraction) -> Fraction:
    return ((r * r - 1) - (s * s - 1) / (u * u)) / 2


def energy_dimension(u: Fraction) -> Fraction:
    return Fraction(3, 2) / u - 1


def energy_dx_du(u: Fraction) -> Fraction:
    return -Fraction(3, 2) / (u * u)


def q_velocity_coefficient(dx_du: Fraction) -> Fraction:
    """At u=2/3, du/dQ=sqrt(3)/(6*pi)."""
    return dx_du / 6


def loop_field(identifier: str, r: Fraction, s: Fraction, u: Fraction) -> dict:
    x = loop_dimension(r, s, u)
    dx_du = loop_dx_du(r, s, u)
    velocity = q_velocity_coefficient(dx_du)
    return {
        "id": identifier,
        "family": "generic_loop_primary",
        "r": fraction_record(r),
        "s": fraction_record(s),
        "leg_count": int(2 * r),
        "conformal_spin": fraction_record(-r * s),
        "scaling_dimension": fraction_record(x),
        "dx_du_at_q1": fraction_record(dx_du),
        "dx_dQ_at_q1": radical_velocity_record(velocity),
        "dlog_transfer_dQ_per_log_area": radical_velocity_record(-velocity / 2),
        "potts_multiplicity": "unresolved",
        "field_normalization_derivative": "unresolved",
    }


def thermal_descendant(identifier: str, u: Fraction, level: int) -> dict:
    primary_x = energy_dimension(u)
    dx_du = energy_dx_du(u)
    velocity = q_velocity_coefficient(dx_du)
    return {
        "id": identifier,
        "family": "thermal_energy_descendant",
        "primary": "energy_epsilon",
        "left_descendant_level": level,
        "leg_count": 0,
        "conformal_spin": fraction_record(Fraction(level)),
        "scaling_dimension": fraction_record(primary_x + level),
        "dx_du_at_q1": fraction_record(dx_du),
        "dx_dQ_at_q1": radical_velocity_record(velocity),
        "dlog_transfer_dQ_per_log_area": radical_velocity_record(-velocity / 2),
        "potts_multiplicity": "singlet_family_expected_but_lattice_overlap_unresolved",
        "field_normalization_derivative": "unresolved",
    }


def analyze(config: dict) -> dict:
    u = Fraction(*config["u_at_q1"])
    if u != Fraction(2, 3):
        raise ValueError("v1 radical normalization is frozen at u=2/3")
    fields = []
    for row in config["generic_loop_fields"]:
        fields.append(loop_field(
            row["id"], Fraction(*row["r"]), Fraction(*row["s"]), u
        ))
    fields.append(thermal_descendant(
        config["thermal_descendant"]["id"], u,
        int(config["thermal_descendant"]["left_level"]),
    ))
    by_id = {field["id"]: field for field in fields}
    four_leg = by_id["loop_V_2_2"]
    thermal = by_id["thermal_Q4_epsilon"]
    x_gap = Fraction(four_leg["scaling_dimension"]["numerator"], four_leg["scaling_dimension"]["denominator"]) - Fraction(
        thermal["scaling_dimension"]["numerator"], thermal["scaling_dimension"]["denominator"]
    )
    v_four = Fraction(**{
        "numerator": four_leg["dx_dQ_at_q1"]["coefficient_of_sqrt3_over_pi"]["numerator"],
        "denominator": four_leg["dx_dQ_at_q1"]["coefficient_of_sqrt3_over_pi"]["denominator"],
    })
    v_thermal = Fraction(**{
        "numerator": thermal["dx_dQ_at_q1"]["coefficient_of_sqrt3_over_pi"]["numerator"],
        "denominator": thermal["dx_dQ_at_q1"]["coefficient_of_sqrt3_over_pi"]["denominator"],
    })
    velocity_gap = v_four - v_thermal
    transfer_examples = []
    for area_multiplier in config["area_multipliers"]:
        transfer_examples.append({
            "area_multiplier": area_multiplier,
            "log_area": math.log(area_multiplier),
            "four_leg_dlog_transfer_dQ": -0.5 * radical_velocity_record(v_four)["decimal"] * math.log(area_multiplier),
            "thermal_dlog_transfer_dQ": -0.5 * radical_velocity_record(v_thermal)["decimal"] * math.log(area_multiplier),
            "relative_thermal_minus_four_leg_dlog_transfer_dQ": -0.5 * (
                radical_velocity_record(v_thermal)["decimal"] - radical_velocity_record(v_four)["decimal"]
            ) * math.log(area_multiplier),
        })
    return {
        "schema_version": 1,
        "issue": 261,
        "u_at_q1": fraction_record(u),
        "dQ_du_at_q1": "2*pi*sqrt(3)",
        "du_dQ_at_q1": "sqrt(3)/(6*pi)",
        "fields": fields,
        "primary_pair": {
            "dimension_four_leg_minus_thermal": fraction_record(x_gap),
            "velocity_four_leg_minus_thermal": radical_velocity_record(velocity_gap),
        },
        "transfer_examples": transfer_examples,
        "angular_aliases_without_q_family": config["angular_aliases_without_q_family"],
        "selection_boundary": config["selection_boundary"],
        "sources": config["sources"],
    }


def render_markdown(result: dict) -> str:
    lines = [
        "# Q=1 spin-sector velocity oracle",
        "",
        "Exact continuum fingerprints; no lattice target or fitted field normalization is used.",
        "",
        "| field | family | legs | spin | x(1) | dx/dQ at 1 |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for field in result["fields"]:
        lines.append(
            f"| {field['id']} | {field['family']} | {field['leg_count']} | "
            f"{field['conformal_spin']['text']} | {field['scaling_dimension']['text']} | "
            f"`{field['dx_dQ_at_q1']['expression']}` |"
        )
    pair = result["primary_pair"]
    lines += [
        "",
        "## Primary discriminator",
        "",
        f"- dimension gap `x_22-x_Q4 = {pair['dimension_four_leg_minus_thermal']['text']}`",
        f"- velocity gap `x'_22-x'_Q4 = {pair['velocity_four_leg_minus_thermal']['expression']}`",
        "",
        "The spin-8/spin-12 generic-loop rows are separately declared `V_(2,4)` and `V_(2,6)` controls. They are not assignments for the experiment-design H8/H12 angular aliases.",
        "",
        "Potts multiplicities, explicit field-definition derivatives, and lattice overlaps remain unresolved inputs to any Q-score measurement.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--json-output", type=Path)
    parser.add_argument("--markdown-output", type=Path)
    args = parser.parse_args()
    result = analyze(json.loads(args.manifest.read_text()))
    if args.json_output:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    else:
        print(json.dumps(result, indent=2, sort_keys=True))
    if args.markdown_output:
        args.markdown_output.parent.mkdir(parents=True, exist_ok=True)
        args.markdown_output.write_text(render_markdown(result))


if __name__ == "__main__":
    main()
