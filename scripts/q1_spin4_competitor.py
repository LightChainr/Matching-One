#!/usr/bin/env python3
"""Exact Q=1 spectrum preflight for the x=17/4 spin-4 competitor."""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path


def delta(r: Fraction, s: Fraction, beta_squared: Fraction) -> Fraction:
    """Coulomb-gas weight without introducing an irrational beta."""
    momentum_square = r * r * beta_squared - 2 * r * s + s * s / beta_squared
    vacuum_square = beta_squared - 2 + 1 / beta_squared
    return (momentum_square - vacuum_square) / 4


def primary(r: Fraction, s: Fraction, beta_squared: Fraction) -> dict:
    left = delta(r, s, beta_squared)
    right = delta(r, -s, beta_squared)
    dimension = left + right
    spin = left - right
    return {
        "r": fraction_record(r),
        "s": fraction_record(s),
        "leg_count": int(2 * r),
        "left_weight": fraction_record(left),
        "right_weight": fraction_record(right),
        "scaling_dimension": fraction_record(dimension),
        "conformal_spin": fraction_record(spin),
        "absolute_spin": fraction_record(abs(spin)),
    }


def fraction_record(value: Fraction) -> dict:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "text": str(value),
        "decimal": float(value),
    }


def analyze(config: dict) -> dict:
    beta_squared = Fraction(*config["beta_squared"])
    r = Fraction(*config["competitor"]["r"])
    s = Fraction(*config["competitor"]["s"])
    competitor = primary(r, s, beta_squared)
    conjugate = primary(r, -s, beta_squared)
    thermal_dimension = Fraction(*config["thermal_q4_dimension"])
    competitor_dimension = Fraction(
        competitor["scaling_dimension"]["numerator"],
        competitor["scaling_dimension"]["denominator"],
    )
    gap = thermal_dimension - competitor_dimension
    dilations = []
    for q in config["gaussian_area_multipliers"]:
        length_dilation = math.sqrt(q)
        dilations.append({
            "area_multiplier": q,
            "length_dilation": length_dilation,
            "x17_amplitude_factor": length_dilation ** (-float(competitor_dimension)),
            "x21_amplitude_factor": length_dilation ** (-float(thermal_dimension)),
            "thermal_to_competitor_relative_factor": length_dilation ** (-float(gap)),
        })
    central_charge = Fraction(13) - 6 * beta_squared - 6 / beta_squared
    return {
        "schema_version": 1,
        "issue": 257,
        "beta_squared": fraction_record(beta_squared),
        "central_charge": fraction_record(central_charge),
        "competitor": competitor,
        "chirality_conjugate": conjugate,
        "thermal_q4_dimension": fraction_record(thermal_dimension),
        "dimension_gap_thermal_minus_competitor": fraction_record(gap),
        "continuum_dilation_oracle": dilations,
        "claim_boundary": config["claim_boundary"],
        "unresolved_gates": config["unresolved_gates"],
        "sources": config["sources"],
    }


def render_markdown(result: dict) -> str:
    c = result["competitor"]
    lines = [
        "# Q=1 spin-4 competitor preflight",
        "",
        "This is exact generic-loop spectrum arithmetic, not a Potts multiplicity or lattice-overlap claim.",
        "",
        "## Exact field data",
        "",
        f"- beta^2 = {result['beta_squared']['text']}; c = {result['central_charge']['text']}",
        f"- legs = {c['leg_count']}",
        f"- (Delta, DeltaBar) = ({c['left_weight']['text']}, {c['right_weight']['text']})",
        f"- x = {c['scaling_dimension']['text']}; spin = {c['conformal_spin']['text']}",
        f"- x(Q4 epsilon) - x(V_(2,2)) = {result['dimension_gap_thermal_minus_competitor']['text']}",
        "",
        "## Continuum two-field dilation oracle",
        "",
        "| area multiplier Q | length dilation | x=17/4 factor | x=21/4 factor | relative Q4/four-leg factor |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in result["continuum_dilation_oracle"]:
        lines.append(
            f"| {row['area_multiplier']} | {row['length_dilation']:.8g} | "
            f"{row['x17_amplitude_factor']:.8g} | {row['x21_amplitude_factor']:.8g} | "
            f"{row['thermal_to_competitor_relative_factor']:.8g} |"
        )
    lines += ["", "## Unresolved gates", ""]
    lines.extend(f"- {gate['id']}: {gate['status']}" for gate in result["unresolved_gates"])
    lines += [
        "",
        "The relative factor follows only from the exact dimension gap. It must not be applied to the normalized local shell observable until that observable's lattice-to-radial normalization is derived.",
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
