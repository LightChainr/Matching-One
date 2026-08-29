#!/usr/bin/env python3
"""Freeze E4-balanced Pell weights without reading root outcomes."""

from __future__ import annotations

import argparse
import cmath
import json
import math
from pathlib import Path


def sigma3(n: int) -> int:
    return sum(d**3 for d in range(1, n + 1) if n % d == 0)


def eisenstein_e4(tau: complex, terms: int = 96) -> complex:
    """Return E4(tau)=1+240 sum sigma_3(n) exp(2 pi i n tau)."""
    q = cmath.exp(2j * math.pi * tau)
    return 1.0 + 240.0 * sum(sigma3(n) * q**n for n in range(1, terms + 1))


def pell_residual(x: int, m: int) -> int:
    return x * x - 3 * m * m


def shape_amplitude(x: int, m: int, terms: int) -> dict:
    linear_scale = 2 * m
    tau = complex(0.5, x / (2 * m))
    e4 = eisenstein_e4(tau, terms)
    if abs(e4.imag) > 1e-12:
        raise ValueError("reflection-symmetric Pell modulus should have real E4")
    return {
        "x": x,
        "m": m,
        "pell_residual": pell_residual(x, m),
        "linear_scale": linear_scale,
        "tau": [tau.real, tau.imag],
        "E4": e4.real,
        "a4": e4.real / linear_scale**4,
    }


def balanced_row(pair: dict, terms: int, scalar_exponent: int) -> dict:
    negative = shape_amplitude(*pair["negative"], terms)
    positive = shape_amplitude(*pair["positive"], terms)
    if negative["pell_residual"] != -2 or positive["pell_residual"] != 1:
        raise ValueError("expected Pell residuals -2 and +1")
    if not negative["E4"] < 0 < positive["E4"]:
        raise ValueError("Pell pair must straddle the E4 zero")
    a_minus, a_plus = negative["a4"], positive["a4"]
    denominator = a_plus - a_minus
    w_plus = -a_minus / denominator
    w_minus = a_plus / denominator
    cancellation = w_plus * a_plus + w_minus * a_minus
    scalar_coefficient = (
        w_plus * positive["linear_scale"] ** (-scalar_exponent)
        + w_minus * negative["linear_scale"] ** (-scalar_exponent)
    )
    return {
        "generation": pair["generation"],
        "negative": negative,
        "positive": positive,
        "weights": {"negative": w_minus, "positive": w_plus},
        "weight_sum": w_minus + w_plus,
        "h4_cancellation_residual": cancellation,
        "h4_cancellation_relative": abs(cancellation) / max(abs(a_minus), abs(a_plus)),
        "scalar_coefficient": scalar_coefficient,
    }


def analyze(config: dict) -> dict:
    terms = int(config["e4_q_series_terms"])
    exponent = int(config["scalar_exponent"])
    rows = [balanced_row(pair, terms, exponent) for pair in config["pell_pairs"]]
    for index in range(len(rows) - 1):
        rows[index]["next_scalar_coefficient_ratio"] = (
            rows[index + 1]["scalar_coefficient"] / rows[index]["scalar_coefficient"]
        )
    rows[-1]["next_scalar_coefficient_ratio"] = None
    unit = 2.0 + math.sqrt(3.0)
    return {
        "schema_version": 1,
        "issue": 160,
        "selection_boundary": config["selection_boundary"],
        "e4_q_series_terms": terms,
        "scalar_exponent": exponent,
        "fundamental_unit": unit,
        "asymptotic_scalar_ratio": unit ** (-exponent),
        "rows": rows,
    }


def render_markdown(result: dict) -> str:
    lines = [
        "# E4-balanced Pell estimator freeze",
        "",
        "Weights use only Pell geometry and the deterministic E4 q-series; no root outcome is an input.",
        "",
        "| generation | negative (x,m) | positive (x,m) | w- | w+ | relative H4 residual | next scalar ratio |",
        "|---:|---|---|---:|---:|---:|---:|",
    ]
    for row in result["rows"]:
        ratio = row["next_scalar_coefficient_ratio"]
        lines.append(
            "| {generation} | ({nx},{nm}) | ({px},{pm}) | {wn:.12f} | {wp:.12f} | {res:.3e} | {ratio} |".format(
                generation=row["generation"],
                nx=row["negative"]["x"], nm=row["negative"]["m"],
                px=row["positive"]["x"], pm=row["positive"]["m"],
                wn=row["weights"]["negative"], wp=row["weights"]["positive"],
                res=row["h4_cancellation_relative"],
                ratio="-" if ratio is None else f"{ratio:.12g}",
            )
        )
    lines += [
        "",
        f"Frozen asymptotic L^-{result['scalar_exponent']} generation ratio: "
        f"`{result['asymptotic_scalar_ratio']:.14g}`.",
        "",
        "The weights must be applied inside every delete-one replicate. H4-null scoring precedes any scalar-law score.",
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
