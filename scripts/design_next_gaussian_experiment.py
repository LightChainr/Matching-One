#!/usr/bin/env python3
"""Rank prospective Gaussian lineages by expected information per CPU-second.

Only source amplitudes and variance/runtime calibration runs are consumed.
Candidate target outcomes are neither accepted by the manifest nor read here.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from itertools import combinations
from pathlib import Path


def norm(g):
    return g[0] * g[0] + g[1] * g[1]


def multiply(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def canonical_d4(g):
    return tuple(sorted((abs(g[0]), abs(g[1])), reverse=True))


def smith(g):
    d = math.gcd(abs(g[0]), abs(g[1]))
    return (d, norm(g) // d)


def cos4(g):
    a, b = g
    n = norm(g)
    return Fraction(a**4 - 6 * a * a * b * b + b**4, n * n)


def harmonic(g, order):
    """Return cos(order*theta) for order in {4,8,12}."""
    x = cos4(g)
    if order == 4:
        return x
    if order == 8:
        return 2 * x * x - 1
    if order == 12:
        return 4 * x * x * x - 3 * x
    raise ValueError("supported harmonics are 4, 8, and 12")


def angular_ratio(first, second, multiplier, order):
    source = harmonic(first, order) - harmonic(second, order)
    if source == 0:
        return None
    target = harmonic(multiply(first, multiplier), order) - harmonic(multiply(second, multiplier), order)
    return target / source


def _power_law(points, value_key):
    if len(points) != 2:
        raise ValueError("v1 requires exactly two calibration points")
    x0, x1 = (math.log(float(p["N"])) for p in points)
    y0, y1 = (math.log(float(p[value_key])) for p in points)
    exponent = (y1 - y0) / (x1 - x0)
    coefficient = math.exp(y0 - exponent * x0)
    return coefficient, exponent


def _predict_power_law(model, size):
    coefficient, exponent = model
    return coefficient * size**exponent


def primitive_multipliers(q_max):
    values = []
    limit = int(math.sqrt(q_max)) + 1
    for a in range(1, limit + 1):
        for b in range(-a, a + 1):
            q = a * a + b * b
            if 1 < q <= q_max and math.gcd(a, abs(b)) == 1:
                values.append((a, b))
    return sorted(set(values), key=lambda g: (norm(g), g[0], g[1]))


def _fraction(value):
    if value is None:
        return None
    return {"numerator": value.numerator, "denominator": value.denominator,
            "decimal": float(value)}


def build_candidates(config):
    runtime_model = _power_law(config["calibration"]["runtime"], "cpu_seconds_per_million")
    se_model = _power_law(config["calibration"]["standard_error"], "se_at_1m")
    sample_millions = float(config["sample_millions"])
    rows = []
    for parent in config["parents"]:
        first, second = tuple(parent["first"]), tuple(parent["second"])
        source_delta = float(parent["source_delta_M"])
        for multiplier in primitive_multipliers(int(config["caps"]["multiplier_norm_max"])):
            q = norm(multiplier)
            target_n = int(parent["N"]) * q
            if target_n > int(config["caps"]["target_N_max"]):
                continue
            children = (multiply(first, multiplier), multiply(second, multiplier))
            if canonical_d4(children[0]) == canonical_d4(children[1]):
                continue
            ratios = {f"H{order}": angular_ratio(first, second, multiplier, order)
                      for order in (4, 8, 12)}
            if any(value is None for value in ratios.values()):
                continue
            means = {
                "H4_x21_over_4": source_delta * q**(-13 / 8) * float(ratios["H4"]),
                "H4_x17_over_4": source_delta * q**(-9 / 8) * float(ratios["H4"]),
                "H8_control": source_delta * q**(-13 / 8) * float(ratios["H8"]),
                "H12_alias": source_delta * q**(-13 / 8) * float(ratios["H12"]),
                "zero_effect": 0.0,
            }
            se = _predict_power_law(se_model, target_n) / math.sqrt(sample_millions)
            variance = se * se
            cpu = _predict_power_law(runtime_model, target_n) * sample_millions
            pairwise = {f"{left}__{right}": (means[left] - means[right])**2 / (2 * variance)
                        for left, right in combinations(config["models"], 2)}
            limiting_pair, min_kl = min(pairwise.items(), key=lambda item: item[1])
            h4_h12 = pairwise["H4_x21_over_4__H12_alias"]
            robust_denominator = (float(config["uncertainty"]["variance_inflation"]) *
                                  float(config["uncertainty"]["runtime_inflation"]))
            rows.append({
                "parent_id": parent["id"], "parent_N": parent["N"],
                "parent_pair": [list(first), list(second)],
                "source_delta_M": source_delta,
                "multiplier": list(multiplier), "multiplier_norm": q,
                "target_N": target_n,
                "child_raw": [list(children[0]), list(children[1])],
                "child_canonical": [list(canonical_d4(children[0])), list(canonical_d4(children[1]))],
                "child_smith_invariants": [list(smith(children[0])), list(smith(children[1]))],
                "harmonic_angular_ratios": {key: _fraction(value) for key, value in ratios.items()},
                "predicted_means": means,
                "estimated_se": se, "estimated_cpu_seconds": cpu,
                "limiting_live_model_pair": limiting_pair,
                "minimum_live_pairwise_kl": min_kl,
                "minimum_live_kl_per_cpu_second": min_kl / cpu,
                "h4_vs_h12_kl": h4_h12,
                "h4_vs_h12_kl_per_cpu_second": h4_h12 / cpu,
                "robust_minimum_kl_per_cpu_second": min_kl / cpu / robust_denominator,
            })
    rows.sort(key=lambda row: (-row["h4_vs_h12_kl_per_cpu_second"], row["target_N"]))
    for rank, row in enumerate(rows, 1):
        row["h4_vs_h12_rank"] = rank
    robust_order = sorted(rows, key=lambda row: -row["robust_minimum_kl_per_cpu_second"])
    for rank, row in enumerate(robust_order, 1):
        row["robust_live_model_rank"] = rank
    return rows, runtime_model, se_model


def design(config):
    rows, runtime_model, se_model = build_candidates(config)
    norm5 = [row for row in rows if (row["parent_id"], tuple(row["multiplier"])) in {
        ("N65", (2, -1)), ("N85", (2, 1))
    }]
    return {
        "schema_version": 1,
        "selection_boundary": config["selection_boundary"],
        "calibration_provenance": config["calibration_provenance"],
        "runtime_power_law": {"coefficient": runtime_model[0], "N_exponent": runtime_model[1]},
        "se_at_1m_power_law": {"coefficient": se_model[0], "N_exponent": se_model[1]},
        "candidate_count": len(rows),
        "candidates": rows,
        "norm5_reference_rows": norm5,
        "output_top_n": int(config["output_top_n"]),
    }


def render_markdown(result):
    top_h4 = result["candidates"][:result["output_top_n"]]
    top_robust = sorted(result["candidates"], key=lambda row: -row["robust_minimum_kl_per_cpu_second"])[:result["output_top_n"]]
    lines = ["# Prospective Gaussian information-per-CPU design", "",
             "No candidate target outcome is an input. Rankings use exact harmonic arithmetic and pilot-only cost/variance fits.", "",
             "## H4 versus H12", "", "| rank | parent | multiplier | Q | target N | KL | CPU s | KL/CPU |", "|---:|---|---|---:|---:|---:|---:|---:|"]
    for row in top_h4:
        lines.append("| {h4_vs_h12_rank} | {parent_id} | {multiplier} | {multiplier_norm} | {target_N} | {h4_vs_h12_kl:.4g} | {estimated_cpu_seconds:.4g} | {h4_vs_h12_kl_per_cpu_second:.4g} |".format(**row))
    lines += ["", "## Worst pair across live models", "", "| rank | parent | multiplier | target N | min KL | robust min KL/CPU |", "|---:|---|---|---:|---:|---:|"]
    for row in top_robust:
        lines.append("| {robust_live_model_rank} | {parent_id} | {multiplier} | {target_N} | {minimum_live_pairwise_kl:.4g} | {robust_minimum_kl_per_cpu_second:.4g} |".format(**row))
    lines += ["", "The robust ranking is diagnostic: a near-degenerate live pair can dominate the maximin score. The H4/H12 table remains the regression target for the successful norm-5 design.", ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--json-output", type=Path)
    parser.add_argument("--markdown-output", type=Path)
    args = parser.parse_args()
    result = design(json.loads(args.manifest.read_text()))
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
