#!/usr/bin/env python3
"""Exact strong-coupling series from the saved N25 (K,g,q) populations."""
from __future__ import annotations

import argparse
import csv
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
SOURCE = "359bde9be45d051b961110e03f0bd70f3ff82b91"
N, DELTA = 25, F(1152, 625)
HASHES = {
    "axis": "2d23fecc98d276d9ad15ad1867199cd308f0570cb5040ef94eb6b923b4c53458",
    "tilted": "225031e612929ed922ba75c55e76703d59990f5283e7ac39b94f022841798da5",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def zero(order):
    return [F(0)] * (order+1)


def add(a, b):
    return [x+y for x, y in zip(a, b)]


def scale(a, c):
    return [c*x for x in a]


def sub(a, b):
    return add(a, scale(b, -1))


def mul(a, b):
    result = zero(len(a)-1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b[:len(a)-i]):
                if y:
                    result[i+j] += x*y
    return result


def divide(a, b):
    if not b[0]:
        raise ArithmeticError("formal division requires a nonzero constant")
    result = zero(len(a)-1)
    for k in range(len(a)):
        result[k] = (a[k]-sum(b[j]*result[k-j] for j in range(1, k+1)))/b[0]
    return result


def read_rows(path, expected):
    if sha(path) != expected:
        raise ValueError(f"saved histogram hash differs: {path}")
    with path.open(newline="") as handle:
        rows = [{key: int(value) for key, value in row.items()}
                for row in csv.DictReader(handle)]
    totals = [0]*(N+1)
    for row in rows:
        totals[row["k"]] += row["count"]
    if totals != [math.comb(N, k) for k in range(N+1)]:
        raise ValueError("not the complete saved N25 population")
    return rows


def evaluate(rows, h):
    """Partial h derivatives, then substitute the formal root h(lambda)."""
    order = len(h)-1
    powers = [zero(order)]
    powers[0][0] = F(1)
    for _ in range(N):
        powers.append(mul(powers[-1], h))
    result = {key: zero(order) for key in ("z", "q", "r1", "zh", "qh", "r1h")}
    for row in rows:
        k, g, q, count = (row[key] for key in ("k", "g", "q", "count"))
        if g > order:
            continue
        for key, weight in (("z", count), ("q", q*count), ("r1", count*(q == 0))):
            if weight:
                for j, value in enumerate(powers[k][:order-g+1]):
                    if value:
                        result[key][g+j] += weight*value
        if k:
            for key, weight in (("zh", k*count), ("qh", k*q*count),
                                ("r1h", k*count*(q == 0))):
                if weight:
                    for j, value in enumerate(powers[k-1][:order-g+1]):
                        if value:
                            result[key][g+j] += weight*value
    return result


def pooled_numerator(pair):
    a, b = pair
    value = add(mul(a["q"], b["z"]), mul(b["q"], a["z"]))
    derivative = add(add(mul(a["qh"], b["z"]), mul(a["q"], b["zh"])),
                     add(mul(b["qh"], a["z"]), mul(b["q"], a["zh"])))
    return value, derivative


def formal_root(rows_pair, order):
    # At lambda=0, the numerator is 2(h^50-1), with derivative100 at h=1.
    h, precision, steps = [F(1)], 1, []
    while precision < order+1:
        target = min(2*precision, order+1)
        h += [F(0)]*(target-len(h))
        value, derivative = pooled_numerator([evaluate(rows, h) for rows in rows_pair])
        h = sub(h, divide(value, derivative))
        steps.append(target-1)
        precision = target
    value, _ = pooled_numerator([evaluate(rows, h) for rows in rows_pair])
    if any(value):
        raise ArithmeticError("root coefficients do not annihilate the pooled numerator")
    return h, steps


def observables(rows_pair, h):
    packets = [evaluate(rows, h) for rows in rows_pair]
    result = {}
    slopes = []
    p1h = []
    for name, packet in zip(("axis", "tilted"), packets):
        z = packet["z"]
        result[f"P1_{name}"] = divide(packet["r1"], z)
        derivative = divide(sub(mul(packet["r1h"], z),
                                mul(packet["r1"], packet["zh"])), mul(z, z))
        result[f"P1_h_{name}"] = derivative
        p1h.append(derivative)
        slopes.append(divide(sub(mul(packet["qh"], z),
                                 mul(packet["q"], packet["zh"])), mul(z, z)))
        result[f"raw_R1_h_{name}"] = packet["r1h"]
    result["D_Q_h"] = scale(add(*slopes), F(1, 2))
    result["Y_h"] = scale(sub(p1h[1], p1h[0]), 1/DELTA)
    result["U_over_A"] = divide(result["Y_h"], result["D_Q_h"])
    return result


def first(series, begin=0):
    for power, coefficient in enumerate(series[begin:], begin):
        if coefficient:
            return {"power": power, "coefficient": str(coefficient)}
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--order", type=int, default=32, choices=range(1, 33))
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    out = args.output_dir.resolve()
    if out.exists():
        raise FileExistsError("refusing to overwrite a prior symbolic result")
    paths = [ROOT/f"results/p337-closed-source-finite-coupling/{name}.csv"
             for name in ("axis", "tilted")]
    pair = [read_rows(path, HASHES[name]) for path, name in zip(paths, HASHES)]
    h, newton_orders = formal_root(pair, args.order)
    moving = observables(pair, h)
    fixed_h = zero(args.order)
    fixed_h[0] = F(1)
    fixed = observables(pair, fixed_h)
    series = {"h0": h, **moving,
              "fixed_h1_Y_h": fixed["Y_h"],
              "fixed_h1_U_over_A": fixed["U_over_A"],
              "root_motion_Y_h": sub(moving["Y_h"], fixed["Y_h"])}
    raw = {}
    for name, rows in zip(HASHES, pair):
        minimum = min(row["g"] for row in rows if row["q"] == 0)
        polynomial = {str(row["k"]): row["count"] for row in rows
                      if row["q"] == 0 and row["g"] == minimum}
        raw[name] = {"minimum_g": minimum, "rank1_h_polynomial": polynomial,
                     "configurations_at_minimum": sum(polynomial.values())}
    lead = first(moving["U_over_A"])
    result = {
        "schema": "matching-one.p337.closed-source-asymptotic.v1",
        "status": "leading_term_found" if lead else f"no_leading_term_through_{args.order}",
        "source_commit": SOURCE, "code_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "script_sha256": sha(Path(__file__)),
        "inputs": [{"path": str(path.relative_to(ROOT)), "sha256": sha(path)} for path in paths],
        "N": N, "delta_cos4": str(DELTA), "lambda": "exp(-t)",
        "branch": "positive pooled root h0(0)=1; h=exp(logit(p)-t)",
        "maximum_degree": args.order, "arithmetic": "fractions.Fraction, exact rational",
        "newton_truncation_degrees": newton_orders,
        "pooled_root_residual": f"identically zero through lambda^{args.order}",
        "raw_rank1_leading_shells": raw,
        "first_nonzero": {key: first(values, 1 if key == "h0" else 0)
                          for key, values in series.items()},
        "series": {key: [str(value) for value in values] for key, values in series.items()},
        "leading_Yh_decomposition": None if not lead else {
            "power": lead["power"],
            "fixed_h1": str(fixed["Y_h"][lead["power"]]),
            "root_motion": str(series["root_motion_Y_h"][lead["power"]]),
            "total": str(moving["Y_h"][lead["power"]]),
            "D_constant": str(moving["D_Q_h"][0])},
        "classification": "deterministic post-reveal symbolic derivation; not prospective or independent evidence",
        "dependency_group": "p337-N25-axis-tilted-exhaustive-configuration-populations",
        "new_enumerations": 0, "new_coupling_points": 0, "new_random_samples": 0,
        "runtime": {"seconds": time.perf_counter()-started,
                    "python": sys.executable, "version": platform.python_version()},
    }
    out.mkdir(parents=True)
    (out/"series.json").write_text(json.dumps(result, indent=2)+"\n")
    with (out/"coefficients.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["lambda_power", *series])
        for power in range(args.order+1):
            writer.writerow([power, *(str(values[power]) for values in series.values())])
    print(json.dumps({"status": result["status"], "leading": result["first_nonzero"],
                      "raw_rank1": raw, "Yh_decomposition": result["leading_Yh_decomposition"],
                      "seconds": result["runtime"]["seconds"]}, indent=2))


if __name__ == "__main__":
    main()
