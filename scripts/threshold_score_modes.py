#!/usr/bin/env python3
"""Project threshold-rank curves onto binomial Krawtchouk score modes.

For X~Bin(N,p), the orthonormal mode H_r(X) is chosen with H_1=(X-Np)/
sqrt(Np(1-p)).  A threshold histogram supplies the microcanonical response
q(x), so c_r=E[q(X) H_r(X)] is available without new simulation.

The script recomputes the intrinsic center inside every delete-one replicate
and reports angular-normalized matching S/D modes with jackknife covariance.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple

import mpmath as mp

from analyze_matching_parity_derivatives_fast import H, combine, cos4, obs, read, remove


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def krawtchouk_mode(n: int, x: int, order: int, p: mp.mpf) -> mp.mpf:
    """Positive-score convention for an orthonormal Bin(N,p) basis."""
    if not 0 < p < 1:
        raise ValueError("p must lie strictly between zero and one")
    if not 0 <= order <= n:
        raise ValueError("mode order must lie in 0..N")
    q = 1 - p
    coefficient = mp.mpf(0)
    lower = max(0, order - (n - x))
    upper = min(order, x)
    for j in range(lower, upper + 1):
        coefficient += (
            math.comb(x, j)
            * math.comb(n - x, order - j)
            * mp.power(-q / p, j)
        )
    polynomial = coefficient / math.comb(n, order)
    return mp.power(-1, order) * polynomial * mp.sqrt(
        math.comb(n, order) * mp.power(p / q, order)
    )


def binomial_weights(n: int, p: mp.mpf) -> List[mp.mpf]:
    q = 1 - p
    weights = [q**n]
    for x in range(n):
        weights.append(weights[-1] * (n - x) * p / ((x + 1) * q))
    total = mp.fsum(weights)
    return [weight / total for weight in weights]


def microcanonical_curves(row: H) -> Tuple[List[mp.mpf], List[mp.mpf]]:
    plus: List[mp.mpf] = []
    partner: List[mp.mpf] = []
    plus_sum = 0
    minus_sum = 0
    for x in range(row.n + 1):
        plus_sum += row.plus[x]
        minus_sum += row.minus[x]
        plus.append(mp.mpf(plus_sum) / row.samples)
        partner.append(1 - mp.mpf(minus_sum) / row.samples)
    return plus, partner


def coefficients(row: H, p: mp.mpf, max_order: int) -> Dict[str, List[mp.mpf]]:
    plus, partner = microcanonical_curves(row)
    weights = binomial_weights(row.n, p)
    result = {"plus": [], "partner": [], "S": [], "D": []}
    for order in range(max_order + 1):
        mode = [krawtchouk_mode(row.n, x, order, p) for x in range(row.n + 1)]
        plus_value = mp.fsum(weights[x] * plus[x] * mode[x] for x in range(row.n + 1))
        partner_value = mp.fsum(
            weights[x] * partner[x] * mode[x] for x in range(row.n + 1)
        )
        result["plus"].append(plus_value)
        result["partner"].append(partner_value)
        result["S"].append((plus_value + partner_value) / 2)
        result["D"].append((plus_value - partner_value) / 2)
    return result


def center(first: H, second: H, iterations: int = 120) -> mp.mpf:
    lo, hi = mp.mpf(0), mp.mpf(1)
    for _ in range(iterations):
        p = (lo + hi) / 2
        mean_matching = (obs(first, p)["M"] + obs(second, p)["M"]) / 2
        if mean_matching < 0:
            lo = p
        else:
            hi = p
    return (lo + hi) / 2


def project(first: H, second: H, max_order: int) -> dict:
    p = center(first, second)
    first_modes = coefficients(first, p, max_order)
    second_modes = coefficients(second, p, max_order)
    delta_cos4 = cos4(first.a, first.b) - cos4(second.a, second.b)
    if delta_cos4 == 0:
        raise ValueError("orientation pair has zero DeltaCos4")
    projected = {
        sector: [
            (first_modes[sector][order] - second_modes[sector][order]) / delta_cos4
            for order in range(max_order + 1)
        ]
        for sector in ("S", "D")
    }
    derivative_factor = mp.sqrt(first.n / (p * (1 - p)))
    return {
        "p0": p,
        "delta_cos4": delta_cos4,
        "P4_S_modes": projected["S"],
        "P4_D_modes": projected["D"],
        "P4_S_prime_from_mode1": derivative_factor * projected["S"][1],
        "P4_D_prime_from_mode1": derivative_factor * projected["D"][1],
    }


def jackknife_covariance(rows: Sequence[Sequence[mp.mpf]]) -> List[List[mp.mpf]]:
    count = len(rows)
    means = [mp.fsum(row[j] for row in rows) / count for j in range(len(rows[0]))]
    factor = mp.mpf(count - 1) / count
    return [[
        factor * mp.fsum(
            (row[i] - means[i]) * (row[j] - means[j]) for row in rows
        )
        for j in range(len(means))] for i in range(len(means))]


def flatten(result: Mapping[str, object], max_order: int) -> List[mp.mpf]:
    return list(result["P4_S_modes"][: max_order + 1]) + list(
        result["P4_D_modes"][: max_order + 1]
    )


def analyze_size(data: Mapping[tuple[int, str, int], H], n: int, max_order: int) -> dict:
    grouped = {
        orientation: [
            data[key] for key in sorted(data)
            if key[0] == n and key[1] == orientation
        ]
        for orientation in ("first", "second")
    }
    if len(grouped["first"]) != len(grouped["second"]) or len(grouped["first"]) < 2:
        raise ValueError(f"N={n}: orientations require aligned nontrivial batches")
    first = combine(grouped["first"])
    second = combine(grouped["second"])
    point = project(first, second, max_order)
    deleted = [
        project(remove(first, left), remove(second, right), max_order)
        for left, right in zip(grouped["first"], grouped["second"])
    ]
    covariance = jackknife_covariance([flatten(row, max_order) for row in deleted])
    labels = [f"P4_S_mode_{r}" for r in range(max_order + 1)] + [
        f"P4_D_mode_{r}" for r in range(max_order + 1)
    ]
    values = flatten(point, max_order)
    standard_errors = [mp.sqrt(max(mp.mpf(0), covariance[i][i])) for i in range(len(labels))]
    direct = {
        "first": obs(first, point["p0"]),
        "second": obs(second, point["p0"]),
    }
    direct_p4_s_prime = (direct["first"]["Sp"] - direct["second"]["Sp"]) / point["delta_cos4"]
    direct_p4_d_prime = (direct["first"]["Dp"] - direct["second"]["Dp"]) / point["delta_cos4"]
    tower_scaled: Dict[str, str] = {}
    for order in range(max_order + 1):
        if order % 2 == 0:
            i_sector, i_value, i_base = "S", point["P4_S_modes"][order], mp.mpf(1)
            t_sector, t_value, t_base = "D", point["P4_D_modes"][order], mp.mpf(13) / 8
        else:
            i_sector, i_value, i_base = "D", point["P4_D_modes"][order], mp.mpf(1)
            t_sector, t_value, t_base = "S", point["P4_S_modes"][order], mp.mpf(13) / 8
        i_exponent = i_base + mp.mpf(order) / 8
        t_exponent = t_base + mp.mpf(order) / 8
        tower_scaled[f"I_{i_sector}_mode_{order}_N_power_{mp.nstr(i_exponent, 8)}"] = mp.nstr(
            mp.power(n, i_exponent) * i_value, 30
        )
        tower_scaled[f"T_{t_sector}_mode_{order}_N_power_{mp.nstr(t_exponent, 8)}"] = mp.nstr(
            mp.power(n, t_exponent) * t_value, 30
        )
    return {
        "N": n,
        "representations": {"first": [first.a, first.b], "second": [second.a, second.b]},
        "samples": first.samples,
        "batches": len(grouped["first"]),
        "p0": mp.nstr(point["p0"], 30),
        "delta_cos4": mp.nstr(point["delta_cos4"], 30),
        "mode_order": labels,
        "point": {label: mp.nstr(value, 30) for label, value in zip(labels, values)},
        "standard_error": {
            label: mp.nstr(value, 15) for label, value in zip(labels, standard_errors)
        },
        "covariance": [[mp.nstr(value, 15) for value in row] for row in covariance],
        "exact_view_identities": {
            "P4_S_mode_0_minus_direct_P4_S": mp.nstr(
                point["P4_S_modes"][0]
                - (direct["first"]["S"] - direct["second"]["S"]) / point["delta_cos4"], 10
            ),
            "P4_D_mode_0_minus_direct_P4_D": mp.nstr(
                point["P4_D_modes"][0]
                - (direct["first"]["D"] - direct["second"]["D"]) / point["delta_cos4"], 10
            ),
            "P4_S_prime_from_mode1_minus_direct": mp.nstr(
                point["P4_S_prime_from_mode1"] - direct_p4_s_prime, 10
            ),
            "P4_D_prime_from_mode1_minus_direct": mp.nstr(
                point["P4_D_prime_from_mode1"] - direct_p4_d_prime, 10
            ),
        },
        "scaled_live_views": {
            "N_times_P4_S_mode_0": mp.nstr(n * point["P4_S_modes"][0], 30),
            "N13_8_times_P4_D_mode_0": mp.nstr(
                mp.power(n, mp.mpf(13) / 8) * point["P4_D_modes"][0], 30
            ),
            "N7_4_times_P4_S_mode_1": mp.nstr(
                mp.power(n, mp.mpf(7) / 4) * point["P4_S_modes"][1], 30
            ),
            "N9_8_times_P4_D_mode_1": mp.nstr(
                mp.power(n, mp.mpf(9) / 8) * point["P4_D_modes"][1], 30
            ),
        },
        "parity_tower_scaled": tower_scaled,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("histograms", nargs="+", type=Path)
    parser.add_argument("--max-order", type=int, default=6)
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if not 1 <= args.max_order <= 12:
        raise SystemExit("--max-order must lie in 1..12")
    mp.mp.dps = args.dps
    data: Dict[tuple[int, str, int], H] = {}
    for path in args.histograms:
        block = read(path)
        overlap = set(data).intersection(block)
        if overlap:
            raise SystemExit(f"duplicate histogram keys from {path}")
        data.update(block)
    payload = {
        "schema": "matching-one/threshold-krawtchouk-score-modes/v1",
        "coordinate": "eta=log(p/(1-p))",
        "basis": "orthonormal Bin(N,p0) Krawtchouk; H1=(K-Np)/sqrt(Np(1-p))",
        "max_order": args.max_order,
        "inputs": [
            {"path": str(path), "sha256": sha256(path)} for path in args.histograms
        ],
        "by_N": {
            str(n): analyze_size(data, n, args.max_order)
            for n in sorted({key[0] for key in data})
        },
        "evidence_guard": (
            "Mode 0 and mode 1 are exact coordinate views of the existing value and "
            "first derivative. They are not independent evidence blocks."
        ),
    }
    rendered = json.dumps(payload, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(args.output)
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
