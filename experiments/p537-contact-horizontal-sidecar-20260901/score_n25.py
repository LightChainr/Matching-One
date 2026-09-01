#!/usr/bin/env python3
"""Exact-population midpoint of the N25 horizontal contact representative."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


N = 25
P = 0.5926655393282267
DELTA = 1152 / 625
ORBITS = {
    "nn_other": {(-1, 0), (0, -1), (0, 1)},
    "diag1": {(-1, -1), (-1, 1), (1, -1), (1, 1)},
    "axial2": {(-2, 0), (0, -2), (2, 0), (0, 2)},
    "knight_a": {(-2, -1), (1, -2), (2, 1), (-1, 2)},
    "knight_b": {(-2, 1), (-1, -2), (2, -1), (1, 2)},
}
FAR = {
    "axis": {(-2, -2), (-2, 2), (2, -2), (2, 2)},
    "tilted": {(-3, 0), (0, -3), (0, 3), (3, 0)},
}


def component(geometry: str, dx: int, dy: int) -> str:
    point = (dx, dy)
    for name, points in ORBITS.items():
        if point in points:
            return name
    if point in FAR[geometry]:
        return "far"
    raise ValueError((geometry, point))


def baseline(path: Path):
    q = e = qt = et = 0.0
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            k = int(row["k"])
            value = int(row["q"])
            count = int(row["count"])
            weight = P ** k * (1 - P) ** (N - k)
            score = k - N * P
            q += count * value * weight
            e += count * value * value * weight
            qt += count * value * weight * score
            et += count * value * value * weight * score
    return {"q": q, "e": e, "qt": qt, "et": et}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--axis", required=True, type=Path)
    parser.add_argument("--tilted", required=True, type=Path)
    parser.add_argument("--aggregates", required=True, type=Path)
    parser.add_argument("--baseline-axis", required=True, type=Path)
    parser.add_argument("--baseline-tilted", required=True, type=Path)
    parser.add_argument("--kernel", required=True, type=Path)
    parser.add_argument("--full-result", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)

    base = {"axis": baseline(args.baseline_axis), "tilted": baseline(args.baseline_tilted)}
    mt = (base["axis"]["qt"] + base["tilted"]["qt"]) / 2
    yt = (base["axis"]["et"] - base["tilted"]["et"]) / DELTA
    ratio = yt / mt
    c = {"axis": 1 / DELTA, "tilted": -1 / DELTA}
    mu_h = {
        g: 2 * c[g] * base[g]["e"] - ratio * base[g]["q"]
        for g in ("axis", "tilted")
    }

    global_rows = []
    with args.aggregates.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["tau"] == "__GLOBAL__":
                global_rows.append(row)
    components = sorted({row["source_component"] for row in global_rows})
    mu_a = {}
    beta = {}
    for g in ("axis", "tilted"):
        for name in components:
            a = qa = 0.0
            for row in global_rows:
                if row["geometry"] != g or row["source_component"] != name:
                    continue
                k = int(row["k_minus"])
                weight = P ** k * (1 - P) ** (N - 1 - k)
                a += weight * (
                    (1 - P) * int(row["sum_a0"]) + P * int(row["sum_a1"])
                ) / (N * 16)
                qa += weight * (
                    (1 - P) * int(row["sum_q0a0"]) + P * int(row["sum_q1a1"])
                ) / (N * 16)
            mu_a[(g, name)] = a
            mu_a[(g, name, "cov")] = qa - base[g]["q"] * a
    for name in components:
        beta[name] = (
            mu_a[("axis", name, "cov")] + mu_a[("tilted", name, "cov")]
        ) / (2 * mt)

    kernel = defaultdict(int)
    with args.kernel.open(newline="") as handle:
        for row in csv.DictReader(
            (line for line in handle if not line.startswith("#")), delimiter="\t"
        ):
            kernel[int(row.get("key", row.get("packed_key")))] = int(row["g16"])

    by_component = defaultdict(lambda: [0.0, 0.0, 0.0])
    v = P * (1 - P)
    selected = 0
    for g, path in (("axis", args.axis), ("tilted", args.tilted)):
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                if not (
                    row["alternating_four_arm"] == "1"
                    and row["rank0"] != row["rank1"]
                    and row["source_absent"] == "0"
                    and row["bell0"] != row["bell1"]
                    and kernel[int(row["bell0"])] != kernel[int(row["bell1"])]
                ):
                    continue
                selected += 1
                name = component(g, int(row["y_dx"]), int(row["y_dy"]))
                k = int(row["k_minus"])
                weight = P ** k * (1 - P) ** (N - 1 - k)
                mu = mu_a[(g, name)]
                b = beta[name]
                s_minus = k - (N - 1) * P
                for state in (0, 1):
                    suffix = str(state)
                    w = (1 - P, P)[state]
                    u = state - P
                    score = s_minus + u
                    bterm = u * score - v
                    count = int(row["count"])
                    q = int(row[f"sum_q{suffix}"])
                    e = int(row[f"sum_E{suffix}"])
                    a = int(row[f"sum_a16_{suffix}"]) / (N * 16)
                    qa = int(row[f"sum_q{suffix}_a16_{suffix}"]) / (N * 16)
                    ea = int(row[f"sum_E{suffix}_a16_{suffix}"]) / (N * 16)
                    sh = 2 * c[g] * e - ratio * q - mu_h[g] * count
                    sha = (
                        2 * c[g] * (ea - mu * e)
                        - ratio * (qa - mu * q)
                        - mu_h[g] * (a - mu * count)
                    )
                    original = 2 * weight * w * (u * sha - b * bterm * sh)
                    cc = 2 * weight * w * v * sh
                    by_component[name][0] += original
                    by_component[name][1] += cc
                    by_component[name][2] += original - b * cc

    contact = [sum(values[i] for values in by_component.values()) for i in range(3)]
    full = json.loads(args.full_result.read_text())
    full_over_mt = float(full["controls"]["all_mode_reproduced_full_J2_over_A"]["midpoint"])
    t_full = full_over_mt * mt
    payload = {
        "schema": "matching-one/p537-contact-horizontal-n25/v1",
        "status": "COMPLETED_EXACT_POPULATION_MIDPOINT",
        "N25": {
            "p": P, "M_t": mt, "T_full": t_full,
            "T_contact_canonical": contact[0], "C_contact": contact[1],
            "beta_C_contact": contact[0] - contact[2],
            "T_contact_horizontal": contact[2],
            "T_remainder_canonical": t_full - contact[0],
            "T_remainder_horizontal": t_full - contact[2],
            "contact_identity_residual": contact[0] - contact[2]
            - sum(beta[name] * by_component[name][1] for name in components),
        },
        "component_order": components,
        "by_component": {
            name: {
                "beta": beta[name], "T_contact_canonical": by_component[name][0],
                "C_contact": by_component[name][1],
                "beta_C_contact": beta[name] * by_component[name][1],
                "T_contact_horizontal": by_component[name][2],
            }
            for name in components
        },
        "selected_exact_row_classes": selected,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    print(json.dumps(payload["N25"], indent=2))


if __name__ == "__main__":
    main()
