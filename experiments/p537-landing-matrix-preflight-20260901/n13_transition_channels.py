#!/usr/bin/env python3
"""Exact N13 rank-transition preflight for the two covariance channels."""
from __future__ import annotations

import argparse
import importlib.util
import json
from collections import Counter
from fractions import Fraction as F
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / "experiments/p337-thermal-gate-audit-20260901/thermal_gate.py"


def exact(value: F) -> dict[str, str | float]:
    return {"fraction": str(value), "decimal": float(value)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=HERE / "n13-transition-result.json")
    args = parser.parse_args()

    spec = importlib.util.spec_from_file_location("thermal_gate", SOURCE)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)

    torus = module.Torus(3, 2)
    n = torus.n
    ranks = []
    sources = []
    for mask in range(1 << n):
        ranks.append(torus.rank(mask))
        sources.append(F(2 * sum(torus.pair_values(mask)), 16 * n * n))

    p = F(3, 5)
    full = [p**k * (1 - p) ** (n - k) for k in range(n + 1)]
    conditional = [p**k * (1 - p) ** (n - 1 - k) for k in range(n)]
    transition_counts: Counter[str] = Counter()
    source_change_counts: Counter[str] = Counter()
    for z in range(n):
        bit = 1 << z
        for mask in range(1 << n):
            if mask & bit:
                continue
            filled = mask | bit
            transition = f"{ranks[mask]}->{ranks[filled]}"
            transition_counts[transition] += 1
            source_change_counts[transition] += sources[mask] != sources[filled]

    observables = {}
    for name, values in (
        ("q", [rank - 1 for rank in ranks]),
        ("E", [(rank - 1) ** 2 for rank in ranks]),
    ):
        mean_source = sum(full[mask.bit_count()] * sources[mask] for mask in range(1 << n))
        mean_observable = sum(full[mask.bit_count()] * values[mask] for mask in range(1 << n))
        cells = {transition: [F(0), F(0)] for transition in ("0->1", "1->2", "0->2")}
        for z in range(n):
            bit = 1 << z
            for mask in range(1 << n):
                if mask & bit:
                    continue
                filled = mask | bit
                transition = f"{ranks[mask]}->{ranks[filled]}"
                if transition not in cells:
                    continue
                weight = conditional[mask.bit_count()]
                o0, o1 = F(values[mask]), F(values[filled])
                a0, a1 = sources[mask], sources[filled]
                cells[transition][0] += weight * ((o0 + o1) / 2 - mean_observable) * (a1 - a0)
                cells[transition][1] += weight * ((a0 + a1) / 2 - mean_source) * (o1 - o0)
        first, second = cells["0->1"], cells["1->2"]
        minor = first[0] * second[1] - second[0] * first[1]
        observables[name] = {
            "cells": {
                key: {"kernel": exact(value[0]), "readout": exact(value[1])}
                for key, value in cells.items()
            },
            "minor_rows_kernel_readout_columns_0to1_1to2": exact(minor),
        }

    payload = {
        "schema": "matching-one/p537-n13-transition-preflight/v1",
        "N": n,
        "geometry": [3, 2],
        "p": str(p),
        "source": "a=2*sum_(x<y) g16/(16*N^2)",
        "transition_counts": {
            key: {
                "records": transition_counts[key],
                "source_change_nonzero": source_change_counts[key],
            }
            for key in sorted(transition_counts)
        },
        "observables": observables,
        "scope_warning": "Transitions aggregate every landing type and have no C4 or root-Schur projection.",
    }
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
