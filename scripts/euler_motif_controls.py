#!/usr/bin/env python3
"""Euler/motif identities and control-variate observables for square tori.

The configuration-level matching identity verified on the implemented
quotients is

    C_black - C_white = q + V - E + F0

where ``q`` is the common wrapping-difference topology variable in
``{-1,0,+1}``, ``V`` is the occupied-site count, ``E`` is the occupied
nearest-neighbour edge count (``2N`` covering NN edges), and ``F0`` is the
occupied elementary-face count (``N`` square plaquettes).  Subtracting the
known Bernoulli means ``N chi(p)`` with ``chi(p)=p-2p^2+p^4`` yields the
equal-mean form in the P34 specification:

    D_cluster = q + (V - N p) - (E - 2 N p^2) + (F0 - N p^4).

Wrapping-only matching-odd channels are configuration-identical on these
quotients and must be rejected rather than GLS-combined.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

from integer_period_torus import (
    IntegerTorusGeometry,
    axis_integer_torus,
    classify_configuration,
    diamond_integer_torus,
    gaussian_integer_torus,
)
from matched_torus_reference import cluster_stats


WRAPPING_CHANNELS = ("cross", "both", "either", "direction_0", "direction_1")
EULER_CONTROLS = ("V", "E", "F0")
EXTRA_MOTIFS = ("nnn_pos", "nnn_neg", "path3_x", "path3_y", "corners")
ALL_MOTIFS = EULER_CONTROLS + EXTRA_MOTIFS


def falling(n: int, k: int) -> int:
    if k < 0:
        raise ValueError("falling factorial degree must be nonnegative")
    if n < k:
        return 0
    value = 1
    for i in range(k):
        value *= n - i
    return value


def chi(p: float) -> float:
    return p - 2.0 * p * p + p**4


def analytic_motif_mean(name: str, n: int, p: float) -> float:
    if name == "V":
        return n * p
    if name == "E":
        return 2.0 * n * p * p
    if name == "F0":
        return n * (p**4)
    if name in ("nnn_pos", "nnn_neg"):
        return n * p * p
    if name in ("path3_x", "path3_y"):
        return n * (p**3)
    if name == "corners":
        return 4.0 * n * (p**3)
    raise ValueError("unknown motif " + name)


def microcanonical_motif_mean(name: str, n: int, k: int) -> float:
    if name == "V":
        return float(k)
    if name == "E":
        denom = falling(n, 2)
        return 0.0 if denom == 0 else 2.0 * n * falling(k, 2) / denom
    if name == "F0":
        denom = falling(n, 4)
        return 0.0 if denom == 0 else n * falling(k, 4) / denom
    if name in ("nnn_pos", "nnn_neg"):
        denom = falling(n, 2)
        return 0.0 if denom == 0 else n * falling(k, 2) / denom
    if name in ("path3_x", "path3_y"):
        denom = falling(n, 3)
        return 0.0 if denom == 0 else n * falling(k, 3) / denom
    if name == "corners":
        denom = falling(n, 3)
        return 0.0 if denom == 0 else 4.0 * n * falling(k, 3) / denom
    raise ValueError("unknown motif " + name)


def _three_distinct(geometry: IntegerTorusGeometry, origin, step) -> bool:
    points = (
        origin,
        (origin[0] + step[0], origin[1] + step[1]),
        (origin[0] + 2 * step[0], origin[1] + 2 * step[1]),
    )
    vertices = [geometry.vertex(point) for point in points]
    return len(set(vertices)) == 3


def count_motifs(geometry: IntegerTorusGeometry, active: Sequence[bool]) -> dict[str, int]:
    if len(active) != geometry.n:
        raise ValueError("active mask length does not match geometry")
    n = geometry.n
    occupied = sum(1 for value in active if value)
    edges = sum(
        1 for edge in geometry.primal_edges if active[edge.i] and active[edge.j]
    )
    faces = 0
    nnn_pos = 0
    nnn_neg = 0
    path3_x = 0
    path3_y = 0
    corners = 0
    for x, y in geometry.coordinates:
        v00 = geometry.vertex((x, y))
        v10 = geometry.vertex((x + 1, y))
        v01 = geometry.vertex((x, y + 1))
        v11 = geometry.vertex((x + 1, y + 1))
        v1m = geometry.vertex((x + 1, y - 1))
        s00 = active[v00]
        s10 = active[v10]
        s01 = active[v01]
        s11 = active[v11]
        if s00 and s10 and s01 and s11:
            faces += 1
        if s00 and s11:
            nnn_pos += 1
        if s00 and active[v1m]:
            nnn_neg += 1
        if s00 and s10 and s01:
            corners += 1
        if s00 and s10 and s11:
            corners += 1
        if s00 and s01 and s11:
            corners += 1
        if s10 and s01 and s11:
            corners += 1
        if _three_distinct(geometry, (x, y), (1, 0)):
            v20 = geometry.vertex((x + 2, y))
            if s00 and s10 and active[v20]:
                path3_x += 1
        if _three_distinct(geometry, (x, y), (0, 1)):
            v02 = geometry.vertex((x, y + 2))
            if s00 and s01 and active[v02]:
                path3_y += 1
    if faces > n or edges > 2 * n:
        raise AssertionError("motif counts exceed covering-space maxima")
    return {
        "V": occupied,
        "E": edges,
        "F0": faces,
        "nnn_pos": nnn_pos,
        "nnn_neg": nnn_neg,
        "path3_x": path3_x,
        "path3_y": path3_y,
        "corners": corners,
    }


def wrapping_differences(
    geometry: IntegerTorusGeometry, active: Sequence[bool]
) -> dict[str, int]:
    white = [not value for value in active]
    black, _ = classify_configuration(geometry, active)
    matching, _ = classify_configuration(geometry, white, matching=True)
    return {
        channel: int(getattr(black, channel)) - int(getattr(matching, channel))
        for channel in WRAPPING_CHANNELS
    }


def cluster_difference(geometry: IntegerTorusGeometry, active: Sequence[bool]) -> tuple[int, int]:
    white = [not value for value in active]
    black_clusters, _ = cluster_stats(list(active), geometry.primal_edges)
    white_clusters, _ = cluster_stats(white, geometry.matching_edges)
    return black_clusters, white_clusters


@dataclass(frozen=True)
class ConfigurationIdentity:
    mask: int
    q: int
    cluster_difference: int
    motifs: dict[str, int]
    wrapping: dict[str, int]
    residual: int


def configuration_identity(
    geometry: IntegerTorusGeometry, active: Sequence[bool], mask: int = -1
) -> ConfigurationIdentity:
    wrapping = wrapping_differences(geometry, active)
    if len(set(wrapping.values())) != 1:
        raise AssertionError("wrapping-difference channels are not identical")
    q = wrapping["either"]
    black_clusters, white_clusters = cluster_difference(geometry, active)
    motifs = count_motifs(geometry, active)
    residual = (black_clusters - white_clusters) - (
        q + motifs["V"] - motifs["E"] + motifs["F0"]
    )
    return ConfigurationIdentity(
        mask, q, black_clusters - white_clusters, motifs, wrapping, residual
    )


def exhaustive_identity(geometry: IntegerTorusGeometry) -> dict[str, object]:
    if geometry.n > 16:
        raise ValueError("exhaustive identity is limited to N<=16")
    residuals = {}
    q_values = set()
    wrapping_not_identical = 0
    identity_failures = 0
    samples = 1 << geometry.n
    for mask in range(samples):
        active = [bool((mask >> vertex) & 1) for vertex in range(geometry.n)]
        wrapping = wrapping_differences(geometry, active)
        if len(set(wrapping.values())) != 1:
            wrapping_not_identical += 1
        record = configuration_identity(geometry, active, mask)
        residuals[record.residual] = residuals.get(record.residual, 0) + 1
        q_values.add(record.q)
        if record.residual != 0:
            identity_failures += 1
    return {
        "name": geometry.name,
        "N": geometry.n,
        "L": geometry.L,
        "physical_period": geometry.physical_period,
        "configurations": samples,
        "q_values": sorted(q_values),
        "residual_histogram": {str(key): value for key, value in sorted(residuals.items())},
        "identity_failures": identity_failures,
        "wrapping_not_identical": wrapping_not_identical,
        "passed": identity_failures == 0 and wrapping_not_identical == 0,
    }


def exhaustive_conditional_means(geometry: IntegerTorusGeometry) -> dict[str, object]:
    if geometry.n > 16:
        raise ValueError("exhaustive conditional means are limited to N<=16")
    n = geometry.n
    totals = {
        k: {name: 0 for name in ALL_MOTIFS} | {"count": 0} for k in range(n + 1)
    }
    path3_available_x = all(
        _three_distinct(geometry, point, (1, 0)) for point in geometry.coordinates
    )
    path3_available_y = all(
        _three_distinct(geometry, point, (0, 1)) for point in geometry.coordinates
    )
    for mask in range(1 << n):
        active = [bool((mask >> vertex) & 1) for vertex in range(n)]
        motifs = count_motifs(geometry, active)
        k = motifs["V"]
        totals[k]["count"] += 1
        for name in ALL_MOTIFS:
            totals[k][name] += motifs[name]

    failures = []
    by_k = []
    for k in range(n + 1):
        count = totals[k]["count"]
        if count == 0:
            continue
        expected_count = math.comb(n, k)
        if count != expected_count:
            failures.append({"K": k, "reason": "binomial_count", "got": count, "expected": expected_count})
        row = {"K": k, "configurations": count, "motifs": {}}
        for name in ALL_MOTIFS:
            if name == "path3_x" and not path3_available_x:
                continue
            if name == "path3_y" and not path3_available_y:
                continue
            observed = totals[k][name]
            expected = microcanonical_motif_mean(name, n, k) * count
            ok = math.isclose(observed, expected, rel_tol=0.0, abs_tol=1e-9)
            row["motifs"][name] = {
                "sum": observed,
                "mean": observed / count,
                "expected_mean": microcanonical_motif_mean(name, n, k),
                "zero_conditional_mean_of_centered": math.isclose(
                    observed / count - microcanonical_motif_mean(name, n, k),
                    0.0,
                    abs_tol=1e-9,
                ),
            }
            if not ok:
                failures.append(
                    {
                        "K": k,
                        "motif": name,
                        "sum": observed,
                        "expected_sum": expected,
                    }
                )
        by_k.append(row)
    return {
        "name": geometry.name,
        "N": n,
        "path3_x_three_distinct": path3_available_x,
        "path3_y_three_distinct": path3_available_y,
        "failures": failures,
        "passed": not failures,
        "by_K": by_k,
    }


def named_tiny_geometries() -> list[IntegerTorusGeometry]:
    return [
        axis_integer_torus(2),
        axis_integer_torus(3),
        gaussian_integer_torus(2, 1),
        diamond_integer_torus(2),
    ]


def run_identity_suite() -> dict[str, object]:
    identities = [exhaustive_identity(geometry) for geometry in named_tiny_geometries()]
    conditionals = [
        exhaustive_conditional_means(geometry) for geometry in named_tiny_geometries()
    ]
    passed = all(item["passed"] for item in identities + conditionals)
    return {
        "schema": "P34 Euler/motif exact identities v1",
        "identity": "C_black - C_white = q + V - E + F0",
        "equal_mean_form": "D_cluster = q + (V-Np) - (E-2Np^2) + (F0-Np^4)",
        "passed": passed,
        "exhaustive_identities": identities,
        "conditional_means": conditionals,
    }


def json_ready(value: object) -> object:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact-json", type=Path, required=True)
    args = parser.parse_args()
    payload = run_identity_suite()
    args.exact_json.parent.mkdir(parents=True, exist_ok=True)
    args.exact_json.write_text(json.dumps(json_ready(payload), indent=2) + "\n", encoding="utf-8")
    status = "PASS" if payload["passed"] else "FAIL"
    print("identity suite " + status)
    for row in payload["exhaustive_identities"]:
        print(
            "{name} N={N}: identity_failures={identity_failures} "
            "wrapping_not_identical={wrapping_not_identical}".format(**row)
        )
    for row in payload["conditional_means"]:
        print(
            "{name} N={N}: conditional_mean_failures={count}".format(
                count=len(row["failures"]), **row
            )
        )
    print("wrote " + str(args.exact_json))
    return 0 if payload["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
