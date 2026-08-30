#!/usr/bin/env python3
"""Reconstruct matching-even/odd spin-4 derivative sectors from threshold ranks.

Input is the long-form batch histogram emitted by the C++ same-N
threshold-rank engine.  The key finite-permutation identities are

    R_G(p)              = E[ I_p(K_plus) ],
    R_hat(1-p)          = 1 - E[ I_p(K_minus) ],

where I_p(k)=P(Bin(N,p)>=k).  Hence

    dR_G/dp             = E[ beta_pdf(K_plus) ],
    dR_hat(1-p)/dp      = -E[ beta_pdf(K_minus) ].

For each same-N orientation pair this script finds the intrinsic center p0
where the direction-average matching function vanishes, then forms

    S=(R_G+R_hat)/2,  D=(R_G-R_hat)/2,
    P4[X]=(X(theta1)-X(theta2))/Delta cos(4 theta).

It reports P4[S], P4[D], P4[S'], P4[D'] and a delete-one-batch jackknife
covariance matrix.  The frozen two-field model predicts powers

    P4[S]   ~ N^-1,
    P4[D]   ~ N^-13/8,
    P4[D']  ~ N^-5/8,
    P4[S']  ~ N^-5/4.

This script is intended for retrospective development on P33 and prospective
scoring on independently frozen geometries/statistics.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

import mpmath as mp


Key = Tuple[int, str, int]
METRICS = ("P4_S", "P4_D", "P4_S_prime", "P4_D_prime")


@dataclass
class BatchHistogram:
    n: int
    a: int
    b: int
    orientation: str
    batch: int
    samples: int
    minus: List[int]
    plus: List[int]


def cos4(a: int, b: int) -> mp.mpf:
    n = a * a + b * b
    return mp.mpf(a**4 - 6 * a * a * b * b + b**4) / (n * n)


def read_histograms(path: Path) -> Dict[Key, BatchHistogram]:
    required = {
        "n", "a", "b", "orientation", "batch", "samples", "kind", "k", "count"
    }
    records: Dict[Key, BatchHistogram] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError("histogram CSV missing: " + ", ".join(sorted(missing)))
        for raw in reader:
            n = int(raw["n"])
            a = int(raw["a"])
            b = int(raw["b"])
            orientation = raw["orientation"]
            batch = int(raw["batch"])
            samples = int(raw["samples"])
            kind = raw["kind"]
            rank = int(raw["k"])
            count = int(raw["count"])
            if orientation not in ("first", "second"):
                raise ValueError(f"unknown orientation {orientation!r}")
            if kind not in ("minus", "plus"):
                raise ValueError(f"unknown histogram kind {kind!r}")
            if n <= 0 or batch < 0 or samples <= 0 or not 1 <= rank <= n or count <= 0:
                raise ValueError("invalid histogram row")
            key = (n, orientation, batch)
            row = records.get(key)
            if row is None:
                row = BatchHistogram(
                    n=n, a=a, b=b, orientation=orientation, batch=batch,
                    samples=samples, minus=[0] * (n + 1), plus=[0] * (n + 1),
                )
                records[key] = row
            elif (row.a, row.b, row.samples) != (a, b, samples):
                raise ValueError("inconsistent metadata within a histogram batch")
            getattr(row, kind)[rank] += count

    if not records:
        raise ValueError("histogram CSV is empty")
    for row in records.values():
        if sum(row.minus) != row.samples or sum(row.plus) != row.samples:
            raise ValueError("histogram total differs from samples")
    return records


def add_histograms(rows: Sequence[BatchHistogram]) -> BatchHistogram:
    if not rows:
        raise ValueError("cannot aggregate an empty batch set")
    first = rows[0]
    total = BatchHistogram(
        n=first.n,
        a=first.a,
        b=first.b,
        orientation=first.orientation,
        batch=-1,
        samples=0,
        minus=[0] * (first.n + 1),
        plus=[0] * (first.n + 1),
    )
    for row in rows:
        if (row.n, row.a, row.b, row.orientation) != (
            first.n, first.a, first.b, first.orientation
        ):
            raise ValueError("attempted to aggregate incompatible histogram batches")
        total.samples += row.samples
        for rank in range(1, first.n + 1):
            total.minus[rank] += row.minus[rank]
            total.plus[rank] += row.plus[rank]
    return total


def subtract_batch(total: BatchHistogram, row: BatchHistogram) -> BatchHistogram:
    if (total.n, total.a, total.b, total.orientation) != (
        row.n, row.a, row.b, row.orientation
    ):
        raise ValueError("cannot subtract incompatible batch")
    result = BatchHistogram(
        n=total.n,
        a=total.a,
        b=total.b,
        orientation=total.orientation,
        batch=-1,
        samples=total.samples - row.samples,
        minus=[0] * (total.n + 1),
        plus=[0] * (total.n + 1),
    )
    if result.samples <= 0:
        raise ValueError("delete-one sample is empty")
    for rank in range(1, total.n + 1):
        result.minus[rank] = total.minus[rank] - row.minus[rank]
        result.plus[rank] = total.plus[rank] - row.plus[rank]
    return result


def threshold_cdf(n: int, rank: int, p: mp.mpf) -> mp.mpf:
    if p <= 0:
        return mp.mpf(0)
    if p >= 1:
        return mp.mpf(1)
    return mp.betainc(rank, n + 1 - rank, 0, p, regularized=True)


def threshold_pdf(n: int, rank: int, p: mp.mpf) -> mp.mpf:
    if not 0 < p < 1:
        return mp.mpf(0)
    return p ** (rank - 1) * (1 - p) ** (n - rank) / mp.beta(
        rank, n + 1 - rank
    )


def orientation_observables(row: BatchHistogram, p: mp.mpf) -> Dict[str, mp.mpf]:
    plus_cdf = mp.mpf(0)
    minus_cdf = mp.mpf(0)
    plus_pdf = mp.mpf(0)
    minus_pdf = mp.mpf(0)
    for rank in range(1, row.n + 1):
        if row.plus[rank] or row.minus[rank]:
            cdf = threshold_cdf(row.n, rank, p)
            pdf = threshold_pdf(row.n, rank, p)
            plus_cdf += row.plus[rank] * cdf
            minus_cdf += row.minus[rank] * cdf
            plus_pdf += row.plus[rank] * pdf
            minus_pdf += row.minus[rank] * pdf
    scale = mp.mpf(row.samples)
    r_g = plus_cdf / scale
    # White/matching cross survives until the black occupation count reaches K_minus.
    r_hat = 1 - minus_cdf / scale
    r_g_prime = plus_pdf / scale
    r_hat_prime = -minus_pdf / scale
    s = (r_g + r_hat) / 2
    d = (r_g - r_hat) / 2
    s_prime = (r_g_prime + r_hat_prime) / 2
    d_prime = (r_g_prime - r_hat_prime) / 2

    # Internal identity checks against the canonical matching representation.
    matching = r_g - r_hat
    matching_prime = r_g_prime - r_hat_prime
    if not mp.almosteq(2 * d, matching):
        raise ArithmeticError("2D != M")
    if not mp.almosteq(2 * d_prime, matching_prime):
        raise ArithmeticError("2D' != M'")

    return {
        "R_G": r_g,
        "R_hat": r_hat,
        "R_G_prime": r_g_prime,
        "R_hat_prime": r_hat_prime,
        "S": s,
        "D": d,
        "S_prime": s_prime,
        "D_prime": d_prime,
        "M": matching,
        "M_prime": matching_prime,
    }


def intrinsic_center(first: BatchHistogram, second: BatchHistogram, iterations: int = 100) -> mp.mpf:
    def mean_matching(p: mp.mpf) -> mp.mpf:
        return (
            orientation_observables(first, p)["M"]
            + orientation_observables(second, p)["M"]
        ) / 2

    lower = mp.mpf(0)
    upper = mp.mpf(1)
    for _ in range(iterations):
        midpoint = (lower + upper) / 2
        if mean_matching(midpoint) < 0:
            lower = midpoint
        else:
            upper = midpoint
    return (lower + upper) / 2


def projected_metrics(first: BatchHistogram, second: BatchHistogram) -> Tuple[mp.mpf, Dict[str, mp.mpf]]:
    if first.n != second.n:
        raise ValueError("orientation pair must have the same N")
    p0 = intrinsic_center(first, second)
    one = orientation_observables(first, p0)
    two = orientation_observables(second, p0)
    delta_cos4 = cos4(first.a, first.b) - cos4(second.a, second.b)
    if delta_cos4 == 0:
        raise ValueError("orientation pair has zero Delta cos4")
    result = {
        "P4_S": (one["S"] - two["S"]) / delta_cos4,
        "P4_D": (one["D"] - two["D"]) / delta_cos4,
        "P4_S_prime": (one["S_prime"] - two["S_prime"]) / delta_cos4,
        "P4_D_prime": (one["D_prime"] - two["D_prime"]) / delta_cos4,
        "delta_cos4": delta_cos4,
    }
    # At the intrinsic average center, the orientation-average D is exactly zero.
    if not mp.almosteq((one["D"] + two["D"]) / 2, 0, rel_eps=mp.mpf("1e-20"), abs_eps=mp.mpf("1e-20")):
        raise ArithmeticError("intrinsic-center condition failed")
    return p0, result


def jackknife_covariance(values: Sequence[Mapping[str, mp.mpf]]) -> List[List[mp.mpf]]:
    if len(values) < 2:
        raise ValueError("at least two delete-one values are required")
    means = {
        name: mp.fsum(row[name] for row in values) / len(values)
        for name in METRICS
    }
    factor = mp.mpf(len(values) - 1) / len(values)
    matrix: List[List[mp.mpf]] = []
    for left in METRICS:
        line = []
        for right in METRICS:
            line.append(
                factor
                * mp.fsum(
                    (row[left] - means[left]) * (row[right] - means[right])
                    for row in values
                )
            )
        matrix.append(line)
    return matrix


def analyze(records: Mapping[Key, BatchHistogram]) -> Dict[str, object]:
    output: Dict[str, object] = {"format_version": 1, "metrics": list(METRICS), "by_N": {}}
    by_n = output["by_N"]
    assert isinstance(by_n, dict)

    for n in sorted({key[0] for key in records}):
        groups: Dict[str, List[BatchHistogram]] = {}
        for orientation in ("first", "second"):
            rows = [records[key] for key in sorted(records) if key[0] == n and key[1] == orientation]
            if len(rows) < 2:
                raise ValueError(f"N={n}: need at least two batches per orientation")
            batch_ids = [row.batch for row in rows]
            if batch_ids != list(range(len(rows))):
                raise ValueError(f"N={n}: batch ids are not complete zero-based sequence")
            groups[orientation] = rows
        if len(groups["first"]) != len(groups["second"]):
            raise ValueError(f"N={n}: orientation batch counts differ")

        first_total = add_histograms(groups["first"])
        second_total = add_histograms(groups["second"])
        p0, point = projected_metrics(first_total, second_total)

        delete_one: List[Dict[str, mp.mpf]] = []
        delete_p0: List[mp.mpf] = []
        for first_batch, second_batch in zip(groups["first"], groups["second"]):
            p_i, metrics_i = projected_metrics(
                subtract_batch(first_total, first_batch),
                subtract_batch(second_total, second_batch),
            )
            delete_p0.append(p_i)
            delete_one.append({name: metrics_i[name] for name in METRICS})

        covariance = jackknife_covariance(delete_one)
        standard_errors = {
            name: mp.sqrt(covariance[index][index])
            for index, name in enumerate(METRICS)
        }
        p0_mean = mp.fsum(delete_p0) / len(delete_p0)
        p0_se = mp.sqrt(
            mp.mpf(len(delete_p0) - 1) / len(delete_p0)
            * mp.fsum((value - p0_mean) ** 2 for value in delete_p0)
        )

        scales = {
            "A_S_N1": n * point["P4_S"],
            "A_D_N13_8": n ** (mp.mpf(13) / 8) * point["P4_D"],
            "A_Dprime_N5_8": n ** (mp.mpf(5) / 8) * point["P4_D_prime"],
            "A_Sprime_N5_4": n ** (mp.mpf(5) / 4) * point["P4_S_prime"],
        }

        by_n[str(n)] = {
            "first_rep": [first_total.a, first_total.b],
            "second_rep": [second_total.a, second_total.b],
            "sample_count_per_orientation": first_total.samples,
            "batch_count": len(groups["first"]),
            "p0": mp.nstr(p0, mp.mp.dps),
            "p0_jackknife_se": mp.nstr(p0_se, 20),
            "delta_cos4": mp.nstr(point["delta_cos4"], 30),
            "point": {name: mp.nstr(point[name], 30) for name in METRICS},
            "standard_error": {name: mp.nstr(standard_errors[name], 20) for name in METRICS},
            "jackknife_covariance": [
                [mp.nstr(value, 20) for value in row] for row in covariance
            ],
            "scaled_diagnostics": {name: mp.nstr(value, 30) for name, value in scales.items()},
        }
    return output


def write_csv(path: Path, payload: Mapping[str, object]) -> None:
    rows = []
    by_n = payload["by_N"]
    assert isinstance(by_n, dict)
    for n_text, raw in by_n.items():
        assert isinstance(raw, dict)
        row = {
            "N": n_text,
            "first_rep": ",".join(map(str, raw["first_rep"])),
            "second_rep": ",".join(map(str, raw["second_rep"])),
            "samples": raw["sample_count_per_orientation"],
            "p0": raw["p0"],
            "p0_se": raw["p0_jackknife_se"],
            "delta_cos4": raw["delta_cos4"],
        }
        point = raw["point"]
        se = raw["standard_error"]
        scaled = raw["scaled_diagnostics"]
        assert isinstance(point, dict) and isinstance(se, dict) and isinstance(scaled, dict)
        for name in METRICS:
            row[name] = point[name]
            row[name + "_se"] = se[name]
        row.update(scaled)
        rows.append(row)
    fields = list(rows[0]) if rows else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--histograms", type=Path, required=True)
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--csv", type=Path, required=True)
    args = parser.parse_args()
    if args.dps < 30:
        raise SystemExit("--dps must be at least 30")
    mp.mp.dps = args.dps
    try:
        payload = analyze(read_histograms(args.histograms))
    except (ValueError, ArithmeticError) as exc:
        raise SystemExit(str(exc)) from exc
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    write_csv(args.csv, payload)
    print(args.json)
    print(args.csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
