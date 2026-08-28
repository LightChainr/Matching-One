#!/usr/bin/env python3
"""Reconstruct orientation-resolved curves from C++ threshold-rank batches."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import mpmath as mp

from analyze_threshold_ranks import (
    matching_derivative,
    matching_root,
    matching_value,
)


Key = Tuple[int, str, int]


def mean_se(values: Sequence[float]) -> Tuple[float, float]:
    if len(values) < 2:
        raise ValueError("at least two batches are required for a batch SE")
    mean = math.fsum(values) / len(values)
    variance = math.fsum((value - mean) ** 2 for value in values)
    return mean, math.sqrt(variance / (len(values) * (len(values) - 1)))


def jackknife_se(values: Sequence[float]) -> float:
    if len(values) < 2:
        raise ValueError("at least two delete-one values are required")
    mean = math.fsum(values) / len(values)
    return math.sqrt(
        (len(values) - 1) / len(values)
        * math.fsum((value - mean) ** 2 for value in values)
    )


def cos4(a: int, b: int) -> float:
    n = a * a + b * b
    return (a**4 - 6 * a * a * b * b + b**4) / float(n * n)


def read_histograms(path: Path) -> Dict[Key, Dict[str, object]]:
    records: Dict[Key, Dict[str, object]] = {}
    required = {"n", "a", "b", "orientation", "batch", "samples", "kind", "k", "count"}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError("histogram CSV missing: " + ", ".join(sorted(missing)))
        for raw in reader:
            n = int(raw["n"])
            orientation = raw["orientation"]
            batch = int(raw["batch"])
            samples = int(raw["samples"])
            rank = int(raw["k"])
            count = int(raw["count"])
            kind = raw["kind"]
            if orientation not in ("first", "second") or kind not in ("minus", "plus"):
                raise ValueError("unknown orientation or histogram kind")
            if n <= 0 or batch < 0 or samples <= 0 or not 1 <= rank <= n or count <= 0:
                raise ValueError("invalid histogram row")
            key = (n, orientation, batch)
            if key not in records:
                records[key] = {
                    "n": n,
                    "a": int(raw["a"]),
                    "b": int(raw["b"]),
                    "orientation": orientation,
                    "batch": batch,
                    "samples": samples,
                    "minus": [0] * (n + 1),
                    "plus": [0] * (n + 1),
                }
            record = records[key]
            if (
                record["samples"] != samples
                or record["a"] != int(raw["a"])
                or record["b"] != int(raw["b"])
            ):
                raise ValueError("inconsistent rows within a histogram batch")
            values = record[kind]
            assert isinstance(values, list)
            values[rank] += count
    if not records:
        raise ValueError("histogram CSV is empty")
    for record in records.values():
        for kind in ("minus", "plus"):
            values = record[kind]
            assert isinstance(values, list)
            if sum(values) != record["samples"]:
                raise ValueError("histogram total differs from batch samples")
    return records


def validate_moments(path: Path, records: Dict[Key, Dict[str, object]]) -> None:
    seen = set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            key = (int(raw["n"]), raw["orientation"], int(raw["batch"]))
            if key not in records or key in seen:
                raise ValueError("moment row has missing or duplicate histogram batch")
            seen.add(key)
            record = records[key]
            minus = record["minus"]
            plus = record["plus"]
            assert isinstance(minus, list) and isinstance(plus, list)
            calculated = {
                "sum_kminus": sum(rank * count for rank, count in enumerate(minus)),
                "sum_kplus": sum(rank * count for rank, count in enumerate(plus)),
                "sum_kminus2": sum(rank * rank * count for rank, count in enumerate(minus)),
                "sum_kplus2": sum(rank * rank * count for rank, count in enumerate(plus)),
            }
            if any(int(raw[name]) != value for name, value in calculated.items()):
                raise ValueError("moment row disagrees with marginal histograms")
            samples = int(record["samples"])
            sum_minus = calculated["sum_kminus"]
            sum_plus = calculated["sum_kplus"]
            sum_product = int(raw["sum_product"])
            if sum_product < 0 or int(raw["sum_gap"]) != sum_plus - sum_minus:
                raise ValueError("invalid joint rank moments")
            if int(raw["samples"]) != samples:
                raise ValueError("moment sample count mismatch")
    if seen != set(records):
        raise ValueError("some histogram batches have no joint moments")


def add_histograms(records: Sequence[Dict[str, object]], kind: str) -> List[int]:
    n = int(records[0]["n"])
    total = [0] * (n + 1)
    for record in records:
        values = record[kind]
        assert isinstance(values, list)
        for rank, count in enumerate(values):
            total[rank] += count
    return total


def evaluate_histogram(record: Dict[str, object], p: mp.mpf) -> Tuple[mp.mpf, mp.mpf]:
    n = int(record["n"])
    samples = int(record["samples"])
    minus = record["minus"]
    plus = record["plus"]
    assert isinstance(minus, list) and isinstance(plus, list)
    return (
        matching_value(n, samples, minus, plus, p),
        matching_derivative(n, samples, minus, plus, p),
    )


def analyze(
    records: Dict[Key, Dict[str, object]], probabilities: Sequence[mp.mpf],
    closure_p: mp.mpf,
) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    by_n: Dict[str, object] = {}
    for n in sorted({key[0] for key in records}):
        orientations: Dict[str, List[Dict[str, object]]] = {}
        for orientation in ("first", "second"):
            selected = [
                records[key] for key in sorted(records)
                if key[0] == n and key[1] == orientation
            ]
            if len(selected) < 2:
                raise ValueError("each orientation requires at least two batches")
            if [int(row["batch"]) for row in selected] != list(range(len(selected))):
                raise ValueError("batch ids must be a complete zero-based range")
            orientations[orientation] = selected
        if len(orientations["first"]) != len(orientations["second"]):
            raise ValueError("orientation batch counts differ")
        sample_sizes = {
            int(record["samples"])
            for orientation in orientations.values()
            for record in orientation
        }
        if len(sample_sizes) != 1:
            raise ValueError("all batches must have equal sample size")

        aggregate = {}
        for orientation, selected in orientations.items():
            minus = add_histograms(selected, "minus")
            plus = add_histograms(selected, "plus")
            samples = sum(int(row["samples"]) for row in selected)
            aggregate[orientation] = {
                "a": selected[0]["a"],
                "b": selected[0]["b"],
                "samples": samples,
                "minus": minus,
                "plus": plus,
            }

        first = aggregate["first"]
        second = aggregate["second"]
        first_root = matching_root(n, first["samples"], first["minus"], first["plus"])
        second_root = matching_root(n, second["samples"], second["minus"], second["plus"])
        batch_root_differences = []
        for first_batch, second_batch in zip(
            orientations["first"], orientations["second"]
        ):
            first_batch_root = matching_root(
                n, int(first_batch["samples"]), first_batch["minus"], first_batch["plus"]
            )
            second_batch_root = matching_root(
                n, int(second_batch["samples"]), second_batch["minus"], second_batch["plus"]
            )
            batch_root_differences.append(float(first_batch_root - second_batch_root))
        _batch_root_mean, root_se = mean_se(batch_root_differences)
        rows.append({
            "N": n,
            "metric": "root",
            "p": "",
            "first": float(first_root),
            "second": float(second_root),
            "difference_first_minus_second": float(first_root - second_root),
            "difference_batch_se": root_se,
        })

        first_closure_m = matching_value(
            n, first["samples"], first["minus"], first["plus"], closure_p
        )
        second_closure_m = matching_value(
            n, second["samples"], second["minus"], second["plus"], closure_p
        )
        first_closure_d = matching_derivative(
            n, first["samples"], first["minus"], first["plus"], closure_p
        )
        second_closure_d = matching_derivative(
            n, second["samples"], second["minus"], second["plus"], closure_p
        )
        delta_m_closure = first_closure_m - second_closure_m
        mean_slope = (first_closure_d + second_closure_d) / 2
        root_gap = first_root - second_root
        closure = -root_gap * mean_slope / delta_m_closure

        delete_one = {
            "closure": [], "delta_M": [], "mean_slope": [], "root_gap": []
        }
        for first_batch, second_batch in zip(
            orientations["first"], orientations["second"]
        ):
            reduced = {}
            for name, aggregate_row, batch_row in (
                ("first", first, first_batch), ("second", second, second_batch)
            ):
                samples = int(aggregate_row["samples"]) - int(batch_row["samples"])
                reduced[name] = {
                    "samples": samples,
                    "minus": [
                        total - removed
                        for total, removed in zip(
                            aggregate_row["minus"], batch_row["minus"]
                        )
                    ],
                    "plus": [
                        total - removed
                        for total, removed in zip(
                            aggregate_row["plus"], batch_row["plus"]
                        )
                    ],
                }
            fm = matching_value(
                n, reduced["first"]["samples"], reduced["first"]["minus"],
                reduced["first"]["plus"], closure_p,
            )
            sm = matching_value(
                n, reduced["second"]["samples"], reduced["second"]["minus"],
                reduced["second"]["plus"], closure_p,
            )
            fd = matching_derivative(
                n, reduced["first"]["samples"], reduced["first"]["minus"],
                reduced["first"]["plus"], closure_p,
            )
            sd = matching_derivative(
                n, reduced["second"]["samples"], reduced["second"]["minus"],
                reduced["second"]["plus"], closure_p,
            )
            fr = matching_root(
                n, reduced["first"]["samples"], reduced["first"]["minus"],
                reduced["first"]["plus"],
            )
            sr = matching_root(
                n, reduced["second"]["samples"], reduced["second"]["minus"],
                reduced["second"]["plus"],
            )
            dm = fm - sm
            slope = (fd + sd) / 2
            gap = fr - sr
            delete_one["delta_M"].append(float(dm))
            delete_one["mean_slope"].append(float(slope))
            delete_one["root_gap"].append(float(gap))
            delete_one["closure"].append(float(-gap * slope / dm))

        delta_cos4 = cos4(int(first["a"]), int(first["b"])) - cos4(
            int(second["a"]), int(second["b"])
        )
        closure_payload = {
            "p": mp.nstr(closure_p, mp.mp.dps),
            "delta_M": mp.nstr(delta_m_closure, mp.mp.dps),
            "delta_M_jackknife_se": jackknife_se(delete_one["delta_M"]),
            "mean_M_prime": mp.nstr(mean_slope, mp.mp.dps),
            "mean_M_prime_jackknife_se": jackknife_se(delete_one["mean_slope"]),
            "delta_M_prime": mp.nstr(first_closure_d - second_closure_d, mp.mp.dps),
            "root_gap": mp.nstr(root_gap, mp.mp.dps),
            "root_gap_jackknife_se": jackknife_se(delete_one["root_gap"]),
            "linearized_root_gap": mp.nstr(-delta_m_closure / mean_slope, mp.mp.dps),
            "C": mp.nstr(closure, mp.mp.dps),
            "C_jackknife_se": jackknife_se(delete_one["closure"]),
            "A_M": mp.nstr(n ** (mp.mpf(13) / 8) * delta_m_closure / delta_cos4, mp.mp.dps),
            "B": mp.nstr(n ** (-mp.mpf(3) / 8) * mean_slope, mp.mp.dps),
            "A_p": mp.nstr(-n * n * root_gap / delta_cos4, mp.mp.dps),
            "delta_cos4": delta_cos4,
            "delete_one_batches": len(delete_one["closure"]),
        }

        evaluations = []
        for p in probabilities:
            first_m = matching_value(n, first["samples"], first["minus"], first["plus"], p)
            second_m = matching_value(n, second["samples"], second["minus"], second["plus"], p)
            first_d = matching_derivative(n, first["samples"], first["minus"], first["plus"], p)
            second_d = matching_derivative(n, second["samples"], second["minus"], second["plus"], p)
            batch_m_differences = []
            batch_d_differences = []
            for first_batch, second_batch in zip(
                orientations["first"], orientations["second"]
            ):
                fm, fd = evaluate_histogram(first_batch, p)
                sm, sd = evaluate_histogram(second_batch, p)
                batch_m_differences.append(float(fm - sm))
                batch_d_differences.append(float(fd - sd))
            _mean_m, se_m = mean_se(batch_m_differences)
            _mean_d, se_d = mean_se(batch_d_differences)
            p_text = mp.nstr(p, mp.mp.dps)
            evaluations.append({
                "p": p_text,
                "first_M": mp.nstr(first_m, mp.mp.dps),
                "second_M": mp.nstr(second_m, mp.mp.dps),
                "delta_M_first_minus_second": mp.nstr(first_m - second_m, mp.mp.dps),
                "delta_M_batch_se": se_m,
                "first_M_prime": mp.nstr(first_d, mp.mp.dps),
                "second_M_prime": mp.nstr(second_d, mp.mp.dps),
                "delta_M_prime_first_minus_second": mp.nstr(first_d - second_d, mp.mp.dps),
                "delta_M_prime_batch_se": se_d,
            })
            for metric, first_value, second_value, standard_error in (
                ("M", first_m, second_m, se_m),
                ("M_prime", first_d, second_d, se_d),
            ):
                rows.append({
                    "N": n,
                    "metric": metric,
                    "p": p_text,
                    "first": float(first_value),
                    "second": float(second_value),
                    "difference_first_minus_second": float(first_value - second_value),
                    "difference_batch_se": standard_error,
                })
        by_n[str(n)] = {
            "first_rep": [first["a"], first["b"]],
            "second_rep": [second["a"], second["b"]],
            "sample_count_per_orientation": first["samples"],
            "root_first": mp.nstr(first_root, mp.mp.dps),
            "root_second": mp.nstr(second_root, mp.mp.dps),
            "root_gap_first_minus_second": mp.nstr(first_root - second_root, mp.mp.dps),
            "root_gap_batch_se": root_se,
            "amplitude_closure": closure_payload,
            "evaluations": evaluations,
        }
    return rows, {"format_version": 1, "by_N": by_n}


def write_csv(path: Path, rows: Sequence[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "N", "metric", "p", "first", "second",
        "difference_first_minus_second", "difference_batch_se",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--histograms", type=Path, required=True)
    parser.add_argument("--moments", type=Path, required=True)
    parser.add_argument("--p", action="append", default=[])
    parser.add_argument(
        "--closure-p",
        help="coordinate for the same-batch amplitude-closure statistic",
    )
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--csv", type=Path, required=True)
    args = parser.parse_args()
    if args.dps < 30:
        raise SystemExit("--dps must be at least 30")
    mp.mp.dps = args.dps
    probabilities = [mp.mpf(value) for value in args.p]
    if not probabilities:
        probabilities = [mp.mpf("0.592746050790")]
    if any(not 0 < p < 1 for p in probabilities):
        raise SystemExit("each --p must lie strictly between zero and one")

    closure_p = (
        mp.mpf(args.closure_p) if args.closure_p is not None
        else probabilities[len(probabilities) // 2]
    )
    if not 0 < closure_p < 1:
        raise SystemExit("--closure-p must lie strictly between zero and one")
    records = read_histograms(args.histograms)
    validate_moments(args.moments, records)
    rows, payload = analyze(records, probabilities, closure_p)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    write_csv(args.csv, rows)
    print(f"wrote {args.json}\nwrote {args.csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
