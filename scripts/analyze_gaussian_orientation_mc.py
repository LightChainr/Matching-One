#!/usr/bin/env python3
"""Summarize batch output from ``src/gaussian_orientation_mc.cpp``.

The four sectors are formed after pairing, so every reported orientation
difference retains the common-random-number covariance.  ``even`` is
``(R_G+R_hat)/2`` and ``odd`` is ``(R_G-R_hat)/2``; ``matching_function`` is
twice ``odd`` and is included to avoid hidden factor-of-two conventions.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, List, Tuple

from control_variate_estimator import minimum_variance_weights, sample_covariance


BatchKey = Tuple[int, str, int]
CHANNELS = ("cross", "both", "either", "direction_0", "direction_1")
SECTORS = ("primal", "matching", "even", "odd", "matching_function")


def mean_se(values: List[float]) -> Tuple[float, float]:
    if len(values) < 2:
        raise ValueError("at least two batches are required")
    mean = math.fsum(values) / len(values)
    variance = math.fsum((value - mean) ** 2 for value in values)
    return mean, math.sqrt(variance / (len(values) * (len(values) - 1)))


def cos4(a: int, b: int) -> float:
    a2 = a * a
    b2 = b * b
    return (a2 * a2 - 6 * a2 * b2 + b2 * b2) / float((a2 + b2) ** 2)


def read_rows(path: Path) -> Tuple[Dict[BatchKey, Dict[str, int]], Dict[int, Tuple[int, int, int, int]]]:
    rows: Dict[BatchKey, Dict[str, int]] = {}
    designs: Dict[int, Tuple[int, int, int, int]] = {}
    required = {
        "n", "batch", "samples", "p_ref", "channel", "a1", "b1", "a2", "b2",
        "first_primal_sum", "first_matching_sum", "second_primal_sum", "second_matching_sum",
    }
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        if not required.issubset(fields):
            raise ValueError("missing CSV fields: " + ", ".join(sorted(required - fields)))
        for raw in reader:
            n = int(raw["n"])
            batch = int(raw["batch"])
            channel = raw["channel"]
            if channel not in CHANNELS:
                raise ValueError("unknown channel: " + channel)
            key = (n, channel, batch)
            if key in rows:
                raise ValueError("duplicate N/channel/batch row: " + repr(key))
            design = (int(raw["a1"]), int(raw["b1"]), int(raw["a2"]), int(raw["b2"]))
            if n in designs and designs[n] != design:
                raise ValueError("inconsistent representations for N=" + str(n))
            if design[0] ** 2 + design[1] ** 2 != n or design[2] ** 2 + design[3] ** 2 != n:
                raise ValueError("representation does not have advertised N")
            designs[n] = design
            samples = int(raw["samples"])
            if samples <= 0:
                raise ValueError("batch sample count must be positive")
            values = {
                "samples": samples,
                "first_primal": int(raw["first_primal_sum"]),
                "first_matching": int(raw["first_matching_sum"]),
                "second_primal": int(raw["second_primal_sum"]),
                "second_matching": int(raw["second_matching_sum"]),
            }
            if any(values[name] < 0 or values[name] > samples for name in values if name != "samples"):
                raise ValueError("indicator sum lies outside [0,samples]")
            rows[key] = values
    if not rows:
        raise ValueError("batch CSV is empty")
    for n in designs:
        expected_batches = None
        for channel in CHANNELS:
            batches = sorted(key[2] for key in rows if key[:2] == (n, channel))
            if len(batches) < 2 or batches != list(range(len(batches))):
                raise ValueError("batches must be a complete zero-based range")
            if expected_batches is None:
                expected_batches = batches
            elif batches != expected_batches:
                raise ValueError("channels do not contain identical batch ids")
        sample_sizes = {
            rows[(n, channel, batch)]["samples"]
            for channel in CHANNELS
            for batch in expected_batches or []
        }
        if len(sample_sizes) != 1:
            raise ValueError("GLS/covariance analysis requires equal batch sizes")
    return rows, designs


def sector_values(row: Dict[str, int], prefix: str) -> Dict[str, float]:
    samples = row["samples"]
    primal = row[prefix + "_primal"] / samples
    matching = row[prefix + "_matching"] / samples
    return {
        "primal": primal,
        "matching": matching,
        "even": 0.5 * (primal + matching),
        "odd": 0.5 * (primal - matching),
        "matching_function": primal - matching,
    }


def analyze(rows: Dict[BatchKey, Dict[str, int]], designs: Dict[int, Tuple[int, int, int, int]]) -> List[Dict[str, object]]:
    output: List[Dict[str, object]] = []
    for n in sorted(designs):
        a1, b1, a2, b2 = designs[n]
        delta_cos4 = cos4(a1, b1) - cos4(a2, b2)
        if delta_cos4 == 0:
            raise ValueError("orientation pair has zero delta cos(4 theta)")
        for channel in CHANNELS:
            keys = sorted((key for key in rows if key[:2] == (n, channel)), key=lambda key: key[2])
            for sector in SECTORS:
                first_batches: List[float] = []
                second_batches: List[float] = []
                difference_batches: List[float] = []
                for key in keys:
                    first = sector_values(rows[key], "first")[sector]
                    second = sector_values(rows[key], "second")[sector]
                    first_batches.append(first)
                    second_batches.append(second)
                    difference_batches.append(first - second)
                first_mean, first_se = mean_se(first_batches)
                second_mean, second_se = mean_se(second_batches)
                difference, difference_se = mean_se(difference_batches)
                normalized = difference / delta_cos4
                normalized_se = difference_se / abs(delta_cos4)
                exponent = 13.0 / 8.0 if sector in ("odd", "matching_function") else 1.0
                scale = n ** exponent
                output.append({
                    "N": n,
                    "channel": channel,
                    "sector": sector,
                    "first_rep": [a1, b1],
                    "second_rep": [a2, b2],
                    "cos4_first": cos4(a1, b1),
                    "cos4_second": cos4(a2, b2),
                    "delta_cos4_first_minus_second": delta_cos4,
                    "first_estimate": first_mean,
                    "first_batch_se": first_se,
                    "second_estimate": second_mean,
                    "second_batch_se": second_se,
                    "difference_first_minus_second": difference,
                    "difference_batch_se": difference_se,
                    "difference_z": difference / difference_se if difference_se else None,
                    "normalized_by_delta_cos4": normalized,
                    "normalized_batch_se": normalized_se,
                    "hypothesis_N_exponent": exponent,
                    "hypothesis_scaled_amplitude": scale * normalized,
                    "hypothesis_scaled_batch_se": scale * normalized_se,
                })
    return output


def write_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    scalar_fields = [name for name in rows[0] if name not in ("first_rep", "second_rep")]
    fields = ["a1", "b1", "a2", "b2"] + scalar_fields
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            flat = dict(row)
            first = flat.pop("first_rep")
            second = flat.pop("second_rep")
            flat.update({"a1": first[0], "b1": first[1], "a2": second[0], "b2": second[1]})
            writer.writerow(flat)


def batch_ids(rows: Dict[BatchKey, Dict[str, int]], n: int) -> List[int]:
    return sorted(key[2] for key in rows if key[:2] == (n, CHANNELS[0]))


def correlation_matrix(covariance: List[List[float]]) -> List[List[object]]:
    output: List[List[object]] = []
    for i, row in enumerate(covariance):
        output_row: List[object] = []
        for j, value in enumerate(row):
            denominator = math.sqrt(max(covariance[i][i], 0.0) * max(covariance[j][j], 0.0))
            if denominator:
                output_row.append(max(-1.0, min(1.0, value / denominator)))
            else:
                output_row.append(None)
        output.append(output_row)
    return output


def matrix_payload(labels: List[str], vectors: List[List[float]]) -> Dict[str, object]:
    covariance = sample_covariance(vectors)
    return {
        "labels": labels,
        "batch_count": len(vectors),
        "unit": "covariance of equal-size batch means",
        "covariance": covariance,
        "correlation": correlation_matrix(covariance),
    }


def covariance_payload(
    rows: Dict[BatchKey, Dict[str, int]],
    designs: Dict[int, Tuple[int, int, int, int]],
) -> Dict[str, object]:
    by_n: Dict[str, object] = {}
    for n in sorted(designs):
        batches = batch_ids(rows, n)
        raw_labels = [
            orientation + "_" + estimator + "_" + channel
            for channel in CHANNELS
            for orientation in ("first", "second")
            for estimator in ("primal", "matching")
        ]
        raw_vectors: List[List[float]] = []
        effect_labels = [
            "delta_" + sector + "_" + channel
            for channel in CHANNELS
            for sector in SECTORS
        ]
        effect_vectors: List[List[float]] = []
        for batch in batches:
            raw_vector: List[float] = []
            effect_vector: List[float] = []
            for channel in CHANNELS:
                row = rows[(n, channel, batch)]
                first = sector_values(row, "first")
                second = sector_values(row, "second")
                raw_vector.extend([
                    first["primal"], first["matching"],
                    second["primal"], second["matching"],
                ])
                effect_vector.extend([first[sector] - second[sector] for sector in SECTORS])
            raw_vectors.append(raw_vector)
            effect_vectors.append(effect_vector)
        by_n[str(n)] = {
            "raw_orientation_channel_matrix": matrix_payload(raw_labels, raw_vectors),
            "orientation_sector_effect_matrix": matrix_payload(effect_labels, effect_vectors),
        }
    return {
        "analysis": "joint covariance/correlation of paired batch means",
        "channels": list(CHANNELS),
        "by_N": by_n,
    }


def gls_rows(
    rows: Dict[BatchKey, Dict[str, int]], n: int, target: str
) -> List[List[float]]:
    output: List[List[float]] = []
    for batch in batch_ids(rows, n):
        vector: List[float] = []
        for channel in CHANNELS:
            row = rows[(n, channel, batch)]
            first = sector_values(row, "first")["odd"]
            second = sector_values(row, "second")["odd"]
            if target == "first":
                vector.append(first)
            elif target == "second":
                vector.append(second)
            else:
                vector.append(first - second)
        output.append(vector)
    return output


def freeze_gls(
    rows: Dict[BatchKey, Dict[str, int]],
    designs: Dict[int, Tuple[int, int, int, int]],
    metadata: Dict[str, object],
    source: Path,
    target: str,
) -> Dict[str, object]:
    by_n: Dict[str, object] = {}
    for n in sorted(designs):
        pilot = gls_rows(rows, n, target)
        covariance = sample_covariance(pilot)
        max_spread = max(max(row) - min(row) for row in pilot)
        scale = max(1.0, max(abs(value) for row in pilot for value in row))
        identical = max_spread <= 1e-14 * scale
        if identical:
            # The exact matching-channel identity makes every sum-one vector
            # equally efficient.  Freeze the symmetric representative rather
            # than report solver noise as an empirical optimization.
            weights = [1.0 / len(CHANNELS)] * len(CHANNELS)
            ridge = 0.0
            best_index = 0
        else:
            weights, ridge = minimum_variance_weights(covariance)
            best_index = min(range(len(CHANNELS)), key=lambda index: covariance[index][index])
        by_n[str(n)] = {
            "weights": weights,
            "applied_diagonal_ridge": ridge,
            "pilot_best_single_channel": CHANNELS[best_index],
            "all_D_channels_identical_batchwise": identical,
            "maximum_pilot_channel_spread": max_spread,
            "pilot_covariance": covariance,
            "pilot_correlation": correlation_matrix(covariance),
            "pilot_batch_count": len(pilot),
        }
    return {
        "schema": "pilot-frozen equal-mean D-channel GLS weights v1",
        "target": target,
        "channel_names": list(CHANNELS),
        "source_batches": str(source),
        "pilot_rng": {
            "seed": metadata.get("seed"),
            "replica_counter_first": metadata.get("replica_counter_first"),
            "replica_counter_last_exclusive": metadata.get("replica_counter_last_exclusive"),
        },
        "by_N": by_n,
    }


def sample_variance(values: List[float]) -> float:
    mean = math.fsum(values) / len(values)
    return math.fsum((value - mean) ** 2 for value in values) / (len(values) - 1)


def combine(vectors: List[List[float]], weights: List[float]) -> List[float]:
    return [math.fsum(value * weight for value, weight in zip(row, weights)) for row in vectors]


def validate_independent_evaluation(
    frozen: Dict[str, object], metadata: Dict[str, object]
) -> None:
    pilot_rng = frozen.get("pilot_rng")
    if not isinstance(pilot_rng, dict):
        raise ValueError("frozen weights lack pilot RNG provenance")
    if pilot_rng.get("seed") != metadata.get("seed"):
        return
    pilot_first = pilot_rng.get("replica_counter_first")
    pilot_last = pilot_rng.get("replica_counter_last_exclusive")
    evaluation_first = metadata.get("replica_counter_first")
    evaluation_last = metadata.get("replica_counter_last_exclusive")
    if not all(isinstance(value, int) for value in (
        pilot_first, pilot_last, evaluation_first, evaluation_last
    )):
        raise ValueError("cannot verify pilot/evaluation RNG counter separation")
    if max(pilot_first, evaluation_first) < min(pilot_last, evaluation_last):
        raise ValueError("evaluation RNG counters overlap the pilot used to freeze GLS weights")


def evaluate_gls(
    rows: Dict[BatchKey, Dict[str, int]],
    designs: Dict[int, Tuple[int, int, int, int]],
    metadata: Dict[str, object],
    frozen: Dict[str, object],
) -> Tuple[Dict[str, object], List[Dict[str, object]]]:
    validate_independent_evaluation(frozen, metadata)
    if frozen.get("channel_names") != list(CHANNELS):
        raise ValueError("frozen GLS channel order does not match analyzer")
    target = frozen.get("target")
    if target not in ("first", "second", "orientation_difference"):
        raise ValueError("frozen GLS target is invalid")
    frozen_by_n = frozen.get("by_N")
    if not isinstance(frozen_by_n, dict):
        raise ValueError("frozen weights lack by_N entries")
    output_rows: List[Dict[str, object]] = []
    by_n: Dict[str, object] = {}
    for n in sorted(designs):
        entry = frozen_by_n.get(str(n))
        if not isinstance(entry, dict):
            raise ValueError("frozen weights do not contain N=" + str(n))
        weights = entry.get("weights")
        if not isinstance(weights, list) or len(weights) != len(CHANNELS):
            raise ValueError("invalid frozen weight vector for N=" + str(n))
        weights = [float(value) for value in weights]
        if not math.isclose(math.fsum(weights), 1.0, rel_tol=1e-10, abs_tol=1e-10):
            raise ValueError("frozen weights do not sum to one")
        vectors = gls_rows(rows, n, target)
        optimized = combine(vectors, weights)
        equal = combine(vectors, [1.0 / len(CHANNELS)] * len(CHANNELS))
        best_name = entry.get("pilot_best_single_channel")
        if best_name not in CHANNELS:
            raise ValueError("invalid pilot-frozen best single channel")
        series: List[Tuple[str, List[float]]] = [
            ("optimized_gls", optimized),
            ("equal_weight", equal),
            ("pilot_best_single_" + str(best_name), [row[CHANNELS.index(best_name)] for row in vectors]),
        ]
        series.extend(("single_" + channel, [row[index] for row in vectors])
                      for index, channel in enumerate(CHANNELS))
        optimized_variance = sample_variance(optimized)
        n_rows: Dict[str, object] = {
            "target": target,
            "weights": weights,
            "evaluation_batch_count": len(vectors),
            "evaluation_covariance": sample_covariance(vectors),
        }
        estimates: Dict[str, object] = {}
        for name, values in series:
            mean = math.fsum(values) / len(values)
            variance = sample_variance(values)
            ratio = variance / optimized_variance if optimized_variance > 0 else None
            result = {
                "N": n,
                "target": target,
                "estimator": name,
                "mean": mean,
                "batch_mean_variance": variance,
                "variance_of_overall_mean": variance / len(values),
                "variance_reduction_vs_optimized": ratio,
            }
            output_rows.append(result)
            estimates[name] = result
        n_rows["estimators"] = estimates
        by_n[str(n)] = n_rows
    return {
        "analysis": "independent evaluation of pilot-frozen equal-mean D-channel GLS weights",
        "frozen_source_batches": frozen.get("source_batches"),
        "evaluation_rng": {
            "seed": metadata.get("seed"),
            "replica_counter_first": metadata.get("replica_counter_first"),
            "replica_counter_last_exclusive": metadata.get("replica_counter_last_exclusive"),
        },
        "by_N": by_n,
    }, output_rows


def write_evaluation_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batches", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--csv", type=Path, required=True, help="long-form sector_effects CSV")
    parser.add_argument("--covariance-json", type=Path)
    gls = parser.add_mutually_exclusive_group()
    gls.add_argument("--freeze-gls", type=Path, help="write pilot-frozen D-channel weights")
    gls.add_argument("--frozen-gls", type=Path, help="evaluate weights frozen on independent data")
    parser.add_argument(
        "--gls-target", choices=("first", "second", "orientation_difference"),
        default="orientation_difference",
    )
    parser.add_argument("--evaluation-json", type=Path)
    parser.add_argument("--evaluation-csv", type=Path)
    args = parser.parse_args()

    with args.metadata.open(encoding="utf-8") as handle:
        metadata = json.load(handle)
    if not isinstance(metadata, dict):
        raise ValueError("metadata must contain an object")
    rows, designs = read_rows(args.batches)
    summaries = analyze(rows, designs)
    payload = {
        "analysis": "paired batch means for same-N Gaussian orientations",
        "source_batches": str(args.batches),
        "source_metadata": str(args.metadata),
        "metadata": metadata,
        "definitions": {
            "even": "(R_primal + R_white_matching)/2",
            "odd": "(R_primal - R_white_matching)/2",
            "matching_function": "R_primal - R_white_matching = 2*odd",
            "orientation_difference": "first representation minus second representation",
        },
        "scaling_diagnostics": {
            "primal_matching_even": "N * difference / delta_cos4 (L^-2 hypothesis)",
            "odd_matching_function": "N^(13/8) * difference / delta_cos4 (L^-13/4 hypothesis)",
        },
        "caveats": [
            "batch standard errors quantify Monte Carlo noise, not finite-size systematic error",
            "a resolved matching-even harmonic is required before interpreting a null matching-odd difference",
            "the five wrapping channels are simultaneous observables, not independent replicas",
        ],
        "summaries": summaries,
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_csv(args.csv, summaries)
    if args.covariance_json:
        args.covariance_json.parent.mkdir(parents=True, exist_ok=True)
        args.covariance_json.write_text(
            json.dumps(covariance_payload(rows, designs), indent=2) + "\n",
            encoding="utf-8",
        )
    if args.freeze_gls:
        frozen = freeze_gls(rows, designs, metadata, args.batches, args.gls_target)
        args.freeze_gls.parent.mkdir(parents=True, exist_ok=True)
        args.freeze_gls.write_text(json.dumps(frozen, indent=2) + "\n", encoding="utf-8")
    if args.frozen_gls:
        if not args.evaluation_json or not args.evaluation_csv:
            parser.error("--frozen-gls requires --evaluation-json and --evaluation-csv")
        with args.frozen_gls.open(encoding="utf-8") as handle:
            frozen = json.load(handle)
        if not isinstance(frozen, dict):
            raise ValueError("frozen weights must contain an object")
        evaluation, evaluation_rows = evaluate_gls(rows, designs, metadata, frozen)
        args.evaluation_json.parent.mkdir(parents=True, exist_ok=True)
        args.evaluation_json.write_text(json.dumps(evaluation, indent=2) + "\n", encoding="utf-8")
        write_evaluation_csv(args.evaluation_csv, evaluation_rows)

    for row in summaries:
        if row["sector"] not in ("even", "matching_function"):
            continue
        z = row["difference_z"]
        z_text = "undefined" if z is None else "{:.2f}".format(z)
        print(
            "N={N} {channel} {sector}: diff={difference_first_minus_second:.6g} "
            "+/- {difference_batch_se:.3g}, z={z}, scaled={hypothesis_scaled_amplitude:.6g}"
            .format(z=z_text, **row)
        )
    print("wrote " + str(args.json))
    print("wrote " + str(args.csv))
    if args.covariance_json:
        print("wrote " + str(args.covariance_json))
    if args.freeze_gls:
        print("wrote " + str(args.freeze_gls))
    if args.frozen_gls:
        print("wrote " + str(args.evaluation_json))
        print("wrote " + str(args.evaluation_csv))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
