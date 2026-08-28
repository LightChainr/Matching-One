#!/usr/bin/env python3
"""Audit the Issue #43 matching-even source/target wrapping-channel contract.

The frozen even-sector amplitude was estimated in the ``either`` wrapping
channel, whereas the threshold-rank production engine records only rank-2
``cross`` wrapping.  On the torus, primal/matching complementarity gives

    R_G^cross + R_hat^either = 1,
    R_G^either + R_hat^cross = 1,

so S_cross = 1 - S_either and orientation differences acquire a minus sign.
This script checks that relation in the pre-target P31 batch artifact and then
reports the deterministic, zero-target-fit channel transport of the already
published Issue #43 score.  It does not replace the literal frozen score.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple


EXPECTED_PRODUCTION_CHANNEL = "rank-2 cross wrapping"


def quadratic_2(vector: Sequence[float], covariance: Sequence[Sequence[float]]) -> float:
    a, c, b = covariance[0][0], covariance[0][1], covariance[1][1]
    determinant = a * b - c * c
    if determinant <= 0.0:
        raise ValueError("score covariance is not positive definite")
    x, y = vector
    return (b * x * x - 2.0 * c * x * y + a * y * y) / determinant


def weighted_constant(rows: Sequence[Tuple[float, float]]) -> Dict[str, float]:
    weights = [1.0 / standard_error**2 for _value, standard_error in rows]
    mean = math.fsum(
        weight * value for weight, (value, _standard_error) in zip(weights, rows)
    ) / math.fsum(weights)
    standard_error = math.sqrt(1.0 / math.fsum(weights))
    chi_square = math.fsum(
        ((value - mean) / row_standard_error) ** 2
        for value, row_standard_error in rows
    )
    return {
        "mean": mean,
        "standard_error": standard_error,
        "chi_square": chi_square,
        "degrees_of_freedom": len(rows) - 1,
    }


def read_source_rows(path: Path) -> Dict[str, List[Dict[str, float]]]:
    selected: Dict[str, List[Dict[str, float]]] = {"cross": [], "either": []}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "N", "channel", "sector", "difference_first_minus_second",
            "difference_batch_se", "delta_cos4_first_minus_second",
            "hypothesis_scaled_amplitude", "hypothesis_scaled_batch_se",
        }
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError("source CSV missing: " + ", ".join(sorted(missing)))
        for raw in reader:
            if raw["sector"] != "even" or raw["channel"] not in selected:
                continue
            selected[raw["channel"]].append({
                "N": int(raw["N"]),
                "difference": float(raw["difference_first_minus_second"]),
                "standard_error": float(raw["difference_batch_se"]),
                "delta_cos4": float(raw["delta_cos4_first_minus_second"]),
                "scaled_amplitude": float(raw["hypothesis_scaled_amplitude"]),
                "scaled_standard_error": float(raw["hypothesis_scaled_batch_se"]),
            })
    for channel in selected:
        selected[channel].sort(key=lambda row: int(row["N"]))
    if not selected["cross"] or [row["N"] for row in selected["cross"]] != [
        row["N"] for row in selected["either"]
    ]:
        raise ValueError("source cross/either size grids differ or are empty")
    for cross, either in zip(selected["cross"], selected["either"]):
        for field in ("difference", "scaled_amplitude"):
            if not math.isclose(cross[field], -either[field], rel_tol=0.0, abs_tol=2e-15):
                raise ValueError("P31 cross/either sign identity failed at N={}".format(cross["N"]))
        for field in ("standard_error", "scaled_standard_error", "delta_cos4"):
            if not math.isclose(cross[field], either[field], rel_tol=2e-12, abs_tol=2e-15):
                raise ValueError("P31 cross/either uncertainty/geometry mismatch at N={}".format(cross["N"]))
    return selected


def verify_batch_identities(path: Path) -> Dict[str, int]:
    required = {
        "n", "batch", "samples", "channel", "first_primal_sum",
        "first_matching_sum", "second_primal_sum", "second_matching_sum",
    }
    records: Dict[Tuple[int, int, str], Dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError("source batch CSV missing: " + ", ".join(sorted(missing)))
        for raw in reader:
            if raw["channel"] not in ("cross", "either"):
                continue
            key = (int(raw["n"]), int(raw["batch"]), raw["channel"])
            if key in records:
                raise ValueError("duplicate source batch row: {}".format(key))
            records[key] = raw
    pairs = sorted({(n, batch) for n, batch, _channel in records})
    if not pairs:
        raise ValueError("source batch CSV has no cross/either rows")
    for n, batch in pairs:
        cross = records.get((n, batch, "cross"))
        either = records.get((n, batch, "either"))
        if cross is None or either is None:
            raise ValueError("source batch cross/either grid is incomplete")
        samples = int(cross["samples"])
        if int(either["samples"]) != samples:
            raise ValueError("cross/either source batch sample counts differ")
        for prefix in ("first", "second"):
            if (
                int(cross[prefix + "_primal_sum"])
                + int(either[prefix + "_matching_sum"])
                != samples
                or int(either[prefix + "_primal_sum"])
                + int(cross[prefix + "_matching_sum"])
                != samples
            ):
                raise ValueError(
                    "configuration complement identity failed at N={}, batch={}, {}"
                    .format(n, batch, prefix)
                )
    return {
        "size_count": len({n for n, _batch in pairs}),
        "batch_count": len(pairs),
        "integer_identity_max_abs_residual": 0,
    }


def read_metadata(paths: Sequence[Path]) -> List[Dict[str, object]]:
    output = []
    for path in paths:
        with path.open(encoding="utf-8") as handle:
            metadata = json.load(handle)
        if metadata.get("channel") != EXPECTED_PRODUCTION_CHANNEL:
            raise ValueError("Issue #43 production is not the expected cross channel")
        designs = metadata.get("designs")
        if not isinstance(designs, list) or len(designs) != 1:
            raise ValueError("each production metadata file must contain one design")
        output.append({
            "path": str(path),
            "N": int(designs[0]["N"]),
            "channel": metadata["channel"],
        })
    output.sort(key=lambda row: int(row["N"]))
    if [row["N"] for row in output] != [185, 265]:
        raise ValueError("metadata must describe N=185 and N=265")
    return output


def audit(
    source_csv: Path,
    source_batches: Path,
    primary_score_json: Path,
    metadata_paths: Sequence[Path],
) -> Dict[str, object]:
    source = read_source_rows(source_csv)
    batch_check = verify_batch_identities(source_batches)
    metadata = read_metadata(metadata_paths)
    with primary_score_json.open(encoding="utf-8") as handle:
        primary = json.load(handle)
    score = primary["scores"]["DeltaS"]
    observed = [float(value) for value in score["observed"]]
    frozen_either_mean = [float(value) for value in score["frozen_mean"]]
    covariance = [[float(value) for value in row] for row in score["target_covariance"]]
    if not all(value > 0.0 for value in frozen_either_mean):
        raise ValueError("literal frozen either-channel means are not positive")

    # Exact channel transport only: no exponent, amplitude, sign, or covariance
    # is estimated from N=185/265.  The literal frozen score remains untouched.
    transported_cross_mean = [-value for value in frozen_either_mean]
    residual = [
        observed[index] - transported_cross_mean[index] for index in range(2)
    ]
    transported_chi_square = quadratic_2(residual, covariance)
    marginal_z = [
        residual[index] / math.sqrt(covariance[index][index]) for index in range(2)
    ]

    fits = {}
    for channel in ("either", "cross"):
        fits[channel] = weighted_constant([
            (float(row["scaled_amplitude"]), float(row["scaled_standard_error"]))
            for row in source[channel]
        ])

    return {
        "audit": "Issue #43 matching-even wrapping-channel sign contract",
        "classification": "frozen prediction generator/protocol channel mismatch",
        "not_classified_as": [
            "physical sign reversal",
            "orientation-order mismatch",
            "delta-cos4 sign mismatch",
            "threshold-rank reconstruction error",
        ],
        "evidence_boundary": (
            "Post-reveal deterministic audit. The literal frozen positive score and its "
            "reported failure remain immutable; the transported score is not a replacement "
            "preregistered endpoint and uses zero target-fit parameters."
        ),
        "source": {
            "path": str(source_csv),
            "batch_path": str(source_batches),
            "frozen_channel": "either",
            "sector": "even",
            "rows": source,
            "constant_amplitude_fits": fits,
            "configuration_level_batch_check": batch_check,
        },
        "target": {
            "score_path": str(primary_score_json),
            "metadata": metadata,
            "production_channel": "cross",
            "orientation_order": "first_minus_second",
        },
        "exact_transport": {
            "configuration_identities": [
                "R_G_cross + R_hat_either = 1",
                "R_G_either + R_hat_cross = 1",
            ],
            "sector_identity": "S_cross = 1 - S_either",
            "orientation_difference_identity": "DeltaS_cross = -DeltaS_either",
            "frozen_either_mean": frozen_either_mean,
            "transported_cross_mean": transported_cross_mean,
            "observed_cross": observed,
            "residual": residual,
            "marginal_z": marginal_z,
            "transported_chi_square": transported_chi_square,
            "degrees_of_freedom": 2,
            "literal_frozen_positive_chi_square": float(score["target_chi_square"]),
            "zero_chi_square": float(score["zero_chi_square"]),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-csv", type=Path, required=True)
    parser.add_argument("--source-batches", type=Path, required=True)
    parser.add_argument("--primary-score", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, action="append", required=True)
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.source_csv, args.source_batches, args.primary_score, args.metadata)
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
