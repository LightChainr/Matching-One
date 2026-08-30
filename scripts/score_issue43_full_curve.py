#!/usr/bin/env python3
"""Score the frozen Issue #43 N=185/265 full-curve predictions.

This standard-library-only scorer consumes two threshold-rank production
triples (histogram, moments, metadata), reconstructs the matching-odd M and
matching-even S sectors at the frozen p_ref, and scores the two-size residuals
with target sampling uncertainty plus the fully correlated source-amplitude
uncertainty frozen in the prediction artifact.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import shlex
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple


P_REF = 0.592746050790
BATCHES = 100
COMMIT_RE = re.compile(r"^[0-9a-fA-F]{40}$")
PREDICTION_SHA256 = "a370e79a10854341fac3ee75e8c518dbf3533e8c077cba2c2ec1018178144f44"
DESIGNS = {
    185: ((13, 4), (11, 8)),
    265: ((16, 3), (12, 11)),
}
PREDICTIONS = {
    "DeltaM": {
        185: (0.00019223156676869253, 0.000008581548700390585),
        265: (0.0001567390401230751, 0.000006997101093636327),
    },
    "DeltaS": {
        185: (0.00006752163745881449, 0.000005964854295017754),
        265: (0.00006891944697034459, 0.000006088336639081946),
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def command_option(command: str, name: str) -> str:
    tokens = shlex.split(command)
    values = [tokens[i + 1] for i, token in enumerate(tokens[:-1]) if token == name]
    if len(values) != 1:
        raise ValueError("metadata command must contain exactly one {}".format(name))
    return values[0]


def validate_metadata(path: Path) -> Dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        metadata = json.load(handle)
    if metadata.get("engine") != "same-N Gaussian threshold-rank Newman-Ziff":
        raise ValueError("unexpected threshold-rank engine")
    designs = metadata.get("designs")
    if not isinstance(designs, list) or len(designs) != 1:
        raise ValueError("each Issue #43 run must contain exactly one design")
    design = designs[0]
    n = int(design.get("N", -1))
    if n not in DESIGNS:
        raise ValueError("Issue #43 primary runs must be N=185 or N=265")
    first, second = DESIGNS[n]
    if tuple(design.get("first", ())) != first or tuple(design.get("second", ())) != second:
        raise ValueError("metadata orientation order differs from the frozen artifact")
    commit = str(metadata.get("git_commit", ""))
    if not COMMIT_RE.fullmatch(commit):
        raise ValueError("source commit must be a full 40-hex id")
    try:
        samples = int(metadata["samples_per_pair"])
        batches = int(metadata["batches"])
        seed = int(metadata["seed"])
        counter_first = int(metadata["replica_counter_first"])
        counter_last = int(metadata["replica_counter_last_exclusive"])
        command = str(metadata["command"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("missing or invalid production metadata") from error
    if batches != BATCHES or samples <= 0 or samples % batches:
        raise ValueError("Issue #43 scoring requires exactly 100 equal batches")
    if seed < 0 or counter_first < 0 or counter_last - counter_first != samples:
        raise ValueError("seed/counter range is inconsistent with samples")
    expected_integer_options = {
        "--n": n, "--samples": samples, "--batches": batches,
        "--seed": seed, "--replica-offset": counter_first,
    }
    for option, expected in expected_integer_options.items():
        try:
            actual = int(command_option(command, option))
        except ValueError as error:
            raise ValueError("command provenance mismatch for {}".format(option)) from error
        if actual != expected:
            raise ValueError("command provenance mismatch for {}".format(option))
    if command_option(command, "--git-commit").lower() != commit.lower():
        raise ValueError("command/source commit mismatch")
    return {
        "N": n, "first": first, "second": second, "commit": commit.lower(),
        "samples": samples, "batches": batches, "seed": seed,
        "counter_first": counter_first, "counter_last": counter_last,
        # A threshold-rank run is microcanonical and has no p_ref.  This is
        # the scorer's frozen reconstruction coordinate, not run metadata.
        "scorer_p_ref": P_REF,
    }


def read_histograms(path: Path, run: Mapping[str, object]) -> Dict[Tuple[str, int], Dict[str, object]]:
    required = {"n", "a", "b", "orientation", "batch", "samples", "kind", "k", "count"}
    records: Dict[Tuple[str, int], Dict[str, object]] = {}
    n = int(run["N"])
    representations = {"first": tuple(run["first"]), "second": tuple(run["second"])}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError("histogram CSV missing: " + ", ".join(sorted(missing)))
        for raw in reader:
            orientation = raw["orientation"]
            batch = int(raw["batch"])
            if int(raw["n"]) != n or orientation not in representations:
                raise ValueError("histogram N/orientation disagrees with metadata")
            if (int(raw["a"]), int(raw["b"])) != representations[orientation]:
                raise ValueError("histogram representation disagrees with metadata")
            samples = int(raw["samples"])
            if samples != int(run["samples"]) // BATCHES:
                raise ValueError("histogram batch samples disagree with metadata")
            kind = raw["kind"]
            rank = int(raw["k"])
            count = int(raw["count"])
            if kind not in ("minus", "plus") or not 1 <= rank <= n or count <= 0:
                raise ValueError("invalid histogram row")
            key = (orientation, batch)
            record = records.setdefault(key, {
                "samples": samples, "minus": [0] * (n + 1), "plus": [0] * (n + 1)
            })
            if record["samples"] != samples:
                raise ValueError("inconsistent samples within histogram")
            record[kind][rank] += count
    expected = {(orientation, batch) for orientation in representations for batch in range(BATCHES)}
    if set(records) != expected:
        raise ValueError("histogram lacks the complete orientation/batch grid")
    for record in records.values():
        if sum(record["minus"]) != record["samples"] or sum(record["plus"]) != record["samples"]:
            raise ValueError("histogram marginal total differs from samples")
    return records


def validate_moments(
    path: Path, run: Mapping[str, object], records: Mapping[Tuple[str, int], Mapping[str, object]]
) -> None:
    required = {
        "n", "a", "b", "orientation", "batch", "samples", "sum_kminus", "sum_kplus",
        "sum_kminus2", "sum_kplus2", "sum_product", "sum_gap", "sum_gap2",
    }
    seen = set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError("moments CSV missing: " + ", ".join(sorted(missing)))
        for raw in reader:
            key = (raw["orientation"], int(raw["batch"]))
            if int(raw["n"]) != run["N"] or key not in records or key in seen:
                raise ValueError("moments N/orientation/batch mismatch")
            seen.add(key)
            first, second = run["first"], run["second"]
            expected_rep = first if key[0] == "first" else second
            if (int(raw["a"]), int(raw["b"])) != tuple(expected_rep):
                raise ValueError("moment representation mismatch")
            record = records[key]
            minus, plus = record["minus"], record["plus"]
            calculated = {
                "sum_kminus": sum(k * count for k, count in enumerate(minus)),
                "sum_kplus": sum(k * count for k, count in enumerate(plus)),
                "sum_kminus2": sum(k * k * count for k, count in enumerate(minus)),
                "sum_kplus2": sum(k * k * count for k, count in enumerate(plus)),
            }
            if int(raw["samples"]) != record["samples"] or any(
                int(raw[name]) != value for name, value in calculated.items()
            ):
                raise ValueError("moments disagree with histogram marginals")
            if int(raw["sum_gap"]) != calculated["sum_kplus"] - calculated["sum_kminus"]:
                raise ValueError("moment gap identity failed")
    if seen != set(records):
        raise ValueError("moments lack the complete histogram batch grid")


def tail(histogram: Sequence[int], samples: int, p: float) -> float:
    n = len(histogram) - 1
    q = 1.0 - p
    probability = q ** n
    cumulative = 0
    value = 0.0
    for occupied in range(n + 1):
        if occupied:
            cumulative += histogram[occupied]
        value += cumulative * probability
        if occupied < n:
            probability *= (n - occupied) * p / ((occupied + 1) * q)
    return value / samples


def reconstruct(records: Mapping[Tuple[str, int], Mapping[str, object]]) -> Dict[str, List[float]]:
    output = {"DeltaM": [], "DeltaS": []}
    for batch in range(BATCHES):
        sector = {}
        for orientation in ("first", "second"):
            record = records[(orientation, batch)]
            r_g = tail(record["plus"], record["samples"], P_REF)
            r_hat = 1.0 - tail(record["minus"], record["samples"], P_REF)
            sector[orientation] = {"M": r_g - r_hat, "S": (r_g + r_hat) / 2.0}
        output["DeltaM"].append(sector["first"]["M"] - sector["second"]["M"])
        output["DeltaS"].append(sector["first"]["S"] - sector["second"]["S"])
    return output


def mean_se(values: Sequence[float]) -> Tuple[float, float]:
    mean = math.fsum(values) / len(values)
    variance = math.fsum((value - mean) ** 2 for value in values)
    return mean, math.sqrt(variance / (len(values) * (len(values) - 1)))


def quadratic_2(vector: Sequence[float], covariance: Sequence[Sequence[float]]) -> float:
    a, c, b = covariance[0][0], covariance[0][1], covariance[1][1]
    determinant = a * b - c * c
    if determinant <= 0.0:
        raise ValueError("score covariance is not positive definite")
    x, y = vector
    return (b * x * x - 2.0 * c * x * y + a * y * y) / determinant


def analyze(runs: Mapping[int, Mapping[str, object]], prediction_path: Path) -> Dict[str, object]:
    if sha256(prediction_path) != PREDICTION_SHA256:
        raise ValueError("frozen prediction artifact hash mismatch")
    if set(runs) != set(DESIGNS):
        raise ValueError("both and only N=185,265 runs are required")
    commits = {str(run["metadata"]["commit"]) for run in runs.values()}
    if len(commits) != 1:
        raise ValueError("N=185 and N=265 must use the same source commit")
    estimates = {}
    scores = {}
    for sector in ("DeltaM", "DeltaS"):
        estimates[sector] = {}
        observed, sampling_se, predicted, source_se = [], [], [], []
        for n in (185, 265):
            mean, se = mean_se(runs[n]["sectors"][sector])
            target, target_se = PREDICTIONS[sector][n]
            estimates[sector][n] = {"estimate": mean, "sampling_se": se}
            observed.append(mean)
            sampling_se.append(se)
            predicted.append(target)
            source_se.append(target_se)
        residual = [observed[i] - predicted[i] for i in range(2)]
        target_covariance = [[
            (sampling_se[i] ** 2 if i == j else 0.0) + source_se[i] * source_se[j]
            for j in range(2)
        ] for i in range(2)]
        zero_covariance = [[sampling_se[i] ** 2 if i == j else 0.0 for j in range(2)] for i in range(2)]
        scores[sector] = {
            "observed": observed, "sampling_se": sampling_se,
            "frozen_mean": predicted, "source_coefficient_se": source_se,
            "source_error_correlation": [[1.0, 1.0], [1.0, 1.0]],
            "target_residual": residual, "target_covariance": target_covariance,
            "target_z_marginal": [residual[i] / math.sqrt(target_covariance[i][i]) for i in range(2)],
            "target_chi_square": quadratic_2(residual, target_covariance), "target_df": 2,
            "zero_chi_square": quadratic_2(observed, zero_covariance), "zero_df": 2,
            "zero_z_marginal": [observed[i] / sampling_se[i] for i in range(2)],
        }
    return {
        "protocol": "Issue #43 prospective N=185/265 two-spin4 full-curve score",
        "status": "frozen primary score; no refit",
        "p_ref": P_REF,
        "prediction_artifact": str(prediction_path),
        "prediction_artifact_sha256": PREDICTION_SHA256,
        "source_commit": next(iter(commits)),
        "runs": {str(n): runs[n]["metadata"] for n in (185, 265)},
        "estimates": estimates,
        "scores": scores,
    }


def write_outputs(result: Mapping[str, object], json_path: Path, csv_path: Path) -> None:
    json_path.parent.mkdir(parents=True, exist_ok=True)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    fields = [
        "sector", "N", "estimate", "sampling_se", "frozen_mean", "source_coefficient_se",
        "residual", "combined_marginal_se", "residual_z", "zero_z",
        "two_size_target_chi_square", "two_size_target_df", "two_size_zero_chi_square",
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for sector in ("DeltaM", "DeltaS"):
            score = result["scores"][sector]
            for i, n in enumerate((185, 265)):
                writer.writerow({
                    "sector": sector, "N": n, "estimate": score["observed"][i],
                    "sampling_se": score["sampling_se"][i], "frozen_mean": score["frozen_mean"][i],
                    "source_coefficient_se": score["source_coefficient_se"][i],
                    "residual": score["target_residual"][i],
                    "combined_marginal_se": math.sqrt(score["target_covariance"][i][i]),
                    "residual_z": score["target_z_marginal"][i], "zero_z": score["zero_z_marginal"][i],
                    "two_size_target_chi_square": score["target_chi_square"],
                    "two_size_target_df": score["target_df"],
                    "two_size_zero_chi_square": score["zero_chi_square"],
                })


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run", nargs=3, action="append", metavar=("HIST", "MOMENTS", "METADATA"), required=True,
        help="repeat once for N=185 and once for N=265",
    )
    parser.add_argument("--predictions", type=Path, default=Path("predictions/two_spin4_heldout_20260828.yaml"))
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--csv", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    runs = {}
    for hist_name, moments_name, metadata_name in args.run:
        metadata = validate_metadata(Path(metadata_name))
        n = int(metadata["N"])
        if n in runs:
            raise ValueError("duplicate N run")
        records = read_histograms(Path(hist_name), metadata)
        validate_moments(Path(moments_name), metadata, records)
        runs[n] = {"metadata": metadata, "sectors": reconstruct(records)}
    result = analyze(runs, args.predictions)
    write_outputs(result, args.json, args.csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
