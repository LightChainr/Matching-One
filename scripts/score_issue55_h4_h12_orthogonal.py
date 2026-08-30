#!/usr/bin/env python3
"""Variance-only pilot and frozen score for Issue #55 H4/H12 rows.

Pilot mode is deliberately blind to target means: it emits only centered
variance estimates, prospective signal-to-noise/power, and the deterministic
sample-count choice.  Final mode applies the preregistered scoring order and
reports the exact two-row coordinates that separately identify A4 and A12.
"""

from __future__ import annotations

import argparse
import csv
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import re
import shlex
from statistics import NormalDist
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "experiments" / "issue55_h4_h12_orthogonal_acquisition_20260830.json"
COMMIT_RE = re.compile(r"^[0-9a-fA-F]{40}$")
NORMAL = NormalDist()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def command_option(command: str, name: str, width: int = 1) -> Tuple[str, ...]:
    tokens = shlex.split(command)
    hits = [index for index, token in enumerate(tokens) if token == name]
    _require(len(hits) == 1, "metadata command must contain exactly one {}".format(name))
    index = hits[0]
    _require(index + width < len(tokens), "metadata command truncates {}".format(name))
    return tuple(tokens[index + 1:index + 1 + width])


def load_manifest(path: Path) -> Dict[str, Any]:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    _require(
        manifest.get("schema") == "matching-one/issue55-h4-h12-orthogonal-acquisition/v1",
        "unknown acquisition schema",
    )
    frozen = manifest["frozen_design"]
    _require(sha256(ROOT / frozen["path"]) == frozen["sha256"], "frozen design hash drift")
    _require(manifest["engine"]["reconstruction_coordinate"] == "0.592746050790", "p_ref drift")
    _require([int(item["N"]) for item in manifest["designs"]] == [305, 325], "design order drift")
    aliases = [Fraction(item["alias_ratio"]) for item in manifest["designs"]]
    _require(aliases[0] < 0 < aliases[1], "opposite alias signs lost")
    return manifest


def validate_metadata(path: Path, manifest: Mapping[str, Any], mode: str) -> Dict[str, Any]:
    metadata = json.loads(path.read_text(encoding="utf-8"))
    _require(
        metadata.get("engine") == "general integer-period threshold-rank Newman-Ziff",
        "unexpected engine",
    )
    designs = metadata.get("designs")
    _require(isinstance(designs, list) and len(designs) == 1, "one design per run required")
    observed = designs[0]
    n = int(observed.get("N", -1))
    expected_by_n = {int(item["N"]): item for item in manifest["designs"]}
    _require(n in expected_by_n, "run N is not frozen")
    expected = expected_by_n[n]
    for field in ("first", "second", "first_period_matrix", "second_period_matrix"):
        _require(observed.get(field) == expected[field], "metadata {} drift at N{}".format(field, n))
    commit = str(metadata.get("git_commit", ""))
    _require(COMMIT_RE.fullmatch(commit) is not None, "source commit must be full 40-hex")
    try:
        samples = int(metadata["samples_per_pair"])
        batches = int(metadata["batches"])
        seed = int(metadata["seed"])
        first = int(metadata["replica_counter_first"])
        last = int(metadata["replica_counter_last_exclusive"])
        command = str(metadata["command"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("invalid run metadata") from error
    _require(samples > 0 and batches >= 2 and samples % batches == 0, "invalid batch partition")
    _require(last - first == samples, "counter interval/sample mismatch")
    expected_seed = int(manifest["rng"][mode]["effective_seed_by_N"][str(n)])
    _require(seed == expected_seed, "RNG domain mismatch")
    expected_options = {
        "--samples": samples,
        "--batches": batches,
        "--seed": seed,
        "--replica-offset": first,
    }
    for option, value in expected_options.items():
        _require(int(command_option(command, option)[0]) == value, "command {} mismatch".format(option))
    _require(command_option(command, "--git-commit")[0].lower() == commit.lower(), "command commit mismatch")
    for option, matrix in (
        ("--first-matrix", expected["first_period_matrix"]),
        ("--second-matrix", expected["second_period_matrix"]),
    ):
        flattened = tuple(str(value) for row in matrix for value in row)
        _require(command_option(command, option, 4) == flattened, "command {} drift".format(option))
    for option, rep in (("--first-rep", expected["first"]), ("--second-rep", expected["second"])):
        _require(command_option(command, option, 2) == tuple(map(str, rep)), "command {} drift".format(option))
    return {
        "N": n,
        "first": tuple(expected["first"]),
        "second": tuple(expected["second"]),
        "samples": samples,
        "batches": batches,
        "seed": seed,
        "counter_first": first,
        "counter_last": last,
        "commit": commit.lower(),
        "metadata_sha256": sha256(path),
    }


def read_histograms(path: Path, run: Mapping[str, Any]) -> Dict[Tuple[str, int], Dict[str, Any]]:
    required = {"n", "a", "b", "orientation", "batch", "samples", "kind", "k", "count"}
    records: Dict[Tuple[str, int], Dict[str, Any]] = {}
    reps = {"first": run["first"], "second": run["second"]}
    n = int(run["N"])
    per_batch = int(run["samples"]) // int(run["batches"])
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        _require(not (required - set(reader.fieldnames or ())), "histogram columns missing")
        for row in reader:
            orientation = row["orientation"]
            batch = int(row["batch"])
            _require(int(row["n"]) == n and orientation in reps, "histogram design mismatch")
            _require((int(row["a"]), int(row["b"])) == tuple(reps[orientation]), "histogram rep mismatch")
            _require(0 <= batch < int(run["batches"]), "histogram batch out of range")
            _require(int(row["samples"]) == per_batch, "histogram batch sample mismatch")
            kind, rank, count = row["kind"], int(row["k"]), int(row["count"])
            _require(kind in ("minus", "plus") and 1 <= rank <= n and count > 0, "bad histogram row")
            record = records.setdefault(
                (orientation, batch),
                {"samples": per_batch, "minus": [0] * (n + 1), "plus": [0] * (n + 1)},
            )
            record[kind][rank] += count
    expected = {(orientation, batch) for orientation in reps for batch in range(int(run["batches"]))}
    _require(set(records) == expected, "histogram batch grid incomplete")
    for record in records.values():
        _require(sum(record["minus"]) == per_batch, "minus histogram total mismatch")
        _require(sum(record["plus"]) == per_batch, "plus histogram total mismatch")
    return records


def validate_moments(
    path: Path, run: Mapping[str, Any], records: Mapping[Tuple[str, int], Mapping[str, Any]]
) -> None:
    required = {
        "n", "a", "b", "orientation", "batch", "samples", "sum_kminus", "sum_kplus",
        "sum_kminus2", "sum_kplus2", "sum_product", "sum_gap", "sum_gap2",
    }
    seen = set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        _require(not (required - set(reader.fieldnames or ())), "moment columns missing")
        for row in reader:
            key = (row["orientation"], int(row["batch"]))
            _require(key in records and key not in seen and int(row["n"]) == run["N"], "moment key mismatch")
            seen.add(key)
            record = records[key]
            minus, plus = record["minus"], record["plus"]
            calculated = {
                "sum_kminus": sum(k * count for k, count in enumerate(minus)),
                "sum_kplus": sum(k * count for k, count in enumerate(plus)),
                "sum_kminus2": sum(k * k * count for k, count in enumerate(minus)),
                "sum_kplus2": sum(k * k * count for k, count in enumerate(plus)),
            }
            _require(int(row["samples"]) == record["samples"], "moment samples mismatch")
            _require(all(int(row[name]) == value for name, value in calculated.items()), "moment marginal mismatch")
            _require(
                int(row["sum_gap"]) == calculated["sum_kplus"] - calculated["sum_kminus"],
                "moment gap identity failed",
            )
    _require(seen == set(records), "moment grid incomplete")


def tail(histogram: Sequence[int], samples: int, p: float) -> float:
    n = len(histogram) - 1
    q = 1.0 - p
    probability = q ** n
    cumulative = 0
    total = 0.0
    for occupied in range(n + 1):
        if occupied:
            cumulative += histogram[occupied]
        total += cumulative * probability
        if occupied < n:
            probability *= (n - occupied) * p / ((occupied + 1) * q)
    return total / samples


def delta_m_batches(records: Mapping[Tuple[str, int], Mapping[str, Any]], batches: int, p: float) -> List[float]:
    values = []
    for batch in range(batches):
        sector: Dict[str, float] = {}
        for orientation in ("first", "second"):
            record = records[(orientation, batch)]
            sector[orientation] = (
                tail(record["plus"], record["samples"], p)
                + tail(record["minus"], record["samples"], p) - 1.0
            )
        values.append(sector["first"] - sector["second"])
    return values


def mean_se(values: Sequence[float]) -> Tuple[float, float]:
    _require(len(values) >= 2, "at least two batches required")
    mean = math.fsum(values) / len(values)
    centered = math.fsum((value - mean) ** 2 for value in values)
    return mean, math.sqrt(centered / (len(values) * (len(values) - 1)))


def two_sided_power(effect: float, se: float, alpha: float = 0.05) -> float:
    if se <= 0:
        return 1.0
    zcrit = NORMAL.inv_cdf(1.0 - alpha / 2.0)
    shift = abs(effect) / se
    return (1.0 - NORMAL.cdf(zcrit - shift)) + NORMAL.cdf(-zcrit - shift)


def noncentral_chi2_df2_power(noncentrality: float, alpha: float = 0.05) -> float:
    """Power of a df=2 chi-square score via its Poisson mixture."""
    _require(noncentrality >= 0.0, "negative noncentrality")
    critical = -2.0 * math.log(alpha)
    poisson_mean = noncentrality / 2.0
    weight = math.exp(-poisson_mean)
    power = 0.0
    total_weight = 0.0
    half_critical = critical / 2.0
    central_term = math.exp(-half_critical)
    central_sum = central_term
    for k in range(10000):
        if k:
            central_term *= half_critical / k
            central_sum += central_term
        power += weight * central_sum
        total_weight += weight
        if total_weight > 1.0 - 1e-14 and k > poisson_mean:
            break
        weight *= poisson_mean / (k + 1) if poisson_mean else 0.0
    return min(1.0, max(0.0, power))


def quadratic_2(vector: Sequence[float], covariance: Sequence[Sequence[float]]) -> float:
    a, c, b = covariance[0][0], covariance[0][1], covariance[1][1]
    determinant = a * b - c * c
    _require(determinant > 0.0, "score covariance is not positive definite")
    x, y = vector
    return (b * x * x - 2.0 * c * x * y + a * y * y) / determinant


def linear_covariance(weights: Sequence[Sequence[float]], covariance: Sequence[Sequence[float]]) -> List[List[float]]:
    output = []
    for left in weights:
        row = []
        for right in weights:
            row.append(sum(left[i] * covariance[i][j] * right[j] for i in range(2) for j in range(2)))
        output.append(row)
    return output


def pilot_result(
    manifest: Mapping[str, Any], campaigns: Mapping[int, Mapping[str, Any]]
) -> Dict[str, Any]:
    smoke = manifest["variance_only_smoke"]
    n_smoke = int(smoke["samples_per_design"])
    frozen_payload = json.loads(
        (ROOT / manifest["frozen_design"]["path"]).read_text(encoding="utf-8")
    )
    frozen_by_n = {int(item["N"]): item for item in frozen_payload["designs"]}
    a4 = float(frozen_payload["frozen_model"]["source_amplitude"])
    estimates = {}
    delta_se = []
    frozen_means = []
    for design in manifest["designs"]:
        n = int(design["N"])
        runs = campaigns[n]["runs"]
        _require(len(runs) == 1, "smoke needs exactly one run per N")
        run = runs[0]
        _require(run["metadata"]["samples"] == n_smoke, "smoke sample count mismatch")
        _require(run["metadata"]["batches"] == int(smoke["batches"]), "smoke batches mismatch")
        _require(run["metadata"]["counter_first"] == int(smoke["counter_first"]), "smoke counter start mismatch")
        _, se = mean_se(run["delta_m"])
        delta_se.append(se)
        frozen_mean = float(frozen_by_n[n]["h4_only_target_mean"])
        frozen_means.append(frozen_mean)
        estimates[str(n)] = {
            "samples": n_smoke,
            "batches": int(smoke["batches"]),
            "centered_sampling_se": se,
            "variance_constant_se_times_sqrt_samples": se * math.sqrt(n_smoke),
            "observed_target_mean": "withheld_by_protocol",
        }
    designs = manifest["designs"]
    aliases = [float(Fraction(item["alias_ratio"])) for item in designs]
    scales = [
        int(item["N"]) ** (13.0 / 8.0) / float(Fraction(item["delta_cos4"]))
        for item in designs
    ]
    denominator = aliases[1] - aliases[0]
    a12_se_smoke = math.sqrt(sum((scales[i] * delta_se[i]) ** 2 for i in range(2))) / abs(denominator)
    h12_shifts = {
        str(int(item["N"])): a4 * float(Fraction(item["delta_cos12"]))
        * int(item["N"]) ** (-13.0 / 8.0)
        for item in designs
    }
    _require(h12_shifts["305"] < 0.0 < h12_shifts["325"], "H12 alternative lost sign flip")
    grids = [int(value) for value in smoke["sample_count_freeze_rule"]["operational_grid_per_design"]]
    projections = []
    chosen = grids[-1]
    found = False
    for samples in grids:
        projected_se = [se * math.sqrt(n_smoke / samples) for se in delta_se]
        snr = [frozen_means[i] / projected_se[i] for i in range(2)]
        projected_a12_se = a12_se_smoke * math.sqrt(n_smoke / samples)
        h4_h12_distance = abs(a4) / projected_a12_se
        zero_noncentrality = sum((frozen_means[i] / projected_se[i]) ** 2 for i in range(2))
        projections.append({
            "samples_per_design": samples,
            "projected_sampling_se_by_N": {"305": projected_se[0], "325": projected_se[1]},
            "projected_h4_signal_to_noise_by_N": {"305": snr[0], "325": snr[1]},
            "conditional_h4_only_vs_equal_A12_mahalanobis": h4_h12_distance,
            "projected_A12_sampling_se": projected_a12_se,
            "zero_effect_h4_vector_mahalanobis": math.sqrt(zero_noncentrality),
            "zero_effect_df2_alpha_0.05_power": noncentral_chi2_df2_power(zero_noncentrality),
        })
        if not found and h4_h12_distance >= 3.0:
            chosen, found = samples, True
    h12_power = []
    for samples in grids:
        projected = a12_se_smoke * math.sqrt(n_smoke / samples)
        h12_power.append({
            "samples_per_design": samples,
            "projected_A12_sampling_se": projected,
            "two_sided_alpha_0.05_power_by_abs_A12_over_A4": {
                str(ratio): two_sided_power(ratio * a4, projected)
                for ratio in (0.10, 0.25, 0.50, 1.00)
            },
        })
    return {
        "schema": "matching-one/issue55-h4-h12-variance-pilot/v1",
        "status": "variance_only_target_means_withheld",
        "p_ref": manifest["engine"]["reconstruction_coordinate"],
        "variance_estimates": estimates,
        "production_grid_projection": projections,
        "sample_count_freeze_rule": smoke["sample_count_freeze_rule"],
        "recommended_samples_per_design": chosen,
        "grid_threshold_reached": found,
        "equal_amplitude_h12_sign_flip_shift_by_N": h12_shifts,
        "orthogonal_A12_power_sensitivity": h12_power,
        "scientific_boundary": "The pilot chooses cost from centered variance only; it contains no target estimate, target residual, or H4/H12 score.",
    }


def final_result(
    manifest: Mapping[str, Any], campaigns: Mapping[int, Mapping[str, Any]]
) -> Dict[str, Any]:
    production = manifest["production"]
    total = production.get("samples_per_design")
    _require(isinstance(total, int) and total > 0, "production sample count is not frozen")
    shard_count = int(production["shard_count"])
    expected_per_shard = total // shard_count
    observed, sampling_se, run_summary = [], [], {}
    commits = set()
    for n in (305, 325):
        runs = campaigns[n]["runs"]
        _require(len(runs) == shard_count, "all three shards required at N{}".format(n))
        runs = sorted(runs, key=lambda item: item["metadata"]["counter_first"])
        expected_first = int(production["counter_first"])
        values: List[float] = []
        run_summary[str(n)] = []
        for run in runs:
            meta = run["metadata"]
            _require(meta["samples"] == expected_per_shard, "production shard sample mismatch")
            _require(meta["batches"] == int(production["batches_per_shard"]), "production batches mismatch")
            _require(meta["counter_first"] == expected_first, "production counter gap/overlap")
            expected_first = meta["counter_last"]
            values.extend(run["delta_m"])
            commits.add(meta["commit"])
            run_summary[str(n)].append(meta)
        _require(expected_first == int(production["counter_last_exclusive"]), "production counter end mismatch")
        mean, se = mean_se(values)
        observed.append(mean)
        sampling_se.append(se)
    _require(len(commits) == 1, "all production shards must use one source commit")
    frozen = json.loads((ROOT / manifest["frozen_design"]["path"]).read_text(encoding="utf-8"))
    frozen_by_n = {int(item["N"]): item for item in frozen["designs"]}
    target = [float(frozen_by_n[n]["h4_only_target_mean"]) for n in (305, 325)]
    source_se = [float(frozen_by_n[n]["source_coefficient_only_se"]) for n in (305, 325)]
    sampling_cov = [[sampling_se[i] ** 2 if i == j else 0.0 for j in range(2)] for i in range(2)]
    target_cov = [[sampling_cov[i][j] + source_se[i] * source_se[j] for j in range(2)] for i in range(2)]
    residual = [observed[i] - target[i] for i in range(2)]
    target_chi2 = quadratic_2(residual, target_cov)
    zero_chi2 = quadratic_2(observed, sampling_cov)

    rows = manifest["designs"]
    aliases = [float(Fraction(row["alias_ratio"])) for row in rows]
    scales = [int(row["N"]) ** (13.0 / 8.0) / float(Fraction(row["delta_cos4"])) for row in rows]
    denominator = aliases[1] - aliases[0]
    weights = [
        [aliases[1] * scales[0] / denominator, -aliases[0] * scales[1] / denominator],
        [-scales[0] / denominator, scales[1] / denominator],
    ]
    amplitudes = [sum(weights[row][i] * observed[i] for i in range(2)) for row in range(2)]
    amplitude_cov = linear_covariance(weights, sampling_cov)
    amplitude_se = [math.sqrt(amplitude_cov[i][i]) for i in range(2)]
    correlation = amplitude_cov[0][1] / (amplitude_se[0] * amplitude_se[1])
    a12_z = amplitudes[1] / amplitude_se[1]
    return {
        "schema": "matching-one/issue55-h4-h12-orthogonal-score/v1",
        "status": "frozen score; no refit and no heldout-third-row claim",
        "p_ref": manifest["engine"]["reconstruction_coordinate"],
        "source_commit": next(iter(commits)),
        "runs": run_summary,
        "observed_DeltaM": {"305": observed[0], "325": observed[1]},
        "sampling_se": {"305": sampling_se[0], "325": sampling_se[1]},
        "sampling_covariance_N305_N325": sampling_cov,
        "frozen_h4_only": {
            "mean": {"305": target[0], "325": target[1]},
            "source_amplitude_se": {"305": source_se[0], "325": source_se[1]},
            "source_error_correlation": [[1.0, 1.0], [1.0, 1.0]],
            "residual": {"305": residual[0], "325": residual[1]},
            "total_covariance": target_cov,
            "chi_square": target_chi2,
            "df": 2,
            "p_value": math.exp(-target_chi2 / 2.0),
        },
        "zero_effect": {
            "chi_square": zero_chi2,
            "df": 2,
            "p_value": math.exp(-zero_chi2 / 2.0),
        },
        "two_column_h4_h12": {
            "normalized_rows": {
                "305": scales[0] * observed[0],
                "325": scales[1] * observed[1],
            },
            "alias_ratios": {"305": aliases[0], "325": aliases[1]},
            "A4": amplitudes[0],
            "A12": amplitudes[1],
            "covariance_A4_A12": amplitude_cov,
            "sampling_se_A4": amplitude_se[0],
            "sampling_se_A12": amplitude_se[1],
            "correlation_A4_A12": correlation,
            "A12_z": a12_z,
            "A12_two_sided_p": math.erfc(abs(a12_z) / math.sqrt(2.0)),
            "fit_df": 0,
            "fit_warning": "Two rows and two columns are saturated; A12 is an exact opposite-alias contrast, not an omnibus goodness-of-fit test.",
        },
        "heldout_third_alias": "not_run; allowed only if A12 is resolved",
    }


def load_campaigns(
    run_args: Iterable[Sequence[str]], manifest: Mapping[str, Any], mode: str
) -> Dict[int, MutableMapping[str, Any]]:
    output: Dict[int, MutableMapping[str, Any]] = {305: {"runs": []}, 325: {"runs": []}}
    for hist_name, moments_name, metadata_name in run_args:
        meta = validate_metadata(Path(metadata_name), manifest, mode)
        records = read_histograms(Path(hist_name), meta)
        validate_moments(Path(moments_name), meta, records)
        values = delta_m_batches(
            records, int(meta["batches"]), float(manifest["engine"]["reconstruction_coordinate"])
        )
        output[int(meta["N"])]["runs"].append({"metadata": meta, "delta_m": values})
    _require(all(output[n]["runs"] for n in output), "both N305 and N325 runs are required")
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("pilot", "final"), required=True)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument(
        "--run", action="append", nargs=3, required=True,
        metavar=("HIST", "MOMENTS", "METADATA"),
    )
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()
    manifest = load_manifest(args.manifest)
    mode = "smoke" if args.mode == "pilot" else "production"
    campaigns = load_campaigns(args.run, manifest, mode)
    result = pilot_result(manifest, campaigns) if args.mode == "pilot" else final_result(manifest, campaigns)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
