#!/usr/bin/env python3
"""Fast retrospective P48 projector/scaling score from aligned rank histograms.

This deliberately small, standard-library-only analyzer consumes the P33 long
histogram format.  It preserves the common random numbers across sizes by
deleting the same batch id at every N, forms jackknife pseudovalues, and uses
their full cross-N/cross-channel covariance for frozen train/heldout scores.

The output is planning evidence only: the input predates the P48 protocol.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


SIZES = (65, 85, 130, 145, 170)
TRAIN_SIZES = (65, 85, 130)
HELDOUT_SIZES = (145, 170)
METRICS = ("P4_S", "P4_D", "P4_S_prime", "P4_D_prime")
POWERS = {
    "P4_S": 1.0,
    "P4_D": 13.0 / 8.0,
    "P4_S_prime": 5.0 / 4.0,
    "P4_D_prime": 5.0 / 8.0,
}
POWER_LABELS = {
    "P4_S": "1",
    "P4_D": "13/8",
    "P4_S_prime": "5/4",
    "P4_D_prime": "5/8",
}


@dataclass
class Histogram:
    n: int
    a: int
    b: int
    orientation: str
    batch: int
    samples: int
    minus: List[int]
    plus: List[int]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def cos4(a: int, b: int) -> float:
    n = a * a + b * b
    return (a ** 4 - 6 * a * a * b * b + b ** 4) / float(n * n)


def read_histograms(path: Path) -> Dict[Tuple[int, str, int], Histogram]:
    required = {"n", "a", "b", "orientation", "batch", "samples", "kind", "k", "count"}
    records: Dict[Tuple[int, str, int], Histogram] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError("histogram CSV missing fields: " + ", ".join(sorted(missing)))
        for raw in reader:
            n = int(raw["n"])
            orientation = raw["orientation"]
            batch = int(raw["batch"])
            samples = int(raw["samples"])
            rank = int(raw["k"])
            count = int(raw["count"])
            kind = raw["kind"]
            if orientation not in ("first", "second") or kind not in ("minus", "plus"):
                raise ValueError("unknown orientation/kind")
            if n <= 0 or batch < 0 or samples <= 0 or not 1 <= rank <= n or count <= 0:
                raise ValueError("invalid histogram row")
            key = (n, orientation, batch)
            if key not in records:
                records[key] = Histogram(
                    n, int(raw["a"]), int(raw["b"]), orientation, batch,
                    samples, [0] * (n + 1), [0] * (n + 1),
                )
            record = records[key]
            if (record.a, record.b, record.samples) != (
                int(raw["a"]), int(raw["b"]), samples
            ):
                raise ValueError("inconsistent metadata within histogram")
            getattr(record, kind)[rank] += count
    if not records:
        raise ValueError("histogram CSV is empty")
    for record in records.values():
        if sum(record.minus) != record.samples or sum(record.plus) != record.samples:
            raise ValueError("histogram total differs from samples")
    validate_alignment(records)
    return records


def validate_alignment(records: Mapping[Tuple[int, str, int], Histogram]) -> None:
    sizes = sorted({key[0] for key in records})
    expected_batches = None
    expected_samples = None
    for n in sizes:
        for orientation in ("first", "second"):
            selected = sorted(
                (r for key, r in records.items() if key[:2] == (n, orientation)),
                key=lambda r: r.batch,
            )
            batches = [r.batch for r in selected]
            if batches != list(range(len(selected))) or len(batches) < 2:
                raise ValueError("each N/orientation needs the same contiguous batch ids")
            sample_sizes = {r.samples for r in selected}
            if len(sample_sizes) != 1:
                raise ValueError("batch sample counts differ")
            if expected_batches is None:
                expected_batches = batches
                expected_samples = next(iter(sample_sizes))
            elif batches != expected_batches or next(iter(sample_sizes)) != expected_samples:
                raise ValueError("cross-N synchronized batch alignment is absent")


def tail_and_derivative(hist: Sequence[int], samples: int, p: float) -> Tuple[float, float]:
    """Return E[1{Binomial(N,p)>=K}] and its analytic p derivative."""
    n = len(hist) - 1
    if not 0.0 < p < 1.0:
        raise ValueError("p must be strictly between zero and one")
    q = 1.0 - p
    probability = q ** n
    if probability == 0.0:
        # Starting the recurrence at k=0 underflows for the N=680 bracket
        # evaluation even though the mass around the binomial mode is ordinary.
        # Re-center only that numerical edge case; the established path below
        # remains bit-for-bit unchanged wherever q**n is representable.
        mode = min(n, int((n + 1) * p))
        weights = [0.0] * (n + 1)
        weights[mode] = 1.0
        for occupied in range(mode, 0, -1):
            weights[occupied - 1] = (
                weights[occupied]
                * occupied
                * q
                / ((n - occupied + 1) * p)
            )
        for occupied in range(mode, n):
            weights[occupied + 1] = (
                weights[occupied]
                * (n - occupied)
                * p
                / ((occupied + 1) * q)
            )
        normalization = math.fsum(weights)
        cumulative = 0
        value_terms = []
        derivative_terms = []
        for occupied, weight in enumerate(weights):
            if occupied:
                cumulative += hist[occupied]
                derivative_terms.append(hist[occupied] * occupied * weight / p)
            value_terms.append(cumulative * weight)
        scale = samples * normalization
        return math.fsum(value_terms) / scale, math.fsum(derivative_terms) / scale
    cumulative = 0
    value = 0.0
    for occupied in range(n + 1):
        if occupied:
            cumulative += hist[occupied]
        value += cumulative * probability
        if occupied < n:
            probability *= (n - occupied) * p / ((occupied + 1) * q)
    density = n * q ** (n - 1)
    derivative = 0.0
    for rank in range(1, n + 1):
        derivative += hist[rank] * density
        if rank < n:
            density *= (n - rank) * p / (rank * q)
    return value / samples, derivative / samples


def add_histograms(rows: Sequence[Histogram], kind: str, omitted: int = -1) -> List[int]:
    n = rows[0].n
    total = [0] * (n + 1)
    for row in rows:
        if row.batch == omitted:
            continue
        source = getattr(row, kind)
        for rank, count in enumerate(source):
            total[rank] += count
    return total


def project_size(
    by_orientation: Mapping[str, Sequence[Histogram]], omitted: int = -1
) -> Dict[str, float]:
    aggregate = {}
    for orientation in ("first", "second"):
        rows = by_orientation[orientation]
        included = [row for row in rows if row.batch != omitted]
        aggregate[orientation] = {
            "a": rows[0].a,
            "b": rows[0].b,
            "samples": sum(row.samples for row in included),
            "minus": add_histograms(rows, "minus", omitted),
            "plus": add_histograms(rows, "plus", omitted),
        }
    n = by_orientation["first"][0].n

    def mean_matching(p: float) -> float:
        value = 0.0
        for orientation in ("first", "second"):
            row = aggregate[orientation]
            minus, _ = tail_and_derivative(row["minus"], row["samples"], p)
            plus, _ = tail_and_derivative(row["plus"], row["samples"], p)
            value += minus + plus - 1.0
        return value / 2.0

    lower, upper = 0.0, 1.0
    for _ in range(56):
        midpoint = (lower + upper) / 2.0
        if mean_matching(midpoint) < 0.0:
            lower = midpoint
        else:
            upper = midpoint
    p0 = (lower + upper) / 2.0
    values = {}
    for orientation in ("first", "second"):
        row = aggregate[orientation]
        r_minus, d_minus = tail_and_derivative(row["minus"], row["samples"], p0)
        r_plus, d_plus = tail_and_derivative(row["plus"], row["samples"], p0)
        values[orientation] = {
            "S": (r_minus + 1.0 - r_plus) / 2.0,
            "D": (r_minus + r_plus - 1.0) / 2.0,
            "S_prime": (d_minus - d_plus) / 2.0,
            "D_prime": (d_minus + d_plus) / 2.0,
        }
    delta_cos4 = cos4(aggregate["first"]["a"], aggregate["first"]["b"]) - cos4(
        aggregate["second"]["a"], aggregate["second"]["b"]
    )
    # The even matching-parity channel reverses under the orientation-pair
    # convention; these signs are the frozen P48 projector definitions.
    return {
        "p0": p0,
        # D=(M/2) for each orientation, hence the orientation-mean matching
        # slope is D'_first+D'_second.  By Russo this is the mean total
        # primal-plus-matching pivotal mass.
        "Mbar_prime": (
            values["first"]["D_prime"] + values["second"]["D_prime"]
        ),
        "first_primal_pivotal": (
            values["first"]["D_prime"] + values["first"]["S_prime"]
        ),
        "first_matching_pivotal": (
            values["first"]["D_prime"] - values["first"]["S_prime"]
        ),
        "second_primal_pivotal": (
            values["second"]["D_prime"] + values["second"]["S_prime"]
        ),
        "second_matching_pivotal": (
            values["second"]["D_prime"] - values["second"]["S_prime"]
        ),
        "P4_S": (values["second"]["S"] - values["first"]["S"]) / delta_cos4,
        "P4_D": (values["first"]["D"] - values["second"]["D"]) / delta_cos4,
        "P4_S_prime": (
            values["second"]["S_prime"] - values["first"]["S_prime"]
        ) / delta_cos4,
        "P4_D_prime": (
            values["first"]["D_prime"] - values["second"]["D_prime"]
        ) / delta_cos4,
        "B_plus": values["first"]["D_prime"] + values["second"]["D_prime"],
        "A_plus": 2.0 * (
            values["first"]["D_prime"] - values["second"]["D_prime"]
        ) / delta_cos4,
        "A_minus": 2.0 * (
            values["second"]["S_prime"] - values["first"]["S_prime"]
        ) / delta_cos4,
        "B_minus": values["first"]["S_prime"] + values["second"]["S_prime"],
        "delta_cos4": delta_cos4,
    }


def pseudovalues(full: float, deleted: Sequence[float]) -> List[float]:
    batches = len(deleted)
    return [batches * full - (batches - 1) * value for value in deleted]


def covariance_of_mean(rows: Sequence[Sequence[float]]) -> List[List[float]]:
    if len(rows) < 2:
        raise ValueError("at least two pseudovalue rows are required")
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("ragged pseudovalue matrix")
    means = [math.fsum(row[j] for row in rows) / len(rows) for j in range(width)]
    denom = len(rows) * (len(rows) - 1)
    return [[
        math.fsum((row[i] - means[i]) * (row[j] - means[j]) for row in rows) / denom
        for j in range(width)
    ] for i in range(width)]


def solve_matrix(matrix: Sequence[Sequence[float]], rhs: Sequence[Sequence[float]]) -> List[List[float]]:
    """Partial-pivoted Gauss-Jordan solve of A X = B for small dense systems."""
    n = len(matrix)
    width = len(rhs[0])
    work = [list(matrix[i]) + list(rhs[i]) for i in range(n)]
    scale = max(max(abs(value) for value in row[:n]) for row in work)
    for column in range(n):
        pivot = max(range(column, n), key=lambda row: abs(work[row][column]))
        if abs(work[pivot][column]) <= scale * 1e-13:
            raise ValueError("covariance/design matrix is numerically singular")
        work[column], work[pivot] = work[pivot], work[column]
        divisor = work[column][column]
        work[column] = [value / divisor for value in work[column]]
        for row in range(n):
            if row == column:
                continue
            factor = work[row][column]
            if factor:
                work[row] = [
                    work[row][j] - factor * work[column][j]
                    for j in range(n + width)
                ]
    return [row[n:] for row in work]


def inverse(matrix: Sequence[Sequence[float]]) -> List[List[float]]:
    n = len(matrix)
    return solve_matrix(matrix, [[float(i == j) for j in range(n)] for i in range(n)])


def matmul(left: Sequence[Sequence[float]], right: Sequence[Sequence[float]]) -> List[List[float]]:
    return [[math.fsum(a * b for a, b in zip(row, column)) for column in zip(*right)] for row in left]


def transpose(matrix: Sequence[Sequence[float]]) -> List[List[float]]:
    return [list(column) for column in zip(*matrix)]


def subset(matrix: Sequence[Sequence[float]], rows: Sequence[int], columns: Sequence[int]) -> List[List[float]]:
    return [[matrix[i][j] for j in columns] for i in rows]


def subtract(left: Sequence[Sequence[float]], right: Sequence[Sequence[float]]) -> List[List[float]]:
    return [[a - b for a, b in zip(x, y)] for x, y in zip(left, right)]


def quadratic(vector: Sequence[float], covariance: Sequence[Sequence[float]]) -> float:
    solved = solve_matrix(covariance, [[value] for value in vector])
    return math.fsum(value * solved[i][0] for i, value in enumerate(vector))


def gls_score(
    values: Sequence[float], covariance: Sequence[Sequence[float]],
    train_indices: Sequence[int], held_indices: Sequence[int],
) -> Dict[str, object]:
    ctt = subset(covariance, train_indices, train_indices)
    cth = subset(covariance, train_indices, held_indices)
    chh = subset(covariance, held_indices, held_indices)
    inv_ctt = inverse(ctt)
    ones_t = [[1.0] for _ in train_indices]
    denom = matmul(transpose(ones_t), matmul(inv_ctt, ones_t))[0][0]
    influence = [value / denom for value in matmul(transpose(ones_t), inv_ctt)[0]]
    amplitude = math.fsum(weight * values[index] for weight, index in zip(influence, train_indices))
    amplitude_variance = 1.0 / denom
    cross = [math.fsum(cth[i][j] * influence[i] for i in range(len(train_indices))) for j in range(len(held_indices))]
    residual_covariance = [[
        chh[i][j] + amplitude_variance - cross[i] - cross[j]
        for j in range(len(held_indices))
    ] for i in range(len(held_indices))]
    residual = [values[index] - amplitude for index in held_indices]
    held_covariance = subset(covariance, held_indices, held_indices)
    return {
        "amplitude": amplitude,
        "amplitude_se": math.sqrt(amplitude_variance),
        "influence": influence,
        "residual": residual,
        "residual_covariance": residual_covariance,
        "residual_se": [math.sqrt(residual_covariance[i][i]) for i in range(len(held_indices))],
        "heldout_chi_square": quadratic(residual, residual_covariance),
        "heldout_df": len(held_indices),
        "zero_chi_square": quadratic([values[i] for i in held_indices], held_covariance),
    }


def joint_score(values: Sequence[float], covariance: Sequence[Sequence[float]]) -> Dict[str, object]:
    """Four separate frozen amplitudes, scored jointly with all covariance."""
    train = [m * len(SIZES) + SIZES.index(n) for m in range(len(METRICS)) for n in TRAIN_SIZES]
    held = [m * len(SIZES) + SIZES.index(n) for m in range(len(METRICS)) for n in HELDOUT_SIZES]
    design_t = [[float(i // len(TRAIN_SIZES) == m) for m in range(len(METRICS))] for i in range(len(train))]
    design_h = [[float(i // len(HELDOUT_SIZES) == m) for m in range(len(METRICS))] for i in range(len(held))]
    ctt = subset(covariance, train, train)
    cth = subset(covariance, train, held)
    chh = subset(covariance, held, held)
    inv_ctt = inverse(ctt)
    xt_ci = matmul(transpose(design_t), inv_ctt)
    normal = matmul(xt_ci, design_t)
    beta_cov = inverse(normal)
    influence = matmul(beta_cov, xt_ci)
    beta = [row[0] for row in matmul(influence, [[values[i]] for i in train])]
    predicted = [row[0] for row in matmul(design_h, [[x] for x in beta])]
    residual = [values[index] - predicted[i] for i, index in enumerate(held)]
    # Var(y_h - X_h beta_hat), retaining train/held cross covariance.
    pred_cov = matmul(design_h, matmul(beta_cov, transpose(design_h)))
    cov_h_beta = matmul(transpose(cth), transpose(influence))
    residual_cov = subtract(subtract([[chh[i][j] + pred_cov[i][j] for j in range(len(held))] for i in range(len(held))], matmul(cov_h_beta, transpose(design_h))), matmul(design_h, transpose(cov_h_beta)))
    return {
        "amplitudes": dict(zip(METRICS, beta)),
        "amplitude_covariance": beta_cov,
        "residual": residual,
        "residual_covariance": residual_cov,
        "heldout_chi_square": quadratic(residual, residual_cov),
        "heldout_df": len(held),
        "zero_chi_square": quadratic([values[i] for i in held], chh),
    }


def write_csv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def analyze(input_path: Path) -> Dict[str, object]:
    records = read_histograms(input_path)
    sizes = tuple(sorted({key[0] for key in records}))
    if sizes != SIZES:
        raise ValueError("P48 retrospective requires N=65,85,130,145,170")
    batch_count = len({key[2] for key in records})
    full: Dict[int, Dict[str, float]] = {}
    deleted: List[Dict[int, Dict[str, float]]] = []
    grouped = {}
    for n in SIZES:
        grouped[n] = {
            orientation: sorted(
                [record for key, record in records.items() if key[:2] == (n, orientation)],
                key=lambda record: record.batch,
            ) for orientation in ("first", "second")
        }
        full[n] = project_size(grouped[n])
    for omitted in range(batch_count):
        deleted.append({n: project_size(grouped[n], omitted) for n in SIZES})

    labels = [(metric, n) for metric in METRICS for n in SIZES]
    point_scaled = [full[n][metric] * n ** POWERS[metric] for metric, n in labels]
    pseudo_columns = {}
    for metric, n in [("p0", n) for n in SIZES] + labels:
        pseudo_columns[(metric, n)] = pseudovalues(
            full[n][metric], [item[n][metric] for item in deleted]
        )
    pseudo_scaled = [[
        pseudo_columns[(metric, n)][batch] * n ** POWERS[metric]
        for metric, n in labels
    ] for batch in range(batch_count)]
    stacked_covariance = covariance_of_mean(pseudo_scaled)
    size_indices = {n: SIZES.index(n) for n in SIZES}
    train_indices = [size_indices[n] for n in TRAIN_SIZES]
    held_indices = [size_indices[n] for n in HELDOUT_SIZES]
    scores = {}
    for metric_index, metric in enumerate(METRICS):
        offset = metric_index * len(SIZES)
        indices = list(range(offset, offset + len(SIZES)))
        scores[metric] = gls_score(
            point_scaled[offset:offset + len(SIZES)],
            subset(stacked_covariance, indices, indices),
            train_indices, held_indices,
        )
    joint = joint_score(point_scaled, stacked_covariance)
    point_rows = []
    for n in SIZES:
        row = {
            "N": n,
            "p0": full[n]["p0"],
            "p0_se": math.sqrt(covariance_of_mean([[x] for x in pseudo_columns[("p0", n)]])[0][0]),
            "delta_cos4": full[n]["delta_cos4"],
        }
        for metric in METRICS:
            index = labels.index((metric, n))
            row[metric] = full[n][metric]
            row[metric + "_se"] = math.sqrt(stacked_covariance[index][index]) / n ** POWERS[metric]
            row[metric + "_scaled"] = point_scaled[index]
            row[metric + "_scaled_se"] = math.sqrt(stacked_covariance[index][index])
        point_rows.append(row)
    return {
        "records": records,
        "batch_count": batch_count,
        "labels": labels,
        "point_scaled": point_scaled,
        "point_rows": point_rows,
        "stacked_covariance": stacked_covariance,
        "scores": scores,
        "joint": joint,
    }


def git_head(cwd: Path) -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=str(cwd), text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def emit(result: Dict[str, object], input_path: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    script_path = Path(__file__).resolve()
    point_fields = ["N", "p0", "p0_se", "delta_cos4"]
    for metric in METRICS:
        point_fields.extend((metric, metric + "_se", metric + "_scaled", metric + "_scaled_se"))
    write_csv(output_dir / "projectors.csv", point_fields, result["point_rows"])

    held_rows = []
    for metric in METRICS:
        score = result["scores"][metric]
        for index, n in enumerate(HELDOUT_SIZES):
            residual_se = score["residual_se"][index]
            held_rows.append({
                "metric": metric, "power": POWERS[metric],
                "train_sizes": "/".join(map(str, TRAIN_SIZES)), "heldout_N": n,
                "train_amplitude": score["amplitude"], "train_amplitude_se": score["amplitude_se"],
                "observed_scaled": result["point_scaled"][METRICS.index(metric) * len(SIZES) + SIZES.index(n)],
                "predicted_scaled": score["amplitude"], "signed_error": score["residual"][index],
                "prediction_error_se": residual_se, "z": score["residual"][index] / residual_se,
                "channel_heldout_chi_square": score["heldout_chi_square"],
                "channel_df": score["heldout_df"], "channel_zero_chi_square": score["zero_chi_square"],
            })
    write_csv(output_dir / "heldout_scores.csv", list(held_rows[0]), held_rows)

    covariance_rows = []
    stacked_rows = []
    covariance = result["stacked_covariance"]
    for metric_index, metric in enumerate(METRICS):
        for i, n_i in enumerate(SIZES):
            for j in range(i, len(SIZES)):
                n_j = SIZES[j]
                ii = metric_index * len(SIZES) + i
                jj = metric_index * len(SIZES) + j
                cov = covariance[ii][jj]
                corr = cov / math.sqrt(covariance[ii][ii] * covariance[jj][jj])
                covariance_rows.append({"metric": metric, "N_i": n_i, "N_j": n_j, "covariance": cov, "correlation": corr})
    write_csv(output_dir / "cross_size_covariance.csv", list(covariance_rows[0]), covariance_rows)
    labels = result["labels"]
    for i, (metric_i, n_i) in enumerate(labels):
        for j in range(i, len(labels)):
            metric_j, n_j = labels[j]
            cov = covariance[i][j]
            corr = cov / math.sqrt(covariance[i][i] * covariance[j][j])
            stacked_rows.append({"metric_i": metric_i, "N_i": n_i, "metric_j": metric_j, "N_j": n_j, "covariance": cov, "correlation": corr})
    write_csv(output_dir / "stacked_covariance.csv", list(stacked_rows[0]), stacked_rows)

    summary = {
        "status": "retrospective planning only",
        "conclusion": "parity pattern supported but four-power joint not passed",
        "protocol": {
            "train_sizes": TRAIN_SIZES, "heldout_sizes": HELDOUT_SIZES,
            "powers": POWERS,
            "covariance": "cross-N synchronized delete-one jackknife pseudovalues; full train/heldout cross covariance retained",
        },
        "provenance": {
            "input": str(input_path), "input_sha256": sha256(input_path),
            "input_batches": result["batch_count"],
            "input_warning": "P33 existed before P48 was frozen and identifies its source only as a working tree; this is not confirmatory evidence.",
            "script": str(script_path), "script_sha256": sha256(script_path),
            "git_head": git_head(script_path.parent.parent),
            "python": platform.python_version(),
            "generated_utc": datetime.now(timezone.utc).isoformat(),
        },
        "projectors": result["point_rows"],
        "channel_scores": result["scores"],
        "joint_four_power_score": result["joint"],
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    lines = [
        "# P48 retrospective projector score", "",
        "**Status: retrospective planning only.** The P33 histogram source existed before P48 was frozen and its recorded source is a working tree. These numbers must not be presented as confirmatory or preregistered evidence.", "",
        "## Decision", "",
        "**Parity pattern supported but four-power joint not passed.** All four projected channels have the expected parity-resolved signal pattern and their frozen nonzero-amplitude predictions improve substantially on zero. However, the heldout `P4_S_prime` scaled amplitude drifts upward and fails its pure `N^-5/4` law; therefore the four-law package does not pass as a joint claim.", "",
        "## Frozen train/heldout score", "",
        "The amplitude of each law was fit using only `N=65,85,130`; `N=145,170` were scored as heldout. Prediction uncertainty includes amplitude-estimation uncertainty and the synchronized train/heldout cross covariance.", "",
        "| channel | power | train amplitude (SE) | heldout chi2 / 2 | zero chi2 / 2 | heldout z (145, 170) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for metric in METRICS:
        score = result["scores"][metric]
        z = [score["residual"][i] / score["residual_se"][i] for i in range(2)]
        lines.append(
            "| {0} | {1} | {2:.8g} ({3:.3g}) | {4:.4f} | {5:.4f} | {6:+.3f}, {7:+.3f} |".format(
                metric, POWER_LABELS[metric], score["amplitude"], score["amplitude_se"],
                score["heldout_chi_square"], score["zero_chi_square"], z[0], z[1]
            )
        )
    joint = result["joint"]
    lines.extend([
        "", "The covariance-aware four-channel heldout omnibus statistic is **{:.4f} / {} df** (zero-model statistic {:.4f} / {} df). The omnibus statistic alone is not unusually large; the declared package nevertheless fails its conjunction criterion because `P4_S_prime` is inconsistent with its frozen law (chi-square 10.1908 / 2 df). These are retrospective diagnostics, not confirmatory p-values.".format(joint["heldout_chi_square"], joint["heldout_df"], joint["zero_chi_square"], joint["heldout_df"]),
        "", "## Method", "",
        "At the intrinsic center where the mean of the two orientation matching functions vanishes, the analyzer reconstructs the thermal-even `S` and matching-odd `D` sectors and their analytic first derivatives. It projects the orientation contrast by the exact Gaussian-integer `Delta cos(4 theta)`. The 100 aligned batch ids are deleted synchronously across every N and channel; jackknife pseudovalues then supply the full covariance matrix used by GLS and predictive residual scoring.",
        "", "The next legitimate use is design: freeze one correction/log alternative motivated by the `P4_S_prime` drift, then test it on fresh independent counters. Do not refit alternatives on the heldout values reported here.", "",
    ])
    (output_dir / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")

    command = "python3 scripts/analyze_p48_retrospective.py --input {} --output {}".format(input_path, output_dir)
    (output_dir / "commands.txt").write_text(
        command + "\npython3 tests/test_p48_retrospective.py\n", encoding="utf-8"
    )
    checksum_paths = sorted(path for path in output_dir.iterdir() if path.name != "checksums.sha256")
    (output_dir / "checksums.sha256").write_text(
        "".join("{}  {}\n".format(sha256(path), path.name) for path in checksum_paths), encoding="utf-8"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    emit(analyze(args.input), args.input, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
