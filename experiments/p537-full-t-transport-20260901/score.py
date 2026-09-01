#!/usr/bin/env python3
"""Score the complete N65 canonical-pair thermal transport from frozen tables."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.stats import binom


N = 65
BATCHES = 100
DELTA = 1152 / 845
P_REF = 0.592746050790
N25 = 25
FIELDS = (
    "count", "sum_q0", "sum_E0", "sum_a16_0", "sum_q0_a16_0",
    "sum_E0_a16_0", "sum_q1", "sum_E1", "sum_a16_1",
    "sum_q1_a16_1", "sum_E1_a16_1",
)
OUTPUT_FIELDS = (
    "p", "T_t", "T_t_over_M_t", "J", "M_t", "R", "R_t",
    "jM", "jM_t", "jY", "jY_t", "J_retained_63", "J_omitted_nn_fill",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_baseline(path: Path) -> tuple[np.ndarray, np.ndarray]:
    raw = np.zeros((BATCHES, 2, 2, N + 1), dtype=np.float64)
    samples = np.zeros(BATCHES, dtype=np.float64)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if int(row["n"]) != N:
                continue
            batch = int(row["batch"])
            geometry = ("first", "second").index(row["orientation"])
            kind = ("minus", "plus").index(row["kind"])
            raw[batch, geometry, kind, int(row["k"])] += int(row["count"])
            samples[batch] = int(row["samples"])
    if np.any(samples <= 0):
        raise ValueError("baseline does not contain all 100 N65 batches")
    return raw, samples


def shard_metadata(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    with path.open() as handle:
        for line in handle:
            if not line.startswith("# "):
                break
            key, value = line[2:].rstrip("\n").split("=", 1)
            result[key] = value
    return result


def load_source(paths: list[Path]) -> tuple[np.ndarray, list[tuple[int, int]], dict]:
    if len(paths) != 4:
        raise ValueError("the frozen contract requires exactly four shards")
    metadata = [shard_metadata(path) for path in paths]
    shards = int(metadata[0]["shard_count"])
    samples = int(metadata[0]["samples"])
    seed = metadata[0]["seed"]
    proposal_p = float(metadata[0]["proposal_p"])
    if shards != 4 or {int(item["shard_index"]) for item in metadata} != set(range(4)):
        raise ValueError("incomplete shard index set")
    if any(
        (int(item["samples"]), item["seed"], float(item["proposal_p"]))
        != (samples, seed, proposal_p)
        for item in metadata
    ):
        raise ValueError("mixed shard contracts")

    displacements: set[tuple[int, int]] = set()
    for path in paths:
        with path.open(newline="") as handle:
            rows = (line for line in handle if not line.startswith("#"))
            for row in csv.DictReader(rows, delimiter="\t"):
                if row["kind"] == "global":
                    displacements.add((int(row["dx"]), int(row["dy"])))
    ordered = sorted(displacements)
    if len(ordered) != N - 2 or (1, 0) in displacements:
        raise ValueError("unexpected common displacement catalogue")
    direction_index = {value: index for index, value in enumerate(ordered)}
    data = np.zeros((BATCHES, 2, len(ordered), N, len(FIELDS)), dtype=np.int64)
    for path in paths:
        with path.open(newline="") as handle:
            rows = (line for line in handle if not line.startswith("#"))
            for row in csv.DictReader(rows, delimiter="\t"):
                if row["kind"] != "global":
                    continue
                batch = int(row["batch"])
                geometry = ("axis", "tilted").index(row["geometry"])
                displacement = direction_index[(int(row["dx"]), int(row["dy"]))]
                k = int(row["k"])
                data[batch, geometry, displacement, k] += np.array(
                    [int(row[field]) for field in FIELDS], dtype=np.int64
                )
    batch_counts = data[..., 0].sum(axis=-1)
    if np.any(batch_counts <= 0):
        raise ValueError("empty source batch")
    if np.any(batch_counts != batch_counts[:, :1, :1]):
        raise ValueError("global sample counts differ by geometry or displacement")
    if int(batch_counts[:, 0, 0].sum()) != samples:
        raise ValueError("source rows do not reconstruct the frozen sample total")
    run = {"samples": samples, "shards": shards, "seed": seed, "proposal_p": proposal_p}
    return data, ordered, run


def interval_map(center: np.ndarray, se: np.ndarray) -> dict[str, dict[str, float]]:
    return {
        field: {
            "value": float(value),
            "se": float(error),
            "ci95": [float(value - 1.96 * error), float(value + 1.96 * error)],
        }
        for field, value, error in zip(OUTPUT_FIELDS, center, se)
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--tables", required=True, nargs=4, type=Path)
    parser.add_argument("--n25-result", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)

    baseline_batches, baseline_samples = load_baseline(args.baseline)
    source_batches, displacements, run = load_source(args.tables)
    baseline_total = baseline_batches.sum(axis=0)
    source_total = source_batches.sum(axis=0)
    baseline_k = np.arange(N + 1)
    rest_k = np.arange(N)
    nn_directions = ((-1, 0), (0, -1), (0, 1))
    if not set(nn_directions).issubset(displacements):
        raise ValueError("the three retained C4-equivalent NN directions are missing")
    nn_indices = [displacements.index(direction) for direction in nn_directions]

    def evaluate(
        fixed_p: float | None = None,
        baseline_omit: int | None = None,
        source_omit: int | None = None,
    ) -> np.ndarray:
        hist = baseline_total - (
            baseline_batches[baseline_omit] if baseline_omit is not None else 0
        )
        sample_count = baseline_samples.sum() - (
            baseline_samples[baseline_omit] if baseline_omit is not None else 0
        )
        cumulative = np.cumsum(hist, axis=-1)
        q_curve = (-sample_count + cumulative[:, 0] + cumulative[:, 1]) / sample_count
        e_curve = (sample_count - cumulative[:, 0] + cumulative[:, 1]) / sample_count

        def baseline_at(p: float) -> tuple[np.ndarray, ...]:
            weight = binom.pmf(baseline_k, N, p)
            score = baseline_k - N * p
            second_score = score * score - N * p * (1 - p)
            return (
                q_curve @ weight,
                e_curve @ weight,
                q_curve @ (weight * score),
                e_curve @ (weight * score),
                q_curve @ (weight * second_score),
                e_curve @ (weight * second_score),
            )

        p = fixed_p
        if p is None:
            p = brentq(lambda value: baseline_at(value)[0].mean(), 0.58, 0.61)
        q, e, q_t, e_t, q_tt, e_tt = baseline_at(p)
        M_t = q_t.mean()
        M_tt = q_tt.mean()
        Y_t = (e_t[0] - e_t[1]) / DELTA
        Y_tt = (e_tt[0] - e_tt[1]) / DELTA
        R = Y_t / M_t
        R_t = (Y_tt - R * M_tt) / M_t

        source = source_total - (source_batches[source_omit] if source_omit is not None else 0)
        denominator = source[0, 0, :, 0].sum()
        proposal_p = run["proposal_p"]
        rest_ratio = (p / proposal_p) ** rest_k * (
            (1 - p) / (1 - proposal_p)
        ) ** (N - 2 - rest_k)
        packets = np.zeros((2, len(displacements), 6), dtype=np.float64)
        for state, columns in enumerate(((3, 4, 5), (8, 9, 10))):
            state_weight = (1 - p, p)[state]
            score = rest_k + state - N * p
            weight = (1 - p) * state_weight * rest_ratio / denominator
            values = source[..., list(columns)] / (16 * N)
            packets[..., :3] += np.einsum("gdkc,k->gdc", values, weight)
            packets[..., 3:] += np.einsum("gdkc,k->gdc", values, weight * score)

        retained = packets.sum(axis=1)
        # The producer omits y=z=+e1.  After Bernoulli integration over z,
        # the four NN pair directions are C4-equivalent at every p.  The
        # average of the retained -e1,+e2,-e2 directions is therefore an
        # unbiased reconstruction of the omitted full-population column.
        omitted_nn = packets[:, nn_indices].mean(axis=1)

        def response(source_moments: np.ndarray) -> tuple[float, float, float]:
            a, qa, ea, a_t, qa_t, ea_t = source_moments.T
            jM_geometry = qa - q * a
            jM_t_geometry = qa_t - q_t * a - q * a_t
            jY_geometry = ea - e * a
            jY_t_geometry = ea_t - e_t * a - e * a_t
            jM = jM_geometry.mean()
            jM_t = jM_t_geometry.mean()
            jY = (jY_geometry[0] - jY_geometry[1]) / DELTA
            jY_t = (jY_t_geometry[0] - jY_t_geometry[1]) / DELTA
            T_t = jY_t - R * jM_t - R_t * jM
            J = (N ** (13 / 8) / 2) * T_t / M_t
            return T_t, T_t / M_t, J

        full = retained + omitted_nn
        T_t, T_t_over_M_t, J = response(full)
        retained_J = response(retained)[2]
        omitted_J = response(omitted_nn)[2]
        a, qa, ea, a_t, qa_t, ea_t = full.T
        jM_geometry = qa - q * a
        jM_t_geometry = qa_t - q_t * a - q * a_t
        jY_geometry = ea - e * a
        jY_t_geometry = ea_t - e_t * a - e * a_t
        return np.array(
            [
                p, T_t, T_t_over_M_t, J, M_t, R, R_t,
                jM_geometry.mean(), jM_t_geometry.mean(),
                (jY_geometry[0] - jY_geometry[1]) / DELTA,
                (jY_t_geometry[0] - jY_t_geometry[1]) / DELTA,
                retained_J, omitted_J,
            ],
            dtype=np.float64,
        )

    root = evaluate()
    reference = evaluate(P_REF)
    root_source_jk = np.array([evaluate(source_omit=batch) for batch in range(BATCHES)])
    root_baseline_jk = np.array([evaluate(baseline_omit=batch) for batch in range(BATCHES)])
    ref_source_jk = np.array(
        [evaluate(P_REF, source_omit=batch) for batch in range(BATCHES)]
    )
    ref_baseline_jk = np.array(
        [evaluate(P_REF, baseline_omit=batch) for batch in range(BATCHES)]
    )

    factor = (BATCHES - 1) / BATCHES

    def independent_se(source_values: np.ndarray, baseline_values: np.ndarray) -> np.ndarray:
        source_delta = source_values - source_values.mean(axis=0)
        baseline_delta = baseline_values - baseline_values.mean(axis=0)
        return np.sqrt(
            factor * np.square(source_delta).sum(axis=0)
            + factor * np.square(baseline_delta).sum(axis=0)
        )

    root_se = independent_se(root_source_jk, root_baseline_jk)
    reference_se = independent_se(ref_source_jk, ref_baseline_jk)
    reference_se[OUTPUT_FIELDS.index("p")] = 0.0
    paired_source = root_source_jk - ref_source_jk
    paired_baseline = root_baseline_jk - ref_baseline_jk
    paired = root - reference
    paired_se = independent_se(paired_source, paired_baseline)

    n25 = json.loads(args.n25_result.read_text())
    n25_J_over_A = float(
        n25["controls"]["all_mode_reproduced_full_J2_over_A"]["midpoint"]
    )
    n25_A = N25 ** (13 / 8) / 2
    n25_J = n25_A * n25_J_over_A
    ratio = root[OUTPUT_FIELDS.index("J")] / n25_J
    ratio_se = root_se[OUTPUT_FIELDS.index("J")] / abs(n25_J)
    exponent = -math.log(abs(ratio)) / math.log(N / N25)
    exponent_se = ratio_se / (abs(ratio) * math.log(N / N25))
    paired_J = paired[OUTPUT_FIELDS.index("J")]
    paired_J_se = paired_se[OUTPUT_FIELDS.index("J")]

    payload = {
        "schema": "matching-one/p537-full-t-transport/v1",
        "status": "COMPLETED_EXISTING_SUFFICIENT_STATISTICS",
        "definition": "complete canonical s2 thermal jet with C4 reconstruction of the one omitted NN column",
        "N65_pooled_root": interval_map(root, root_se),
        "N65_p_ref": interval_map(reference, reference_se),
        "paired_root_minus_p_ref": interval_map(paired, paired_se),
        "N25_to_N65": {
            "N25_J_over_A": n25_J_over_A,
            "N25_J": n25_J,
            "J65_over_J25": {
                "value": ratio,
                "se": ratio_se,
                "ci95": [ratio - 1.96 * ratio_se, ratio + 1.96 * ratio_se],
            },
            "effective_power_for_abs_J": {
                "value": exponent,
                "se_delta_method": exponent_se,
                "ci95_delta_method": [
                    exponent - 1.96 * exponent_se,
                    exponent + 1.96 * exponent_se,
                ],
                "definition": "-log(abs(J65/J25))/log(65/25)",
            },
        },
        "transport": {
            "p_ref": P_REF,
            "p_ref_scope": "prescribed square-site critical reference; not a rigorous exact-pc enclosure",
            "root_minus_p_ref_J": paired_J,
            "root_minus_p_ref_J_se": paired_J_se,
            "absolute_fraction_of_root_J": abs(paired_J / root[OUTPUT_FIELDS.index("J")]),
        },
        "C4_omitted_column": {
            "omitted_displacement": [1, 0],
            "retained_reconstruction_directions": [list(item) for item in nn_directions],
            "formula": "F_full=sum_retained_63 F_d + (F_-e1+F_+e2+F_-e2)/3",
            "moments": ["a", "qa", "Ea", "aS", "qaS", "EaS"],
        },
        "jackknife": {
            "groups": [
                {"id": "new_N65_source_MC", "batches": BATCHES, "shared_across_four_shards": True},
                {"id": "P45_N65_baseline", "batches": BATCHES, "independent_of_source_MC": True},
            ],
            "variance": "sum over groups of (99/100)*sum_b(theta_-b-mean(theta_-b))^2",
            "root_minus_p_ref": "paired within each group before adding the two independent variance contributions",
        },
        "run": run,
        "inputs": {
            "baseline": {"path": str(args.baseline), "sha256": sha256(args.baseline)},
            "tables": [
                {"path": str(path), "sha256": sha256(path)} for path in args.tables
            ],
            "N25_result": {"path": str(args.n25_result), "sha256": sha256(args.n25_result)},
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    print(
        json.dumps(
            {
                "status": payload["status"],
                "J65": payload["N65_pooled_root"]["J"],
                "J65_over_J25": payload["N25_to_N65"]["J65_over_J25"],
                "root_minus_p_ref_J": payload["paired_root_minus_p_ref"]["J"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
