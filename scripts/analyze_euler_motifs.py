#!/usr/bin/env python3
"""Pilot-frozen Euler/motif control variates from C++ batch moments.

Reads ``*.motifs.jsonl`` written by ``src/gaussian_orientation_mc.cpp --euler-motifs``.
Coefficients are frozen on a leading pilot fraction of independent batches and
evaluated on the remaining replicas.  Duplicate wrapping channels are rejected
rather than GLS-combined.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Iterable, Sequence

from control_variate_estimator import (
    DuplicateChannelError,
    FrozenEstimator,
    _solve,
    invertibility_test,
)
from euler_motif_controls import (
    ALL_MOTIFS,
    EULER_CONTROLS,
    EXTRA_MOTIFS,
    WRAPPING_CHANNELS,
    analytic_motif_mean,
    chi,
    microcanonical_motif_mean,
)


def sample_variance(values: Sequence[float]) -> float:
    if len(values) < 2:
        raise ValueError("at least two values are required")
    mean = math.fsum(values) / len(values)
    return math.fsum((value - mean) ** 2 for value in values) / (len(values) - 1)


def read_jsonl(path: Path) -> list[dict]:
    rows = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    if not rows:
        raise ValueError(path.as_posix() + " is empty")
    return rows


def pool_geometry(rows: Sequence[dict], which: str) -> dict[str, object]:
    names = list(rows[0]["names"])
    width = len(names)
    samples = 0
    total = [0.0] * width
    gram = [[0.0] * width for _ in range(width)]
    wrapping_l1 = 0.0
    identity_l1 = 0.0
    a = b = None
    n = None
    p_ref = None
    for row in rows:
        geom = row[which]
        if a is None:
            a, b = geom["a"], geom["b"]
            n = row["n"]
            p_ref = float(row["p_ref"])
        elif (geom["a"], geom["b"], row["n"]) != (a, b, n):
            raise ValueError("inconsistent geometry in motif JSONL")
        if list(row["names"]) != names:
            raise ValueError("inconsistent observable names")
        samples += int(row["samples"])
        wrapping_l1 += float(geom["wrapping_l1"])
        identity_l1 += float(geom["identity_l1"])
        for i, value in enumerate(geom["sum"]):
            total[i] += float(value)
        for i, gram_row in enumerate(geom["gram"]):
            for j, value in enumerate(gram_row):
                gram[i][j] += float(value)
    means = [value / samples for value in total]
    covariance = [[0.0] * width for _ in range(width)]
    for i in range(width):
        for j in range(width):
            covariance[i][j] = (gram[i][j] - samples * means[i] * means[j]) / (samples - 1)
    return {
        "n": n,
        "a": a,
        "b": b,
        "p_ref": p_ref,
        "samples": samples,
        "names": names,
        "sum": total,
        "means": means,
        "gram": gram,
        "covariance": covariance,
        "wrapping_l1": wrapping_l1,
        "identity_l1": identity_l1,
    }


def index_map(names: Sequence[str]) -> dict[str, int]:
    return {name: i for i, name in enumerate(names)}


def greedy_beta(
    covariance: Sequence[Sequence[float]], names: Sequence[str], target_index: int,
    control_indices: Sequence[int],
) -> tuple[list[int], list[float]]:
    selected: list[int] = []
    dropped: list[int] = []
    for index in control_indices:
        trial = selected + [index]
        sub = [[covariance[i][j] for j in trial] for i in trial]
        try:
            invertibility_test(sub)
        except ArithmeticError:
            dropped.append(index)
            continue
        selected = trial
    if not selected:
        raise ArithmeticError("no full-rank motif control subset")
    sub = [[covariance[i][j] for j in selected] for i in selected]
    cross = [covariance[i][target_index] for i in selected]
    return selected, _solve(sub, cross)


def quadratic_form(covariance: Sequence[Sequence[float]], weights: Sequence[float]) -> float:
    total = 0.0
    for i, wi in enumerate(weights):
        for j, wj in enumerate(weights):
            total += wi * covariance[i][j] * wj
    return total


def estimator_variance(
    covariance: Sequence[Sequence[float]],
    target_index: int,
    control_indices: Sequence[int],
    beta: Sequence[float],
) -> float:
    width = len(covariance)
    weights = [0.0] * width
    weights[target_index] = 1.0
    for index, value in zip(control_indices, beta):
        weights[index] -= value
    return quadratic_form(covariance, weights)


def motif_analytic_means(n: int, p: float, names: Sequence[str]) -> dict[str, float]:
    means = {}
    for name in names:
        if name in ALL_MOTIFS:
            means[name] = analytic_motif_mean(name, n, p)
    return means


def build_control_plan(names: Sequence[str], n: int, p: float) -> list[dict[str, object]]:
    lookup = index_map(names)
    plans = [
        {"id": "euler_canonical", "controls": list(EULER_CONTROLS), "centering": "canonical"},
        {
            "id": "euler_plus_motifs_canonical",
            "controls": list(EULER_CONTROLS) + list(EXTRA_MOTIFS),
            "centering": "canonical",
        },
    ]
    return plans


def microcanonical_rows_from_moments_not_available() -> None:
    return None


def evaluate_plan(
    pooled_eval: dict[str, object],
    pooled_pilot: dict[str, object],
    control_names: Sequence[str],
    centering: str,
) -> dict[str, object]:
    names = pooled_pilot["names"]
    lookup = index_map(names)
    n = int(pooled_pilot["n"])
    p = float(pooled_pilot["p_ref"])
    target_index = lookup["q"]
    control_indices = [lookup[name] for name in control_names]
    selected, beta = greedy_beta(
        pooled_pilot["covariance"], names, target_index, control_indices
    )
    selected_names = [names[i] for i in selected]
    dropped = [name for name in control_names if name not in selected_names]
    analytic = [
        0.0 if name.endswith("_mc") else analytic_motif_mean(name, n, p)
        for name in selected_names
    ]
    var_q = pooled_eval["covariance"][target_index][target_index]
    var_hat = estimator_variance(pooled_eval["covariance"], target_index, selected, beta)
    # Euler identity estimator D_cluster = (C_black - C_white) - N chi(p)
    i_cb, i_cw = lookup["C_black"], lookup["C_white"]
    d_weights = [0.0] * len(names)
    d_weights[i_cb] = 1.0
    d_weights[i_cw] = -1.0
    var_d = quadratic_form(pooled_eval["covariance"], d_weights)
    means = pooled_eval["means"]
    mean_q = means[target_index]
    mean_d = means[i_cb] - means[i_cw] - n * chi(p)
    mean_hat = mean_q
    for index, weight, mu in zip(selected, beta, analytic):
        mean_hat -= weight * (means[index] - mu)
    best_single_var = min(var_q, var_d)
    best_single = "q" if var_q <= var_d else "D_cluster"
    return {
        "centering": centering,
        "control_names": selected_names,
        "dropped_controls": dropped,
        "weights": beta,
        "analytic_control_means": analytic,
        "eval_mean_q": mean_q,
        "eval_mean_D_cluster": mean_d,
        "eval_mean_adjusted": mean_hat,
        "eval_var_q": var_q,
        "eval_var_D_cluster": var_d,
        "eval_var_adjusted": var_hat,
        "best_single": best_single,
        "best_single_variance": best_single_var,
        "variance_reduction_vs_q": (var_q / var_hat) if var_hat > 0 else None,
        "variance_reduction_vs_D_cluster": (var_d / var_hat) if var_hat > 0 else None,
        "variance_reduction_vs_best_single": (best_single_var / var_hat) if var_hat > 0 else None,
        "identity_l1_eval": pooled_eval["identity_l1"],
        "wrapping_l1_eval": pooled_eval["wrapping_l1"],
    }


def orientation_difference_pool(rows: Sequence[dict]) -> dict[str, object]:
    names = list(rows[0]["names"])
    width = len(names)
    samples = 0
    total = [0.0] * width
    gram = [[0.0] * width for _ in range(width)]
    wrapping_l1 = 0.0
    identity_l1 = 0.0
    n = rows[0]["n"]
    p_ref = float(rows[0]["p_ref"])
    for row in rows:
        first = row["first"]
        second = row["second"]
        samples += int(row["samples"])
        wrapping_l1 += float(first["wrapping_l1"]) + float(second["wrapping_l1"])
        identity_l1 += float(first["identity_l1"]) + float(second["identity_l1"])
        delta_sum = [float(a) - float(b) for a, b in zip(first["sum"], second["sum"])]
        for i in range(width):
            total[i] += delta_sum[i]
        # Gram of (X-Y): Gx + Gy - Cxy - Cyx. Cross gram is not stored separately,
        # so orientation-difference covariance is formed from concatenated
        # first-second replica vectors reconstructed via the stored per-geometry
        # grams only for the marginals.  For paired Δ we need the joint gram.
        # The JSONL stores per-geometry grams, not the 22x22 joint.  Use batch
        # means of Δq, ΔX as a conservative fallback; the replica-level Δ
        # analysis is computed below from per-batch sums (equal-size batches).
    batch_vectors = []
    for row in rows:
        count = float(row["samples"])
        first = row["first"]["sum"]
        second = row["second"]["sum"]
        batch_vectors.append([(float(a) - float(b)) / count for a, b in zip(first, second)])
    means = [math.fsum(vector[i] for vector in batch_vectors) / len(batch_vectors)
             for i in range(width)]
    covariance = [[0.0] * width for _ in range(width)]
    for i in range(width):
        for j in range(i, width):
            value = math.fsum(
                (vector[i] - means[i]) * (vector[j] - means[j]) for vector in batch_vectors
            ) / (len(batch_vectors) - 1)
            covariance[i][j] = covariance[j][i] = value
    replica_means = [value / samples for value in total]
    return {
        "n": n,
        "p_ref": p_ref,
        "samples": samples,
        "batch_count": len(rows),
        "names": names,
        "means": replica_means,
        "batch_mean_covariance": covariance,
        "wrapping_l1": wrapping_l1,
        "identity_l1": identity_l1,
        "unit": "covariance of equal-size batch means of first-minus-second",
    }


def reject_wrapping_gls(batch_csv: Path, n: int) -> dict[str, object]:
    rows_by_batch: dict[int, dict[str, float]] = {}
    with batch_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            if int(raw["n"]) != n:
                continue
            batch = int(raw["batch"])
            samples = float(raw["samples"])
            q = (float(raw["first_primal_sum"]) - float(raw["first_matching_sum"])) / samples
            rows_by_batch.setdefault(batch, {})[raw["channel"]] = q
    ordered = []
    for batch in sorted(rows_by_batch):
        payload = rows_by_batch[batch]
        ordered.append([payload[name] for name in WRAPPING_CHANNELS])
    try:
        FrozenEstimator.fit(WRAPPING_CHANNELS, ordered)
        rejected = False
        message = "GLS accepted duplicate wrapping channels"
    except DuplicateChannelError as exc:
        rejected = True
        message = str(exc)
    return {
        "N": n,
        "duplicate_wrapping_rejected": rejected,
        "message": message,
        "batch_count": len(ordered),
    }


def analyze(
    motif_rows: Sequence[dict],
    batch_csv: Path,
    pilot_fraction: float,
) -> dict[str, object]:
    by_n: dict[int, list[dict]] = {}
    for row in motif_rows:
        by_n.setdefault(int(row["n"]), []).append(row)
    report = {"schema": "P34 Euler/motif control-variate evaluation v1", "by_N": {}}
    for n, rows in sorted(by_n.items()):
        rows = sorted(rows, key=lambda item: int(item["batch"]))
        if len(rows) < 4:
            raise ValueError("need at least four batches to split pilot/evaluation")
        split = max(2, int(round(pilot_fraction * len(rows))))
        if split >= len(rows) - 1:
            raise ValueError("pilot_fraction leaves too few evaluation batches")
        pilot_rows = rows[:split]
        eval_rows = rows[split:]
        entry = {
            "N": n,
            "pilot_batches": [int(row["batch"]) for row in pilot_rows],
            "evaluation_batches": [int(row["batch"]) for row in eval_rows],
            "pilot_fraction": len(pilot_rows) / len(rows),
            "wrapping_gls_rejection": reject_wrapping_gls(batch_csv, n),
            "geometries": {},
            "orientation_difference": {},
        }
        for which in ("first", "second"):
            pilot = pool_geometry(pilot_rows, which)
            evaluation = pool_geometry(eval_rows, which)
            if evaluation["identity_l1"] != 0.0 or pilot["identity_l1"] != 0.0:
                identity_status = "FAIL"
            else:
                identity_status = "PASS"
            if evaluation["wrapping_l1"] != 0.0 or pilot["wrapping_l1"] != 0.0:
                wrap_status = "FAIL"
            else:
                wrap_status = "PASS"
            plans = {}
            for control_names, plan_id, centering in (
                (EULER_CONTROLS, "euler_canonical", "canonical"),
                (("V", "E_mc", "F0_mc"), "euler_microcanonical", "microcanonical"),
                (EULER_CONTROLS + EXTRA_MOTIFS, "euler_plus_motifs_canonical", "canonical"),
                (("V", "E_mc", "F0_mc", "nnn_pos_mc", "nnn_neg_mc",
                  "path3_x_mc", "path3_y_mc", "corners_mc"),
                 "euler_plus_motifs_microcanonical", "microcanonical"),
            ):
                plans[plan_id] = evaluate_plan(evaluation, pilot, control_names, centering)
            entry["geometries"][which] = {
                "a": evaluation["a"],
                "b": evaluation["b"],
                "identity_status": identity_status,
                "wrapping_duplicate_status": wrap_status,
                "pilot_samples": pilot["samples"],
                "evaluation_samples": evaluation["samples"],
                "eval_mean_q": evaluation["means"][index_map(evaluation["names"])["q"]],
                "plans": plans,
            }
        # Orientation difference uses batch-mean covariance (paired).
        pilot_delta = orientation_difference_pool(pilot_rows)
        eval_delta = orientation_difference_pool(eval_rows)
        names = pilot_delta["names"]
        lookup = index_map(names)
        delta_plans = {}
        for control_names, plan_id in (
            (EULER_CONTROLS, "euler_canonical"),
            (("V", "E_mc", "F0_mc"), "euler_microcanonical"),
            (EULER_CONTROLS + EXTRA_MOTIFS, "euler_plus_motifs_canonical"),
            (("V", "E_mc", "F0_mc", "nnn_pos_mc", "nnn_neg_mc",
              "path3_x_mc", "path3_y_mc", "corners_mc"),
             "euler_plus_motifs_microcanonical"),
        ):
            selected, beta = greedy_beta(
                pilot_delta["batch_mean_covariance"],
                names,
                lookup["q"],
                [lookup[name] for name in control_names],
            )
            var_q = eval_delta["batch_mean_covariance"][lookup["q"]][lookup["q"]]
            var_hat = estimator_variance(
                eval_delta["batch_mean_covariance"], lookup["q"], selected, beta
            )
            i_cb, i_cw = lookup["C_black"], lookup["C_white"]
            d_weights = [0.0] * len(names)
            d_weights[i_cb] = 1.0
            d_weights[i_cw] = -1.0
            var_d = quadratic_form(eval_delta["batch_mean_covariance"], d_weights)
            best = min(var_q, var_d)
            delta_plans[plan_id] = {
                "control_names": [names[i] for i in selected],
                "dropped_controls": [name for name in control_names if name not in {names[i] for i in selected}],
                "weights": beta,
                "eval_batch_var_delta_q": var_q,
                "eval_batch_var_delta_D_cluster": var_d,
                "eval_batch_var_adjusted": var_hat,
                "variance_reduction_vs_best_single": (best / var_hat) if var_hat > 0 else None,
                "best_single": "delta_q" if var_q <= var_d else "delta_D_cluster",
            }
        entry["orientation_difference"] = {
            "unit": eval_delta["unit"],
            "pilot_batches": len(pilot_rows),
            "evaluation_batches": len(eval_rows),
            "plans": delta_plans,
        }
        report["by_N"][str(n)] = entry
    return report


def gpu_gate(report: dict[str, object]) -> dict[str, object]:
    ratios = []
    for n, entry in report["by_N"].items():
        for which, geom in entry["geometries"].items():
            plan = geom["plans"]["euler_plus_motifs_microcanonical"]
            ratio = plan["variance_reduction_vs_best_single"]
            ratios.append(
                {
                    "N": int(n),
                    "geometry": which,
                    "a": geom["a"],
                    "b": geom["b"],
                    "vr_vs_best_single": ratio,
                    "passes_2x": ratio is not None and ratio >= 2.0,
                }
            )
    passing_sizes = sorted({
        item["N"] for item in ratios if item["passes_2x"]
    })
    return {
        "threshold": 2.0,
        "started_gpu": False,
        "per_geometry": ratios,
        "sizes_with_vr_at_least_2x": passing_sizes,
        "gate_pass": len(passing_sizes) >= 2,
    }


def write_variance_csv(path: Path, report: dict[str, object]) -> None:
    rows = []
    for n, entry in report["by_N"].items():
        for which, geom in entry["geometries"].items():
            for plan_id, plan in geom["plans"].items():
                rows.append({
                    "N": n,
                    "geometry": which,
                    "a": geom["a"],
                    "b": geom["b"],
                    "plan": plan_id,
                    "best_single": plan["best_single"],
                    "var_q": plan["eval_var_q"],
                    "var_D_cluster": plan["eval_var_D_cluster"],
                    "var_adjusted": plan["eval_var_adjusted"],
                    "vr_vs_q": plan["variance_reduction_vs_q"],
                    "vr_vs_D_cluster": plan["variance_reduction_vs_D_cluster"],
                    "vr_vs_best_single": plan["variance_reduction_vs_best_single"],
                })
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--motifs", type=Path, required=True)
    parser.add_argument("--batches", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--weights-json", type=Path, required=True)
    parser.add_argument("--variance-csv", type=Path, required=True)
    parser.add_argument("--pilot-fraction", type=float, default=0.2)
    args = parser.parse_args()
    if not 0.0 < args.pilot_fraction < 0.5:
        raise SystemExit("pilot-fraction must lie in (0, 0.5)")
    motif_rows = read_jsonl(args.motifs)
    report = analyze(motif_rows, args.batches, args.pilot_fraction)
    report["gpu_gate"] = gpu_gate(report)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    weights = {
        "schema": "P34 frozen Euler/motif OLS weights v1",
        "pilot_fraction": args.pilot_fraction,
        "by_N": {},
    }
    for n, entry in report["by_N"].items():
        weights["by_N"][n] = {
            which: geom["plans"] for which, geom in entry["geometries"].items()
        }
        weights["by_N"][n]["orientation_difference"] = entry["orientation_difference"]["plans"]
    args.weights_json.write_text(json.dumps(weights, indent=2) + "\n", encoding="utf-8")
    write_variance_csv(args.variance_csv, report)
    print("wrote " + str(args.json))
    print("wrote " + str(args.weights_json))
    print("wrote " + str(args.variance_csv))
    gate = report["gpu_gate"]
    print("GPU 2x gate: " + ("PASS" if gate["gate_pass"] else "FAIL") +
          " sizes=" + ",".join(str(v) for v in gate["sizes_with_vr_at_least_2x"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
