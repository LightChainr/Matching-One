#!/usr/bin/env python3
"""Score the frozen P337 ambient-H1 source-column pilot."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import score_marked_birth_path as base  # noqa: E402


METRICS = (
    "C_far_W_re", "C_far_W_im", "C_far_S_re", "C_far_S_im",
    "C_sep4_W_re", "C_sep4_W_im", "C_sep4_S_re", "C_sep4_S_im",
    "det_re", "det_im", "normalized_wedge",
    "W_residual_fraction_first", "W_residual_fraction_second",
    "conditional_mark12_h4", "conditional_completion",
    "observed_mark12_h4", "risk_composition_remainder",
)


def mean(row: base.PathRow, column: str) -> mp.mpf:
    return row.values[column] / row.samples


def orientation_lane(row: base.PathRow) -> dict[str, mp.mpf | mp.mpc]:
    w = mp.mpc(mean(row, "sum_W_line_re"), mean(row, "sum_W_line_im"))
    source = mp.mpc(mean(row, "sum_J_S_re"), mean(row, "sum_J_S_im"))

    def observer(name: str) -> mp.mpf:
        return mean(row, f"sum_{name}")

    def coupling(observer_name: str, source_name: str) -> mp.mpc:
        product = mp.mpc(
            mean(row, f"sum_{observer_name}_{source_name}_re"),
            mean(row, f"sum_{observer_name}_{source_name}_im"),
        )
        source_mean = w if source_name == "W_line" else source
        return product - observer(observer_name) * source_mean

    far_w = coupling("O_ext", "W_line") - coupling("O_near", "W_line")
    far_s = coupling("O_ext", "J_S") - coupling("O_near", "J_S")
    sep_w = coupling("O_sep4", "W_line")
    sep_s = coupling("O_sep4", "J_S")

    w2 = mean(row, "sum_abs_W_line2") - abs(w) ** 2
    s2 = mean(row, "sum_abs_J_S2") - abs(source) ** 2
    cross = mp.mpc(
        mean(row, "sum_W_line_conj_J_S_re"),
        mean(row, "sum_W_line_conj_J_S_im"),
    ) - mp.conj(w) * source
    if w2 <= 0 or s2 <= 0:
        residual = mp.nan
    else:
        residual = 1 - abs(cross) ** 2 / (w2 * s2)
        if residual < 0 and abs(residual) < mp.mpf("1e-30"):
            residual = mp.mpf(0)
    return {
        "far_W": far_w, "far_S": far_s,
        "sep4_W": sep_w, "sep4_S": sep_s,
        "W_mean": w, "JS_mean": source,
        "W_centered_norm2": w2, "JS_centered_norm2": s2,
        "W_conj_JS_centered": cross,
        "W_residual_fraction": residual,
    }


def path_metrics(
    first: Sequence[base.PathRow], second: Sequence[base.PathRow], k0: int
) -> tuple[dict[str, mp.mpf], dict[str, Any]]:
    left = orientation_lane(first[k0])
    right = orientation_lane(second[k0])
    leverage = base.cos4(first[0].a, first[0].b) - base.cos4(second[0].a, second[0].b)
    if leverage == 0:
        raise ValueError("zero frozen H4 leverage")

    def p4(name: str) -> mp.mpc:
        return (left[name] - right[name]) / leverage

    far_w, far_s = p4("far_W"), p4("far_S")
    sep_w, sep_s = p4("sep4_W"), p4("sep4_S")
    determinant = far_w * sep_s - far_s * sep_w
    row_norm = mp.sqrt(abs(far_w) ** 2 + abs(far_s) ** 2)
    other_norm = mp.sqrt(abs(sep_w) ** 2 + abs(sep_s) ** 2)
    wedge = abs(determinant) / (row_norm * other_norm) if row_norm and other_norm else mp.nan
    metrics = {
        "C_far_W_re": mp.re(far_w), "C_far_W_im": mp.im(far_w),
        "C_far_S_re": mp.re(far_s), "C_far_S_im": mp.im(far_s),
        "C_sep4_W_re": mp.re(sep_w), "C_sep4_W_im": mp.im(sep_w),
        "C_sep4_S_re": mp.re(sep_s), "C_sep4_S_im": mp.im(sep_s),
        "det_re": mp.re(determinant), "det_im": mp.im(determinant),
        "normalized_wedge": wedge,
        "W_residual_fraction_first": left["W_residual_fraction"],
        "W_residual_fraction_second": right["W_residual_fraction"],
    }
    detail = {
        "k0": k0,
        "leverage_cos4_first_minus_second": leverage,
        "first": left,
        "second": right,
    }
    return metrics, detail


MarkCell = list[int]


def read_marks(path: Path, n: int, k0: int) -> tuple[dict[str, dict[int, dict[tuple[int, int, int], MarkCell]]], list[int]]:
    output: dict[str, dict[int, dict[tuple[int, int, int], MarkCell]]] = {
        side: defaultdict(lambda: defaultdict(lambda: [0, 0, 0]))
        for side in ("first", "second")
    }
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            if int(raw["n"]) != n:
                raise ValueError("marked-birth size changed")
            side = raw["orientation"]
            k1, k2 = int(raw["k1"]), int(raw["k2"])
            if side not in output or not (k1 <= k0 < k2):
                continue
            if raw["direct_0_to_2"] == "1" or raw["line_null"] == "1":
                raise ValueError("direct/null line entered strict rank-one risk set")
            if raw["mark12_valid"] != "1":
                raise ValueError("rank-one risk row has invalid completion mark")
            key = (int(raw["ell_u"]), int(raw["ell_v"]), k0 - k1)
            cell = output[side][int(raw["batch"])][key]
            count = int(raw["count"])
            event = int(k2 == k0 + 1)
            cell[0] += count
            cell[1] += count * event
            cell[2] += count * event * int(raw["mark12_h4"])
    batches = sorted(set(output["first"]) & set(output["second"]))
    if len(batches) < 2:
        raise ValueError("marked score needs aligned batches")
    return output, batches


def merge_marks(
    source: Mapping[int, Mapping[tuple[int, int, int], MarkCell]],
    kept: Sequence[int],
) -> dict[tuple[int, int, int], MarkCell]:
    merged: dict[tuple[int, int, int], MarkCell] = defaultdict(lambda: [0, 0, 0])
    for batch in kept:
        for key, values in source[batch].items():
            for index, value in enumerate(values):
                merged[key][index] += value
    return merged


def mark_metrics(
    first: Mapping[tuple[int, int, int], MarkCell],
    second: Mapping[tuple[int, int, int], MarkCell],
    leverage: mp.mpf,
) -> tuple[dict[str, mp.mpf], dict[str, Any]]:
    tables = (first, second)
    totals = [sum(values[0] for values in table.values()) for table in tables]
    common = sorted(set(first) & set(second))
    if not common or min(totals) <= 0:
        raise ValueError("empty common rank-one line-age support")
    coverage = [mp.mpf(sum(table[key][0] for key in common)) / total
                for table, total in zip(tables, totals)]
    target = {
        key: (mp.mpf(first[key][0]) / totals[0] + mp.mpf(second[key][0]) / totals[1]) / 2
        for key in common
    }
    target_norm = mp.fsum(target.values())
    target = {key: value / target_norm for key, value in target.items()}
    standardized_mark = []
    standardized_exit = []
    observed_mark = []
    for table, total in zip(tables, totals):
        standardized_mark.append(mp.fsum(
            target[key] * mp.mpf(table[key][2]) / table[key][0] for key in common
        ))
        standardized_exit.append(mp.fsum(
            target[key] * mp.mpf(table[key][1]) / table[key][0] for key in common
        ))
        observed_mark.append(mp.mpf(sum(row[2] for row in table.values())) / total)
    conditional = (standardized_mark[0] - standardized_mark[1]) / leverage
    conditional_exit = (standardized_exit[0] - standardized_exit[1]) / leverage
    observed = (observed_mark[0] - observed_mark[1]) / leverage
    metrics = {
        "conditional_mark12_h4": conditional,
        "conditional_completion": conditional_exit,
        "observed_mark12_h4": observed,
        "risk_composition_remainder": observed - conditional,
    }
    detail = {
        "common_line_age_strata": len(common),
        "common_support_coverage": {"first": coverage[0], "second": coverage[1]},
        "risk_totals": {"first": totals[0], "second": totals[1]},
        "target": "equal-orientation mixture of risk distributions on common (ell_u,ell_v,age) support",
        "standardized_mark12_h4": {"first": standardized_mark[0], "second": standardized_mark[1]},
        "standardized_completion": {"first": standardized_exit[0], "second": standardized_exit[1]},
        "observed_mark12_h4": {"first": observed_mark[0], "second": observed_mark[1]},
    }
    return metrics, detail


def covariance(rows: Sequence[dict[str, mp.mpf]]) -> list[list[mp.mpf]]:
    count = len(rows)
    means = {name: mp.fsum(row[name] for row in rows) / count for name in METRICS}
    factor = mp.mpf(count - 1) / count
    return [[
        factor * mp.fsum(
            (row[left] - means[left]) * (row[right] - means[right]) for row in rows
        )
        for right in METRICS
    ] for left in METRICS]


def determinant_chi2(point: Mapping[str, mp.mpf], cov: Sequence[Sequence[mp.mpf]]) -> mp.mpf:
    i, j = METRICS.index("det_re"), METRICS.index("det_im")
    a, b, c = cov[i][i], cov[i][j], cov[j][j]
    denominator = a * c - b * b
    if denominator <= 0:
        return mp.nan
    x, y = point["det_re"], point["det_im"]
    return (c * x * x - 2 * b * x * y + a * y * y) / denominator


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_audit(path: Path) -> dict[str, Any]:
    fields = (
        "endpoint_failures", "site_failures", "line_failures",
        "local_mark_failures", "index_mismatches", "separated_mark_failures",
    )
    totals = {field: 0 for field in fields}
    rows = 0
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            rows += 1
            for field in fields:
                totals[field] += int(raw[field])
    if any(totals.values()):
        raise ValueError(f"nonzero exact audit: {totals}")
    return {"rows": rows, "totals": totals}


def serializable(value: Any) -> Any:
    if isinstance(value, (mp.mpf, mp.mpc)):
        return base._text(value)
    if isinstance(value, dict):
        return {key: serializable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serializable(item) for item in value]
    return value


def score_prefix(prefix: Path, spec: Mapping[str, Any], actual_host: str) -> dict[str, Any]:
    path_file = Path(str(prefix) + ".path.csv")
    marks_file = Path(str(prefix) + ".marked_births.csv")
    audit_file = Path(str(prefix) + ".complement_audit.csv")
    metadata_file = Path(str(prefix) + ".metadata.json")
    groups = base.read_path(path_file)
    n = int(spec["first_matrix"][0][0] ** 2 + spec["first_matrix"][1][0] ** 2)
    # General determinant rather than the Gaussian shortcut above.
    matrix = spec["first_matrix"]
    n = abs(int(matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]))
    path_batches = sorted(
        {key[2] for key in groups if key[:2] == (n, "first")}
        & {key[2] for key in groups if key[:2] == (n, "second")}
    )
    mark_batches, sparse_batches = read_marks(marks_file, n, int(spec["k0"]))
    batches = sorted(set(path_batches) & set(sparse_batches))
    if len(batches) != int(spec["batches"]):
        raise ValueError("path/sparse batch alignment changed")

    def evaluate(kept: Sequence[int]) -> tuple[dict[str, mp.mpf], dict[str, Any]]:
        first_path = base.combine([groups[(n, "first", batch)] for batch in kept])
        second_path = base.combine([groups[(n, "second", batch)] for batch in kept])
        point, path_detail = path_metrics(first_path, second_path, int(spec["k0"]))
        mark_point, mark_detail = mark_metrics(
            merge_marks(mark_batches["first"], kept),
            merge_marks(mark_batches["second"], kept),
            path_detail["leverage_cos4_first_minus_second"],
        )
        point.update(mark_point)
        return point, {"path": path_detail, "completion": mark_detail}

    point, detail = evaluate(batches)
    delete_rows = [evaluate([batch for batch in batches if batch != omitted])[0]
                   for omitted in batches]
    cov = covariance(delete_rows)
    standard_error = {
        name: mp.sqrt(max(mp.mpf(0), cov[index][index]))
        for index, name in enumerate(METRICS)
    }
    chi2 = determinant_chi2(point, cov)
    mark_t = point["conditional_mark12_h4"] / standard_error["conditional_mark12_h4"]
    metadata = json.loads(metadata_file.read_text(encoding="utf-8"))
    expected = {
        "git_commit": "7d89171", "samples_per_pair": int(spec["samples"]),
        "batches": int(spec["batches"]), "seed": int(spec["seed"]),
        "replica_counter_first": int(spec["replica_offset"]),
    }
    for key, value in expected.items():
        if metadata.get(key) != value:
            raise ValueError(f"metadata mismatch {key}: {metadata.get(key)!r} != {value!r}")
    return serializable({
        "schema": "matching-one/p337-ambient-line-source-pilot-score/v1",
        "N": n,
        "k0": int(spec["k0"]),
        "samples": int(spec["samples"]),
        "batches": len(batches),
        "actual_host": actual_host,
        "generated_utc": metadata.get("generated_utc"),
        "elapsed_seconds": metadata.get("elapsed_seconds"),
        "metric_order": list(METRICS),
        "point": point,
        "standard_error": standard_error,
        "delete_one_covariance": cov,
        "primary": {
            "determinant_chi2_2d": chi2,
            "chi2_2d_upper_tail_p": mp.e ** (-chi2 / 2) if mp.isfinite(chi2) else mp.nan,
            "resolved_at_alpha_0_01": bool(chi2 > mp.mpf("9.21034037198")),
        },
        "completion_control": {
            "conditional_mark12_h4_student_t": mark_t,
            "degrees_of_freedom": len(batches) - 1,
            "resolved_at_two_sided_alpha_0_01": bool(abs(mark_t) > mp.mpf("2.860934606")),
        },
        "detail": detail,
        "audit": validate_audit(audit_file),
        "provenance": {
            "metadata": str(metadata_file),
            "metadata_sha256": sha256(metadata_file),
            "path_sha256": sha256(path_file),
            "marked_births_sha256": sha256(marks_file),
            "counter_range": [int(spec["replica_offset"]), int(spec["replica_offset"]) + int(spec["samples"])],
        },
    })


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze", type=Path, default=ROOT / "analysis/p337_ambient_line_source_pilot_freeze.json")
    parser.add_argument("--raw-root", type=Path, default=ROOT / "results/server-20260830/P337-ambient-line-source-pilot/raw")
    parser.add_argument("--output", type=Path, default=ROOT / "results/server-20260830/P337-ambient-line-source-pilot/p337_ambient_line_source_pilot_score.json")
    parser.add_argument("--dps", type=int, default=50)
    args = parser.parse_args()
    mp.mp.dps = args.dps
    freeze = json.loads(args.freeze.read_text(encoding="utf-8"))
    hosts = {"N325": "Huawei-CodeBuddy-TgFr7R", "N425": "Huawei-CodeBuddy-XPk2PZ"}
    scores = {}
    for run_id in ("N325", "N425"):
        prefix = args.raw_root / run_id / f"{run_id}_20k"
        scores[run_id] = score_prefix(prefix, freeze["runs"][run_id], hosts[run_id])
    result = {
        "schema": "matching-one/p337-ambient-line-source-pilot-result/v1",
        "status": "completed_frozen_pilot",
        "freeze": str(args.freeze),
        "new_samples": {"N325": 20000, "N425": 20000},
        "scores": scores,
        "decision": {
            "source_rank_lift_by_size": {
                run_id: score["primary"]["resolved_at_alpha_0_01"]
                for run_id, score in scores.items()
            },
            "completion_cause_by_size": {
                run_id: score["completion_control"]["resolved_at_two_sided_alpha_0_01"]
                for run_id, score in scores.items()
            },
            "pooling": "none; N325 and N425 are independent evidence blocks",
        },
        "claim_boundary": [
            "This pilot tests whether W_line adds an observer-space source direction beyond JS at fixed occupancy; it does not identify a field.",
            "The completion control standardizes on common primitive-line and plateau-age support and does not make generations or sizes into replicas.",
            "A 20k null does not prove exact source dependence; a resolved determinant only authorizes a larger frozen production lane.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result["decision"], indent=2))


if __name__ == "__main__":
    main()
