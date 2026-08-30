#!/usr/bin/env python3
"""Score the frozen P337 projective-F3 line-observer rank pilot."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import score_marked_birth_path as base  # noqa: E402


COUPLINGS = ("H_W", "H_S", "D_W", "D_S")
METRICS = tuple(
    f"{side}_{name}_{part}"
    for side in ("first", "second")
    for name in COUPLINGS
    for part in ("re", "im")
) + tuple(
    name
    for side in ("first", "second")
    for name in (f"{side}_det_re", f"{side}_det_im", f"{side}_normalized_wedge")
) + (
    "first_W_residual_fraction", "second_W_residual_fraction",
)

PRIMARY = ("first_det_re", "first_det_im", "second_det_re", "second_det_im")


def mean(row: base.PathRow, column: str) -> mp.mpf:
    return row.values[column] / row.samples


def orientation(row: base.PathRow) -> tuple[dict[str, mp.mpf], dict[str, Any]]:
    h, d = mean(row, "sum_F3_H"), mean(row, "sum_F3_D")
    w = mp.mpc(mean(row, "sum_W_line_re"), mean(row, "sum_W_line_im"))
    source = mp.mpc(mean(row, "sum_J_S_re"), mean(row, "sum_J_S_im"))

    def complex_mean(stem: str) -> mp.mpc:
        return mp.mpc(mean(row, stem + "_re"), mean(row, stem + "_im"))

    couplings = {
        "H_W": complex_mean("sum_F3_H_W_line") - h * w,
        "H_S": complex_mean("sum_F3_H_J_S") - h * source,
        "D_W": complex_mean("sum_F3_D_W_line") - d * w,
        "D_S": complex_mean("sum_F3_D_J_S") - d * source,
    }
    determinant = couplings["H_W"] * couplings["D_S"] - couplings["H_S"] * couplings["D_W"]
    h_row = mp.sqrt(abs(couplings["H_W"]) ** 2 + abs(couplings["H_S"]) ** 2)
    d_row = mp.sqrt(abs(couplings["D_W"]) ** 2 + abs(couplings["D_S"]) ** 2)
    wedge = abs(determinant) / (h_row * d_row) if h_row and d_row else mp.nan

    w2 = mean(row, "sum_abs_W_line2") - abs(w) ** 2
    s2 = mean(row, "sum_abs_J_S2") - abs(source) ** 2
    cross = mp.mpc(
        mean(row, "sum_W_line_conj_J_S_re"),
        mean(row, "sum_W_line_conj_J_S_im"),
    ) - mp.conj(w) * source
    residual = 1 - abs(cross) ** 2 / (w2 * s2) if w2 > 0 and s2 > 0 else mp.nan
    point = {
        **{f"{name}_re": mp.re(value) for name, value in couplings.items()},
        **{f"{name}_im": mp.im(value) for name, value in couplings.items()},
        "det_re": mp.re(determinant),
        "det_im": mp.im(determinant),
        "normalized_wedge": wedge,
        "W_residual_fraction": residual,
    }
    detail = {
        "mean": {"F3_H": h, "F3_D": d, "W_line": w, "JS": source},
        "observer_gram": {
            "var_F3_H": mean(row, "sum_F3_H2") - h * h,
            "var_F3_D": mean(row, "sum_F3_D2") - d * d,
            "cov_F3_H_D": mean(row, "sum_F3_H_D") - h * d,
        },
        "source_gram": {
            "var_W_line": w2,
            "var_JS": s2,
            "conj_W_line_JS": cross,
            "W_line_residual_fraction_after_JS": residual,
        },
    }
    return point, detail


def evaluate(
    groups: Mapping[tuple[int, str, int], list[base.PathRow]],
    n: int,
    batches: Sequence[int],
    k0: int,
) -> tuple[dict[str, mp.mpf], dict[str, Any]]:
    point: dict[str, mp.mpf] = {}
    detail = {}
    for side in ("first", "second"):
        rows = base.combine([groups[(n, side, batch)] for batch in batches])
        one, one_detail = orientation(rows[k0])
        for name in COUPLINGS:
            for part in ("re", "im"):
                point[f"{side}_{name}_{part}"] = one[f"{name}_{part}"]
        for name in ("det_re", "det_im", "normalized_wedge"):
            point[f"{side}_{name}"] = one[name]
        point[f"{side}_W_residual_fraction"] = one["W_residual_fraction"]
        detail[side] = one_detail
    return point, detail


def covariance(rows: Sequence[Mapping[str, mp.mpf]]) -> list[list[mp.mpf]]:
    count = len(rows)
    center = {name: mp.fsum(row[name] for row in rows) / count for name in METRICS}
    factor = mp.mpf(count - 1) / count
    return [[
        factor * mp.fsum(
            (row[left] - center[left]) * (row[right] - center[right]) for row in rows
        )
        for right in METRICS
    ] for left in METRICS]


def chi2(point: Mapping[str, mp.mpf], cov: Sequence[Sequence[mp.mpf]], names: Sequence[str]) -> mp.mpf:
    indices = [METRICS.index(name) for name in names]
    matrix = mp.matrix([[cov[i][j] for j in indices] for i in indices])
    vector = mp.matrix([point[name] for name in names])
    try:
        return (vector.T * matrix ** -1 * vector)[0]
    except (ZeroDivisionError, ValueError):
        return mp.nan


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_audit(path: Path) -> dict[str, Any]:
    fields = (
        "endpoint_failures", "site_failures", "line_failures", "local_mark_failures",
        "index_mismatches", "separated_mark_failures",
    )
    totals = {name: 0 for name in fields}
    rows = 0
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            rows += 1
            for name in fields:
                totals[name] += int(raw[name])
    if any(totals.values()):
        raise ValueError(f"nonzero exact audit: {totals}")
    return {"rows": rows, "totals": totals}


def render(value: Any) -> Any:
    if isinstance(value, (mp.mpf, mp.mpc)):
        return base._text(value)
    if isinstance(value, dict):
        return {key: render(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [render(item) for item in value]
    return value


def score(prefix: Path, spec: Mapping[str, Any], host: str) -> dict[str, Any]:
    path_file = Path(str(prefix) + ".path.csv")
    audit_file = Path(str(prefix) + ".complement_audit.csv")
    metadata_file = Path(str(prefix) + ".metadata.json")
    groups = base.read_path(path_file)
    matrix = spec["first_matrix"]
    n = abs(int(matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]))
    batches = sorted(
        {key[2] for key in groups if key[:2] == (n, "first")}
        & {key[2] for key in groups if key[:2] == (n, "second")}
    )
    if len(batches) != int(spec["batches"]):
        raise ValueError("batch count differs from freeze")
    point, detail = evaluate(groups, n, batches, int(spec["k0"]))
    delete_rows = [
        evaluate(groups, n, [batch for batch in batches if batch != omitted], int(spec["k0"]))[0]
        for omitted in batches
    ]
    cov = covariance(delete_rows)
    standard_error = {
        name: mp.sqrt(max(mp.mpf(0), cov[index][index]))
        for index, name in enumerate(METRICS)
    }
    joint = chi2(point, cov, PRIMARY)
    side_tests = {
        side: chi2(point, cov, (f"{side}_det_re", f"{side}_det_im"))
        for side in ("first", "second")
    }
    metadata = json.loads(metadata_file.read_text(encoding="utf-8"))
    expected = {
        "git_commit": "c7eef55", "samples_per_pair": int(spec["samples"]),
        "batches": int(spec["batches"]), "seed": int(spec["seed"]),
        "replica_counter_first": int(spec["replica_offset"]),
    }
    for name, value in expected.items():
        if metadata.get(name) != value:
            raise ValueError(f"metadata mismatch {name}")
    return render({
        "schema": "matching-one/p337-f3-line-observer-pilot-score/v1",
        "N": n, "k0": int(spec["k0"]), "samples": int(spec["samples"]),
        "batches": len(batches), "actual_host": host,
        "generated_utc": metadata.get("generated_utc"),
        "elapsed_seconds": metadata.get("elapsed_seconds"),
        "metric_order": list(METRICS),
        "point": point,
        "standard_error": standard_error,
        "delete_one_covariance": cov,
        "primary": {
            "joint_determinant_chi2_4d": joint,
            "resolved_at_alpha_0_01": bool(joint > mp.mpf("13.27670413599")),
            "per_orientation_chi2_2d": side_tests,
        },
        "detail": detail,
        "audit": validate_audit(audit_file),
        "provenance": {
            "path_sha256": sha256(path_file),
            "metadata_sha256": sha256(metadata_file),
            "counter_range": [int(spec["replica_offset"]), int(spec["replica_offset"]) + int(spec["samples"])],
        },
    })


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze", type=Path, default=ROOT / "analysis/p337_f3_line_observer_pilot_freeze.json")
    parser.add_argument("--raw-root", type=Path, default=ROOT / "results/server-20260830/P337-f3-line-observer-pilot/raw")
    parser.add_argument("--output", type=Path, default=ROOT / "results/server-20260830/P337-f3-line-observer-pilot/p337_f3_line_observer_pilot_score.json")
    parser.add_argument("--dps", type=int, default=50)
    args = parser.parse_args()
    mp.mp.dps = args.dps
    freeze = json.loads(args.freeze.read_text(encoding="utf-8"))
    hosts = {"N325": "Huawei-CodeBuddy-TgFr7R", "N425": "Huawei-CodeBuddy-XPk2PZ"}
    scores = {
        run_id: score(args.raw_root / run_id / f"{run_id}_20k", freeze["runs"][run_id], hosts[run_id])
        for run_id in ("N325", "N425")
    }
    result = {
        "schema": "matching-one/p337-f3-line-observer-pilot-result/v1",
        "status": "completed_frozen_pilot",
        "new_samples": {"N325": 20000, "N425": 20000},
        "scores": scores,
        "decision": {
            "coupling_rank_lift_by_size": {
                run_id: value["primary"]["resolved_at_alpha_0_01"]
                for run_id, value in scores.items()
            },
            "pooling": "none; N325 and N425 are independent evidence blocks",
        },
        "claim_boundary": [
            "The observer rows are exact finite-abelian functions of the current projective winding line.",
            "A coupling-rank lift is finite-size evidence that W_line and JS address distinct projective-line responses; it is not a field or exponent identity.",
            "The same-state construction is intentionally line-addressed and is not independent of the ambient-H1 filtration itself.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result["decision"], indent=2))


if __name__ == "__main__":
    main()
