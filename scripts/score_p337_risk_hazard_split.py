#!/usr/bin/env python3
"""One-off zero-sample P337 completion-current mechanism split."""

from __future__ import annotations

import csv
from collections import defaultdict
from fractions import Fraction
import json
import math
from pathlib import Path
import sys
from typing import Any, Mapping

import numpy as np
from scipy.stats import t
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import analyze_p337_birth_state_current as state  # noqa: E402


SCHEMA = "matching-one/p337-risk-hazard-split/v1"
P334_RAW_ROOT = Path(
    "/Volumes/Mac Data/Research/论文项目总库/Matching-One-large-artifacts/"
    "P267-two-observer-source-rank-2m"
)
MARK_HEADER = (
    "n,a,b,orientation,batch,samples,k1,k2,direct_0_to_2,site01,site12,"
    "line_null,ell_u,ell_v,iota01,iota12,physical_x,physical_y,chi4_re,chi4_im,"
    "mark01_valid,mark01_axis,mark01_diagonal,mark01_landed,mark01_h4,"
    "mark12_valid,mark12_axis,mark12_diagonal,mark12_landed,mark12_h4,count"
)


def covariance(rows: list[list[float]]) -> np.ndarray:
    values = np.asarray(rows, dtype=float)
    centered = values - values.mean(axis=0)
    return (len(values) - 1.0) / len(values) * centered.T @ centered


def split_completion(
    first: Mapping[str, Any],
    second: Mapping[str, Any],
    n: int,
    p: float,
    delta_cos4: float,
) -> dict[str, Any]:
    """Apply the symmetric product identity cellwise on (line,current layer)."""
    phi = [0.0] + [n * value for value in state.binomial_pmf(n - 1, p)]
    first_samples = float(first["samples"])
    second_samples = float(second["samples"])
    composition = 0.0
    hazard = 0.0
    unmatched = 0.0
    unmatched_cells = 0
    for line in range(len(state.LINE_ORDER)):
        for current in range(n):
            weight = phi[current + 1]
            r1 = first["risk"][line][current] / first_samples
            r2 = second["risk"][line][current] / second_samples
            y1 = first["exit_y"][line][current] / first_samples
            y2 = second["exit_y"][line][current] / second_samples
            if r1 <= 0.0 and r2 <= 0.0:
                continue
            if r1 > 0.0 and r2 > 0.0:
                h1, h2 = y1 / r1, y2 / r2
            else:
                # A conditional hazard is undefined off support.  The unique observed
                # hazard is extended to the absent side, assigning the entire support
                # entry/exit contribution to risk composition rather than fabricating
                # a conditional-hazard contrast.
                present_h = y1 / r1 if r1 > 0.0 else y2 / r2
                h1 = h2 = present_h
                unmatched_cells += 1
                unmatched += abs(weight * (y2 - y1))
            r1 *= weight
            r2 *= weight
            composition += 0.5 * (h1 + h2) * (r2 - r1)
            hazard += 0.5 * (r1 + r2) * (h2 - h1)
    line_exit_first = [
        sum(first["exit_y"][line][k] for line in range(len(state.LINE_ORDER)))
        / first_samples
        for k in range(n + 1)
    ]
    line_exit_second = [
        sum(second["exit_y"][line][k] for line in range(len(state.LINE_ORDER)))
        / second_samples
        for k in range(n + 1)
    ]
    # exit_y is indexed by current=k2-1.
    delta_j12 = math.fsum(
        phi[k + 1] * (line_exit_second[k] - line_exit_first[k]) for k in range(n)
    )
    scale = p * (1.0 - p) / delta_cos4
    return {
        "vector": [scale * delta_j12, scale * composition, scale * hazard],
        "unscaled": {
            "delta_J12": delta_j12,
            "risk_composition": composition,
            "conditional_hazard": hazard,
            "closure_residual": delta_j12 - composition - hazard,
        },
        "support_extension": {
            "rule": "absent_side_hazard_equals_present_side_hazard; unmatched support assigned to composition",
            "unmatched_cells": unmatched_cells,
            "absolute_weighted_J12_mass": unmatched,
            "fraction_of_absolute_delta_J12": unmatched / max(abs(delta_j12), 1e-300),
        },
    }


def score_state_run(root: Path, run: Mapping[str, Any]) -> dict[str, Any]:
    loaded = state.read_run(root, run)
    n = int(run["N"])
    batch_ids = loaded["batch_ids"]
    rows = {
        orientation: [loaded["batches"][(orientation, batch)] for batch in batch_ids]
        for orientation in state.ORIENTATIONS
    }
    totals = {orientation: state.totals(rows[orientation], n) for orientation in state.ORIENTATIONS}
    delta_cos4 = float(
        Fraction(str(run["cos4"]["second"])) - Fraction(str(run["cos4"]["first"]))
    )

    def evaluate(omitted: int | None) -> tuple[float, dict[str, Any]]:
        omit = {
            orientation: rows[orientation][omitted] if omitted is not None else None
            for orientation in state.ORIENTATIONS
        }
        reduced = {
            orientation: state.totals(
                [row for index, row in enumerate(rows[orientation]) if index != omitted], n
            ) if omitted is not None else totals[orientation]
            for orientation in state.ORIENTATIONS
        }
        p = state.matching_root(totals, omit, n)
        return p, split_completion(reduced["first"], reduced["second"], n, p, delta_cos4)

    point_p, point = evaluate(None)
    leave = [evaluate(index)[1]["vector"] for index in range(len(batch_ids))]
    cov = covariance(leave)
    vector = np.asarray(point["vector"])
    se = np.sqrt(np.maximum(np.diag(cov), 0.0))
    difference = vector[1] - vector[2]
    difference_variance = cov[1, 1] + cov[2, 2] - 2.0 * cov[1, 2]
    difference_se = math.sqrt(max(float(difference_variance), 0.0))
    dominant = "risk_composition" if abs(vector[1]) > abs(vector[2]) else "conditional_hazard"
    coordinate_tests = {}
    for index, label in enumerate(("J12", "risk_composition", "conditional_hazard")):
        statistic = vector[index] / se[index] if se[index] else math.nan
        coordinate_tests[label] = {
            "student_t": float(statistic),
            "degrees_of_freedom": len(batch_ids) - 1,
            "two_sided_p": float(2.0 * t.sf(abs(statistic), len(batch_ids) - 1)),
            "resolved_at_alpha_0_01": bool(
                2.0 * t.sf(abs(statistic), len(batch_ids) - 1) < 0.01
            ),
        }
    return {
        "id": run["id"],
        "N": n,
        "role": run["role"],
        "dependency_group": run["dependency_group"],
        "p": point_p,
        "vector_order": [
            "angular_J12_log", "angular_risk_composition_log", "angular_conditional_hazard_log"
        ],
        "vector": vector.tolist(),
        "standard_error": se.tolist(),
        "delete_one_covariance": cov.tolist(),
        "delete_one_unit": "same batch removed from both orientations",
        "coordinate_tests": coordinate_tests,
        "exact_midpoint_split": point["unscaled"],
        "support_extension": point["support_extension"],
        "composition_minus_hazard": {
            "estimate": float(difference),
            "standard_error": difference_se,
            "student_t": float(difference / difference_se) if difference_se else None,
            "two_sided_p": float(2.0 * t.sf(abs(difference / difference_se), len(batch_ids) - 1))
            if difference_se else None,
        },
        "dominant_absolute_component": dominant,
        "mark12_h4": {
            "status": "not_scoreable",
            "reason": "the primary state-current birth archive has no mark12 columns",
        } if run["role"] == "four_generation_primary" else {
            "status": "available_in_source_and_scored_separately_below"
        },
    }


def empty_mark_batches(batch_count: int) -> dict[str, list[dict[tuple[int, int, int], list[int]]]]:
    return {
        orientation: [defaultdict(lambda: [0, 0, 0]) for _ in range(batch_count)]
        for orientation in state.ORIENTATIONS
    }


def parse_mark_archive(path: Path, n: int, k0: int, batches: int) -> tuple[dict[str, Any], dict[str, float]]:
    output = empty_mark_batches(batches)
    chi4: dict[str, float] = {}
    with path.open(encoding="utf-8") as handle:
        if handle.readline().rstrip("\n") != MARK_HEADER:
            raise ValueError(f"unexpected marked-birth header: {path}")
        for line in handle:
            fields = line.rstrip("\n").split(",")
            orientation = fields[3]
            k1, k2 = int(fields[6]), int(fields[7])
            direct, line_null = int(fields[8]), int(fields[11])
            if int(fields[0]) != n or orientation not in output:
                raise ValueError("marked archive geometry changed")
            chi4.setdefault(orientation, float(fields[18]))
            if not (k1 <= k0 < k2) or direct or line_null:
                continue
            if int(fields[25]) != 1:
                raise ValueError("rank-one risk row lacks valid completion mark")
            key = (int(fields[12]), int(fields[13]), k0 - k1)
            row = output[orientation][int(fields[4])][key]
            count = int(fields[30])
            event = k2 == k0 + 1
            row[0] += count
            row[1] += count * event
            row[2] += count * event * int(fields[29])
    return output, chi4


def merge_mark(
    batches: list[dict[tuple[int, int, int], list[int]]], omitted: int | None
) -> dict[tuple[int, int, int], list[int]]:
    result: dict[tuple[int, int, int], list[int]] = defaultdict(lambda: [0, 0, 0])
    for index, batch in enumerate(batches):
        if index == omitted:
            continue
        for key, values in batch.items():
            for j, value in enumerate(values):
                result[key][j] += value
    return result


def standardized_mark_vector(
    first: Mapping[tuple[int, int, int], list[int]],
    second: Mapping[tuple[int, int, int], list[int]],
    delta_cos4: float,
) -> tuple[list[float], dict[str, Any]]:
    totals = [sum(row[0] for row in table.values()) for table in (first, second)]
    common = sorted(set(first) & set(second))
    coverage = [sum(table[key][0] for key in common) / total for table, total in zip((first, second), totals)]
    qbar = {
        key: 0.5 * (first[key][0] / totals[0] + second[key][0] / totals[1])
        for key in common
    }
    qnorm = math.fsum(qbar.values())
    qbar = {key: value / qnorm for key, value in qbar.items()}
    std_mark = []
    std_exit = []
    observed_mark = []
    for table, total in zip((first, second), totals):
        std_mark.append(math.fsum(qbar[key] * table[key][2] / table[key][0] for key in common))
        std_exit.append(math.fsum(qbar[key] * table[key][1] / table[key][0] for key in common))
        observed_mark.append(math.fsum(row[2] for row in table.values()) / total)
    observed_delta = (observed_mark[1] - observed_mark[0]) / delta_cos4
    conditional_delta = (std_mark[1] - std_mark[0]) / delta_cos4
    vector = [
        std_mark[0], std_mark[1], conditional_delta,
        (std_exit[1] - std_exit[0]) / delta_cos4,
        observed_delta, observed_delta - conditional_delta,
    ]
    return vector, {
        "common_age_line_strata": len(common),
        "common_support_coverage": {"first": coverage[0], "second": coverage[1]},
        "target_distribution": "equal-orientation mixture of risk distributions on common (ell,age) support",
        "observed_marked_hazard": {"first": observed_mark[0], "second": observed_mark[1]},
        "standardized_marked_hazard": {"first": std_mark[0], "second": std_mark[1]},
        "standardized_unmarked_hazard": {"first": std_exit[0], "second": std_exit[1]},
    }


def score_mark_source(
    source_id: str,
    n: int,
    k0: int,
    batches: dict[str, list[dict[tuple[int, int, int], list[int]]]],
    chi4: Mapping[str, float],
    role: str,
) -> dict[str, Any]:
    delta_cos4 = chi4["second"] - chi4["first"]

    def evaluate(omitted: int | None) -> tuple[list[float], dict[str, Any]]:
        tables = {orientation: merge_mark(batches[orientation], omitted) for orientation in state.ORIENTATIONS}
        return standardized_mark_vector(tables["first"], tables["second"], delta_cos4)

    point, detail = evaluate(None)
    leave = [evaluate(index)[0] for index in range(len(batches["first"]))]
    cov = covariance(leave)
    se = np.sqrt(np.maximum(np.diag(cov), 0.0))
    statistic = point[2] / se[2]
    return {
        "id": source_id,
        "N": n,
        "k0": k0,
        "role": role,
        "delta_cos4": delta_cos4,
        "vector_order": [
            "first_age_standardized_mark12_h4_hazard",
            "second_age_standardized_mark12_h4_hazard",
            "angular_age_standardized_mark12_h4_conditional_hazard",
            "angular_age_standardized_unmarked_completion_hazard",
            "angular_observed_mark12_h4_hazard",
            "angular_mark12_h4_risk_composition_remainder",
        ],
        "vector": point,
        "standard_error": se.tolist(),
        "delete_one_covariance": cov.tolist(),
        "delete_one_unit": "same batch removed from both orientations",
        "conditional_mark12_test": {
            "student_t": float(statistic),
            "degrees_of_freedom": len(batches["first"]) - 1,
            "two_sided_p": float(2.0 * t.sf(abs(statistic), len(batches["first"]) - 1)),
            "resolved_at_alpha_0_01": bool(2.0 * t.sf(abs(statistic), len(batches["first"]) - 1) < 0.01),
        },
        "standardization": detail,
    }


def main() -> None:
    manifest = yaml.safe_load((ROOT / "analysis/p337_birth_state_current_manifest.yaml").read_text())
    state_runs = [score_state_run(ROOT, run) for run in manifest["runs"]]

    mark_scores = []
    freeze = json.loads((ROOT / "analysis/p334_birth_age_production_freeze.json").read_text())
    for source_id in ("N325", "N425"):
        spec = freeze["inputs"][source_id]
        path = P334_RAW_ROOT / f"{source_id}_2m/{source_id}_2m.marked_births.csv"
        batches, chi4 = parse_mark_archive(path, int(spec["N"]), int(spec["k0"]), int(spec["batches"]))
        score = score_mark_source(
            source_id, int(spec["N"]), int(spec["k0"]), batches, chi4,
            "independent_P334_production_age_block_not_four_generation_replica",
        )
        score["source"] = {
            "path": str(path),
            "sha256": spec["sha256"],
            "hash_verification": "reused_from_p334_frozen_production_score",
        }
        mark_scores.append(score)

    n130 = next(run for run in manifest["runs"] if run["id"] == "N130-control")
    n130_state = next(run for run in state_runs if run["id"] == "N130-control")
    with state.git_blob(ROOT, n130["archive_commit"], n130["path"], n130["sha256"]) as path:
        batches, chi4 = parse_mark_archive(
            path, int(n130["N"]), int(math.floor(n130["N"] * n130_state["p"] + 0.5)), int(n130["batches"])
        )
    mark_scores.append(score_mark_source(
        "N130-control", int(n130["N"]), int(math.floor(n130["N"] * n130_state["p"] + 0.5)),
        batches, chi4, "independent_low_statistics_cross_lineage_control",
    ))

    primary = [run for run in state_runs if run["role"] == "four_generation_primary"]
    result = {
        "schema": SCHEMA,
        "status": "completed_zero_new_sample_mechanism_split",
        "new_samples": False,
        "issues": [337, 334],
        "state_current_midpoint_splits": state_runs,
        "age_standardized_mark12_h4": mark_scores,
        "decision": {
            "primary_absolute_dominance_by_size": {
                run["id"]: run["dominant_absolute_component"] for run in primary
            },
            "primary_mark12_status": "not_scoreable_in_all_four_primary_archives",
            "point_direction": (
                "all_four_primary_generations_have_positive_angular_conditional_hazard_and_"
                "negative_angular_risk_composition_so_composition_partially_cancels_the_completion_H4"
            ),
            "interpretation_rule": (
                "Each generation is scored as its own covariance block; no cross-generation vote, pooled replicate, "
                "shared-amplitude fit or cross-lineage substitution is made."
            ),
        },
        "claim_boundary": [
            "The midpoint identity is an exact algebraic attribution on observed (primitive line,current-age layer) cells.",
            "Off-support hazard is not identified; assigning unmatched support entirely to composition avoids fabricating a hazard contrast and its weighted mass is reported.",
            "N85/N170/N340/N680 contain no mark12 field; N130 and P334 N325/N425 mark scores are controls, not replacements or extra generation votes.",
            "No field, exponent, intrinsic memory variable or continuum mechanism is identified.",
        ],
    }
    output = ROOT / "results/p337-risk-hazard-split/latest.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["decision"], indent=2))


if __name__ == "__main__":
    main()
