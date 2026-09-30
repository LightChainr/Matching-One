#!/usr/bin/env python3
"""One retrospective directional Markov-necessity check; no new samples.

Reuse the saved nu_b, unchanged cohorts/cutoffs, exact direction and overlap
weights from the completed block. Do not mutate its scorer or results.
"""
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import platform
import sys
import time


HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / "completion-hazard-production-20260929"
SPEC = importlib.util.spec_from_file_location("completion_parent", PARENT / "analyze.py")
BASE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = BASE
SPEC.loader.exec_module(BASE)
NAMES = ["pooled_nu_early_minus_late", "within_D_nu_early_minus_late", "within_D_survival_early_minus_late"]


def score(stats):
    early, late = stats.groups
    if not early.n or not late.n:
        raise ValueError("Both risk cohorts are required")
    denominator, numerator, survival = 0.0, 0.0, 0.0
    support = [0, 0]
    rows = []
    for direction, (e, l) in sorted(stats.directions.items()):
        weight = e.n*l.n/(e.n+l.n) if e.n and l.n else 0.0
        if weight:
            denominator += weight
            numerator += weight*(e.nu/e.n-l.nu/l.n)
            survival += weight*(e.z/e.n-l.z/l.n)
            support[0] += e.n
            support[1] += l.n
        rows.append({"direction": list(direction), "early": e.report(), "late": l.report(),
                     "weight": weight,
                     "mean_nu_difference": e.nu/e.n-l.nu/l.n if weight else None})
    if not denominator:
        raise ValueError("No common direction support")
    return [early.nu/early.n-late.nu/late.n, numerator/denominator, survival/denominator], {
        "supported_early": support[0], "supported_late": support[1],
        "early_total": early.n, "late_total": late.n, "overlap_weight_sum": denominator,
        "direction_rows": rows,
    }


def main():
    started = time.perf_counter()
    contract = json.loads((PARENT / "contract.json").read_text())
    run = json.loads((PARENT / "data/run.json").read_text())
    assert run["status"] == "completed"
    assert len(run["batches"]) == contract["batches"]
    assert all(b["samples"] == contract["samples_per_batch"] for b in run["batches"])
    blocks, provenance = [], []
    for entry in run["batches"]:
        stats, ref = BASE.read_batch((PARENT / "data").resolve(), entry, contract)
        blocks.append(stats); provenance.append(ref)
    pooled, detail = score(BASE.aggregate(blocks))
    deletes = [score(BASE.aggregate(blocks, omit=i))[0] for i in range(len(blocks))]
    B = len(blocks)
    center = [math.fsum(row[j] for row in deletes)/B for j in range(len(NAMES))]
    covariance = [[(B-1)/B*math.fsum((row[j]-center[j])*(row[k]-center[k]) for row in deletes)
                   for k in range(len(NAMES))] for j in range(len(NAMES))]
    se = [math.sqrt(covariance[j][j]) for j in range(len(NAMES))]
    # The already-published overlap survival statistic must retain its value:
    # we are replacing only its outcome by nu_b for the new necessity check.
    original = json.loads((PARENT / "results/result.json").read_text())
    assert abs(pooled[2]-original["point"]["scalars"]["delta_overlap_D"]["estimate"]) < 1e-14
    assert abs(se[2]-original["jackknife"]["uncertainty"]["delta_overlap_D"]["standard_error"]) < 1e-14
    step_factor = 1/(contract["N"]-contract["b"])
    window_factor = (contract["c"]-contract["b"])*step_factor
    derived = {
        "within_D_next_insertion_exit_difference": {"estimate": pooled[1]*step_factor, "se": se[1]*step_factor},
        "within_D_window_scaled_exit_difference": {"estimate": pooled[1]*window_factor, "se": se[1]*window_factor},
    }
    result = {
        "schema": "matching-one.directional-hazard-contrast.v1", "date": "2026-09-29",
        "status": "retrospective_single_necessary_condition_not_new_validation",
        "source_commit": "87f6644d3b57c7dbd67c44c63bcceffb6b7693ac",
        "source_manifest": "analysis/completion-hazard-production-20260929/data/run.json",
        "source_manifest_sha256": hashlib.sha256((PARENT / "data/run.json").read_bytes()).hexdigest(),
        "source_reader_sha256": hashlib.sha256((PARENT / "analyze.py").read_bytes()).hexdigest(),
        "samples": sum(block.samples for block in blocks),
        "cutoffs": [contract[k] for k in ("a", "b", "c")], "N": contract["N"],
        "question": "Does direction-only Markov closure pass its exact next-insertion hazard necessity?",
        "null": "For every direction d, E[nu_b|early,d]=E[nu_b|late,d] under marked-rank Markovness",
        "selection": "Question chosen after seeing the published finite-lag readout; no time, threshold, direction subset or feature scan",
        "weights": "same exact-D overlap weights e*l/(e+l), recomputed in every whole-batch deletion",
        "scalar_order": NAMES, "point": dict(zip(NAMES, pooled)), "se": dict(zip(NAMES, se)),
        "covariance": covariance, "delete_one_vectors": deletes, "derived": derived,
        "support": detail, "batches": provenance,
        "boundaries": [
            "Raw counts, next-insertion probability and window-scaled intensity have different units.",
            "Within-(D,nu_b) immediate nu equality is tautological and is not tested as closure evidence.",
            "Same 140k block, not independent corroboration; original prospective estimands and results unchanged.",
            "One aggregate zero does not prove sectorwise equality; nonzero population weighted contrast contradicts the necessary null.",
            "Batch SE is approximate; no multiplicity-adjusted certification or continuum extrapolation.",
        ],
        "python": platform.python_version(), "wall_seconds": time.perf_counter()-started,
    }
    (HERE / "result.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# Directional next-insertion hazard: one retrospective check", "",
             "Same completed 140k block, unchanged cutoffs, cohorts and exact-D overlap weights.", "",
             "| Quantity | Estimate | Aligned batch SE |", "|---|---:|---:|"]
    lines += [f"| {name} | {pooled[j]:.10g} | {se[j]:.10g} |" for j, name in enumerate(NAMES)]
    lines += [f"| {name} | {value['estimate']:.10g} | {value['se']:.10g} |" for name, value in derived.items()]
    lines += ["", "Positive nu contrast means greater early-cohort next-insertion exit probability.",
              "Window scaling is not a finite-window exit estimate. The last survival contrast is the unchanged original result.",
              "", "This necessity check was selected after the original report; it is not independent validation.",
              "No new observations, subgroup searches, thresholds, cutoffs or samples are added.", ""]
    (HERE / "RESULT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
