#!/usr/bin/env python3
"""Score the exact conditional mean on the SAME fixed 140k source prefixes.

Reuse the committed original scorer for validation, cells and overlap weights;
only actual original probes enter its Moments. Exact edge sums stay integers.
"""
import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import importlib.util
from itertools import zip_longest
import json
import math
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent.parent
OLD = ROOT.parent / "safe-insertion-production-20260929"
OLD_COMMIT = "2a1677e1efe1137dffc07c3c4ae4ceebf15ecfd3"
METRICS = ("RB_deltaY", "original_probe_deltaY", "paired_RB_minus_probe",
           "RB_deltaE", "RB_deltaZ2", "RB_earlyY", "RB_lateY", "secondary_deltaNu")
WARNINGS = [
    "Estimator refinement after the original readout; same 140k prefixes and batch blocks, zero new samples.",
    "No independent confirmation: reproducing old primary/secondary is source/estimand identity, not scientific replication.",
    "Exact means conditional mean 2e/(m-c); finite-prefix uncertainty remains. Four actual old probes are averaged per prefix.",
    "SE_RB/SE_probe is descriptive; variance reduction does not guarantee empirical delete-batch jackknife SE shrinkage.",
    "No automatic Markov claim: one weighted moment can cancel across cells; finite L512 count clock only, no continuum claim.",
]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def git(*args):
    return subprocess.check_output(["git", "-C", str(REPO), *args])


def reference(path):
    raw = path.read_bytes()
    return {"path": str(path), "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}


def aligned_batch(old, original, replay):
    identity, old_path, old_hash, old_entry = original
    new_id, new_path, new_hash, new_entry = replay
    require(identity == new_id and old_entry["seed"] == new_entry["seed"], "batch/seed mismatch")
    require(new_entry["original_file"] == old_entry["file"] and
            new_entry["dependency_group"] == old_entry["dependency_group"], "source block mismatch")
    refs = [reference(p) for p in (old_path, new_path)]
    require([r["sha256"] for r in refs] == [old_hash, new_hash], f"batch {identity}: hash mismatch")
    stats, edges, candidates = old.BatchStats(), {}, 0
    m = old.FIXED["N"] - old.FIXED["b"]
    with gzip.open(old_path, "rt", newline="", encoding="utf-8-sig") as fo, \
            gzip.open(new_path, "rt", newline="", encoding="utf-8-sig") as fn:
        ro, rn = csv.reader(fo, strict=True), csv.reader(fn, strict=True)
        require(next(ro, None) == list(old.COLUMNS), f"batch {identity}: old header mismatch")
        require(next(rn, None) == list(old.COLUMNS[:5]) + ["synergy_edges", "nonlocal_candidates"],
                f"batch {identity}: replay header mismatch")
        for line, (o, n) in enumerate(zip_longest(ro, rn), 2):
            label = f"batch {identity}, CSV line {line}"
            require(o is not None and n is not None, f"{label}: unequal row counts")
            require(len(o) == 9 and len(n) == 7 and o[:5] == n[:5], f"{label}: prefix mismatch")
            values = [old.integer(v, label) for v in o]
            stats.observe(values, old.FIXED)
            j1, rank, dx, dy, c = values[:5]
            e, candidate = [old.integer(v, label) for v in n[5:]]
            if rank != 1:
                require((e, candidate) == (-1, -1), f"{label}: outside-risk sentinel mismatch")
                continue
            s = m - c
            require(0 <= e <= s * (s - 1) // 2 and candidate >= 0, f"{label}: invalid pair counts")
            require(s > 1 or (e, candidate) == (0, 0), f"{label}: impossible pair counts")
            candidates += candidate
            if not s:
                continue
            if dx < 0 or (dx == 0 and dy < 0):
                dx, dy = -dx, -dy
            edges.setdefault((dx, dy, c), [0, 0])[int(j1 > old.FIXED["a"])] += e
    require(stats.samples == 10000, f"batch {identity}: expected exactly 10000 aligned rows")
    return (stats, edges), {"batch": identity, "rows": stats.samples, "all_first5_identical": True,
                           "seed": old_entry["seed"], "dependency_group": old_entry["dependency_group"],
                           "sum_nonlocal_candidates": candidates, "original": refs[0], "replay": refs[1]}


def estimate(old, blocks, omit=None, detailed=False):
    chosen = [block for i, block in enumerate(blocks) if i != omit]
    total = old.aggregate([block[0] for block in chosen])
    _, primary, secondary = old.estimate(total, old.FIXED, detailed=True)
    require(primary["status"] == secondary["status"] == "scoreable", f"absent support, deletion={omit}")
    edges = {}
    for _, cells in chosen:
        for key, sums in cells.items():
            target = edges.setdefault(key, [0, 0])
            for g in range(2):
                target[g] += sums[g]
    m = old.FIXED["N"] - old.FIXED["b"]
    terms, rows = [[] for _ in range(7)], []
    for cell in primary.pop("cells"):
        key = (*cell["direction"], cell["nu_b"])
        s, means, cohorts = m - key[2], [], {}
        for g, name in enumerate(old.GROUPS):
            moment, sum_e = cell[name], edges[key][g]
            n, sum_probe = moment["prefixes"], moment["sum_probe_Y"]
            mean_e = sum_e / n if n else None
            mean_y = 2 * sum_e / (s * n) if n else None
            means.append((mean_e, mean_y))
            cohorts[name] = dict(n=n, sum_e=sum_e, sum_old_probe_Y=sum_probe,
                                 mean_e=mean_e, meanY=mean_y, mean_old_prefixprobe_Y=moment["mean_Ybar"])
        record = {k: cell[k] for k in ("direction", "nu_b", "common_support", "overlap_weight", "normalized_weight")}
        record.update(safe_sites=s, **cohorts)
        if cell["common_support"]:
            de, dy = (means[0][j] - means[1][j] for j in range(2))
            dp = cell["deltaY"]  # Actual four-probe prefix mean, unchanged original arithmetic.
            vector = [dy, dp, dy - dp, de, -2 * de / (m * (m - 1)), means[0][1], means[1][1]]
            for target, value in zip(terms, vector):
                target.append(cell["overlap_weight"] * value)
            record.update(zip(METRICS[:5], vector[:5]))
        rows.append(record)
    vector = [math.fsum(t) / primary["overlap_weight_sum"] for t in terms] + [secondary["delta"]]
    vector[2], vector[4] = vector[0] - vector[1], -2 * vector[3] / (m * (m - 1))
    primary.update(delta=vector[0], deltaE=vector[3], deltaZ2=vector[4], original_probe_deltaY=vector[1],
                   weighted_early_mean=vector[5], weighted_late_mean=vector[6])
    secondary_rows = secondary.pop("cells")
    if detailed:
        primary["cells"] = rows
        secondary["cells"] = [{**{k: r[k] for k in ("direction", "common_support", "overlap_weight", "normalized_weight", "deltaNu")},
                               **{g: dict(n=r[g]["prefixes"], sumNu=r[g]["sum_nu_b"], meanNu=r[g]["mean_nu_b"])
                                  for g in old.GROUPS}} for r in secondary_rows]
    return vector, primary, secondary, total.totals(4)


def markdown(result):
    lines = ["# Same-prefix exact completion-pair refinement", "",
             "All 140000 rows aligned across 14 original batches; every first-five-field tuple matched.", "",
             "Early minus late; fixed a,b and original overlap weights. No new random prefixes.", "",
             "| Metric | Estimate | Batch jackknife SE |", "|---|---:|---:|"]
    for key in METRICS:
        lines.append(f"| {key} | {result['point_estimates'][key]:.15g} | {result['standard_errors'][key]:.15g} |")
    lines += ["", f"Descriptive SE_RB/SE_probe: {result['refinement']['SE_RB_over_SE_probe']}. "
              f"Variance ratio: {result['refinement']['variance_RB_over_probe']}.", "",
              "| Support / cohort | Rank-one | Eligible | Supported |", "|---|---:|---:|---:|"]
    for name in ("primary", "secondary"):
        for group in ("early", "late", "combined"):
            c = result[name]["coverage"][group]
            lines.append(f"| {name} / {group} | {c['rank1_prefixes']} | {c['eligible_prefixes']} | {c['supported_prefixes']} |")
    lines += ["", "Original committed primary and secondary point/SE reproduced (absolute tolerances 1e-14 / 1e-12).",
              "This is a source/estimand identity check, not scientific replication.", "",
              "JSON retains exact integer cell sums, mean_e, meanY, all 14 aligned deletion vectors/supports, "
              "full 8x8 covariance, and input/source hashes and commits. First seven metrics share exact (D,c) cells; "
              "the secondary uses all rank-one prefixes within D, including c=m.", "",
              "Y=2e/(m-c); original prefix probe mean=(Y0+Y1+Y2+Y3)/4. "
              "RB_deltaE averages cell mean-edge differences; RB_deltaZ2=-2*RB_deltaE/[m(m-1)]. "
              "c=m has no conditional safe mean. Singular covariance is expected; no inverse is used.", ""]
    return "\n".join(lines + [f"- {w}" for w in WARNINGS]) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", type=Path, default=ROOT / "data")
    ap.add_argument("--output-dir", type=Path, default=ROOT / "results")
    args = ap.parse_args()
    output = args.output_dir.resolve()
    require(ROOT in output.parents, "outputs must be inside this new analysis directory")
    destinations = [output / "result.json", output / "RESULT.md"]
    require(not any(p.exists() for p in destinations), "refusing to overwrite outputs; choose a fresh output directory")
    for path in (OLD / "analyze.py", OLD / "results/result.json", OLD / "contract.json", OLD / "data/run.json"):
        require(path.read_bytes() == git("show", f"{OLD_COMMIT}:{path.relative_to(REPO)}"), f"original changed: {path}")
    sys.dont_write_bytecode = True
    spec = importlib.util.spec_from_file_location("original_safe_insertion_score", OLD / "analyze.py")
    old = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = old
    spec.loader.exec_module(old)
    contract, contract_ref = old.read_json(OLD / "contract.json")
    old.check_contract(contract)
    original, original_ref = old.read_json(OLD / "data/run.json")
    run, run_ref = old.read_json(args.data_dir.resolve() / "run.json")
    old.check_contract(run["contract"], "replay.contract")
    require(run["contract"] == original["contract"] == contract, "contract changed")
    require(run["new_random_prefixes"] == 0 and run["original_manifest_sha256"] == original_ref["sha256"], "wrong source manifest")
    original_entries = old.batch_entries(original, OLD / "data", old.FIXED)
    replay_entries = old.batch_entries(run, args.data_dir.resolve(), old.FIXED)
    require([v[0] for v in original_entries] == [v[0] for v in replay_entries] == list(range(14)), "require original 14 batch IDs")
    blocks, manifests = [], []
    for o, n in zip(original_entries, replay_entries):
        block, manifest = aligned_batch(old, o, n)
        blocks.append(block)
        manifests.append(manifest)
    point, primary, secondary, totals = estimate(old, blocks, detailed=True)
    deletions = []
    for i in range(14):
        v, p, s, t = estimate(old, blocks, omit=i)
        deletions.append(dict(deleted_batch=manifests[i]["batch"], vector=v, primary_support=p, secondary_support=s, totals=t))
    vectors = [r["vector"] for r in deletions]
    mean = [math.fsum(v[k] for v in vectors) / 14 for k in range(8)]
    covariance = [[13 / 14 * math.fsum((v[j] - mean[j]) * (v[k] - mean[k]) for v in vectors)
                   for k in range(8)] for j in range(8)]
    se = [math.sqrt(covariance[k][k]) for k in range(8)]
    committed, committed_ref = old.read_json(OLD / "results/result.json")
    require(committed["status"] == "scored" and totals == committed["totals"], "original total identity failed")
    checks = {}
    for index, key in ((1, "deltaY"), (7, "secondary_deltaNu")):
        dp, ds = point[index] - committed["point_estimates"][key], se[index] - committed["standard_errors"][key]
        require(abs(dp) <= 1e-14 and abs(ds) <= 1e-12, f"original {key} point/SE identity failed: {dp}, {ds}")
        checks[key] = dict(committed_point=committed["point_estimates"][key], committed_SE=committed["standard_errors"][key],
                           reproduced_point=point[index], reproduced_SE=se[index], point_error=dp, SE_error=ds)
    ratio = se[0] / se[1] if se[1] else None
    sources = [ROOT / n for n in ("analyze.py", "engine.cpp", "run.py")]
    sources += [OLD / n for n in ("analyze.py", "engine.cpp", "run.py")]
    geometry = ROOT.parent / "completion-hazard-production-20260929/engine.cpp"
    require(reference(geometry)["sha256"] == run["included_geometry_sha256"] == original["included_geometry_sha256"], "geometry changed")
    for name, digest in run["source_sha256"].items():
        require(reference(ROOT / name)["sha256"] == digest, f"replay source changed: {name}")
    result = dict(schema="matching-one.exact-pair-score.v1", status="scored", parameters=old.FIXED,
                  created_utc=datetime.now(timezone.utc).isoformat(), m=old.FIXED["N"] - old.FIXED["b"],
                  metric_order=METRICS, point_vector=point, point_estimates=dict(zip(METRICS, point)),
                  standard_errors=dict(zip(METRICS, se)), primary=primary, secondary=secondary, totals=totals,
                  original_identity=dict(status="passed", point_tolerance=1e-14, SE_tolerance=1e-12, checks=checks),
                  refinement=dict(SE_RB_over_SE_probe=ratio, variance_RB_over_probe=ratio ** 2 if ratio is not None else None),
                  joint_uncertainty=dict(method="aligned whole-batch delete-one; support and weights recomputed",
                      covariance_formula="(B-1)/B * sum_j (theta_(-j)-mean)(theta_(-j)-mean)^T", batch_count=14,
                      metric_order=METRICS, covariance=covariance, delete_one_mean_vector=mean, delete_one=deletions,
                      point_estimator="pooled plug-in; no jackknife bias correction"),
                  input_manifest=dict(contract=contract_ref, original_run=original_ref, replay_run=run_ref,
                                      committed_original_result=committed_ref, aligned_batches=manifests),
                  source_manifest=dict(original_result_commit=OLD_COMMIT, original_acquisition_commit=original["source_commit"],
                      replay_source_commit=run["source_commit"], scoring_head_commit=git("rev-parse", "HEAD").decode().strip(),
                      files=[reference(p) for p in sources + [geometry, REPO / "notes/completion-pair-synergy-20260929.md"]],
                      note="Commit references and working-file hashes are distinct; replay code may be uncommitted."),
                  warnings=WARNINGS)
    encoded = json.dumps(result, indent=2, allow_nan=False) + "\n"
    prose = markdown(result)
    output.mkdir(parents=True, exist_ok=True)
    for path, content in zip(destinations, (encoded, prose)):
        with path.open("x", encoding="utf-8") as handle:
            handle.write(content)
    print(json.dumps({"status": "scored", "outputs": [str(p) for p in destinations]}))


if __name__ == "__main__":
    try:
        main()
    except (OSError, EOFError, ValueError, KeyError, csv.Error, subprocess.CalledProcessError) as error:
        sys.exit(f"exact-pair scoring error: {error}")
