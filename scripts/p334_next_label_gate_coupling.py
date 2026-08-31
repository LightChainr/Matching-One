#!/usr/bin/env python3
"""Fixed 01/10 gate/clock half-difference products from completed quartet CSVs."""
import argparse
import csv
import gzip
from hashlib import sha256
import io
from itertools import groupby
import json
from pathlib import Path
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE_DIRECTORY = "results/p334-nested-next-label-forks"
NORMALIZATION_COMMIT = "bb79fd47"
NORMALIZATION_PATH = "results/p334-nine-layer-complete-ae/batch_vectors.json"
FORMULA_COMMIT = "47d6eb418694fd1941e612bff13cff8f412e34ea"
RUNNER_COMMIT = "a3249a598e7af6db940fd602df31533ab5c13d38"
HEADER = ("N,batch,counter,k0,first_rank,second_rank,quartet,group,replica,"
          "next_label,first_next_rank,second_next_rank,first_k1,first_k2,second_k1,second_k2").split(",")
CELLS = ("01", "10")
FIELDS = (
    "prevalence", "mean_gate_R0_first", "mean_gate_R1_second",
    "same_label_joint_gate", "gate_cov_R0_R0", "gate_cov_R0_R1", "gate_cov_R1_R1",
    "gate_R0_cov_A", "gate_R0_cov_E", "gate_R1_cov_A", "gate_R1_cov_E",
    "paired_mean_A", "paired_mean_E",
)


def quartet_vector(ranks, next_ranks, births, n, delta_cos4):
    """Return (cell, 13-vector), or (None, zeros) outside 01/10.

    ranks: original (first,second) ranks.
    next_ranks: [label U/V, original orientation first/second].
    births: [label U/V, suffix 0/1, original orientation first/second, K1/K2].
    These are the common next-label ranks, not terminal ranks or future marks.
    """
    ranks = np.asarray(ranks, dtype=int)
    cell = "".join(str(x) for x in ranks)
    if cell not in CELLS:
        return None, np.zeros(len(FIELDS))
    next_ranks = np.asarray(next_ranks, dtype=int)
    births = np.asarray(births, dtype=float)
    if next_ranks.shape != (2, 2) or births.shape != (2, 2, 2, 2):
        raise ValueError("A quartet needs two labels by two paired suffixes")
    low, high = (0, 1) if cell == "01" else (1, 0)
    g = np.column_stack((next_ranks[:, low] >= 1, next_ranks[:, high] == 2)).astype(float)
    a = 1-(births[..., 0]+births[..., 1])/(n+1)
    e = 1-(births[..., 1]-births[..., 0])/(n+1)
    observations = np.stack(((a[..., 0]-a[..., 1])/delta_cos4,
                             (e[..., 0]-e[..., 1])/delta_cos4), axis=-1)
    m = observations.mean(axis=1)
    dg, dm = g[0]-g[1], m[0]-m[1]
    gate_cov = .5*np.outer(dg, dg)
    response_cov = .5*np.outer(dg, dm)
    return cell, np.r_[1., g.mean(axis=0), np.prod(g, axis=1).mean(),
                       gate_cov[0, 0], gate_cov[0, 1], gate_cov[1, 1],
                       response_cov.ravel(), m.mean(axis=0)]


def read_git(commit, path):
    return subprocess.check_output(["git", "show", f"{commit}:{path}"], cwd=ROOT)


def batch_vector(blob, delta_by_n):
    reader = csv.DictReader(io.TextIOWrapper(gzip.GzipFile(fileobj=io.BytesIO(blob))))
    if reader.fieldnames != HEADER:
        raise ValueError("Production header differs from the fixed producer schema")
    batches = np.zeros((2, len(FIELDS)))
    prefix_states, prefix_quartets = {}, {}
    identity = None
    for (counter, quartet), group in groupby(reader, key=lambda r: (int(r["counter"]), int(r["quartet"]))):
        rows = [{key: int(value) for key, value in r.items()} for r in group]
        if len(rows) != 4 or {(r["group"], r["replica"]) for r in rows} != {(0,0), (0,1), (1,0), (1,1)}:
            raise ValueError("A quartet must retain exactly two labels by two suffixes")
        row = rows[0]
        current_identity = (row["N"], row["batch"])
        identity = current_identity if identity is None else identity
        if any((r["N"], r["batch"]) != identity for r in rows):
            raise ValueError("A raw batch mixes original N/batch identities")
        prefix = (row["k0"], row["first_rank"], row["second_rank"])
        if any((r["k0"], r["first_rank"], r["second_rank"]) != prefix for r in rows):
            raise ValueError("The quartet does not share one original prefix")
        if counter in prefix_states and prefix_states[counter] != prefix:
            raise ValueError("Original checkpoint rank changed between quartets")
        prefix_states[counter] = prefix
        if quartet in prefix_quartets.setdefault(counter, set()):
            raise ValueError("Repeated original-counter quartet")
        prefix_quartets[counter].add(quartet)
        next_ranks = np.zeros((2, 2), dtype=int)
        births = np.zeros((2, 2, 2, 2))
        for label_group in (0, 1):
            tails = [r for r in rows if r["group"] == label_group]
            if len({(r["next_label"], r["first_next_rank"], r["second_next_rank"]) for r in tails}) != 1:
                raise ValueError("Next label/rank differs between its independent suffixes")
            next_ranks[label_group] = (tails[0]["first_next_rank"], tails[0]["second_next_rank"])
            for tail in tails:
                births[label_group, tail["replica"]] = ((tail["first_k1"], tail["first_k2"]),
                                                        (tail["second_k1"], tail["second_k2"]))
        cell, vector = quartet_vector(prefix[1:], next_ranks, births, identity[0], delta_by_n[identity[0]])
        if cell is not None:
            batches[CELLS.index(cell)] += vector
    if identity is None or len(prefix_states) != 1000 or any(q != set(range(8)) for q in prefix_quartets.values()):
        raise ValueError("The full original 1000-prefix by 8-quartet batch is incomplete")
    return identity, batches.ravel()/8000, {
        cell: sum("".join(str(x) for x in state[1:]) == cell for state in prefix_states.values())
        for cell in CELLS}


def covariance_of_mean(batches):
    centered = batches-batches.mean(axis=0)
    return centered.T@centered/(20*19)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-commit", required=True,
                        help="Producer's completed immutable 40-gzip source commit, never a partial directory")
    parser.add_argument("--output", type=Path, default=ROOT/"results/p334-next-label-gate-coupling")
    args = parser.parse_args()
    source = subprocess.check_output(["git", "rev-parse", args.source_commit+"^{commit}"], cwd=ROOT, text=True).strip()
    norm_blob = read_git(NORMALIZATION_COMMIT, NORMALIZATION_PATH)
    normalization = json.loads(norm_blob)
    deltas = {int(n): row["delta_cos4"] for n, row in normalization["sizes"].items()}
    paths = subprocess.check_output(["git", "ls-tree", "-r", "--name-only", source, SOURCE_DIRECTORY], cwd=ROOT, text=True).splitlines()
    paths = [path for path in paths if path.endswith(".gz")]
    if len(paths) != 40:
        raise ValueError("Exactly forty completed original-batch gzip files are required")
    batches, counts, hashes = {}, {}, {}
    for path in paths:
        blob = read_git(source, path)
        identity, vector, cell_counts = batch_vector(blob, deltas)
        if identity in batches:
            raise ValueError("Duplicate original N/batch")
        batches[identity], counts[identity], hashes[path] = vector, cell_counts, sha256(blob).hexdigest()
    if set(batches) != {(n, b) for n in (325, 425) for b in range(20)}:
        raise ValueError("Source must retain both sizes and the original twenty batch IDs")
    labels = [f"cell.{cell}.{field}" for cell in CELLS for field in FIELDS]
    sizes = {}
    for n in (325, 425):
        raw = np.asarray([batches[n, b] for b in range(20)])
        point, cov = raw.mean(axis=0), covariance_of_mean(raw)
        old_sum = raw.reshape(20, 2, len(FIELDS)).sum(axis=1)
        # Only X-bearing fields change sign under first/second -> low/high.
        x_fields = [FIELDS.index(field) for field in FIELDS if "_cov_A" in field or "_cov_E" in field or field.startswith("paired_mean_")]
        low_high = deltas[n]*(raw[:, :len(FIELDS)]-raw[:, len(FIELDS):])[:, x_fields]
        sizes[str(n)] = {"delta_cos4": deltas[n], "batch_ids": list(range(20)),
            "batch_denominators": [1000]*20, "quartets_per_prefix": 8,
            "cell_counts_by_original_batch": [counts[n, b] for b in range(20)],
            "labels": labels, "joint_20_batch_means": raw.tolist(),
            "joint_LOO_vectors": ((20*point-raw)/19).tolist(),
            "estimate": point.tolist(), "se": np.sqrt(np.diag(cov)).tolist(),
            "joint_covariance": cov.tolist(),
            "pooled_01_10": {"labels": list(FIELDS), "joint_20_batch_means": old_sum.tolist(),
                "estimate": old_sum.mean(axis=0).tolist(),
                "se": np.sqrt(np.diag(covariance_of_mean(old_sum))).tolist(),
                "X_orientation": "first-minus-second divided by the original delta_cos4"},
            "low_minus_high_unscaled_X": {"labels": [FIELDS[i] for i in x_fields],
                "joint_20_batch_means": low_high.tolist(), "estimate": low_high.mean(axis=0).tolist(),
                "se": np.sqrt(np.diag(covariance_of_mean(low_high))).tolist(),
                "definition": "delta_cos4 times (01-minus-10); raw low-minus-high A/E, not a new H4 coefficient"}}
    result = {"schema": "matching-one/p334-next-label-gate-coupling/v1", "source_commit": source,
        "source_directory": SOURCE_DIRECTORY, "source_sha256": hashes,
        "formula_freeze_commit": FORMULA_COMMIT, "producer_code_commit": RUNNER_COMMIT,
        "normalization": {"commit": NORMALIZATION_COMMIT, "path": NORMALIZATION_PATH, "sha256": sha256(norm_blob).hexdigest()},
        "original_full_birth_commit": normalization["full_birth_commit"],
        "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "sizes": sizes, "new_samples_by_this_scorer": 0,
        "boundary": "Mass-weighted original-prefix conditional covariances. Two cells, all readouts and other providers share the same original prefixes and fresh quartet block. No population rate-product replacement, covariance inverse, unknown-pi Doob fraction, or universal location coupling."}
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/"batch_vectors.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# Shared next-label birth gates and complete A/E", ""]
    for n, row in sizes.items():
        lines += [f"## N{n}", "", "| Full-population mass / covariance | 01+10 estimate | shared-batch SE |", "|---|---:|---:|"]
        for label, value, se in zip(FIELDS, row["pooled_01_10"]["estimate"], row["pooled_01_10"]["se"]):
            lines.append(f"| {label} | {value:.10g} | {se:.6g} |")
        lines += ["", "Raw low-minus-high A/E projection (distinct from old first-minus-second/H4):", ""]
        for label, value, se in zip(row["low_minus_high_unscaled_X"]["labels"], row["low_minus_high_unscaled_X"]["estimate"], row["low_minus_high_unscaled_X"]["se"]):
            lines.append(f"- {label}: {value:.10g} +/- {se:.6g}")
        lines.append("")
    lines += [result["boundary"], ""]
    (args.output/"REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
