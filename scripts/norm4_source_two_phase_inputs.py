"""Load the independent q/E complement of the marked norm-4 production subset.

This reader performs no simulation or scoring.  ``marked_profiles`` is the
array returned by ``analyze_norm4_source_thermal.read_raw``: its first two
coordinates are exact integer q/E sums; its three source coordinates may
already have been divided by N and are not used here.  Dependency grouping
remains the caller's responsibility under the source-chain manifest.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import subprocess
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE_COMMIT = "8b26a30a785bc142a9d17bfed99a8d0e98ddc4dc"
ARCHIVE_ROOT = "results/server-20260829/P154-norm4-production/raw"
ORIENTATIONS = ("first", "second")
DESIGNS = {
    65: ((8, 1), (7, 4)),
    85: ((9, 2), (7, 6)),
    130: ((11, 3), (9, 7)),
    170: ((13, 1), (11, 7)),
    260: ((16, 2), (14, 8)),
    340: ((18, 4), (14, 12)),
}
MARKED_BATCHES = 100
MARKED_PER_BATCH = 1000
MARKED_SAMPLES = MARKED_BATCHES * MARKED_PER_BATCH


def _git_bytes(path: str) -> bytes:
    return subprocess.check_output(
        ["git", "show", "--no-textconv", f"{ARCHIVE_COMMIT}:{path}"],
        cwd=ROOT,
    )


def _check_qe(sums: np.ndarray, counts: np.ndarray, label: str) -> None:
    """Reject impossible integer three-state counts, including the endpoints."""
    q, e = sums[..., 0], sums[..., 1]
    limit = counts[:, None, None]
    if (np.any(e < 0) or np.any(e > limit) or np.any(e < np.abs(q))
            or np.any((e + q) % 2)):
        raise ValueError(f"{label}: invalid q/E three-state sums")
    if (not np.all(q[:, :, 0] == -counts[:, None])
            or not np.all(q[:, :, -1] == counts[:, None])
            or not np.all(e[:, :, 0] == counts[:, None])
            or not np.all(e[:, :, -1] == counts[:, None])):
        raise ValueError(f"{label}: incorrect empty/full occupation endpoints")


def load_complement(n: int, marked_profiles: np.ndarray) -> dict:
    """Return integer per-batch q/E sums for old production minus its first 100k.

    Returned ``counts`` has shape (100,) and ``sums`` has shape
    (100, 2, N+1, 2), ordered first/second and q/E.  Batch zero is smaller
    after subtraction; callers must use the returned counts for weighting.
    No source statistics or new random counters are generated.
    """
    if n not in DESIGNS:
        raise ValueError(f"N{n}: no declared original norm-4 production archive")
    marked = np.asarray(marked_profiles)
    if marked.shape != (MARKED_BATCHES, 2, n + 1, 5):
        raise ValueError(f"N{n}: expected marked profile shape (100,2,{n + 1},5)")
    qe = marked[..., :2]
    if (not np.isfinite(qe).all() or np.any(np.abs(qe) > MARKED_PER_BATCH)
            or not np.equal(qe, np.rint(qe)).all()):
        raise ValueError(f"N{n}: marked q/E coordinates must be exact integer sums")
    marked_qe = qe.astype(np.int64)
    _check_qe(marked_qe, np.full(MARKED_BATCHES, MARKED_PER_BATCH), "marked subset")

    suffix = "1b" if n in (260, 340) else "1900m"
    prefix = f"{ARCHIVE_ROOT}/n{n}_{suffix}"
    hist_path = f"{prefix}.hist.csv"
    metadata_path = f"{prefix}.metadata.json"
    metadata_bytes = _git_bytes(metadata_path)
    metadata = json.loads(metadata_bytes)
    original_count = int(metadata["samples_per_pair"])
    batches = int(metadata["batches"])
    designs = metadata["designs"]
    if (batches != 100 or original_count % batches or len(designs) != 1
            or int(designs[0]["N"]) != n
            or tuple(designs[0]["first"]) != DESIGNS[n][0]
            or tuple(designs[0]["second"]) != DESIGNS[n][1]):
        raise ValueError(f"{metadata_path}: unexpected production batching or geometry order")
    per_batch = original_count // batches
    first_counter = int(metadata["replica_counter_first"])
    last_counter = int(metadata["replica_counter_last_exclusive"])
    if last_counter - first_counter != original_count or per_batch <= MARKED_SAMPLES:
        raise ValueError(f"{metadata_path}: inconsistent production counter interval")

    hist_bytes = _git_bytes(hist_path)
    histogram = np.zeros((batches, 2, 2, n + 1), dtype=np.int64)
    seen = set()
    reader = csv.DictReader(io.StringIO(hist_bytes.decode("utf-8")))
    expected_header = ("n", "a", "b", "orientation", "batch", "samples", "kind", "k", "count")
    if tuple(reader.fieldnames or ()) != expected_header:
        raise ValueError(f"{hist_path}: unexpected threshold histogram columns")
    for row in reader:
        g = ORIENTATIONS.index(row["orientation"])
        kind = ("minus", "plus").index(row["kind"])
        batch, k, count = int(row["batch"]), int(row["k"]), int(row["count"])
        if (int(row["n"]) != n or (int(row["a"]), int(row["b"])) != DESIGNS[n][g]
                or not 0 <= batch < batches or not 1 <= k <= n
                or int(row["samples"]) != per_batch or not 0 <= count <= per_batch):
            raise ValueError(f"{hist_path}: inconsistent histogram row {row}")
        key = (batch, g, kind, k)
        if key in seen:
            raise ValueError(f"{hist_path}: duplicate threshold bin {key}")
        seen.add(key)
        histogram[batch, g, kind, k] = count
    if not np.all(histogram.sum(axis=-1) == per_batch):
        raise ValueError(f"{hist_path}: incomplete per-batch threshold histograms")

    cumulative = histogram.cumsum(axis=-1)
    minus, plus = cumulative[:, :, 0], cumulative[:, :, 1]
    sums = np.stack((-per_batch + minus + plus, per_batch - minus + plus), axis=-1)
    counts = np.full(batches, per_batch, dtype=np.int64)
    _check_qe(sums, counts, "full production")
    sums[0] -= marked_qe.sum(axis=0, dtype=np.int64)
    counts[0] -= MARKED_SAMPLES
    _check_qe(sums, counts, "production complement")

    provenance = {
        "commit": ARCHIVE_COMMIT,
        "path": hist_path,
        "sha": hashlib.sha256(hist_bytes).hexdigest(),
        "metadata_path": metadata_path,
        "metadata_sha256": hashlib.sha256(metadata_bytes).hexdigest(),
        "engine_commit": metadata["git_commit"],
        "originalcount": original_count,
        "complementcount": int(counts.sum()),
        "batchcounts": counts.tolist(),
        "original_batchcounts": [per_batch] * batches,
        "removedcount": MARKED_SAMPLES,
        "removed_from_original_batch": 0,
        "original_counter_interval": [first_counter, last_counter],
        "removed_counter_interval": [first_counter, first_counter + MARKED_SAMPLES],
        "seed": metadata["seed"],
        "orientation_order": list(ORIENTATIONS),
        "design": designs[0],
        "fields": ["sum_q", "sum_e"],
        "marked_fields_used": "Only exact q/E sums; density-scaled source columns are ignored",
        "operation": "Reconstruct all-K integer q/E sums from marginal threshold CDFs and subtract all 100 marked batches from original production batch zero",
        "dependency_boundary": "Complement of the first-100k marked subset; preserve the four cyclic shared-counter sizes and the two separate endpoint seed groups from the main manifest",
    }
    return {"counts": counts, "sums": sums, "provenance": provenance}
