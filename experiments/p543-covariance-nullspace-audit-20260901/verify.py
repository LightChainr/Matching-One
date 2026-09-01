#!/usr/bin/env python3
"""Reproduce the archived-score audit and verify its decision boundary."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


HERE = Path(__file__).resolve().parent


def main() -> int:
    with tempfile.TemporaryDirectory() as directory:
        reproduced = Path(directory) / "RESULT.json"
        subprocess.run(
            [
                sys.executable,
                str(HERE / "audit.py"),
                "--input",
                str(HERE / "INPUT_VECTORS.json"),
                "--output",
                str(reproduced),
            ],
            check=True,
            stdout=subprocess.DEVNULL,
        )
        if reproduced.read_bytes() != (HERE / "RESULT.json").read_bytes():
            raise SystemExit("RESULT.json does not reproduce byte for byte")

    result = json.loads((HERE / "RESULT.json").read_text(encoding="utf-8"))
    if result["archived_vector_scores"] != 16:
        raise SystemExit("audit must cover exactly 16 archived vector scores")
    if result["classification_counts"] != {
        "unchanged": 15,
        "numerically_changed": 0,
        "interpretation_changed": 1,
    }:
        raise SystemExit("historical classification changed")
    if result["default_displayed_statistics_changed"] != 0:
        raise SystemExit("an archived default displayed statistic changed")

    changed = [
        row for row in result["results"]
        if row["classification"] == "interpretation_changed"
    ]
    if len(changed) != 1 or not changed[0]["id"].startswith("p50_fullcurve"):
        raise SystemExit("P50 must be the sole interpretation change")
    p50 = changed[0]["rescored"]
    if p50["nullspace_status"] != "estimated_near_null_incompatibility":
        raise SystemExit("P50 discarded projection is no longer flagged")
    rows = {
        item["relative_eigenvalue_cutoff"]: item
        for item in p50["cutoff_sensitivity"]
    }
    expected_ranks = {"1.0e-14": 3, "1.0e-10": 2, "0.000001": 1}
    if any(rows[key]["numerical_rank"] != rank for key, rank in expected_ranks.items()):
        raise SystemExit("P50 cutoff-sensitivity ranks changed")

    for line in (HERE / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
        expected, name = line.split("  ", 1)
        actual = hashlib.sha256((HERE / name).read_bytes()).hexdigest()
        if actual != expected:
            raise SystemExit(f"checksum mismatch: {name}")
    print(json.dumps({"status": "verified", "archived_vector_scores": 16,
                      "new_random_samples": 0}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
