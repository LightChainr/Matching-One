#!/usr/bin/env python3
"""Read fixed-reference positive two-birth shape coordinates, no new samples."""
import csv
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/p334-centered-birth-shape"
SOURCE_COMMIT = "9c495ab13e65f2bc93dc0849ee3b73f88724c4b1"
FIELDS = ["rank_center_energy", "rank_lifetime_energy", "canonical_energy", "J0"]


def main():
    config = json.loads((ROOT / "analysis/p334_full_birth_readout.json").read_text())
    p0 = config["p_ref"]
    payload = {
        "source_commit": SOURCE_COMMIT,
        "analysis_code_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "p_ref": p0,
        "new_random_samples": 0,
        "new_conditional_solves": 0,
        "new_path_replays": 0,
        "labels": [f"{orientation}_{field}" for orientation in ("first", "second")
                   for field in FIELDS],
        "covariance_handoff": "Original 20-batch vectors to the single global coordinator; no separate SE or independence assertion.",
        "identities": {
            "rank_energy": "[(C-(N+1)*p_ref)^2+W^2/4]/(N+1)^2",
            "canonical_energy": "[C^2+C+W^2/4]/[(N+1)(N+2)]-2*p_ref*C/(N+1)+p_ref^2",
            "a_N": "[(N+1)*p_ref+1/2]/(N+2)",
            "contrast": "Delta R_ref=-(N+2)/(N+1)*Delta(J1-a_N*J0)",
        },
        "sizes": {},
    }
    for n in (325, 425):
        path = ROOT / "results/p334-full-birth-archive" / f"N{n}.csv"
        batches = np.zeros((20, 8))
        counts = np.zeros(20, dtype=int)
        n1 = n + 1
        with path.open() as stream:
            for row in csv.DictReader(stream):
                batch = int(row["batch"])
                counts[batch] += 1
                for oi, orientation in enumerate(("first", "second")):
                    k1, k2 = int(row[orientation + "_k1"]), int(row[orientation + "_k2"])
                    c, w = (k1 + k2) / 2, k2 - k1
                    center = ((c - n1*p0) / n1)**2
                    lifetime = (w / (2*n1))**2
                    canonical = (c*c + c + w*w/4) / (n1*(n1+1)) - 2*p0*c/n1 + p0*p0
                    batches[batch, 4*oi:4*oi+4] += [center, lifetime, canonical, 1-2*c/n1]
        if not np.all(counts == 1000):
            raise ValueError("incomplete original batches")
        batches /= counts[:, None]
        payload["sizes"][str(n)] = {
            "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "a_N": (n1*p0 + .5) / (n1+1),
            "reference_birth_center": n1*p0,
            "batch_ids": list(range(20)),
            "batch_denominators": counts.tolist(),
            "joint_20_batch_means": batches.tolist(),
            "mean": batches.mean(axis=0).tolist(),
        }
        print(n, "a_N", payload["sizes"][str(n)]["a_N"], "means", batches.mean(axis=0).tolist())
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "batch_vectors.json").write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
