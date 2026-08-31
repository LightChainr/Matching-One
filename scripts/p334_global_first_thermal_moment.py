#!/usr/bin/env python3
"""Read fixed J0/J1 center and lifetime coordinates on the completed paths."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/p334-global-first-thermal-moment"


def main():
    payload = {"source_commit": "9c495ab13e65f2bc93dc0849ee3b73f88724c4b1",
               "theory_commit": "f0dbc070b826761b097171315ab64750dee90823",
               "new_random_samples": 0, "new_conditional_solves": 0,
               "labels": [f"{orientation}_{field}" for orientation in ("first", "second")
                          for field in ("J0", "J1_center", "J1_width")],
               "identity": "J1=J1_center+J1_width; width=-W^2/[4(N+1)(N+2)]",
               "sizes": {}, "covariance_handoff": "Original20 batch vectors, to shared global covariance coordinator; no independent errors computed here."}
    for n in (325,425):
        path = ROOT / "results/p334-full-birth-archive" / f"N{n}.csv"
        with path.open() as stream:
            rows = list(csv.DictReader(stream))
        batches = np.zeros((20,6)); counts = np.zeros(20, dtype=int)
        for row in rows:
            batch = int(row["batch"]); counts[batch] += 1
            for oi, orientation in enumerate(("first", "second")):
                k1, k2 = int(row[orientation+"_k1"]), int(row[orientation+"_k2"])
                c, w = (k1+k2)/2, k2-k1
                d2 = (n+1)*(n+2)
                batches[batch,3*oi:3*oi+3] += [1-2*c/(n+1),
                    .5-(c*c+c)/d2, -w*w/(4*d2)]
        if not np.all(counts == 1000):
            raise ValueError("incomplete original batches")
        batches /= counts[:,None]
        payload["sizes"][str(n)] = {"source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "batch_ids": list(range(20)), "batch_denominators": counts.tolist(),
            "joint_20_batch_means": batches.tolist(), "mean": batches.mean(axis=0).tolist()}
        print(n, "mean first/second J0,J1_center,J1_width", batches.mean(axis=0).tolist())
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"batch_vectors.json").write_text(json.dumps(payload,indent=2)+"\n")


if __name__ == "__main__":
    main()
