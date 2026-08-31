#!/usr/bin/env python3
"""One frozen primary contrast, computed after both independent sizes finish."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

FORECAST = np.array([4.116, 3.233, 3.300, 3.977])*1e-8
DELTA = 1e-8


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    result = {"schema": "p334.independent-normal-decision.v1", "sizes": {},
              "forecast_four_C": FORECAST.tolist(), "delta": DELTA,
              "source_groups": "fresh N325 and N425 independent; discovery block is not pooled",
              "input_sha256": {}}
    size_means, size_variances, four = [], [], []
    for n in (325, 425):
        directory = a.input / f"N{n}"
        run = json.loads((directory / "run.json").read_text())
        if run["status"] != "completed":
            raise ValueError(f"N{n} incomplete")
        rows = []
        counts = []
        for b in range(20):
            path = directory / f"batch-{b:02}" / "batch.json"
            data = path.read_bytes()
            result["input_sha256"][str(path.relative_to(a.input))] = hashlib.sha256(data).hexdigest()
            r = json.loads(data)
            if (r["N"], r["batch"], r["full_prefix_denominator"], r["paired_reps_per_own_axis"], r["master_seed"], r["technical_smoke"]) != (n, b, 25000, 8, 202608311920334, False):
                raise ValueError("Frozen block identity differs")
            if b == 0:
                labels = r["labels"]
            elif labels != r["labels"]:
                raise ValueError("Batch labels differ")
            rows.append(r["mean_full_prefix"])
            counts.append(r["original00_count"])
        matrix = np.array(rows)
        mean = matrix.mean(axis=0)
        factor = (matrix-mean)/np.sqrt(20*19)
        ix = [labels.index(f"source_{o}__receiver_{o}__C") for o in ("first", "second")]
        primary_batches = matrix[:, ix].mean(axis=1)
        size_means.append(float(primary_batches.mean()))
        size_variances.append(float(primary_batches.var(ddof=1)/20))
        four.extend(mean[ix].tolist())
        result["sizes"][str(n)] = {"labels": labels, "batch_ids": list(range(20)),
            "raw_batch_means": matrix.tolist(), "mean": mean.tolist(),
            "se": np.linalg.norm(factor, axis=0).tolist(), "factor": factor.tolist(),
            "original00_counts": counts, "full_prefix_count": 500000,
            "own_C_mean": mean[ix].tolist(), "primary_size_mean": size_means[-1],
            "primary_size_se": np.sqrt(size_variances[-1])}
    estimate = sum(size_means)/2
    se = np.sqrt(sum(size_variances))/2
    low, high = estimate-3*se, estimate+3*se
    decision = ("stop_complete_two_score_label_closure" if low > DELTA else
                "stop_material_positive_archive_forecast" if high < DELTA else
                "inconclusive_at_fixed_budget_stop")
    result["primary"] = {"estimate": estimate, "se": se, "three_se_interval": [low, high],
                         "frozen_forecast": float(FORECAST.mean()), "decision": decision,
                         "four_coordinate_errors": (np.array(four)-FORECAST).tolist(),
                         "meaning": "fixed-batch diagnostic; not anytime-valid or a proof of a continuum/global mechanism"}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps(result["primary"], indent=2))
    for n in ("325", "425"):
        r = result["sizes"][n]
        for o in ("first", "second"):
            for f in ("C", "A_ref", "W"):
                i = r["labels"].index(f"source_{o}__receiver_{o}__{f}")
                print(n, o, f, r["mean"][i], "+/-", r["se"][i])


if __name__ == "__main__":
    main()
