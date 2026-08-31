#!/usr/bin/env python3
"""Resolve between-prefix birth transport into within-cell and between-cell parts."""
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np

from p334_birth_covariance_hierarchy import ROOT, CELLS, ORIS, MARKS, PAIRS, PAIR_NAMES, shape_coordinates

SOURCE = "44dc9e3396e39105cae85a29d04b39d0afc82d84"
PATH = "results/p334-birth-covariance-hierarchy/score.json"
OUT = ROOT/"results/p334-rankcell-covariance-transport"


def derive(row, labels, p, delta):
    d, physical = dict(zip(labels, row)), {}
    for mark in MARKS:
        for ori in ORIS:
            mu = [d[f"all.B.{ori}.{v}"] for v in ("x", "y")]
            for cell in CELLS[1:]:
                pi = d[f"{cell}.mass"]
                count = int(round(p*pi))
                if count < 2:
                    raise ValueError("Rank-cell mean-product U-statistic needs two prefixes")
                bg = [d[f"{cell}.B.{ori}.{v}"]/pi for v in ("x", "y")]
                hg = [d[f"{cell}.H.{mark}.{ori}.{v}"]/pi for v in ("x", "y")]
                within, between = [], []
                for (i, j), pair in zip(PAIRS, PAIR_NAMES):
                    off = d[f"{cell}.H.{mark}.{ori}.offquartet.{pair}"]
                    diag = d[f"{cell}.H.{mark}.{ori}.prefixdiag.{pair}"]
                    group_product = pi*(count*(bg[i]*hg[j]+bg[j]*hg[i])-diag/pi)/(count-1)
                    global_product = (p*pi*(mu[i]*hg[j]+mu[j]*hg[i])-diag)/(p-1)
                    within.append(off-group_product)
                    between.append(group_product-global_product)
                for part, values in (("within_rankcell_prefixes", within), ("between_rankcells", between)):
                    for field, value in shape_coordinates(*values).items():
                        physical[(cell, mark, ori, field, part)] = value
            for field in shape_coordinates(0, 0, 0):
                for part in ("within_rankcell_prefixes", "between_rankcells"):
                    physical[("all", mark, ori, field, part)] = sum(
                        physical[(c, mark, ori, field, part)] for c in CELLS[1:])
    result = {}
    for cell in CELLS:
        for mark in MARKS:
            for field in shape_coordinates(0, 0, 0):
                for part in ("within_rankcell_prefixes", "between_rankcells"):
                    a, b = [physical[(cell, mark, ori, field, part)] for ori in ORIS]
                    for axis, value in (("first", a), ("second", b), ("S", (a+b)/2), ("D", (a-b)/delta)):
                        result[f"{cell}.{mark}->{axis}.{field}.{part}"] = value
    return result


def main():
    blob = subprocess.check_output(["git", "show", SOURCE+":"+PATH], cwd=ROOT)
    data = json.loads(blob)
    OUT.mkdir(parents=True, exist_ok=False)
    result = {"schema": "matching-one/p334-rankcell-covariance-transport/v1",
              "source_commit": SOURCE, "source_path": PATH, "source_sha256": hashlib.sha256(blob).hexdigest(),
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "identity": "between_prefixes=within_rankcell_prefixes+between_rankcells",
              "interpretation": "Rank-cell probabilities do not change. Between-cell covariance changes when their conditional means respond differently. Five cells exhaust the tangent; the global baseline still includes all nine.",
              "estimator": "Conditional distinct-prefix U-products within each observed rank cell, weighted by its original population mass; global products and cell counts recomputed in each original-batch deletion.",
              "new_raw_reads": 0, "new_samples": 0, "sizes": {}}
    for ns, s in data["sizes"].items():
        x = np.array(s["raw_batch_means"])
        mean = x.mean(axis=0)
        n = len(x); p = sum(s["prefix_counts"])
        point = derive(mean, s["raw_labels"], p, s["delta_cos4"])
        loo = np.array([list(derive((n*mean-row)/(n-1), s["raw_labels"], p-s["prefix_counts"][j], s["delta_cos4"]).values())
                        for j, row in enumerate(x)])
        factor = np.sqrt((n-1)/n)*(loo-loo.mean(axis=0))
        result["sizes"][ns] = {"batch_ids": s["batch_ids"], "labels": list(point), "estimate": list(point.values()),
                              "se": np.linalg.norm(factor, axis=0).tolist(), "LOO": loo.tolist(), "factor": factor.tolist()}
        for mark, axis in (("plus", "S"), ("minus", "D")):
            for field in ("cov_xy", "var_C", "var_W"):
                for part in ("within_rankcell_prefixes", "between_rankcells"):
                    key = f"all.{mark}->{axis}.{field}.{part}"
                    j = list(point).index(key)
                    print(ns, key, f"{point[key]:+.10g} +/- {np.linalg.norm(factor[:,j]):.6g}", flush=True)
    (OUT/"score.json").write_text(json.dumps(result, separators=(",", ":"), allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
