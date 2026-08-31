#!/usr/bin/env python3
"""Exact policy-tangent projections of saved same-degree contact covariances."""
import json
from pathlib import Path
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = "b9f79bfb6e1ba4177ff245f74f7b2e51c3bd2fdc"
source_path = "results/p334-safe-contact-response/score.json"
source = json.loads(subprocess.check_output(["git", "show", SOURCE+":"+source_path], cwd=ROOT))
out = {"source_commit": SOURCE, "p_ref": source["p_ref"], "new_samples": 0,
       "policy": "At each R0 prefix preserve each safe contact-degree class mass pi_e, and tilt the within-class uniform next-label distribution by exp(t*pi_e*loop). Keep all other labels unchanged. Average the two orientation-specific policies equally.",
       "identity": "Derivative at t=0 equals sum_e pi_e^2 Cov(loop,m|Z,safe,e), exactly the saved half-difference contact numerator.",
       "sizes": {}}
for ntext, saved in source["sizes"].items():
    n = int(ntext)
    raw = np.array(saved["raw_joint_20_batch_means"])
    names = saved["raw_labels"]
    base = "R0_safe_equal_contact_degree.GX[contractible_cycles,"
    h = {name: raw[:, names.index(base+name+"]")] for name in ("p_ref.F1", "p_ref.F2", "p_integral.F1", "p_integral.F2")}
    k1, k2 = -(n+1)*h["p_integral.F1"], -(n+1)*h["p_integral.F2"]
    cols = {"K1": k1, "K2": k2, "C": (k1+k2)/2, "W": k2-k1}
    for ep in ("p_ref", "p_integral"):
        cols[ep+".A"] = h[ep+".F1"]+h[ep+".F2"]
        cols[ep+".E"] = -h[ep+".F1"]+h[ep+".F2"]
    x = np.column_stack(list(cols.values()))
    mean = x.mean(0)
    factor = (x-mean)/np.sqrt(20*19)
    se = np.linalg.norm(factor, axis=0)
    out["sizes"][ntext] = {"batch_ids": list(range(20)), "labels": list(cols),
        "joint_20_batch_means": x.tolist(), "estimate": mean.tolist(), "se": se.tolist(),
        "factor": factor.tolist(), "covariance": (factor.T@factor).tolist()}
    for name, value, error in zip(cols, mean, se):
        print(n, name, value, error)
path = ROOT/"results/p334-euler-invisible-tangent"
path.mkdir(parents=True, exist_ok=True)
(path/"score.json").write_text(json.dumps(out, indent=2, allow_nan=False)+"\n")
