#!/usr/bin/env python3
"""Apply the unchanged common-symmetric-kernel moment gate to saved N900.

Forward-evaluate stored LOO three-center coordinates back to their matched
moments3..6; never rerun the source realization or select another branch.
"""
from collections import Counter
import hashlib
import json
from math import factorial
from pathlib import Path
import subprocess

import numpy as np

from etop_two_lobe_moment_closure import SYMMETRIC_LABELS, covariance, symmetric_kernel

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/"results/etop-n900-symmetric-two-lobe-moments"
SOURCE_REV="54430ea7"
SOURCE_PATH="results/p267-max-gaussian-three-center-n900/score.json"


def forward_moments(vector):
    """m3..m6 only; an algebraic readout of the saved matched coordinates."""
    alpha=float(vector[0])
    centers=np.array(vector[1:4]);weights=np.array(vector[4:7])
    return np.array([sum(alpha**j*factorial(r)/(2**j*factorial(j)*factorial(r-2*j))
                         *float(weights@centers**(r-2*j)) for j in range(r//2+1))
                     for r in range(3,7)])


def classify(moments):
    m3,m4,m5,m6=(float(v) for v in moments)
    diagnostics={"m3":m3,"odd_root_equation_at_r1":None if abs(m3)<1e-14 else 8+(m5-10*m3)/m3-m3*m3}
    try:
        value=symmetric_kernel(moments)
    except ValueError as error:
        if abs(m3)<1e-14:
            status="degenerate_zero_skew_identification"
        elif diagnostics["odd_root_equation_at_r1"]<0:
            status="no_admissible_r_le1_branch"
        else:
            status="frozen_solver_failure"
        return {"status":status,"reason":str(error),"diagnostics":diagnostics}
    if not np.all(np.isfinite(value)):
        return {"status":"nonfinite_frozen_branch","diagnostics":diagnostics}
    return {"status":"identified_unique_admissible_branch","estimate":value.tolist(),
            "diagnostics":diagnostics,"kernel_m6_negative":bool(value[4]<0),
            "moment_Hankel_determinant_negative":bool(value[5]<0)}


def main():
    revision=subprocess.check_output(["git","rev-parse",SOURCE_REV],cwd=ROOT,text=True).strip()
    blob=subprocess.check_output(["git","show",revision+":"+SOURCE_PATH],cwd=ROOT)
    source=json.loads(blob)["sources"]["900"]
    point_moments=np.array(source["moments_0_to_8"][3:7])
    point=classify(point_moments)
    replicas=[]
    for index,vector in enumerate(source["leave_one_common_batch_out_vectors"]):
        moments=forward_moments(vector)
        replicas.append({"deleted_common_batch_index":index,"moments3_to6":moments.tolist(),**classify(moments)})
    counts=Counter(r["status"] for r in replicas)
    branch_statuses=("identified_unique_admissible_branch",
                     "no_admissible_r_le1_branch",
                     "degenerate_zero_skew_identification",
                     "frozen_solver_failure", "nonfinite_frozen_branch")
    all_valid=all(r["status"]=="identified_unique_admissible_branch" for r in replicas)
    uncertainty={"status":"not_reported_incomplete_frozen_branches"}
    if all_valid and point["status"]=="identified_unique_admissible_branch":
        values=np.array([r["estimate"] for r in replicas])
        cov=covariance(values)
        se=np.sqrt(np.diag(cov))
        uncertainty={"status":"full_aligned800_delete_one_covariance",
                     "standard_errors":se.tolist(),"covariance":cov.tolist(),
                     "marginal_estimate_over_SE":(np.array(point["estimate"])/se).tolist(),
                     "kernel_m6_negative_LOO_count":sum(r["kernel_m6_negative"] for r in replicas),
                     "hankel_negative_LOO_count":sum(r["moment_Hankel_determinant_negative"] for r in replicas),
                     "LOO_coordinate_minima":values.min(axis=0).tolist(),
                     "LOO_coordinate_maxima":values.max(axis=0).tolist(),
                     "interpretation":"Same full-source delete-one propagation as ddf7d564. Marginal z values are descriptive, not a boundary-calibrated cone likelihood or exact confidence certificate."}
    b=len(replicas)
    output={"necessary_condition_commit":"ddf7d564","condition_function":"scripts/etop_two_lobe_moment_closure.py:symmetric_kernel (unchanged)",
            "shape_source_commit":revision,"shape_source_path":SOURCE_PATH,"shape_source_sha256":hashlib.sha256(blob).hexdigest(),
            "production_source_commit":source["source_commit"],"production_source_directory":source["source_directory"],
            "production_source_hashes":source["source_sha256"],"common_batches":b,
            "labels":SYMMETRIC_LABELS,"point_moments3_to6":point_moments.tolist(),"point":point,
            "LOO_moment_input":"Forward moment evaluation of the already-saved matched three-center LOO coordinates; no realization/refit or raw-profile reconstruction.",
            "source_max_LOO_construction_moment_error":source["LOO_construction"]["max_m0_to_m6_error"],
            "branch_counts":{key:counts[key] for key in branch_statuses},
            "branch_proportions":{key:counts[key]/b for key in branch_statuses},
            "no_admissible_branch_count":counts["no_admissible_r_le1_branch"],
            "degenerate_skew_branch_count":counts["degenerate_zero_skew_identification"],
            "other_failed_branch_count":counts["frozen_solver_failure"]+counts["nonfinite_frozen_branch"],
            "uncertainty":uncertainty,"leave_one_common_batch_out":replicas,
            "new_MC":0,"source_construction_reruns":0,"old_scale_reruns":0,"alternative_branch_or_kernel_search":False,
            "scope":"Existing N900 D_A signed rank-step moments; any positive common symmetric two-translation kernel is the candidate, not a positivity assumption about the source. Same N900 dependency group as width/three-center scores."}
    OUTPUT.mkdir(parents=True,exist_ok=True)
    (OUTPUT/"score.json").write_text(json.dumps(output,indent=2,allow_nan=False)+"\n")
    print("point",point)
    print("branch counts",counts)
    print("uncertainty",uncertainty)


if __name__=="__main__":
    main()
