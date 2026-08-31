#!/usr/bin/env python3
"""Coordinate-only N900 comparisons; never refit the old source realizations."""
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.stats import chi2

from p267_max_gaussian_three_center import affine_shape_chart, jackknife_covariance

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/"results/p267-max-gaussian-three-center-n900"
OLD_COMMIT = "2a824e96"
OLD_PATH = "results/p267-max-gaussian-three-center/score.json"
CHART_LABELS = ["gaussian_variance_fraction", "early_weight", "late_weight", "relative_middle_gap"]


def chart(row):
    mean = affine_shape_chart(np.array(row["estimate"]))
    replicas = np.array([affine_shape_chart(np.array(v))
                         for v in row["leave_one_common_batch_out_vectors"]])
    return mean, jackknife_covariance(replicas)


def main():
    old = json.loads(subprocess.check_output(["git", "show", OLD_COMMIT+":"+OLD_PATH], cwd=ROOT))
    target = json.loads((OUT/"score.json").read_text())
    rows = dict(old["sources"], **target["sources"])
    if any(row["status"] != "positive_three_center_realization" or "full_covariance" not in row
           for row in rows.values()):
        raise ValueError("Only existing successful, full-LOO realizations can be compared")
    stats = {n: chart(row) for n, row in rows.items()}
    sources = {n: {"source_commit": row["source_commit"], "common_batches": row["common_batches"],
                   "chart": stats[n][0].tolist(), "chart_se": np.sqrt(np.diag(stats[n][1])).tolist(),
                   "chart_covariance": stats[n][1].tolist(),
                   "quarter_coordinate_width_decomposition": row["quarter_coordinate_width_decomposition"]}
               for n, row in rows.items()}
    comparisons = {}
    for previous in ("400", "100"):
        difference = stats["900"][0]-stats[previous][0]
        cov = stats["900"][1]+stats[previous][1]
        blocks = {}
        for name, indices in (("pure_affine_low_moment_gate", np.arange(4)),
                              ("affine_plus_common_Gaussian_blur_necessary_gate", np.arange(1, 4))):
            delta, block = difference[indices], cov[np.ix_(indices, indices)]
            statistic = float(delta@np.linalg.solve(block, delta))
            blocks[name] = {"indices": indices.tolist(), "chi2": statistic, "df": len(indices),
                            "nominal_p": float(chi2.sf(statistic, len(indices))),
                            "scope": "Low-moment necessary gate only; not a complete functional or physical-state test"}
        before = rows[previous]["quarter_coordinate_width_decomposition"]
        after = rows["900"]["quarter_coordinate_width_decomposition"]
        width_cov = np.array(before["full_covariance"])+after["full_covariance"]
        comparisons["N900_minus_N"+previous] = {
            "chart_difference": difference.tolist(), "chart_se": np.sqrt(np.diag(cov)).tolist(),
            "chart_covariance": cov.tolist(), "fixed_coordinate_gates": blocks,
            "quarter_width_difference": (np.array(after["estimate"])-before["estimate"]).tolist(),
            "quarter_width_difference_se": np.sqrt(np.diag(width_cov)).tolist(),
            "quarter_width_difference_covariance": width_cov.tolist()}
    result = {"schema": "matching-one/p267-three-center-independent-N900-coordinate-comparison/v1",
              "construction_freeze": "191c20e2", "pre_reveal_chart_and_CLI_commit": "e97010a4",
              "old_realization_commit": OLD_COMMIT, "old_realization_path": OLD_PATH,
              "chart_labels": CHART_LABELS, "quarter_width_labels": ["total", "common_gaussian", "between_centers"],
              "sources": sources, "comparisons": comparisons,
              "consecutive_increment_cross_covariance": (-stats["400"][1]).tolist(),
              "cross_covariance_definition": "Cov(phi400-phi100,phi900-phi400)=-Cov(phi400). Do not treat consecutive increments as independent.",
              "boundary": "Only stored N100/N400 LOO parameter vectors were re-expressed, with zero old refits. N900 is independent of both older streams. The two N900 comparisons share N900 and are not combined into a meta-p. Three centers are truncated-moment coordinates, not fields. The primary width forecast is not changed; this is its source-defined shape auxiliary.",
              "new_MC": 0, "old_construction_reruns": 0}
    (OUT/"comparison.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"sources": sources, "comparisons": comparisons}, indent=2))


if __name__ == "__main__":
    main()
