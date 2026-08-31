#!/usr/bin/env python3
"""Describe a lossless factor append without copying the large inherited block."""
import hashlib
import json

from p334_safe_contact_response import ROOT, blob

BASE = "621dfbd0325c0ccd2e3ed26d9612487d155267f6"
BASE_DIR = "results/p334-mixed-source-curvature-joint"
APPEND = "2c3a5ca2"
APPEND_PATH = "results/p334-source-normal-curvature/score.json"


def main():
    data = blob(BASE, BASE_DIR+"/score.json")
    base = json.loads(data)
    extra_data = blob(APPEND, APPEND_PATH)
    extra = json.loads(extra_data)
    out = {"schema": "p334.source-normal-covariance-link.v1", "base_commit": BASE,
           "append_commit": APPEND, "append_path": APPEND_PATH,
           "base_score_sha256": hashlib.sha256(data).hexdigest(),
           "append_score_sha256": hashlib.sha256(extra_data).hexdigest(),
           "meaning": "Virtual complete factor F=[F_base,F_append[:,indices]]. Both use the same original20 deletion rows with sqrt(19/20)*(LOO-mean_LOO) sign. Covariance between any blocks is their factor cross-product; no independent-block assumption.",
           "storage": "References only. The large inherited factor and new numerical rows are already committed; no duplicate numerical archive is created.",
           "sizes": {}}
    for n in ("325", "425"):
        b, a = base["sizes"][n], extra["sizes"][n]
        if b["batch_ids"] != a["batch_ids"]:
            raise ValueError("Deleted-batch rows differ")
        indices = [i for i, k in enumerate(a["labels"]) if not k.startswith("raw.")]
        out["sizes"][n] = {"batch_ids": b["batch_ids"],
            "base_factor_path": BASE_DIR+"/"+b["complete_covariance_factor_file"],
            "base_factor_sha256": b["complete_covariance_factor_sha256"],
            "base_coordinate_count": b["complete_coordinate_count"],
            "append_indices": indices, "append_labels": ["source_normal_decomposition."+a["labels"][i] for i in indices],
            "complete_coordinate_count": b["complete_coordinate_count"]+len(indices)}
    path = ROOT/"results/p334-source-normal-curvature/covariance-link.json"
    path.write_text(json.dumps(out, indent=2)+"\n")
    print(path)


if __name__ == "__main__":
    main()
