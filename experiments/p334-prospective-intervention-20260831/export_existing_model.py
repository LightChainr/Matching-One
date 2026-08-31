#!/usr/bin/env python3
"""Export already fitted coefficients; no refit, new prefix, or tail sampling."""
from pathlib import Path
import hashlib
import json
import subprocess

ROOT = Path(__file__).resolve().parent
COMMIT = "323de7d5ee4a980b3c77e1a972cb6c812a9f88e5"
PATH = "results/p334-new64-feature-loading/score.json"
EXPECTED_SHA256 = "8116979ddcadbea69c6134bcee8bd41744b54ebb0830aac883016e468cf2c3db"
FEATURES = ("joint_safe_mass", "own_score_energy", "own_safe_degree", "own_safe_loop")
RESPONSES = ("source_first.C", "source_first.W", "source_second.C", "source_second.W")


def export_model(source):
    result = {"status": "draft_existing_parameters_not_prospective_freeze", "training_commit": COMMIT,
              "features": FEATURES, "responses": RESPONSES,
              "rule": "old mean response + (new exact contact features - old mean features) @ old beta",
              "new_fits": 0, "new_samples": 0, "sizes": {}}
    for n, item in source["sizes"].items():
        raw = item["raw_batch_means"]
        labels = item["raw_labels"]
        point_raw = [sum(row[j] for row in raw)/20 for j in range(len(labels))]
        raw_loo = [[(20*point_raw[j]-row[j])/19 for j in range(len(labels))] for row in raw]

        def parameters(raw_row, estimate):
            d = dict(zip(labels, raw_row))
            fitted = dict(zip(item["labels"], estimate))
            pi = d["cell00.mass"]
            out = {}
            for ori in ("first", "second"):
                out[ori] = {
                    "mean_features": [d[f"{ori}.meanX.{f}"]/pi for f in FEATURES],
                    "mean_responses": [d[f"{ori}.new64.meanY.{r}"]/pi for r in RESPONSES],
                    "beta": [[fitted[f"{ori}.new64.contact.beta.{f}|{r}"] for r in RESPONSES] for f in FEATURES],
                    "R_old": fitted[f"{ori}.new64.contact.residual_signed_loading"],
                }
            return out

        result["sizes"][n] = {"point": parameters(point_raw, item["estimate"]),
                              "training_LOO": [parameters(row, est) for row, est in zip(raw_loo, item["LOO"])],
                              "training_batch_ids": item["batch_ids"],
                              "original00_prefixes": item["original00_prefixes"]}
    return result


if __name__ == "__main__":
    target = ROOT / "inputs" / "existing_contact_score.json"
    data = target.read_bytes() if target.exists() else subprocess.check_output(["git", "show", f"{COMMIT}:{PATH}"], cwd=ROOT)
    if hashlib.sha256(data).hexdigest() != EXPECTED_SHA256:
        raise ValueError("fixed training blob hash changed")
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(data)
    model = export_model(json.loads(data))
    model["training_blob_sha256"] = hashlib.sha256(data).hexdigest()
    (ROOT / "existing_model.json").write_text(json.dumps(model, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"status":model["status"], "R_old":{n:{o:p["R_old"] for o,p in v["point"].items()} for n,v in model["sizes"].items()},"new_fits":0,"new_samples":0},indent=2))
