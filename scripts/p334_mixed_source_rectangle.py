#!/usr/bin/env python3
"""Read mixed two-source response from the four already computed corners."""
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = "8ad30617b0a3076a5c01a208eb213096d8879b32"
INPUT = "experiments/p334-finite-source-20260831/output/N{n}.json"
OUT = ROOT / "results/p334-mixed-source-rectangle"
CELLS = ("all", "00", "01", "02", "10", "20")
DELTA = {325: -0.7634556213017751, 425: -0.8928996539792388}
OBS = {"A_ref": ("p_ref", "A", 1.), "E_ref": ("p_ref", "E", 1.),
       "C": ("p_integral", "A", -.5), "W": ("p_integral", "E", -1.)}
CORNERS = {"++": "plus.+1", "--": "plus.-1", "+-": "minus.+1", "-+": "minus.-1"}


def main():
    OUT.mkdir(parents=True, exist_ok=False)
    result = {"schema": "p334.mixed-source-rectangle.v1", "source_commit": SOURCE,
              "reader_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "source_coordinates": "q_t within class a proportional to exp[pi_a(t_first L_first+t_second L_second)], class mass pi_a fixed",
              "corners": {k: [(.5 if c == "+" else -.5) for c in k] for k in CORNERS},
              "estimand": "F(+,+)+F(-,-)-F(+,-)-F(-,+) at physical +/-0.5 source coordinates; all responses are changes from the same zero-source baseline",
              "interpretation": "Finite mixed interaction equals the box integral of d_first d_second F, not exactly its value at zero. It vanishes for additive fixed-coordinate response functions and for every prefix outside00. No path ordering or noncommutativity is implied.",
              "sampling": "Reuse only saved original8 finite-source20-batch vectors. No raw forks, weights, prefixes or suffixes recalculated.",
              "new_samples": 0, "p_ref": .59274605079, "input_sha256": {}, "sizes": {}}
    for n in (325, 425):
        path = INPUT.format(n=n)
        data = subprocess.check_output(["git", "show", f"{SOURCE}:{path}"], cwd=ROOT)
        source = json.loads(data)
        result["input_sha256"][path] = hashlib.sha256(data).hexdigest()
        source_x = np.asarray(source["batch_values"])
        pos = {k: j for j, k in enumerate(source["labels"])}
        columns = {}
        for cell in CELLS:
            for observable, (endpoint, ae, scale) in OBS.items():
                by_corner = {}
                for corner, policy in CORNERS.items():
                    s, d = [scale*source_x[:, pos[f"{cell}.{endpoint}.{sector}.{ae}.{policy}"]]
                            for sector in ("S", "D")]
                    by_corner[corner] = {"S": s, "D": d,
                                         "first": s+DELTA[n]*d/2,
                                         "second": s-DELTA[n]*d/2}
                for receiver in ("first", "second", "S", "D"):
                    stem = f"{cell}.{receiver}.{observable}."
                    z = {corner: by_corner[corner][receiver] for corner in CORNERS}
                    for corner in CORNERS:
                        columns[stem+f"corner{corner}"] = z[corner]
                    columns[stem+"mixed_rectangle"] = z["++"]+z["--"]-z["+-"]-z["-+"]
        x = np.column_stack(list(columns.values()))
        mean = x.mean(axis=0)
        loo = (20*mean-x)/19
        factor = np.sqrt(19/20)*(loo-loo.mean(axis=0))
        errors = np.linalg.norm(factor, axis=0)
        result["sizes"][str(n)] = {"batch_ids": source["batch_ids"],
            "population_per_batch": source["population_per_batch"],
            "labels": list(columns), "batch_values": x.tolist(), "estimate": mean.tolist(),
            "LOO": loo.tolist(), "factor": factor.tolist(), "se": errors.tolist(),
            "factor_convention": "sqrt(19/20)*(LOO-mean_LOO); opposite sign to old finite-source raw-centered factor; same original deletion row"}
        for key, value, error in zip(columns, mean, errors):
            if key.startswith("all.") and key.endswith("mixed_rectangle"):
                print(n, key, f"{value:.11g} +/- {error:.6g}", flush=True)
    (OUT/"score.json").write_text(json.dumps(result, separators=(",", ":"), allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
