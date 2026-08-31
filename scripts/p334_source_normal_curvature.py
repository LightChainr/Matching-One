#!/usr/bin/env python3
"""Separate second source response from its exact first-score tangent part."""
import hashlib
import io
import json
from pathlib import Path
import subprocess
import time

import numpy as np

from p334_safe_contact_response import ROOT, blob, array_csv

CLASS = "1cfa4ae892a2f7f4168e9a71690efd7a5560d4cd"
DATA = "8ad30617b0a3076a5c01a208eb213096d8879b32"
HESSIAN = "c48fa360a37a9887ef32ff6d3ce947c4e4601b53"
OUT = ROOT / "results/p334-source-normal-curvature"
ORIS = ("first", "second")
OBS = ("A_ref", "E_ref", "C", "W")
PAIRS = ((0, 0), (0, 1), (1, 1))
COMP = ("ff", "fs", "ss")


def main():
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=False)
    hashes = {}

    def read(commit, path):
        data = blob(commit, path)
        hashes[f"{commit}:{path}"] = hashlib.sha256(data).hexdigest()
        return data

    result = {"schema": "p334.source-normal-curvature.v1", "descriptor_commit": CLASS,
              "first_response_and_census_commit": DATA, "second_response_commit": HESSIAN,
              "reader_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
              "reader_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "definition": "G=E[s s^T], M_ij,k=E[T_ij s_k], alpha=M G^-1; tangent response alpha H_1 F; normal response H_2 F-alpha H_1 F, pointwise in original00 prefix. G and M use exact full-census moments, not response estimates.",
              "interpretation": "Source-normal density direction is orthogonal to both original first scores. A nonzero mean normal response rejects affine first-score sufficiency for the future conditional label mean; zero at these directions is not a complete closure theorem. Not a continuum field count, or temporal curvature.",
              "new_samples": 0, "new_fork_gzip_reads": 0, "new_DP": 0,
              "basis": {"receiver": ORIS, "observable": OBS, "component": COMP},
              "dependence": "Original00 new64 and same original20 batches only; shared-prefix and shared-suffix dependence retained, no independent replication.",
              "sizes": {}}
    for n in (325, 425):
        with np.load(io.BytesIO(read(CLASS, f"results/p334-exact-prefix-structure/N{n}.npz")), allow_pickle=False) as f:
            ix00 = np.flatnonzero(f["rankcell"] == 0)
            counter, batch, d = f["counter"][ix00], f["batch"][ix00], f["d"][ix00]
            gram = f["score_gram_physical"][ix00]
            count, sums = f["class_count"][ix00], f["class_loop_sum"][ix00]
        header, census = array_csv(read(DATA, f"experiments/p334-finite-source-20260831/census/N{n}/census.csv.gz"))
        ci = {k: j for j, k in enumerate(header)}
        census = census[(census[:, ci["first_rank"]] == 0)&(census[:, ci["second_rank"]] == 0)]
        order = np.argsort(counter)
        ids = order[np.searchsorted(counter[order], census[:, ci["counter"]])]
        if not np.array_equal(counter[ids], census[:, ci["counter"]]):
            raise ValueError("Census and original00 identities differ")
        clas = (ids, census[:, ci["first_e"]], census[:, ci["second_e"]])
        marks = census[:, [ci["L_first"], ci["L_second"]]]
        a = count[clas][:, None]*marks-sums[clas]
        cubic_num = np.zeros((len(counter), 3, 2), dtype=np.int64)
        for k, (i, j) in enumerate(PAIRS):
            for z in range(2):
                np.add.at(cubic_num[:, k, z], ids, census[:, ci["count"]]*a[:, i]*a[:, j]*a[:, z])
        moment = cubic_num / d[:, None, None].astype(float)**4
        alpha = moment @ np.linalg.inv(gram)

        with np.load(io.BytesIO(read(HESSIAN, f"results/p334-mixed-source-curvature/new64_N{n}.npz")), allow_pickle=False) as h:
            if not np.array_equal(counter, h["counter"]) or not np.array_equal(batch, h["batch"]):
                raise ValueError("Second responses use different original00 prefixes")
            second = h["mean_response2"]
        first_path = f"experiments/p334-mechanism-response-20260831/results-extension/prefix_statistics_N{n}.npz"
        with np.load(io.BytesIO(read(DATA, first_path)), allow_pickle=False) as h:
            if not np.array_equal(counter, h["counter"]) or not np.array_equal(batch, h["batch"]):
                raise ValueError("First responses use different original00 prefixes")
            names, values = list(h["labels"]), h["new64"]
            first = np.empty((len(counter), 2, 4, 2))
            for oi, ori in enumerate(ORIS):
                for fi, (endpoint, obs, scale) in enumerate((("p_ref", "A", 1.), ("p_ref", "E", 1.),
                                                           ("p_integral", "A", -.5), ("p_integral", "E", -1.))):
                    for si, source in enumerate(ORIS):
                        first[:, oi, fi, si] = scale*values[:, names.index(f"{endpoint}.{obs}.mean_J[{ori},{source}]")]
        tangent = np.einsum("pjk,pofk->pofj", alpha, first)
        normal = second-tangent
        fields = {}
        for mode, v in (("raw", second), ("source_tangent", tangent), ("source_normal", normal)):
            for oi, ori in enumerate(ORIS):
                for fi, obs in enumerate(OBS):
                    for cj, component in enumerate(COMP):
                        fields[f"{mode}.{ori}.{obs}.{component}"] = np.array([
                            v[batch == b, oi, fi, cj].sum()/1000 for b in range(20)])
        raw = np.column_stack(list(fields.values())); mean = raw.mean(axis=0)
        loo = (20*mean-raw)/19
        factor = np.sqrt(19/20)*(loo-loo.mean(axis=0))
        se = np.linalg.norm(factor, axis=0)
        path = OUT/f"N{n}.npz"
        np.savez_compressed(path, counter=counter, batch=batch, d=d, gram=gram,
                            cubic_numerator=cubic_num, tangent_coefficients=alpha,
                            mean_first_response=first, mean_second_response=second,
                            mean_tangent_response=tangent, mean_normal_response=normal)
        result["sizes"][str(n)] = {"prefix_count": len(counter), "batch_ids": list(range(20)),
            "population_per_batch": 1000, "labels": list(fields), "raw_batch_means": raw.tolist(),
            "estimate": mean.tolist(), "se": se.tolist(), "LOO": loo.tolist(), "factor": factor.tolist(),
            "prefix_file": path.name, "prefix_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        for key, value, error in zip(fields, mean, se):
            if key.startswith("source_normal.") and (".C." in key or ".A_ref." in key):
                print(n, key, f"{value:.11g} +/- {error:.6g}", flush=True)
    result["input_sha256"] = hashes
    result["elapsed_seconds"] = time.perf_counter()-started
    (OUT/"score.json").write_text(json.dumps(result, separators=(",", ":"), allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
