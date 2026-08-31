#!/usr/bin/env python3
"""Unbiased within-prefix common-label response rank from existing eight quartets.

No new prefixes, next labels, tails, census, or fitted source centering.
The half-difference estimator exactly removes each label-class source center.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import gzip
import hashlib
from itertools import combinations
import json
import os
from pathlib import Path
import platform
import time

import numpy as np
import scipy
from scipy.stats import binom

P_REF = 0.59274605079
DELTA = {325: -0.7634556213017751, 425: -0.8928996539792388}
OBS = ["p_ref.A", "p_ref.E", "p_integral.A", "p_integral.E", "input_G"]
PAIR = list(combinations(range(8), 2))
PAIR_INDEX = {x: i for i, x in enumerate(PAIR)}
DISJOINT = []
for a, b, c, d in combinations(range(8), 4):
    DISJOINT.extend([(PAIR_INDEX[(a, b)], PAIR_INDEX[(c, d)]),
                     (PAIR_INDEX[(a, c)], PAIR_INDEX[(b, d)]),
                     (PAIR_INDEX[(a, d)], PAIR_INDEX[(b, c)])])

def array_csv(path):
    with gzip.open(path, "rt") as stream:
        header = stream.readline().strip().split(",")
        data = np.loadtxt(stream, delimiter=",", dtype=np.int64)
    return {x: i for i, x in enumerate(header)}, data

def determinant_polarization(x, y):
    return (x[..., 0, 0] * y[..., 1, 1] + y[..., 0, 0] * x[..., 1, 1]
            - x[..., 0, 1] * y[..., 1, 0] - y[..., 0, 1] * x[..., 1, 0]) / 2

def algebra_checks():
    # Exact population enumeration: iid draws from a finite arbitrary matrix law.
    matrices = np.array([[[1., 2.], [3., 4.]], [[-1., 2.], [0., 3.]], [[2., 0.], [1., -1.]]])
    expected = np.linalg.det(matrices.mean(axis=0))
    d = determinant_polarization(matrices[:, None], matrices[None, :])
    assert abs(d.mean() - expected) < 1e-12
    assert abs((d[:, :, None, None] * d[None, None, :, :]).mean() - expected**2) < 1e-12
    # A fixed common response column or a single active source has local rank <= 1.
    rank_one = np.einsum("i,qj->qij", [1., -2.], [[1., 2.], [3., -1.], [2., 4.]])
    assert np.max(np.abs(determinant_polarization(rank_one[:, None], rank_one[None, :]))) == 0
    return {"finite_iid_population_det_identity_error": float(d.mean() - expected),
            "finite_iid_population_det2_identity_error": float((d[:, :, None, None] * d[None, None, :, :]).mean() - expected**2),
            "rank_one_common_direction_kernel": "exact zero"}

def read_batch(task):
    root, n, batch = task
    root = Path(root)
    hi, raw = array_csv(root / "forks" / f"N{n}.batch{batch:02}.csv.gz")
    ci, contact = array_csv(root / "contact" / f"N{n}.batch{batch:02}.csv.gz")
    raw = raw[np.lexsort(tuple(raw[:, hi[k]] for k in ("replica", "group", "quartet", "counter")))].reshape(1000, 8, 2, 2, -1)
    contact = contact[np.lexsort(tuple(contact[:, ci[k]] for k in ("group", "quartet", "counter")))].reshape(1000, 8, 2, -1)
    for k in ("counter", "quartet", "group", "next_label"):
        if not np.array_equal(raw[..., 0, hi[k]], contact[..., ci[k]]):
            raise ValueError(f"Join mismatch: N{n}/batch{batch}/{k}")
    for k in ("N", "batch"):
        assert np.all(raw[..., hi[k]] == {"N": n, "batch": batch}[k])
    counter = raw[:, 0, 0, 0, hi["counter"]]
    ranks = np.stack([raw[:, 0, 0, 0, hi[f"{o}_rank"]] for o in ("first", "second")], axis=-1)
    cells = 3 * ranks[:, 0] + ranks[:, 1]
    selected = np.ones((1000, 8), dtype=bool)
    loop_marks = []
    responses = []
    tail = binom.sf(np.arange(n + 1) - 1, n, P_REF)
    for io, o in enumerate(("first", "second")):
        e, c = (contact[..., ci[f"{o}_{k}"]] for k in ("e", "c"))
        nr = contact[..., ci[f"{o}_rank_after"]]
        selected &= np.all(nr == ranks[:, io, None, None], axis=2)
        selected &= e[:, :, 0] == e[:, :, 1]
        loop_marks.append((ranks[:, io] == 0)[:, None, None] * (e - c))
        k1, k2 = (raw[..., hi[f"{o}_{k}"]] for k in ("k1", "k2"))
        assert np.all((0 <= k1) & (k1 <= k2) & (k2 <= n))
        # Constants in A/E cancel in the label half-difference.
        responses.append(np.stack((tail[k1] + tail[k2] - 1, 1 - tail[k1] + tail[k2],
                                   1 - (k1 + k2) / (n + 1), (k1 - k2) / (n + 1) + 1), axis=-1))
    marks = np.stack(loop_marks, axis=-1)  # prefix, quartet, U/V, input L_first/L_second
    dm = marks[:, :, 0] - marks[:, :, 1]
    y = np.stack(responses, axis=-1)  # prefix, quartet, U/V, replica, observable, output
    dy = y[:, :, 0].mean(axis=2) - y[:, :, 1].mean(axis=2)
    matrices = selected[..., None, None, None] * dy[..., :, :, None] * dm[..., None, None, :] / 2
    gram = selected[..., None, None] * dm[..., :, None] * dm[..., None, :] / 2
    matrices = np.concatenate((matrices, gram[:, :, None]), axis=2)
    pol = np.stack([determinant_polarization(matrices[:, q], matrices[:, r]) for q, r in PAIR], axis=1)
    udet = pol.mean(axis=1)
    udet2 = sum(pol[:, a] * pol[:, b] for a, b in DISJOINT) / len(DISJOINT)
    energy = sum(np.sum(matrices[:, q] * matrices[:, r], axis=(-2, -1)) for q, r in PAIR) / len(PAIR)
    # Conditional support makes local determinant exactly zero outside cell00.
    assert np.all(udet[cells != 0] == 0)
    assert np.all(udet2[cells != 0] == 0)
    # Input Gram determinants and their disjoint products are nonnegative by PSD geometry.
    assert np.min(pol[..., -1]) >= -1e-15
    active = selected & np.any(dm != 0, axis=-1)
    matrix_mean = matrices.mean(axis=1)
    fields = {}
    for j, obs in enumerate(OBS):
        fields[obs + ".E_det_JZ"] = udet[:, j]
        fields[obs + ".E_det_JZ_squared"] = udet2[:, j]
        fields[obs + ".E_frobenius_JZ_squared"] = energy[:, j]
        for i, o in enumerate(("first", "second")):
            for k, inp in enumerate(("first", "second")):
                fields[obs + f".mean_J[{o},{inp}]"] = matrix_mean[:, j, i, k]
    fields["diagnostic.selected_quartets_per_prefix"] = selected.sum(axis=1)
    fields["diagnostic.active_quartets_per_prefix"] = active.sum(axis=1)
    fields["diagnostic.first_source_active_quartets"] = (selected & (dm[..., 0] != 0)).sum(axis=1)
    fields["diagnostic.second_source_active_quartets"] = (selected & (dm[..., 1] != 0)).sum(axis=1)
    fields["diagnostic.noncollinear_source_quartet_pairs"] = (pol[..., -1] > 0).sum(axis=1)
    fields["diagnostic.usable_four_quartet_pairings"] = sum((pol[:, a, -1] > 0) & (pol[:, b, -1] > 0) for a, b in DISJOINT)
    fields["diagnostic.at_least_one_noncollinear_pair"] = fields["diagnostic.noncollinear_source_quartet_pairs"] > 0
    fields["diagnostic.at_least_one_usable_four_pairing"] = fields["diagnostic.usable_four_quartet_pairings"] > 0
    x = np.column_stack(list(fields.values()))
    return {"N": n, "batch": batch, "counter": counter, "cells": cells,
            "labels": list(fields), "prefix_values": x,
            "active_quartet_histogram_by_cell": np.array([[np.sum((cells == cell) & (active.sum(1) == a)) for a in range(9)] for cell in range(9)])}

def derive(values, labels):
    out = dict(zip(labels, values))
    for group in ["all"] + [f"cell{a}{b}" for a in range(3) for b in range(3)]:
        for obs in OBS:
            base = group + "." + obs
            j = np.array([[out[base + f".mean_J[{a},{b}]"] for b in ("first", "second")] for a in ("first", "second")])
            out[base + ".det_mean_JZ"] = float(np.linalg.det(j))
            out[base + ".ensemble_minus_mean_local_det"] = float(np.linalg.det(j)) - out[base + ".E_det_JZ"]
    return out

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=Path(__file__).parent / "inputs")
    parser.add_argument("--output", type=Path, default=Path(__file__).parent / "results")
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--attempt", type=int, default=1)
    args = parser.parse_args()
    started = time.time()
    args.output.mkdir(parents=True, exist_ok=False)
    manifest = json.loads((args.inputs / "manifest.json").read_text())
    for entry in manifest["files"]:
        data = (args.inputs / entry["local_path"]).read_bytes()
        assert hashlib.sha256(data).hexdigest() == entry["sha256"], entry["local_path"]
    algebra = algebra_checks()
    common = json.loads((args.inputs / "common_label_score.json").read_text())
    tasks = [(str(args.inputs), n, b) for n in (325, 425) for b in range(20)]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        blocks = list(pool.map(read_batch, tasks))
    result = {"estimand": "E_Z det J(Z), E_Z[(det J(Z))^2], versus det(E_Z J(Z)); same common-label class-preserving source",
              "source_basis": "L_first,L_second. g_plus=(L_first+L_second)/2, g_minus=(L_first-L_second)/2; ensemble J matches 73608ba9.",
              "matrix_kernel": "X_q[o,s] = selected*(L_s(U)-L_s(V))*(Ybar_o(U)-Ybar_o(V))/2; E[X_q|Z]=J(Z)",
              "class_law": "Same joint safe and equal (e_first,e_second) label class; fixed class mass pi_a; exp(t*pi_a*L_s) tilt.",
              "independence": "Eight separate quartet RNG stream keys conditional on original prefix; independent U,V labels with replacement and two tails per label; paired geometries share each tail. No quartet-pooled source centering.",
              "det_kernel": "D(q,r)=(a_q*d_r+a_r*d_q-b_q*c_r-b_r*c_q)/2; average all28 distinct-quartet pairs.",
              "det2_kernel": "Average 70 four-subsets and all3 disjoint pairings D(i,j)D(k,l), i,j,k,l all distinct.",
              "denominators": "All original20 batches x1000 prefixes per size; each cell contribution uses full1000 denominator, never selection-normalized.",
              "uncertainty": "Original20 batch delete-one factors; no clipping, Gaussian p-values, rank-one acceptance, or independent-source claim. Nonnegative population det-square can have negative unbiased estimate.",
              "new_prefixes": 0, "new_labels": 0, "new_tails": 0, "algebra_checks": algebra,
              "sizes": {}, "manifest_sha256": hashlib.sha256((args.inputs / "manifest.json").read_bytes()).hexdigest()}
    for n in (325, 425):
        bb = [b for b in blocks if b["N"] == n]
        names = bb[0]["labels"]
        vectors = []
        count_vectors = []
        for b in bb:
            values = b["prefix_values"]
            assert len(values) == 1000
            vectors.append(np.concatenate([values.mean(axis=0)] + [(values * (b["cells"] == cell)[:, None]).mean(axis=0) for cell in range(9)]))
            count_vectors.append([np.sum(b["cells"] == cell) for cell in range(9)])
        labels = [f"{group}.{label}" for group in ["all"] + [f"cell{a}{b}" for a in range(3) for b in range(3)] for label in names]
        x = np.asarray(vectors)
        mean = x.mean(axis=0)
        derived = derive(mean, labels)
        loo = np.array([list(derive((20 * mean - row) / 19, labels).values()) for row in x])
        factors = np.sqrt(19 / 20) * (loo - loo.mean(axis=0))
        # Reproduce every original-batch orientation-basis matrix from old S/D output.
        old = common["sizes"][str(n)]
        oldx = np.asarray(old["joint_20_batch_means"])
        ix = {label: i for i, label in enumerate(old["labels"])}
        error = 0.
        for obs in OBS[:-1]:
            sp, sm, dp, dm = [oldx[:, ix[f"all.{mark}->{out}.{obs}"]] for mark, out in (("plus", "S"), ("minus", "S"), ("plus", "D"), ("minus", "D"))]
            target = np.stack([sp + sm + DELTA[n] * (dp + dm) / 2,
                               sp - sm + DELTA[n] * (dp - dm) / 2,
                               sp + sm - DELTA[n] * (dp + dm) / 2,
                               sp - sm - DELTA[n] * (dp - dm) / 2], axis=-1)
            actual = x[:, [labels.index(f"all.{obs}.mean_J[{o},{s}]") for o in ("first", "second") for s in ("first", "second")]]
            error = max(error, float(np.max(np.abs(actual - target))))
        assert error < 2e-15, error
        np.savez_compressed(args.output / f"prefix_statistics_N{n}.npz", counter=np.concatenate([b["counter"] for b in bb]),
                            batch=np.repeat(np.arange(20), 1000), cell=np.concatenate([b["cells"] for b in bb]),
                            labels=np.array(names), values=np.concatenate([b["prefix_values"] for b in bb]))
        result["sizes"][str(n)] = {"batch_ids": list(range(20)), "labels": list(derived), "estimate": list(derived.values()),
             "se": np.linalg.norm(factors, axis=0).tolist(), "LOO": loo.tolist(), "factor": factors.tolist(),
             "base_labels": labels, "joint_20_batch_means": x.tolist(),
             "cell_prefix_counts_by_batch": np.asarray(count_vectors).tolist(),
             "active_quartet_histogram_by_cell": sum(b["active_quartet_histogram_by_cell"] for b in bb).tolist(),
             "ensemble_J_reproduction_max_abs_batch_error": error,
             "four_quartet_support": {
                 "prefixes_with_any_pairing": int(sum(np.sum(b["prefix_values"][:, names.index("diagnostic.at_least_one_usable_four_pairing")]) for b in bb)),
                 "batches_with_any_pairing": int(sum(np.any(b["prefix_values"][:, names.index("diagnostic.at_least_one_usable_four_pairing")]) for b in bb)),
                 "interpretation_gate": "A zero/low-support estimate and zero batch SE cannot establish local rank one; inspect this support before interpreting det-square."}}
        print("N", n, "ensemble J reproduction max error", error, flush=True)
        for label, val, se in zip(derived, derived.values(), np.linalg.norm(factors, axis=0)):
            if label.startswith("all.") and ("E_det" in label or "det_mean" in label or "usable" in label or "noncollinear" in label):
                print(label, f"{val:.12g} +/- {se:.6g}", flush=True)
    (args.output / "score.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    receipt = {"started_unix": started, "finished_unix": time.time(), "elapsed_seconds": time.time() - started,
               "hostname": platform.node(), "python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
               "cpu_count": os.cpu_count(), "workers": args.workers,
               "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               "input_manifest_sha256": result["manifest_sha256"], "command": "analyze_prefix_local_rank.py --workers " + str(args.workers),
               "new_samples": 0, "successful_analysis_runs": 1, "attempt": args.attempt,
               "prior_attempts": "Attempt1 completed numeric calculations but failed JSON int64 serialization" if args.attempt == 2 else "none"}
    (args.output / "run_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt), flush=True)

if __name__ == "__main__":
    main()
