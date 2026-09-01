#!/usr/bin/env python3
"""Repeat the same local-J estimands with exact census-centered label scores.

This efficiency amendment follows the first estimator's zero four-quartet support.
It reuses every saved label/tail and introduces no new sampling.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import platform
import time
import numpy as np
import scipy
from scipy.stats import binom
from analyze_prefix_local_rank import array_csv, PAIR, DISJOINT, determinant_polarization, algebra_checks, P_REF

OBS = ["p_ref.A", "p_ref.E", "p_integral.A", "p_integral.E", "score_G"]

def read_batch(task):
    root, n, batch, census = task
    root = Path(root)
    hi, raw = array_csv(root / "forks" / f"N{n}.batch{batch:02}.csv.gz")
    ci, contact = array_csv(root / "contact" / f"N{n}.batch{batch:02}.csv.gz")
    raw = raw[np.lexsort(tuple(raw[:, hi[k]] for k in ("replica", "group", "quartet", "counter")))].reshape(1000, 8, 2, 2, -1)
    contact = contact[np.lexsort(tuple(contact[:, ci[k]] for k in ("group", "quartet", "counter")))].reshape(1000, 8, 2, -1)
    for key in ("counter", "quartet", "group", "next_label"):
        assert np.array_equal(raw[..., 0, hi[key]], contact[..., ci[key]])
    counter = raw[:, 0, 0, 0, hi["counter"]]
    ranks = np.stack([raw[:, 0, 0, 0, hi[f"{o}_rank"]] for o in ("first", "second")], axis=-1)
    cells = 3 * ranks[:, 0] + ranks[:, 1]
    k0 = int(raw[0, 0, 0, 0, hi["k0"]])
    vacant = n-k0
    assert np.all(raw[..., hi["k0"]] == k0)
    groups = {}
    for row in census:
        nn, bb, ctr, kk, r0, r1, e0, e1, l0, l1, count = map(int, row)
        assert (nn, bb, kk) == (n, batch, k0) and count > 0
        prefix = groups.setdefault(ctr, {"rank": (r0, r1), "classes": {}})
        assert prefix["rank"] == (r0, r1)
        prefix["classes"].setdefault((e0, e1), []).append((l0, l1, count))
    score_numerator = np.zeros((1000, 8, 2, 2), dtype=np.int64)
    exact_g = np.zeros((1000, 2, 2))
    exact_det = np.zeros(1000)
    exact_rank = np.zeros(1000, dtype=np.int64)
    safe_labels = np.zeros(1000, dtype=np.int64)
    for p, ctr in enumerate(counter):
        prefix = groups[int(ctr)]
        assert prefix["rank"] == tuple(ranks[p])
        gt = np.zeros((2, 2), dtype=np.int64)
        classes = {}
        for key, rows in prefix["classes"].items():
            mass = sum(row[2] for row in rows)
            totals = np.array([sum(row[i]*row[2] for row in rows) for i in (0, 1)], dtype=np.int64)
            safe_labels[p] += mass
            center_check = np.zeros(2, dtype=np.int64)
            support = set()
            for lf, ls, count in rows:
                z = mass*np.array([lf, ls], dtype=np.int64)-totals
                center_check += count*z
                gt += count*np.outer(z, z)
                support.add((lf, ls))
            assert np.all(center_check == 0)
            classes[key] = (mass, totals, support)
        assert safe_labels[p] <= vacant
        gd = int(gt[0, 0])*int(gt[1, 1])-int(gt[0, 1])*int(gt[1, 0])
        assert gd >= 0
        exact_g[p] = gt / float(vacant**3)
        exact_det[p] = gd / float(vacant**6)
        exact_rank[p] = 2 if gd > 0 else int(np.trace(gt) > 0)
        for q in range(8):
            for group in (0, 1):
                row = contact[p, q, group]
                after = tuple(int(row[ci[f"{o}_rank_after"]]) for o in ("first", "second"))
                if after != tuple(ranks[p]):
                    continue
                degrees = tuple(int(row[ci[f"{o}_e"]]) for o in ("first", "second"))
                marks = tuple(int(ranks[p, i] == 0)*int(row[ci[f"{o}_e"]]-row[ci[f"{o}_c"]]) for i, o in enumerate(("first", "second")))
                mass, totals, support = classes[degrees]
                assert marks in support
                score_numerator[p, q, group] = mass*np.asarray(marks)-totals
    # H=(n_class*L-sum_class L)/vacant. E[H|Z]=0 exactly by the integer census.
    dh_numerator = score_numerator[:, :, 0]-score_numerator[:, :, 1]
    dh = dh_numerator / float(vacant)
    tail = binom.sf(np.arange(n+1)-1, n, P_REF)
    responses = []
    for o in ("first", "second"):
        k1, k2 = (raw[..., hi[f"{o}_{k}"]] for k in ("k1", "k2"))
        responses.append(np.stack((tail[k1]+tail[k2]-1, 1-tail[k1]+tail[k2],
                                   1-(k1+k2)/(n+1), 1+(k1-k2)/(n+1)), axis=-1))
    y = np.stack(responses, axis=-1)
    dy = y[:, :, 0].mean(axis=2)-y[:, :, 1].mean(axis=2)
    matrices = dy[..., :, :, None]*dh[..., None, None, :]/2
    gram = dh[..., :, None]*dh[..., None, :]/2
    matrices = np.concatenate((matrices, gram[:, :, None]), axis=2)
    pol = np.stack([determinant_polarization(matrices[:, q], matrices[:, r]) for q, r in PAIR], axis=1)
    udet = pol.mean(axis=1)
    udet2 = sum(pol[:, a]*pol[:, b] for a, b in DISJOINT)/len(DISJOINT)
    energy = sum(np.sum(matrices[:, q]*matrices[:, r], axis=(-2, -1)) for q, r in PAIR)/len(PAIR)
    # Count true source noncollinearity with integers, not floating cancellation.
    noncollinear = np.stack([dh_numerator[:, q, 0]*dh_numerator[:, r, 1]-dh_numerator[:, q, 1]*dh_numerator[:, r, 0] != 0 for q, r in PAIR], axis=1)
    usable = sum(noncollinear[:, a] & noncollinear[:, b] for a, b in DISJOINT)
    assert np.max(np.abs(udet[cells != 0])) == 0
    assert np.max(np.abs(udet2[cells != 0])) == 0
    mean_j = matrices.mean(axis=1)
    fields = {}
    for j, obs in enumerate(OBS):
        fields[obs+".E_det_JZ"] = udet[:, j]
        fields[obs+".E_det_JZ_squared"] = udet2[:, j]
        fields[obs+".E_frobenius_JZ_squared"] = energy[:, j]
        for i, o in enumerate(("first", "second")):
            for k, s in enumerate(("first", "second")):
                fields[obs+f".mean_J[{o},{s}]"] = mean_j[:, j, i, k]
    fields["exact_score_G.det"] = exact_det
    fields["exact_score_G.det_squared"] = exact_det**2
    fields["exact_score_G.rank2_prefix"] = exact_rank == 2
    fields["exact_score_G.rank1_prefix"] = exact_rank == 1
    fields["exact_score_G.rank0_prefix"] = exact_rank == 0
    for i, o in enumerate(("first", "second")):
        for j, s in enumerate(("first", "second")):
            fields[f"exact_score_G.matrix[{o},{s}]"] = exact_g[:, i, j]
    fields["diagnostic.active_quartets"] = np.any(dh_numerator != 0, axis=-1).sum(axis=1)
    fields["diagnostic.noncollinear_source_quartet_pairs"] = noncollinear.sum(axis=1)
    fields["diagnostic.usable_four_quartet_pairings"] = usable
    fields["diagnostic.any_noncollinear_pair"] = noncollinear.any(axis=1)
    fields["diagnostic.any_usable_four_pairing"] = usable > 0
    fields["diagnostic.safe_label_fraction"] = safe_labels/float(vacant)
    return {"N": n, "batch": batch, "counter": counter, "cell": cells,
            "labels": list(fields), "values": np.column_stack(list(fields.values()))}

def derive(v, labels):
    out = dict(zip(labels, v))
    for group in ["all"]+[f"cell{a}{b}" for a in range(3) for b in range(3)]:
        for obs in OBS:
            base = group+"."+obs
            j = np.array([[out[base+f".mean_J[{o},{s}]"] for s in ("first", "second")] for o in ("first", "second")])
            out[base+".det_mean_JZ"] = float(np.linalg.det(j))
            out[base+".ensemble_minus_mean_local_det"] = float(np.linalg.det(j))-out[base+".E_det_JZ"]
    return out

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=Path(__file__).parent/"inputs")
    parser.add_argument("--census", type=Path, default=Path(__file__).parent/"census")
    parser.add_argument("--pair-results", type=Path, default=Path(__file__).parent/"results/score.json")
    parser.add_argument("--output", type=Path, default=Path(__file__).parent/"results-exact-score")
    parser.add_argument("--workers", type=int, default=16)
    args = parser.parse_args()
    started = time.time()
    args.output.mkdir(parents=True, exist_ok=False)
    manifest = json.loads((args.inputs/"manifest.json").read_text())
    for entry in manifest["files"]:
        assert hashlib.sha256((args.inputs/entry["local_path"]).read_bytes()).hexdigest() == entry["sha256"]
    census_manifest = json.loads((args.census/"manifest.json").read_text())
    for entry in census_manifest["files"]:
        assert hashlib.sha256((args.census/entry["local_path"]).read_bytes()).hexdigest() == entry["sha256"]
    census_hashes, tasks = {}, []
    for n in (325, 425):
        path = args.census/f"N{n}/census.csv.gz"
        census_hashes[str(n)] = hashlib.sha256(path.read_bytes()).hexdigest()
        header, data = array_csv(path)
        assert list(header) == ["N", "batch", "counter", "k0", "first_rank", "second_rank", "first_e", "second_e", "L_first", "L_second", "count"]
        assert len(np.unique(data[:, 2])) == 20000
        for b in range(20):
            tasks.append((str(args.inputs), n, b, data[data[:, 1] == b]))
    algebra = algebra_checks()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        blocks = list(pool.map(read_batch, tasks))
    pair = json.loads(args.pair_results.read_text())
    result = {"estimand": "Same J(Z), E det J(Z), E[(det J(Z))^2] as pair-difference plan; exact-census efficiency amendment after no four-quartet support was observed.",
              "H_definition": "H_s(u)=pi_a*(L_s(u)-mean_a L_s) on joint-safe class a; zero outside. Integer numerator=(class_count*L_s-class_sum_L_s), denominator=vacant count.",
              "X_definition": "X_q[o,s]=(H_s(U)-H_s(V))*(Ybar_o(U)-Ybar_o(V))/2. All label pairs, including different classes. E H=0 exactly, so E X=original J.",
              "input_Gram_distinction": "score_G=E[HH^T]=sum_a pi_a^3 Cov(L|a), unlike original masked input_G=sum_a pi_a^2 Cov(L|a). Same source-null directions; numeric values must not be equated.",
              "dependence": "Same original20 paired batches per N; no new prefixes/labels/tails. Census is deterministic enumeration, not an independent sample. Four distinct quartets for determinant-square.",
              "uncertainty": "Delete-one original batch; signed unbiased det-square retained; no normal p-values or rank-one acceptance from nonresolution. Conditional on prefix, exact centering does not couple the eight quartet random streams.",
              "new_samples": 0, "algebra_checks": algebra, "census_sha256": census_hashes, "sizes": {}}
    for n in (325, 425):
        bb = [b for b in blocks if b["N"] == n]
        names = bb[0]["labels"]
        labels = [f"{g}.{s}" for g in ["all"]+[f"cell{a}{b}" for a in range(3) for b in range(3)] for s in names]
        x = np.array([np.concatenate([b["values"].mean(0)]+[(b["values"]*(b["cell"] == c)[:, None]).mean(0) for c in range(9)]) for b in bb])
        mean = x.mean(0)
        derived = derive(mean, labels)
        loo = np.array([list(derive((20*mean-row)/19, labels).values()) for row in x])
        factors = np.sqrt(19/20)*(loo-loo.mean(0))
        ps = pair["sizes"][str(n)]
        px = np.array(ps["joint_20_batch_means"])
        difference = {}
        for label in labels:
            if label.startswith("all.") and ".mean_J[" in label and ".score_G." not in label:
                z = x[:, labels.index(label)]-px[:, ps["base_labels"].index(label)]
                difference[label] = {"estimate": float(z.mean()), "se": float(z.std(ddof=1)/np.sqrt(20)), "batch_differences": z.tolist()}
        for a, b in (("score_G.E_det_JZ", "exact_score_G.det"), ("score_G.E_det_JZ_squared", "exact_score_G.det_squared")):
            z = x[:, labels.index("all."+a)]-x[:, labels.index("all."+b)]
            difference["sample_Gram_vs_exact_census."+a] = {"estimate": float(z.mean()), "se": float(z.std(ddof=1)/np.sqrt(20)), "batch_differences": z.tolist()}
        np.savez_compressed(args.output/f"prefix_statistics_N{n}.npz", counter=np.concatenate([b["counter"] for b in bb]), batch=np.repeat(np.arange(20), 1000),
                            cell=np.concatenate([b["cell"] for b in bb]), labels=np.array(names), values=np.concatenate([b["values"] for b in bb]))
        result["sizes"][str(n)] = {"batch_ids": list(range(20)), "labels": list(derived), "estimate": list(derived.values()), "se": np.linalg.norm(factors, axis=0).tolist(),
            "LOO": loo.tolist(), "factor": factors.tolist(), "base_labels": labels, "joint_20_batch_means": x.tolist(), "paired_comparison": difference,
            "cell_counts": np.bincount(np.concatenate([b["cell"] for b in bb]), minlength=9).tolist(),
            "four_quartet_prefix_support": int(sum(np.sum(b["values"][:, names.index("diagnostic.any_usable_four_pairing")]) for b in bb)),
            "four_quartet_batch_support": int(sum(np.any(b["values"][:, names.index("diagnostic.any_usable_four_pairing")]) for b in bb))}
        for label, val, err in zip(derived, derived.values(), np.linalg.norm(factors, axis=0)):
            if label.startswith("all.") and ("E_det_JZ" in label or "any_usable" in label or "rank2_prefix" in label or "exact_score_G.det" in label):
                print(n, label, f"{val:.12g} +/- {err:.6g}", flush=True)
    (args.output/"score.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    receipt = {"hostname": platform.node(), "python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, "workers": args.workers, "cpu_count": os.cpu_count(),
               "started_unix": started, "finished_unix": time.time(), "elapsed_seconds": time.time()-started,
               "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "shared_script_sha256": hashlib.sha256(Path(__file__).with_name("analyze_prefix_local_rank.py").read_bytes()).hexdigest(),
               "census_sha256": census_hashes, "new_prefixes": 0, "new_labels": 0, "new_tails": 0, "analysis_attempt": 1}
    (args.output/"run_receipt.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print(json.dumps(receipt), flush=True)

if __name__ == "__main__":
    main()
