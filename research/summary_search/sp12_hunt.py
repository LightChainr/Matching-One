#!/usr/bin/env python3
"""Hunt the series-parallel family up to n=12 at radius 1, memory-bounded and parallel.

Differences from run_full_search.hunt_full:
  * never keeps safe_table alive (n=12 -> 4096 bools each; 160k graphs would need ~18 GB)
  * both passes run over a process pool
  * pass 1 computes S/H2/b2/neighbourhood/edges only; pass 2 recomputes tables
    only for graphs sitting in a multi-member summary class
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction
from math import comb
from multiprocessing import get_context
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bounded_summary_search import (  # noqa: E402
    FROZEN_EXPERIMENT_ORDER,
    S_coeffs,
    behavior_tuple,
    edge_count,
    experiments_fast,
    first_split,
    neighborhood_key,
    safe_table,
)

HERE = Path(__file__).resolve().parent


def _summary_worker(chunk):
    out = []
    for idx, g in chunk:
        table = safe_table(g)
        s = S_coeffs(table, g.n)
        h2 = g.n - s[1] if g.n >= 1 else 0
        b2 = (comb(g.n, 2) - s[2]) if g.n >= 2 else 0
        out.append(
            (idx, g.n, s, h2, b2, neighborhood_key(g, 1), neighborhood_key(g, 2),
             edge_count(g), g.name)
        )
    return out


def _exps_worker(chunk):
    out = []
    for idx, g in chunk:
        table = safe_table(g)
        out.append((idx, experiments_fast(table, g.n)))
    return out


def chunks(seq, size):
    for i in range(0, len(seq), size):
        yield seq[i:i + size]


def hunt(nets, procs, chunk_size=400, tag="SP", radius=1):
    """radius selects which neighbourhood (1 or 2) enters the summary key.

    Both radii are always computed in pass 1 (cheap); radius only picks the key.
    """
    nets = list(enumerate(nets))
    t0 = time.time()
    recs = {}
    with get_context("fork").Pool(procs) as pool:
        for done, part in enumerate(pool.imap_unordered(_summary_worker, chunks(nets, chunk_size)), 1):
            for row in part:
                recs[row[0]] = row
            if done % 20 == 0:
                print(f"  pass1 chunks {done} ({len(recs)}/{len(nets)}) {time.time()-t0:.0f}s", flush=True)
    print(f"  pass1 done {len(recs)} graphs in {time.time()-t0:.1f}s", flush=True)

    groups = defaultdict(list)
    groups_S = defaultdict(list)
    for idx, n, s, h2, b2, n1, n2, e, name in recs.values():
        neigh = n1 if radius == 1 else n2
        groups[(n, s, h2, b2, neigh)].append(idx)
        groups_S[(n, s)].append(idx)

    hist = Counter(n for idx, n, *_ in recs.values())
    print(
        f"  S-classes={len(groups_S)} multiS={sum(1 for v in groups_S.values() if len(v) > 1)} "
        f"summary_classes={len(groups)} multi_summary={sum(1 for v in groups.values() if len(v) > 1)}",
        flush=True,
    )

    multi_idx = sorted(i for idxs in groups.values() if len(idxs) > 1 for i in idxs)
    print(f"  graphs needing experiments: {len(multi_idx)}", flush=True)
    exps = {}
    if multi_idx:
        want = set(multi_idx)
        work = [(i, g) for i, g in nets if i in want]
        t1 = time.time()
        with get_context("fork").Pool(procs) as pool:
            for done, part in enumerate(pool.imap_unordered(_exps_worker, chunks(work, 64)), 1):
                for idx, e in part:
                    exps[idx] = e
                if done % 20 == 0:
                    print(f"  pass2 chunks {done} ({len(exps)}/{len(work)}) {time.time()-t1:.0f}s", flush=True)
        print(f"  pass2 done in {time.time()-t1:.1f}s", flush=True)

    splits = []
    for key, idxs in groups.items():
        if len(idxs) < 2:
            continue
        behav = defaultdict(list)
        for i in idxs:
            behav[behavior_tuple({"exps": exps[i]})].append(i)
        if len(behav) > 1:
            reps = sorted(
                (v[0] for v in behav.values()),
                key=lambda i: (recs[i][1], recs[i][7], recs[i][8]),
            )
            a, b = reps[0], reps[1]
            ra = {"n": recs[a][1], "S": recs[a][2], "edges": recs[a][7],
                  "name": recs[a][8], "exps": exps[a]}
            rb = {"n": recs[b][1], "S": recs[b][2], "edges": recs[b][7],
                  "name": recs[b][8], "exps": exps[b]}
            splits.append((ra, rb, first_split(ra, rb), len(idxs), len(behav)))
    splits.sort(key=lambda t: (t[0]["n"], t[0]["edges"] + t[1]["edges"]))
    print(f"{tag} r={radius} N={len(nets)} splits={len(splits)} total_t={time.time()-t0:.1f}s", flush=True)
    return {
        "N": len(nets),
        "S_classes": len(groups_S),
        "summary_classes": len(groups),
        "multi_summary": sum(1 for v in groups.values() if len(v) > 1),
        "n_hist": {str(k): v for k, v in sorted(hist.items())},
        "r1_splits": len(splits),
        "splits": splits,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n", type=int, default=12)
    ap.add_argument("--only-n", type=int, default=0,
                    help="hunt a single layer n in isolation (valid because n is part of the summary key)")
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--out", default=str(HERE / "sp12_hunt.json"))
    args = ap.parse_args()

    nets = []
    layers = [args.only_n] if args.only_n else range(1, args.max_n + 1)
    for n in layers:
        p = HERE / f"cache_sp_n{n}.pkl"
        if not p.exists():
            print(f"missing layer {p}", flush=True)
            return 2
        nets.extend(pickle.loads(p.read_bytes()))
    print(f"loaded {len(nets)} SP graphs, layers={list(layers)}", flush=True)

    tag = f"SP_n{args.only_n}" if args.only_n else f"SP_nle{args.max_n}"
    stats = hunt(nets, args.procs, tag=tag)

    by_n = Counter()
    for a, b, sp, sz, nb in stats["splits"]:
        by_n[a["n"]] += 1
    print("splits by n:", dict(sorted(by_n.items())), flush=True)
    print("split experiments:", Counter(sp for _, _, sp, _, _ in stats["splits"]), flush=True)
    for a, b, sp, sz, nb in stats["splits"][:15]:
        gap = abs(a["exps"][sp] - b["exps"][sp])
        print(f"  n={a['n']} |E|={a['edges']}+{b['edges']} class={sz} behaviours={nb} "
              f"split={sp} {a['exps'][sp]} vs {b['exps'][sp]} gap={gap}", flush=True)

    blob = {
        "max_n": args.max_n,
        "N": stats["N"],
        "n_hist": stats["n_hist"],
        "S_classes": stats["S_classes"],
        "summary_classes": stats["summary_classes"],
        "multi_summary": stats["multi_summary"],
        "r1_splits": stats["r1_splits"],
        "splits_by_n": {str(k): v for k, v in sorted(by_n.items())},
        "splits": [
            {
                "n": a["n"],
                "edges": [a["edges"], b["edges"]],
                "class_size": sz,
                "behaviours": nb,
                "split_experiment": sp,
                "p_A": str(a["exps"][sp]),
                "p_B": str(b["exps"][sp]),
                "gap": str(abs(a["exps"][sp] - b["exps"][sp])),
                "A": a["name"],
                "B": b["name"],
                "S": list(a["S"]),
            }
            for a, b, sp, sz, nb in stats["splits"]
        ],
    }
    Path(args.out).write_text(json.dumps(blob, indent=2) + "\n", encoding="utf-8")
    print("wrote", args.out, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
