#!/usr/bin/env python3
"""Generate the n=12 series-parallel layer in parallel pieces, then merge.

Why: one-shot generation needs ~8-15 minutes of canonical-key work, which exceeds a
single command's lifetime. The n=12 cartesian product is cut into 18 independent
pieces (12 series pairs (n1, 11-n1) + 6 parallel pairs n1+n2=12, n1<=n2). Each piece
canonical-dedups in-slice and writes a pickle; the merge step then canonical-dedups
across pieces with a process pool and writes cache_sp_n12.pkl.

Usage:
  python3 sp12_gen_pieces.py --generate [--procs N]   # parallel piece generation
  python3 sp12_gen_pieces.py --merge    [--procs N]   # across-piece dedup -> layer
"""
from __future__ import annotations

import argparse
import pickle
import sys
import time
from multiprocessing import get_context
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bounded_summary_search import (  # noqa: E402
    Net,
    canonical_key,
    connected_carrier,
    parallel,
    series,
)
from sp_gen12 import layer_path  # noqa: E402

HERE = Path(__file__).resolve().parent
PIECE_DIR = HERE / "pieces_n12"


def piece_path(tag: str) -> Path:
    return PIECE_DIR / f"piece_{tag}.pkl"


def load_layer(n: int) -> list:
    return pickle.loads(layer_path(n).read_bytes())


def gen_one(args) -> dict:
    tag, kind, n1, n2 = args
    A_pool = load_layer(n1) if n1 > 0 else [load_layer(0)[0]]
    B_pool = load_layer(n2) if n2 > 0 else [load_layer(0)[0]]
    seen = {}
    out = []
    cands = 0
    if kind == "series":
        for A in A_pool:
            for B in B_pool:
                cands += 1
                g = series(A, B, f"ser_n12")
                if g.lr_edge or not connected_carrier(g):
                    continue
                k = canonical_key(g)
                if k not in seen:
                    seen[k] = True
                    out.append(g)
    else:  # parallel
        if n1 == n2:
            for i, A in enumerate(A_pool):
                for B in B_pool[i:]:
                    cands += 1
                    g = parallel(A, B, f"par_n12")
                    if g.lr_edge or not connected_carrier(g):
                        continue
                    k = canonical_key(g)
                    if k not in seen:
                        seen[k] = True
                        out.append(g)
        else:
            for A in A_pool:
                for B in B_pool:
                    cands += 1
                    g = parallel(A, B, f"par_n12")
                    if g.lr_edge or not connected_carrier(g):
                        continue
                    k = canonical_key(g)
                    if k not in seen:
                        seen[k] = True
                        out.append(g)
    return {"tag": tag, "cands": cands, "kept": len(out), "nets": out}


def generate(procs: int) -> None:
    PIECE_DIR.mkdir(exist_ok=True)
    tasks = []
    for n1 in range(0, 12):
        tasks.append((f"ser_{n1}_{11 - n1}", "series", n1, 11 - n1))
    for n1 in range(1, 7):
        tasks.append((f"par_{n1}_{12 - n1}", "parallel", n1, 12 - n1))
    print(f"{len(tasks)} pieces, procs={procs}", flush=True)
    t0 = time.time()
    with get_context("fork").Pool(procs) as pool:
        for i, res in enumerate(pool.imap_unordered(gen_one, tasks), 1):
            tag = res["tag"]
            piece_path(tag).write_bytes(pickle.dumps(res["nets"], protocol=5))
            print(f"  [{i}/{len(tasks)}] {tag}: {res['cands']} candidates -> "
                  f"{res['kept']} graphs ({time.time()-t0:.0f}s)", flush=True)
    print(f"generate done in {time.time()-t0:.1f}s", flush=True)


def key_one_chunk(chunk):
    """Module-level worker: canonicalize a chunk of (piece_idx, net) pairs."""
    out = []
    for net in chunk:
        out.append((canonical_key(net), net))
    return out


def merge(procs: int) -> None:
    t0 = time.time()
    files = sorted(PIECE_DIR.glob("piece_*.pkl"))
    print(f"merging {len(files)} pieces", flush=True)

    loaded = [(i, pickle.loads(f.read_bytes())) for i, f in enumerate(files)]
    total = sum(len(nets) for _, nets in loaded)
    print(f"loaded {total} graphs from pieces", flush=True)

    work = []
    for i, nets in loaded:
        work.extend(nets)
    chunk_size = max(2000, len(work) // (procs * 8))
    chunks = [work[i:i + chunk_size] for i in range(0, len(work), chunk_size)]

    seen = {}
    with get_context("fork").Pool(procs) as pool:
        for done, part in enumerate(pool.imap_unordered(key_one_chunk, chunks), 1):
            for canon, net in part:
                if canon not in seen:
                    seen[canon] = net
            if done % 10 == 0:
                print(f"  merged chunks {done}/{len(chunks)} ({len(seen)}) {time.time()-t0:.0f}s",
                      flush=True)
    layer = list(seen.values())
    print(f"merge dedup: {total} -> {len(layer)} in {time.time()-t0:.1f}s", flush=True)
    layer_path(12).write_bytes(pickle.dumps(layer, protocol=5))
    print(f"wrote {layer_path(12)}", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--generate", action="store_true")
    ap.add_argument("--merge", action="store_true")
    ap.add_argument("--procs", type=int, default=9)
    args = ap.parse_args()
    if args.generate:
        generate(args.procs)
    if args.merge:
        merge(args.procs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
