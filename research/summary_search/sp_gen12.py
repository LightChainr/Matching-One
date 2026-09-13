#!/usr/bin/env python3
"""Generate two-terminal series-parallel graphs up to n=12, layer by layer, with per-layer cache.

Same generator contract as bounded_summary_search.generate_spsp, but resumable:
each layer n is written to cache_sp_n{n}.pkl so an interrupted run continues.
"""
from __future__ import annotations

import pickle
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bounded_summary_search import (  # noqa: E402
    canonical_key,
    connected_carrier,
    path_net,
    series,
    parallel,
)

HERE = Path(__file__).resolve().parent


def layer_path(n: int) -> Path:
    return HERE / f"cache_sp_n{n}.pkl"


def build(max_n: int = 12) -> dict:
    by_n: dict[int, list] = {0: [path_net(0)]}
    seen: dict[int, set] = {0: {canonical_key(path_net(0))}}
    p1 = path_net(1)
    by_n[1] = [p1]
    seen[1] = {canonical_key(p1)}
    # n=0/1 are seeded, not produced by the loop, so cache them explicitly.
    for seed_n in (0, 1):
        if not layer_path(seed_n).exists():
            layer_path(seed_n).write_bytes(pickle.dumps(by_n[seed_n], protocol=5))

    def add(net) -> bool:
        n = net.n
        if n > max_n or n <= 0:
            return False
        if net.lr_edge:
            return False
        if not connected_carrier(net):
            return False
        key = canonical_key(net)
        bucket = seen.setdefault(n, set())
        if key in bucket:
            return False
        bucket.add(key)
        by_n.setdefault(n, []).append(net)
        return True

    for n in range(2, max_n + 1):
        if layer_path(n).exists():
            cached = pickle.loads(layer_path(n).read_bytes())
            by_n[n] = cached
            seen[n] = {canonical_key(g) for g in cached}
            print(f"n={n}: loaded {len(cached)} from cache", flush=True)
            continue
        t0 = time.time()
        cands = 0
        for n1 in range(0, n):
            n2 = n - 1 - n1
            if n2 < 0:
                continue
            for A in by_n.get(n1, []):
                for B in by_n.get(n2, []):
                    cands += 1
                    add(series(A, B, f"ser_n{n}"))
        for n1 in range(1, n):
            n2 = n - n1
            if n1 > n2:
                break
            pool1 = by_n.get(n1, [])
            pool2 = by_n.get(n2, [])
            if n1 == n2:
                # parallel is symmetric: only i <= j to halve the cartesian product
                for i, A in enumerate(pool1):
                    for B in pool2[i:]:
                        cands += 1
                        add(parallel(A, B, f"par_n{n}"))
            else:
                for A in pool1:
                    for B in pool2:
                        cands += 1
                        add(parallel(A, B, f"par_n{n}"))
        layer_path(n).write_bytes(pickle.dumps(by_n[n], protocol=5))
        print(
            f"n={n}: {len(by_n[n])} graphs from {cands} candidates "
            f"in {time.time()-t0:.1f}s",
            flush=True,
        )
    return by_n


if __name__ == "__main__":
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    t = time.time()
    out = build(limit)
    total = sum(len(v) for k, v in out.items() if k > 0)
    print(f"TOTAL n<= {limit}: {total} graphs in {time.time()-t:.1f}s", flush=True)
