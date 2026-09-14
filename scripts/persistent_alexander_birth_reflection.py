"""Finite controls for persistent digital-Alexander birth reflection.

For an LxL honest square torus, black sites use NN connectivity and white sites
use the NN+diagonal matching graph.  We check configuration-level rank/count/
rank-one-line identities.  At L=3 we also exhaust all 9! label orders; at
larger L the birth-path control uses fixed-seed random permutations after
precomputing all mask ranks.

This is a finite control, not the proof.  The birth-index reflection follows
algebraically from r4(S)+r8(S^c)=2 along nested occupied sets S_k.
"""
from __future__ import annotations

import argparse
import json
import random
from itertools import permutations
from math import factorial, gcd
from pathlib import Path


def run(L: int, permutation_samples: int, seed: int):
    N = L * L
    NN = [(1, 0), (0, 1)]
    G8 = [(1, 0), (0, 1), (1, 1), (1, -1)]
    full = (1 << N) - 1

    def vid(x, y):
        return (y % L) * L + (x % L)

    def xy(v):
        return (v % L, v // L)

    def canon_line(v):
        a, b = v
        if a == 0 and b == 0:
            return None
        g = gcd(abs(a), abs(b))
        a //= g
        b //= g
        if a < 0 or (a == 0 and b < 0):
            a, b = -a, -b
        return (a, b)

    def rank_and_components(mask, steps):
        occ = [bool(mask >> i & 1) for i in range(N)]
        adj = [[] for _ in range(N)]
        for u in range(N):
            if not occ[u]:
                continue
            x, y = xy(u)
            for dx, dy in steps:
                v = vid(x + dx, y + dy)
                if occ[v]:
                    adj[u].append((v, (dx, dy)))
                    adj[v].append((u, (-dx, -dy)))

        seen = [False] * N
        all_gains = []
        winding_components = 0
        for root in range(N):
            if not occ[root] or seen[root]:
                continue
            seen[root] = True
            lift = {root: (0, 0)}
            stack = [root]
            gains = []
            while stack:
                u = stack.pop()
                ux, uy = lift[u]
                for v, (dx, dy) in adj[u]:
                    cand = (ux + dx, uy + dy)
                    if not seen[v]:
                        seen[v] = True
                        lift[v] = cand
                        stack.append(v)
                    else:
                        vx, vy = lift[v]
                        defect = (cand[0] - vx, cand[1] - vy)
                        if defect != (0, 0):
                            assert defect[0] % L == 0 and defect[1] % L == 0
                            gains.append((defect[0] // L, defect[1] // L))
            nz = [g for g in gains if g != (0, 0)]
            if nz:
                winding_components += 1
            all_gains.extend(nz)

        rank = 0
        line = None
        if all_gains:
            rank = 1
            line = canon_line(all_gains[0])
            a, b = all_gains[0]
            for c, d in all_gains[1:]:
                if a * d - b * c != 0:
                    rank = 2
                    line = None
                    break
        return rank, winding_components, line

    size = 1 << N
    if N > 16:
        raise ValueError("This finite control intentionally caps exhaustive configuration enumeration at N<=16.")

    rank4 = [0] * size
    rank8_direct = [0] * size
    for mask in range(size):
        rank4[mask] = rank_and_components(mask, NN)[0]
        rank8_direct[mask] = rank_and_components(mask, G8)[0]

    violations = []
    rank_counts = [0, 0, 0]
    rank1_count = 0
    pair_counts = {}
    for mask in range(size):
        comp = full ^ mask
        r4, w4, l4 = rank_and_components(mask, NN)
        r8, w8, l8 = rank_and_components(comp, G8)
        pair_counts[(w4, w8)] = pair_counts.get((w4, w8), 0) + 1
        rank_counts[r4] += 1
        if r4 == 1:
            rank1_count += 1
        if r4 + r8 != 2:
            violations.append(("rank", mask, r4, r8))
        if w4 - w8 != r4 - 1:
            violations.append(("component-count", mask, w4, w8, r4))
        if r4 == 1 and l4 != l8:
            violations.append(("rank1-line", mask, l4, l8))
        if violations:
            break
    if violations:
        raise AssertionError(violations[:3])

    def births(perm, rank_table):
        mask = 0
        k1 = k2 = None
        for k, u in enumerate(perm, 1):
            mask |= 1 << u
            r = rank_table[mask]
            if k1 is None and r >= 1:
                k1 = k
            if k2 is None and r >= 2:
                k2 = k
        return k1, k2

    rng = random.Random(seed)
    if N <= 9:
        iterator = permutations(range(N))
        nperm = factorial(N)
        mode = "exhaustive"
    else:
        def samples():
            base = list(range(N))
            for _ in range(permutation_samples):
                a = base.copy()
                rng.shuffle(a)
                yield tuple(a)
        iterator = samples()
        nperm = permutation_samples
        mode = f"random(seed={seed})"

    birth_viol = 0
    for perm in iterator:
        k1, k2 = births(perm, rank4)
        j1, j2 = births(tuple(reversed(perm)), rank8_direct)
        if (j1, j2) != (N + 1 - k2, N + 1 - k1):
            birth_viol += 1
            break
    if birth_viol:
        raise AssertionError("birth reflection violation")

    return {
        "L": L,
        "N": N,
        "configurations_checked": size,
        "permutation_mode": mode,
        "permutations_checked": nperm,
        "digital_alexander_rank_identity": True,
        "winding_component_count_identity": True,
        "rank1_homology_line_matches_dual": True,
        "birth_index_relation": "K1_matching(reverse)=N+1-K2_NN; K2_matching(reverse)=N+1-K1_NN",
        "birth_index_violations": 0,
        "all_checks_passed": True,
        "rank1_configurations": rank1_count,
        "rank_counts_NN": {str(i): rank_counts[i] for i in range(3)},
        "winding_component_pair_counts": {
            f"{a},{b}": n for (a, b), n in sorted(pair_counts.items())
        },
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=int, default=3)
    ap.add_argument("--permutation-samples", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=20260914)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    result = run(args.L, args.permutation_samples, args.seed)
    text = json.dumps(result, indent=2)
    if args.output:
        args.output.write_text(text + "\n")
    print(text)
