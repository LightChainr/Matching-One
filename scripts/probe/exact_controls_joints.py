#!/usr/bin/env python3
"""
C4 + C5.

C4 — joint (n_black, r_b, r_w) and Alexander.
  full joint histograms (n_black, r_b) and (n_black, r_w); confirm r_b+r_w=2
  on every config (failure list if any); six in-repo hand configs plus three
  NEW hand configs declared before running.

C5 — two algorithms, every rank.
  per-rank confusion table (L=3, L=4); max |r_A - r_B| (must be 0); runtime.
"""
import sys
import json
import time
from pathlib import Path
from fractions import Fraction
from itertools import product
from collections import defaultdict

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "homological_balance"))
from exact_torus_enum import ambient_rank, neighbors  # noqa: E402
from verify_independent import ambient_rank_v2  # noqa: E402


# ---- three NEW hand configs, predictions declared BEFORE running ----
# (i,j) in {0..L-1}^2.  L=3.
NEW_HAND = {
    "2x2_block":  {(0, 0), (0, 1), (1, 0), (1, 1)},
    "checkerboard": {(i, j) for i in range(3) for j in range(3) if (i + j) % 2 == 0},
    "two_horizontal_rings": {(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)},
}
# predicted (r_b, r_w) via Alexander r_b+r_w=2 and direct inspection
NEW_HAND_PREDICT = {
    "2x2_block": (0, 2),
    "checkerboard": (1, 1),
    "two_horizontal_rings": (1, 1),
}


def enumerate_full(L):
    """All configs: (n_black, r_b, r_w), dual failures, runtime."""
    N = L * L
    jb = defaultdict(int)
    jw = defaultdict(int)
    dual_fail = []
    t0 = time.time()
    for mask in range(1 << N):
        black = [(k // L, k % L) for k in range(N) if (mask >> k) & 1]
        n = len(black)
        rb = ambient_rank(black, L, False)
        rw = ambient_rank([(k // L, k % L) for k in range(N) if not ((mask >> k) & 1)],
                          L, True)
        jb[(n, rb)] += 1
        jw[(n, rw)] += 1
        if rb + rw != 2:
            dual_fail.append((mask, rb, rw))
    dt = time.time() - t0
    return jb, jw, dual_fail, dt


def confusion_table(L):
    """per-rank confusion between ambient_rank (A) and ambient_rank_v2 (B)."""
    N = L * L
    conf = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    maxdiff = 0
    mismatch = 0
    first_bad = None
    t0 = time.time()
    for mask in range(1 << N):
        black = [(k // L, k % L) for k in range(N) if (mask >> k) & 1]
        rA = ambient_rank(black, L, False)
        rB = ambient_rank_v2(black, L, False)
        conf[rA][rB] += 1
        d = abs(rA - rB)
        maxdiff = max(maxdiff, d)
        if d != 0:
            mismatch += 1
            if first_bad is None:
                first_bad = (mask, black, rA, rB)
    dt = time.time() - t0
    return conf, maxdiff, mismatch, first_bad, dt


def main():
    out = {}
    for L in (3, 4):
        jb, jw, dual_fail, dt_enum = enumerate_full(L)
        conf, maxdiff, mismatch, first_bad, dt_conf = confusion_table(L)
        out[str(L)] = {
            "sites": L * L,
            "configs": 2 ** (L * L),
            "joint_black_count_rank": {f"{n},{r}": jb[(n, r)] for (n, r) in sorted(jb)},
            "joint_white_count_rank": {f"{n},{r}": jw[(n, r)] for (n, r) in sorted(jw)},
            "dual_fail_count": len(dual_fail),
            "dual_fail_examples": dual_fail[:5],
            "confusion_rA_rB": conf,
            "max_abs_rank_diff": maxdiff,
            "mismatch_count": mismatch,
            "first_mismatch": first_bad,
            "enumerate_seconds": round(dt_enum, 3),
            "confusion_seconds": round(dt_conf, 3),
        }
        print(f"===== L={L} =====")
        print(f"configs={2**(L*L)}  dual_fail={len(dual_fail)}")
        print(f"joint (n,r_b) entries: {len(jb)}   joint (n,r_w) entries: {len(jw)}")
        print(f"confusion rA\\rB rows: {conf}")
        print(f"max|rA-rB|={maxdiff}  mismatches={mismatch}  first={first_bad}")
        print(f"enumerate {dt_enum:.2f}s  confusion {dt_conf:.2f}s")
        print()

    # hand configs: six in-repo (re-run here) + three new
    print("===== hand configs (L=3) =====")
    hand = {}
    allpts = set(product(range(3), range(3)))
    inrepo = {
        "full_black": set(allpts),
        "empty": set(),
        "single_point": {(0, 0)},
        "row0_horiz_ring": {(0, 0), (0, 1), (0, 2)},
        "col0_vert_ring": {(0, 0), (1, 0), (2, 0)},
        "diagonal": {(0, 0), (1, 1), (2, 2)},
        "cross_row0_col0": {(0, 0), (0, 1), (0, 2), (1, 0), (2, 0)},
    }
    for name, B in {**inrepo, **NEW_HAND}.items():
        rb = ambient_rank(B, 3, False)
        rw = ambient_rank(allpts - B, 3, True)
        hand[name] = {"r_b": rb, "r_w": rw, "sum": rb + rw}
        pred = NEW_HAND_PREDICT.get(name)
        mark = ""
        if pred is not None:
            mark = "NEW" if pred == (rb, rw) else f"NEW-PREDICT-MISMATCH(expected {pred})"
        print(f"  {name:22s} r_b={rb} r_w={rw} sum={rb+rw}  {mark}")
    out["hand_configs_L3"] = hand

    base = ROOT / "results" / "probe-exact-controls"
    base.mkdir(parents=True, exist_ok=True)
    (base / "joints.json").write_text(json.dumps({
        "schema": "matching-one.probe-exact-controls.joints.v1",
        "by_L": {k: out[k] for k in ("3", "4")},
        "hand_configs_L3": hand,
    }, indent=2))
    print("\nwrote", base / "joints.json")


if __name__ == "__main__":
    main()
