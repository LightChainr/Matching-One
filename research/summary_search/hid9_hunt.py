#!/usr/bin/env python3
"""Hunt the HID family over the FULL contracted hop grid, radius 1 and radius 2.

The frozen contract (PROMPT sec.1.4 and the stopping rule sec.9.2) declares HID as
"corridor hiding on both sides with hop in {1,2,3}, n<=12". The shipped certificate
only covered the two symmetric combinations L1R1 and L2R2. This script sweeps all
nine (Lh,Rh) pairs so the declared class is actually exhausted.

Corridor convention is the one the primary witness uses (run_full_search.hops_hide):
    hops_hide(core, Lh, Rh) = series(path(Lh-1), series(core, path(Rh-1)))
    path(0) = unit edge (n=0)          =>  n = core.n + Lh + Rh
so L1R1 on an n=5 core gives n=7, matching the certified witness.
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bounded_summary_search import (  # noqa: E402
    canonical_key,
    connected_carrier,
    path_net,
    series,
)
from sp12_hunt import hunt  # noqa: E402

HERE = Path(__file__).resolve().parent
HOPS = (1, 2, 3)


def hops_hide(core, Lh: int, Rh: int):
    g = series(path_net(Lh - 1), core)
    g = series(g, path_net(Rh - 1))
    g.name = f"hide(L{Lh},{core.name},R{Rh})"
    return g


def build_layer(cores, Lh: int, Rh: int, max_n: int = 12):
    out = []
    seen = set()
    for core in cores:
        if core.n + Lh + Rh > max_n:
            continue
        g = hops_hide(core, Lh, Rh)
        if g.n > max_n or g.lr_edge or not connected_carrier(g):
            continue
        k = canonical_key(g)
        if k in seen:
            continue
        seen.add(k)
        out.append(g)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=6)
    ap.add_argument("--radii", default="1,2")
    ap.add_argument("--out", default=str(HERE / "hid9_hunt.json"))
    args = ap.parse_args()
    radii = [int(x) for x in args.radii.split(",") if x]

    cache = HERE / "cache_exh.pkl"
    cores = pickle.loads(cache.read_bytes())
    print(f"loaded {len(cores)} exhaustive cores (n<=5)", flush=True)

    result = {}
    for Lh in HOPS:
        for Rh in HOPS:
            t0 = time.time()
            nets = build_layer(cores, Lh, Rh)
            tag = f"HID_L{Lh}R{Rh}"
            print(f"\n=== {tag}: {len(nets)} graphs (n from {min(g.n for g in nets)} "
                  f"to {max(g.n for g in nets)}) ===", flush=True)
            entry = {"N": len(nets),
                     "n_hist": dict(sorted(Counter(g.n for g in nets).items()))}
            for r in radii:
                st = hunt(nets, args.procs, tag=f"{tag}", radius=r)
                entry[f"r{r}"] = {
                    "splits": st["r1_splits"],
                    "summary_classes": st["summary_classes"],
                    "multi_summary": st["multi_summary"],
                    "S_classes": st["S_classes"],
                    "splits_detail": [
                        {
                            "n": s[0]["n"],
                            "edges": [s[0]["edges"], s[1]["edges"]],
                            "class_size": s[3],
                            "behaviours": s[4],
                            "split_experiment": s[2],
                            "p_A": str(s[0]["exps"][s[2]]),
                            "p_B": str(s[1]["exps"][s[2]]),
                            "gap": str(abs(s[0]["exps"][s[2]] - s[1]["exps"][s[2]])) if s[2] else None,
                            "A": s[0]["name"],
                            "B": s[1]["name"],
                            "S": list(s[0]["S"]),
                        }
                        for s in st["splits"]
                    ],
                }
            entry["wall_s"] = round(time.time() - t0, 1)
            result[f"L{Lh}R{Rh}"] = entry
            Path(args.out).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
            print(f"  wrote {args.out} ({tag} done in {entry['wall_s']}s)", flush=True)

    print("\n=== SUMMARY ===", flush=True)
    for key, e in result.items():
        line = f"{key}: N={e['N']}"
        for r in radii:
            line += f"  r{r}_splits={e[f'r{r}']['splits']}"
        print(line, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
