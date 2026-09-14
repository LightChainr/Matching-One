#!/usr/bin/env python3
"""n325rec step 4 -- COST PROBE (spec 2.3B cost gate).

Uses the REFERENCE automaton (PATH B's own `step`/`transformed_edges`) and
reports, per (geometry, sector): the safe-state count if the build finishes,
otherwise the count reached when the state cap or the wall-clock cap trips.

Also measures the size law along two families that pass through the N=325
targets --  u=(1,k), n=1  (memory ~ k+1)  and  u=(2,3), n=1..5 (the (10,15)
orientation realizes as (2,3) with n=5) -- so the target size can be
extrapolated instead of guessed.
"""
from __future__ import annotations

import json
import math
import sys
import time
from collections import deque

sys.path.insert(0, "/workspace/dpfloor/scripts")
import fl_path_oblique as B            # noqa: E402


def build_probe(width, direction, matching, cap, tcap):
    comp, edges, memory = B.transformed_edges(direction, matching)
    start = B.empty_state(width, memory)
    states = [start]
    index = {start: 0}
    queue = deque([start])
    t0 = time.time()
    while queue:
        if time.time() - t0 > tcap:
            return None, len(states), time.time() - t0, "wallclock-cap"
        s = queue.popleft()
        for mask in range(1 << width):
            nxt = B.step(s, mask, width, edges, memory)
            if nxt is None:
                continue
            if nxt not in index:
                if len(states) >= cap:
                    return None, len(states), time.time() - t0, "state-cap"
                index[nxt] = len(states)
                states.append(nxt)
                queue.append(nxt)
    return states, len(states), time.time() - t0, "done"


TARGETS = [
    ("n325_1_18",   (1, 18), 1, 18.027756377319946),
    ("n325_6_17",   (6, 17), 1, 18.027756377319946),
    ("n325_2_3_n5", (2, 3),  5, 18.027756377319946),
    # the same three "directions" at a SMALLER index/ell, as a size law
    ("f_1_18_n1",   (1, 18), 1, None),
    ("f_1_12_n1",   (1, 12), 1, None),
    ("f_1_10_n1",   (1, 10), 1, None),
    ("f_1_8_n1",    (1, 8),  1, None),
    ("f_1_6_n1",    (1, 6),  1, None),
    ("f_2_3_n4",    (2, 3),  4, None),
    ("f_2_3_n3",    (2, 3),  3, None),
    ("f_2_3_n2",    (2, 3),  2, None),
    ("f_3_4_n1",    (3, 4),  1, None),
    ("f_3_4_n2",    (3, 4),  2, None),
    ("f_3_4_n3",    (3, 4),  3, None),
]


def main():
    cap = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000
    tcap = float(sys.argv[2]) if len(sys.argv) > 2 else 180.0
    outp = sys.argv[3] if len(sys.argv) > 3 else \
        "/workspace/n325rec/out/s4_cost.json"
    only = sys.argv[4].split(",") if len(sys.argv) > 4 else None
    res = {"schema": "n325rec.cost.v1", "state_cap": cap, "wall_cap_s": tcap,
           "records": {}}
    for tag, u, n, _ell in TARGETS:
        if only and tag not in only:
            continue
        rec = {"tag": tag, "direction": list(u), "n": n,
               "ell": n * math.hypot(*u)}
        for sect, matching in (("G4", False), ("G8", True)):
            _st, cnt, secs, why = build_probe(n, u, matching, cap, tcap)
            rec[sect] = {"states": cnt, "seconds": round(secs, 2),
                         "status": why}
            print("COST %-13s %s %-13s states=%-9d %6.2fs"
                  % (tag, sect, why, cnt, secs), flush=True)
        res["records"][tag] = rec
        with open(outp, "w") as fh:
            json.dump(res, fh, indent=1)
    print("->", outp)


if __name__ == "__main__":
    main()
