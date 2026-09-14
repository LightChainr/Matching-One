#!/usr/bin/env python3
"""n325rec step 2 -- cheap telemetry: how big are the oblique safe automata?

Builds PATH B's automaton (verbatim copy of rev769's
scripts/oblique_charge_transfer.py, i.e. dpfloor's B_oblique) for a list of
(direction, n) geometries and records, per sector, the safe-state count, the
row memory and the wall time.  NO root finding, NO eigendecomposition.

Cost gate (spec 2.3B): stop if a sector exceeds --state-cap.
"""
from __future__ import annotations
import json
import math
import sys
import time

sys.path.insert(0, "/workspace/dpfloor/scripts")
import fl_path_oblique as B          # the B_oblique reference implementation


# tag, u, n     (ell = n*|u|)
GEOMS = [
    ("diag_n4",      (1, 1),  4),        # ell 5.657   (dpfloor regression)
    ("diag_n5",      (1, 1),  5),        # ell 7.071
    ("slope21_n3",   (2, 1),  3),        # ell 6.708
    ("slope32_n2",   (3, 2),  2),        # ell 7.211
    ("slope31_n3",   (3, 1),  3),        # ell 9.487
    ("n325_1_18",    (1, 18), 1),        # N=325, ell 18.0278
    ("n325_6_17",    (6, 17), 1),        # N=325, ell 18.0278
    ("n325_2_3_n5",  (2, 3),  5),        # N=325, ell 5*sqrt(13)=18.0278
    ("n25_3_4",      (3, 4),  1),        # N=25,  ell 5.0  (regression point)
]


def main():
    cap = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
    out_path = sys.argv[2] if len(sys.argv) > 2 else \
        "/workspace/n325rec/out/s2_scout.json"
    keep = sys.argv[3].split(",") if len(sys.argv) > 3 else None
    res = {"schema": "n325rec.scout.v1", "path": "B_oblique",
           "md5_note": "fl_path_oblique.py c402cdf6a0d92d317864d793b6eece44",
           "state_cap": cap, "records": {}}
    for tag, u, n in GEOMS:
        if keep and tag not in keep:
            continue
        rec = {"tag": tag, "direction": list(u), "n": n,
               "ell": n * math.hypot(*u)}
        try:
            for sector, matching in (("G4", False), ("G8", True)):
                t0 = time.time()
                a = B.ObliqueSafeTransfer(n, u, matching, state_cap=cap)
                rec["%s_n_safe" % sector] = a.n
                rec["%s_memory" % sector] = a.memory
                rec["%s_edges" % sector] = [list(e) for e in a.edges]
                rec["%s_seconds" % sector] = round(time.time() - t0, 2)
                print("SCOUT %-13s %s n_safe=%-7d memory=%-2d edges=%s %.2fs"
                      % (tag, sector, a.n, a.memory, a.edges,
                         rec["%s_seconds" % sector]), flush=True)
            rec["ok"] = True
        except RuntimeError as exc:
            rec["ok"] = False
            rec["error"] = "%s (state cap %d)" % (exc, cap)
            print("SCOUT %-13s ABORT %s" % (tag, exc), flush=True)
        except Exception as exc:                       # noqa: BLE001
            rec["ok"] = False
            rec["error"] = "%s: %s" % (type(exc).__name__, exc)
            print("SCOUT %-13s ERROR %s" % (tag, exc), flush=True)
        res["records"][tag] = rec
        with open(out_path, "w") as fh:
            json.dump(res, fh, indent=1)
    print("->", out_path)


if __name__ == "__main__":
    main()
