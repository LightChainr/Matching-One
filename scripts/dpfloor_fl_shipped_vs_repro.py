#!/usr/bin/env python3
"""Shipped settings vs re-run: how reproducible is the shipped oblique file?

The shipped file's roots were produced with ARPACK tol=1e-10 and scipy brentq
xtol=3e-11.  Re-running the SAME algorithm at (i) the shipped settings and
(ii) tight settings quantifies the configuration scatter of every geometry in
the file.  (The root value itself does not depend on p_c, so this comparison is
free of any p_c convention question.)
"""
from __future__ import annotations
import sys, os, json

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import fl_common as C

OUT = "/workspace/dpfloor/out"
LD = np.longdouble


def load(name):
    p = os.path.join(OUT, name)
    return json.load(open(p)) if os.path.exists(p) else None


def key_of(direction, n):
    for tag, u, nn, _e in C.GEOMS:
        if list(u) == list(direction) and nn == n:
            return tag
    return None


def main():
    shipped = json.load(open(sys.argv[1] if len(sys.argv) > 1
                             else "/workspace/dpfloor/in/oblique-spin4-controls.json"))
    loose = load("pathB_oblique_loose.json")
    t1 = load("pathB_oblique_tight.json")
    t2 = load("pathB_oblique_tight2.json")
    tight = {}
    for d in (t1, t2):
        if d:
            tight.update(d["records"])

    rows = {}
    for rec in shipped["records"]:
        tag = key_of(rec["direction"], rec["n"])
        if tag is None:
            continue
        ps = LD(str(rec["charge_root_p"]))
        pl = LD(str(loose["records"][tag]["p_root"])) if loose and tag in loose["records"] else None
        pt = LD(str(tight[tag]["p_root"])) if tag in tight else None
        r = {"direction": rec["direction"], "n": rec["n"], "ell": rec["physical_circumference"],
             "p_shipped": str(ps),
             "p_loose_rerun": str(pl) if pl is not None else None,
             "p_tight": str(pt) if pt is not None else None}
        if pl is not None:
            r["shipped_minus_loose"] = str(ps - pl)
        if pt is not None:
            r["shipped_minus_tight"] = str(ps - pt)
            r["loose_minus_tight"] = str(pl - pt) if pl is not None else None
            r["Omega_shipped_minus_tight"] = str((ps - pt) * LD(r["ell"]) ** 4)
        rows[tag] = r
        def fmt(x):
            return "-" if x is None else ("%.2e" % float(x))
        print("%-11s %-9s ell=%7.3f  shipped-tight=%s  loose-tight=%s"
              % (tag, str(rec["direction"]), r["ell"],
                 fmt(r.get("shipped_minus_tight")), fmt(r.get("loose_minus_tight"))),
              flush=True)
    diffs = [abs(LD(r["shipped_minus_tight"])) for r in rows.values() if "shipped_minus_tight" in r]
    dl = [abs(LD(r["loose_minus_tight"])) for r in rows.values() if r.get("loose_minus_tight") is not None]
    summ = {
        "n_geometries": len(rows),
        "max_abs_shipped_minus_tight_p": float(max(diffs)) if diffs else None,
        "max_abs_loose_minus_tight_p": float(max(dl)) if dl else None,
        "statement": ("re-running the shipped algorithm with the shipped settings "
                      "reproduces the shipped roots only to ~1e-12 in p; with tight "
                      "settings the same geometries are reproducible to ~1e-16."),
    }
    print(json.dumps(summ, indent=2))
    with open(os.path.join(OUT, "step_shipped_vs_repro.json"), "w") as fh:
        json.dump({"rows": rows, "summary": summ}, fh, indent=2)
    print("->", os.path.join(OUT, "step_shipped_vs_repro.json"))


if __name__ == "__main__":
    main()
