#!/usr/bin/env python3
"""n325rec step 9 -- FINAL two-path comparison on the reference geometry set.

Reference numbers are dpfloor's OWN recorded PATH B values, read from
  /workspace/dpfloor/out/pathB_oblique_tight.json   (tight, p_c=0.5927460507921)
  /workspace/dpfloor/out/pathB_oblique_tight2.json  (tight, same p_c)
(these were produced by fl_path_oblique.py, md5 c402cdf6a0d92d317864d793b6eece4,
i.e. the same code as B_oblique), so the comparison needs no re-run of B's
expensive brentq/dense-eig pipeline.

For every geometry:
  * build the INDEPENDENT automaton (G4, G8) and record states / memory / time;
  * compare n_safe and row memory against the recorded PATH B values;
  * solve Delta(p)=0 with the independent longdouble machinery and compare with
    the recorded PATH B p_root;
  * for n_safe <= 1200 also rebuild PATH B's matrix and compare the sorted real
    spectrum of R(0.6) elementwise -- a permutation-invariant fingerprint.

The two geometries with no recorded reference (N=25's two orientations) get a
full PATH B run (they are tiny).
"""
from __future__ import annotations

import json
import math
import sys
import time

import numpy as np

sys.path.insert(0, "/workspace/dpfloor/scripts")
sys.path.insert(0, "/workspace/n325rec/scripts")

import fl_path_oblique as B            # noqa: E402
import fl_common as C                  # noqa: E402
import oblique_indep as N              # noqa: E402

PC = 0.5927460507921
PROBE_P = 0.6
SPECTRUM_CAP = 1200
REF_FILES = [
    "/workspace/dpfloor/out/pathB_oblique_tight.json",
    "/workspace/dpfloor/out/pathB_oblique_tight2.json",
]

GEOMS = [
    ("axis_n4",    (1, 0), 4),
    ("axis_n6",    (1, 0), 6),
    ("axis_n8",    (1, 0), 8),
    ("diag_n4",    (1, 1), 4),
    ("diag_n5",    (1, 1), 5),
    ("slope21_n3", (2, 1), 3),
    ("slope21_n4", (2, 1), 4),
    ("slope31_n3", (3, 1), 3),
    ("slope32_n2", (3, 2), 2),
    ("slope52_n2", (5, 2), 2),
    ("n25_3_4",    (3, 4), 1),
    ("n25_0_1_n5", (0, 1), 5),
]


def b_tight(u, n):
    """PATH B tight pipeline (exactly fl_path_oblique.run's tight branch)."""
    g4 = B.ObliqueSafeTransfer(n, u, False)
    g8 = B.ObliqueSafeTransfer(n, u, True)
    mode = "dense" if (g4.n <= 2500 and g8.n <= 2500) else "arpack"

    def eq(p):
        l4, _ = g4.lambda0(p, mode=mode, arpack_tol=1e-15)
        l8, _ = g8.lambda0(1.0 - p, mode=mode, arpack_tol=1e-15)
        return math.log(l4) - math.log(l8)

    root, _w, _nf, fres = C.solve_root(eq, 0.5, 0.8, xtol=1e-16)
    return root, fres, g4, g8


def main():
    outp = sys.argv[1] if len(sys.argv) > 1 else \
        "/workspace/n325rec/out/s9_final.json"
    ref = {}
    for f in REF_FILES:
        d = json.load(open(f))
        for k, v in d["records"].items():
            ref[k] = {"p_root": v["p_root"], "n_safe_G4": v["n_safe_G4"],
                      "n_safe_G8": v["n_safe_G8"],
                      "memory_G4": v["row_memory_G4"],
                      "memory_G8": v["row_memory_G8"],
                      "source": f.split("/")[-1]}
    res = {"schema": "n325rec.final-comparison.v1", "p_c": PC,
           "reference": "dpfloor PATH B recorded values (" +
                        ", ".join(f.split("/")[-1] for f in REF_FILES) + ")",
           "records": {}}
    for tag, u, n in GEOMS:
        t0 = time.time()
        rec = {"tag": tag, "direction": list(u), "n": n,
               "ell": n * math.hypot(*u)}
        t4 = N.IndependentOblique(n, u, False)
        t8 = N.IndependentOblique(n, u, True)
        build_s = time.time() - t0
        root, resid, bracket = N.solve_p_root(t4, t8, PC, dtype=np.longdouble)
        ell = rec["ell"]
        a, b = u
        rec.update({
            "edges_G4": [list(e) for e in t4.edges],
            "edges_G8": [list(e) for e in t8.edges],
            "comp": list(t4.comp),
            "nnz_G4": t4.nnz, "nnz_G8": t8.nnz,
            "memory_G4": t4.memory, "memory_G8": t8.memory,
            "n_safe_G4_indep": t4.n, "n_safe_G8_indep": t8.n,
            "build_seconds": round(build_s, 2),
            "p_root_indep": repr(root),
            "p_root_indep_float": float(root),
            "root_minus_pc": float(root - PC),
            "Delta_at_root": float(resid),
            "bracket_width": float(bracket),
            "Omega": float(-(root - PC) * ell ** 4),
            "cos4": (a ** 4 - 6 * a * a * b * b + b ** 4)
                    / (a * a + b * b) ** 2,
            "seconds": round(time.time() - t0, 2)})
        if tag in ref:
            rec["reference_source"] = ref[tag]["source"]
            rec["p_root_ref"] = ref[tag]["p_root"]
            rec["p_root_diff"] = float(root) - ref[tag]["p_root"]
            rec["n_safe_match"] = (t4.n == ref[tag]["n_safe_G4"]
                                   and t8.n == ref[tag]["n_safe_G8"])
            rec["memory_match"] = (t4.memory == ref[tag]["memory_G4"]
                                   and t8.memory == ref[tag]["memory_G8"])
        else:
            rb, fresb, b4, b8 = b_tight(u, n)
            rec["reference_source"] = "computed here (PATH B tight)"
            rec["p_root_ref"] = rb
            rec["p_root_diff"] = float(root) - rb
            rec["Delta_at_root_B"] = float(fresb)
            rec["n_safe_match"] = (t4.n == b4.n and t8.n == b8.n)
            rec["memory_match"] = (t4.memory == b4.memory
                                   and t8.memory == b8.memory)
        if max(t4.n, t8.n) <= SPECTRUM_CAP and tag in ref:
            b4 = B.ObliqueSafeTransfer(n, u, False)
            b8 = B.ObliqueSafeTransfer(n, u, True)
            devs = {}
            for sect, tb, tn in (("G4", b4, t4), ("G8", b8, t8)):
                sb = np.sort(np.linalg.eigvals(tb.dense(PROBE_P)).real)
                sn = np.sort(np.linalg.eigvals(tn.matrix_f64(PROBE_P)).real)
                devs[sect] = {"max_abs_dev": float(np.max(np.abs(sb - sn))),
                              "scale": float(np.max(np.abs(sb)))}
            rec["spectrum_fingerprint"] = devs
            rec["spectrum_max_abs_dev"] = max(
                devs[s]["max_abs_dev"] for s in devs)
        res["records"][tag] = rec
        print("FIN %-12s n_safe %s  p_root_diff=%+.3e  spec=%s  %.1fs"
              % (tag, "OK" if rec["n_safe_match"] else "MISMATCH",
                 rec["p_root_diff"], rec.get("spectrum_max_abs_dev"),
                 rec["seconds"]), flush=True)
        with open(outp, "w") as fh:
            json.dump(res, fh, indent=1)
    ds = [abs(v["p_root_diff"]) for v in res["records"].values()
          if "p_root_diff" in v]
    obl = [abs(v["p_root_diff"]) for k, v in res["records"].items()
           if "p_root_diff" in v and not k.startswith("axis")]
    res["summary"] = {
        "n_geometries": len(res["records"]),
        "all_state_counts_match": all(v.get("n_safe_match")
                                      for v in res["records"].values()),
        "all_memories_match": all(v.get("memory_match")
                                  for v in res["records"].values()),
        "max_abs_p_root_diff_all": max(ds) if ds else None,
        "max_abs_p_root_diff_oblique": max(obl) if obl else None,
        "F_oblique_from_this_comparison": max(obl) if obl else None,
    }
    print("SUMMARY", json.dumps(res["summary"], indent=1))
    with open(outp, "w") as fh:
        json.dump(res, fh, indent=1)
    print("->", outp)


if __name__ == "__main__":
    main()
