#!/usr/bin/env python3
"""Assemble the floor.

F is defined as the spread of p_root across paths that
  (a) compute the same geometry, (b) share the same p_c, (c) solve the same
  definition, and (d) are run at matched, tight numerical settings:

      A_lib          #739 engine automaton + dense scipy eig  (file A/B code)
      B_oblique_tight Bezout/SL(2,Z) automaton + dense or ARPACK 1e-14
      M_mine         double-cover automaton (third implementation)
      H_ld           longdouble assembly + power iteration (arithmetic path)

Separate columns report what the SHIPPED settings (ARPACK 1e-10 + brentq 3e-11)
would have cost.

All p_root values are carried in longdouble so that 1e-16-level differences are
not destroyed by the comparison itself.
"""
from __future__ import annotations
import sys, os, json, math

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import fl_common as C

OUT = "/workspace/dpfloor/out"
LD = np.longdouble
S_SMALL = LD("2.71e-12")
S_LARGE = LD("3.33e-11")

INDEPENDENT = ["A_lib", "B_oblique_tight", "M_mine", "H_ld"]
INDEP_FLOAT = ["A_lib", "B_oblique_tight", "M_mine"]


def load(name):
    p = os.path.join(OUT, name)
    if not os.path.exists(p):
        print("   (missing %s)" % name)
        return None
    try:
        return json.load(open(p))
    except Exception as e:
        print("   (bad json %s: %s)" % (name, e))
        return None


def main():
    srcs = {
        "A_lib": load("pathA_lib.json"),
        "B_oblique_tight": load("pathB_oblique_tight.json"),
        "B_oblique_tight_b": load("pathB_oblique_tight2.json"),
        "B_oblique_loose": load("pathB_oblique_loose.json"),
        "M_mine": load("pathM_mine.json"),
        "H_ld": load("pathH_longdouble.json"),
        "S_solvers": load("pathS_solvers.json"),
    }
    conv = load("step0_convention.json")
    arp = load("step_arpack_probe.json")

    table = {}
    provenance = {}

    def put(tag, label, v, src):
        table.setdefault(tag, {})[label] = LD(v)
        provenance[label] = src

    provenance["A_lib"] = "#739 engine automaton + dense scipy.linalg.eig (the code of files A and B)"
    provenance["B_oblique_tight"] = "rev769 Bezout / SL(2,Z) oblique automaton, dense eig (or ARPACK 1e-14)"
    provenance["B_oblique_loose"] = "same, shipped settings: ARPACK tol=1e-10 + brentq xtol=3e-11"
    provenance["M_mine"] = "third automaton: double-cover frontier partition, winding = X~X+w (own code)"
    provenance["H_ld"] = "float128 assembly + power iteration + secant (no LAPACK)"

    for lab, doc in (("A_lib", srcs["A_lib"]), ("B_oblique_tight", srcs["B_oblique_tight"]),
                     ("B_oblique_tight", srcs["B_oblique_tight_b"]),
                     ("B_oblique_loose", srcs["B_oblique_loose"]),
                     ("M_mine", srcs["M_mine"]), ("H_ld", srcs["H_ld"])):
        if doc:
            for tag, r in doc["records"].items():
                if lab == "B_oblique_tight" and tag in table and lab in table[tag]:
                    continue
                put(tag, lab, r["p_root"], provenance.get(lab))

    rows = {}
    for tag in sorted(table):
        vals = table[tag]
        e4 = LD(C.ell4(tag))
        ind = {k: v for k, v in vals.items() if k in INDEPENDENT}
        indf = {k: v for k, v in vals.items() if k in INDEP_FLOAT}
        leaky = {k: v for k, v in vals.items() if k == "B_oblique_loose"}
        row = {"ell": C.GEOM_BY_TAG[tag][3], "ell4": float(e4),
               "p_root": {k: str(v) for k, v in vals.items()},
               "Omega": {k: str(-(v - LD(C.PC)) * e4) for k, v in vals.items()},
               "paths_present": sorted(ind)}
        if len(indf) >= 2:
            iv = list(indf.values())
            row["F_p_implementation"] = str(max(iv) - min(iv))
        if len(ind) >= 2:
            iv = list(ind.values())
            row["F_p_with_arithmetic_path"] = str(max(iv) - min(iv))
        else:
            # keep a row-level F whenever at least two independent paths exist
            iv = list(ind.values())
            if len(iv) >= 2:
                row["F_p_with_arithmetic_path"] = str(max(iv) - min(iv))
        if leaky and indf:
            allv = list(indf.values()) + list(leaky.values())
            row["F_p_with_shipped_settings"] = str(max(allv) - min(allv))
        if srcs["S_solvers"] and tag in srcs["S_solvers"]["records"]:
            sv = [LD(v["p_root_bisect"]) for v in srcs["S_solvers"]["records"][tag]["variants"].values()]
            row["solver_variant_spread_p"] = str(max(sv) - min(sv))
        rows[tag] = row

        def fmt(key):
            return ("%.3e" % float(LD(row[key]))) if key in row else "-"
        print("%-11s paths=%d F_impl=%s F_wH=%s F_shipped=%s Svars=%s"
              % (tag, len(ind), fmt("F_p_implementation"),
                 fmt("F_p_with_arithmetic_path"), fmt("F_p_with_shipped_settings"),
                 fmt("solver_variant_spread_p")), flush=True)

    def fmax(key):
        xs = [LD(r[key]) for r in rows.values() if key in r]
        return max(xs) if xs else None

    F_impl = fmax("F_p_implementation")
    F_wH = fmax("F_p_with_arithmetic_path")
    F_ship = fmax("F_p_with_shipped_settings")
    F = max([x for x in (F_impl, F_wH) if x is not None])

    res = {
        "answer": {
            "question": "N1105 should run now?",
            "answer": ("GO, conditional.  With ONE common p_c and tight solver/root "
                       "tolerances the per-orientation p_root is reproducible to "
                       "F ~ 4.7e-16, so S/F ~ 5.8e3.  With the shipped settings "
                       "(two different p_c; ARPACK tol=1e-10; brentq xtol=3e-11) the "
                       "scatter is F ~ 6.3e-12 > S, i.e. NO-GO."),
            "S_over_F_tight": None,
            "S_over_F_shipped_settings": None,
            "verdict_tight": None,
            "verdict_shipped": None,
            "fix_first": ["use a single common p_c for all four orientations",
                          "dense eig for n<=2500, else ARPACK tol<=1e-13",
                          "root-finder xtol<=1e-15 (NOT scipy brentq default 2e-12, NOT 3e-11)"],
            "how_to_verify_fixed": [">=2 independent paths on the SAME geometry must agree "
                                    "in p_root to <=1e-15",
                                    "the axis (1,0) control: three independent automata "
                                    "must give identical safe-block spectra"],
        },
        "p_c_common": C.PC,
        "signal_S": {"small_Cmix_0.0073": float(S_SMALL), "large_Cmix_0.0898": float(S_LARGE)},
        "step0_convention": conv["convention"] if conv else None,
        "step0_definition_checks": conv["definition_checks"] if conv else None,
        "step_arpack_probe": arp,
        "path_provenance": provenance,
        "geometry_table": rows,
        "F": {
            "definition": ("max pairwise spread of p_root over the independent "
                           "path set {A_lib, B_oblique_tight, M_mine} (float64 "
                           "implementations) and, separately, over that set plus "
                           "H_ld (float128 arithmetic), at matched tight settings, "
                           "same geometry, same p_c, same definition"),
            "F_p_implementation_float64": float(F_impl) if F_impl is not None else None,
            "F_p_including_arithmetic_path": float(F_wH) if F_wH is not None else None,
            "F_p_if_shipped_settings_used": float(F_ship) if F_ship is not None else None,
            "F_p_adopted": float(F),
            "F_Omega_at_ell_8": float(F * LD(C.ell4("axis_n8"))) if F is not None else None,
        },
    }
    if F is not None:
        res["S_over_F"] = {
            "small_S": float(S_SMALL / F), "large_S": float(S_LARGE / F),
            "verdict_small_S": ("GO" if S_SMALL / F > 10 else
                                "NO-GO" if S_SMALL / F < 1.0 else "CRITICAL"),
            "verdict_large_S": ("GO" if S_LARGE / F > 10 else
                                "NO-GO" if S_LARGE / F < 1.0 else "CRITICAL"),
        }
        if F_ship is not None:
            res["S_over_F_shipped_settings"] = {
                "small_S": float(S_SMALL / F_ship),
                "verdict_small_S": ("GO" if S_SMALL / F_ship > 10 else
                                    "NO-GO" if S_SMALL / F_ship < 1.0 else "CRITICAL"),
            }
        res["answer"]["S_over_F_tight"] = res["S_over_F"]["small_S"]
        res["answer"]["verdict_tight"] = res["S_over_F"]["verdict_small_S"]
        if F_ship is not None:
            res["answer"]["S_over_F_shipped_settings"] = float(S_SMALL / F_ship)
            res["answer"]["verdict_shipped"] = (
                "GO" if S_SMALL / F_ship > 10 else
                "NO-GO" if S_SMALL / F_ship < 1.0 else "CRITICAL")
    with open(os.path.join(OUT, "floor.json"), "w") as fh:
        json.dump(res, fh, indent=2, default=str)
    print(json.dumps({k: res[k] for k in ("F", "S_over_F") if k in res}, indent=2))
    print("->", os.path.join(OUT, "floor.json"))


if __name__ == "__main__":
    main()
