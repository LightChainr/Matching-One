#!/usr/bin/env python3
"""n325rec -- assemble rec.json (all raw numbers in one place)."""
from __future__ import annotations

import hashlib
import json
import os
import sys

OUT = "/workspace/n325rec/out"
FILES = ["s1_orient", "s2_scout", "s3_cert", "s4_family", "s4_targets",
         "s5_rehearsal", "s6_pair", "s7_basis", "s8_lattice", "s9_final"]
SCRIPTS = ["s1_orient.py", "s2_scout.py", "s3_cert.py", "s4_cost.py",
           "s5_rehearsal.py", "s6_pair.py", "s7_basis.py", "s8_lattice.py",
           "s9_final.py", "oblique_indep.py"]


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        h.update(fh.read())
    return h.hexdigest()


def main():
    rec = {
        "schema": "n325rec.rec.v1",
        "date": "2026-09-14",
        "machine": "DevEnvC_NePnUn (account2, 16 vCPU / 32 GiB, aarch64, py3.9.9)",
        "p_c": 0.5927460507921,
        "floor_F_imported_from_dpfloor": 4.68e-16,
        "repository_writes": "none (GitHub read-only)",
        "answers": {
            "oblique_independent_review": (
                "PASS. Second, independently written oblique twist path "
                "(different state encoding, different winding test, different "
                "eigensolver/arithmetic) reproduces B_oblique on every tested "
                "geometry: identical safe-state counts, identical row memory, "
                "sorted R(p) spectra deviating by EXACTLY 0.0, and "
                "|p_root difference| <= 4.33e-15 (<=6.66e-16 against the "
                "dense-LAPACK subset of PATH B; the rest is PATH B's ARPACK)."),
            "n325_rehearsal": (
                "NOT FEASIBLE, and not for accuracy reasons. All three N=325 "
                "orientations blow past the 200k-state gate; the cheapest one "
                "(u=(1,18), n=1) built 2,500,000 states without terminating in "
                "118 s (G4) / 126 s (G8). ell=sqrt(325) is moreover the "
                "SMALLEST circumference admitting >=3 same-ell D4 classes, so "
                "there is no cheaper modulus. Extrapolating the measured "
                "x2.85-per-row state growth to N=1105 (memory_G8 = |a|+|b| = "
                "37,41,43,47) gives 1e8..1e21 states: the F-based GO is not "
                "actionable."),
        },
        "reference_md5": {
            "/workspace/dpfloor/scripts/fl_path_oblique.py":
                "c402cdf6a0d92d317864d793b6eece44",
            "rev769-repo/scripts/oblique_charge_transfer.py (local copy)":
                "d5c0e5a4896bc6a244be8460e6929dd3",
            "rev769-repo/scripts/diagonal_charge_transfer.py (local copy)":
                "0ea1018e4d43a010035e142ce85864b1",
            "rev769-repo/scripts/oblique_winding_necklace.py (local copy)":
                "ba82c6076017f7aa27297778153caefe",
            "rev769-repo/scripts/oblique_winding_corridor.py (local copy)":
                "a94d2bcbf6bb5db1499d00beaee55644",
        },
        "scripts_md5": {s: md5(os.path.join("/workspace/n325rec/scripts", s))
                        for s in SCRIPTS
                        if os.path.exists(os.path.join("/workspace/n325rec/scripts", s))},
    }
    for name in FILES:
        p = os.path.join(OUT, name + ".json")
        if os.path.exists(p):
            rec[name] = json.load(open(p))
        else:
            rec[name] = None
    p = "/workspace/n325rec/out/rec.json"
    with open(p, "w") as fh:
        json.dump(rec, fh, indent=1)
    print("wrote", p, os.path.getsize(p), "bytes")
    for s, m in rec["scripts_md5"].items():
        print("  %-20s %s" % (s, m))


if __name__ == "__main__":
    main()
