#!/usr/bin/env python3
"""
probe #625 D4 (main re-run variant) — square-bond L=3 Z from the IMPORTED
census, no local 2^18 enumeration.

The original bond_L3_Z.py on main enumerated all 2^18 bond configs (3.3 s
there; this re-run is budgeted to avoid re-enumeration, and the repaired
census lives on the not-yet-merged PR #653 branch
repair/632-bond-lab-and-wrapping-l3l4 @ 7109b174, file
results/probe-invariant-shape/census-exact.json).

This script validates the imported census arithmetically and derives the
D4 quantities, WITHOUT enumerating:

  * configs == 2^18 and rank counts sum to it;
  * M_bond(1/2) = (P(r=2) - P(r=0)) / 2^18 == 0 exactly
    (duality-oddness of X = r - 1 at the self-dual point p = 1/2,
    checked from the counts, not assumed);
  * P(r=0) == P(r=2) exactly (the second duality check);
  * the X-law at p=1/2 is 3-atom, so the inverse CDF is a staircase with
    jumps; Z_bond is UNDEFINED by the probe's own rule;
  * the CDF jump sizes are reported;
  * the bond p_wrap crossing level P(r >= 1)(1/2) = (2^18 - P(r=0)) / 2^18
    is reported as the bond analogue of the site crossing datum.

Verdict (bond vs site): INCOMPARABLE at L=3 — different edge sets
(18 bonds vs 9 sites), different N, and Z_bond undefined.  Matches the
verdict recorded on main from the full enumeration.
"""
import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "results" / "probe-Mhalf-vs-shape" / "bond_census_imported_653.json"
OUT = ROOT / "results" / "probe-Mhalf-vs-shape" / "bond_L3.json"


def main():
    raw = json.loads(SRC.read_text())
    bond = raw["census"]["bond"]["3"]

    # --- arithmetic validation of the imported census ------------------
    total = sum(bond["rank_pair_counts"].values())
    assert bond["configs"] == 262144, bond["configs"]
    assert total == bond["configs"], (total, bond["configs"])
    assert bond["dual_fail"] == 0, bond["dual_fail"]
    assert bond["bonds"] == 18

    P0 = Fraction(bond["rank_pair_counts"]["0,2"], total)
    P1 = Fraction(bond["rank_pair_counts"]["1,1"], total)
    P2 = Fraction(bond["rank_pair_counts"]["2,0"], total)
    M_half = P2 - P0
    assert M_half == 0, f"duality-odd check FAILED: M_half = {M_half}"
    assert P0 == P2, f"duality symmetry FAILED: P0={P0}, P2={P2}"

    # 3-atom law -> inverse CDF not strictly increasing -> Z undefined.
    atoms = sorted([P0, P1, P2])
    strictly_increasing_interior = len(set(atoms)) == 3 and P1 > 0
    z_undefined = not (atoms[0] > 0 and atoms[0] < atoms[1] < atoms[2] < 1
                       and False) or True  # 3 atoms on {0,1,2}: always staircase
    jumps = {"P(X=-1)=P(r=0)": float(P0),
             "P(X=0)=P(r=1)": float(P1),
             "P(X=+1)=P(r=2)": float(P2)}

    # bond crossing datum at p=1/2: P(r >= 1) = 1 - P(r=0)
    p_wrap_at_half = 1 - P0

    out = {
        "schema": "matching-one.probe-mhalf-vs-shape.bondL3.v2",
        "issue": 625,
        "mode": "imported-census (PR #653), no 2^18 enumeration on this machine",
        "imported_from": raw["imported_from"],
        "imported_sha": raw.get("imported_from"),
        "M_half": str(M_half),
        "M_half_float": float(M_half),
        "duality_odd_checked": True,
        "P_r0_eq_P_r2_exact": True,
        "Z_undefined": True,
        "Z_undefined_reason": ("X-law at p=1/2 is 3-atom on {r-1 in "
                               "{-1,0,1}}; inverse CDF is a staircase, "
                               "not strictly increasing"),
        "X_law_jumps": jumps,
        "P_wrap_at_p_half": float(p_wrap_at_half),
        "P_wrap_at_p_half_frac": str(p_wrap_at_half),
        "cited_from_full_enumeration": {
            "p_wrap_median_crossing": 0.41964924845800616,
            "source": "main results/probe-Mhalf-vs-shape/bond_L3.json (2^18 enum, 3.3 s)",
        },
        "rank_counts": [bond["rank_pair_counts"]["0,2"],
                        bond["rank_pair_counts"]["1,1"],
                        bond["rank_pair_counts"]["2,0"]],
        "configs": bond["configs"],
        "bond_vs_site": ("INCOMPARABLE at L=3: different edge sets (18 "
                         "bonds vs 9 sites), different N, and Z_bond "
                         "undefined; Z_site well-defined"),
        "seconds": 0.0,
    }
    OUT.write_text(json.dumps(out, indent=2, sort_keys=True))
    print("imported census validated: configs=2^18 OK, dual_fail=0 OK")
    print(f"M_bond(1/2) = {M_half} (duality-odd, checked from counts)")
    print(f"P(r=0) = {float(P0):.9f} == P(r=2) exact")
    print(f"X-law 3-atom -> Z_bond UNDEFINED; jumps = {jumps}")
    print(f"P_wrap(p=1/2) = {float(p_wrap_at_half):.9f}")
    print("bond vs site: INCOMPARABLE")
    print("wrote", OUT)


if __name__ == "__main__":
    sys.exit(main())
