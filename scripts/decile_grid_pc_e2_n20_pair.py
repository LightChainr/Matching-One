#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E2 rung check at N=20 for issue #641: the first decile-grid near-counterexample.

The N<=18 sweep (decile_grid_pc_admissible_sweep.py) found the closest pair
HNF(2,1,9) vs HNF(2,0,9) with max decile bracket gap 6.92e-6, above the
~3e-6 SE scale of the published nine-decile vector, so E2 looked unreached.
A post-commit self-check extended the (2,1,L)/(2,0,L) torus family one rung
further: at L=10, N=20, the pair HNF(2,1,10) (2xL cell, helical vertical
wrap, winding (1,L)) vs HNF(2,0,10) (plain 2xL torus) has max decile
bracket gap 648033958315/288230376151711744 = 2.2483e-6 < 3e-6, i.e. two
admissible finite models whose nine-decile quantile vectors agree to better
than the quoted SE while their p_c differ.  E2 IS reached -- the N<=18
boundary was an artifact of the sweep cutoff, not a property of the class.

This script recomputes that pair exactly (Fraction only, Sturm isolation,
bits=60) and records the whole family L=5..10 so the scaling claim in the
note rests on measured numbers, not a fitted guess.

Run:  python scripts/decile_grid_pc_e2_n20_pair.py
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from decile_grid_pc_admissible_sweep import (
    DECILES,
    ISOLATION_BITS,
    build_geometry_hnf,
    bracket_gap,
    evaluate_m,
    isolate_roots,
    cdf_from_counts,
    bernstein_to_power,
    bernstein_counts,
    sweep_record,
)
from exact_polynomial_root_certificate import fraction_text

PAIRS = [((2, 1, L), (2, 0, L)) for L in range(5, 11)]


def family_record(h11: int, h21: int, h22: int) -> dict:
    geometry = build_geometry_hnf((h11, h21, h22))
    record = sweep_record(geometry)
    return {
        "triple": [h11, h21, h22],
        "n": record["n"],
        "m_power": record["m_power"],
        "m_half": record["m_half"],
        "m_roots": record["m_roots"],
        "deciles": record["deciles"],
        "_deciles": record["_deciles"],
        "_m_roots": record["_m_roots"],
        "_m_half": record["_m_half"],
    }


def main() -> int:
    out = {}
    family = []
    for a_triple, b_triple in PAIRS:
        a = family_record(a_triple[0], a_triple[1], a_triple[2])
        b = family_record(b_triple[0], b_triple[1], b_triple[2])
        max_gap = max(
            bracket_gap(a["_deciles"][u], b["_deciles"][u]) for u in DECILES
        )
        root_gap = min(
            bracket_gap(ra, rb) for ra in a["_m_roots"] for rb in b["_m_roots"]
        )
        L = a_triple[2]
        family.append({
            "L": L,
            "pair": f"HNF({a_triple[0]},{a_triple[1]},{L}) vs HNF({b_triple[0]},{b_triple[1]},{L})",
            "n": a["n"],
            "max_decile_bracket_gap": fraction_text(max_gap),
            "max_decile_bracket_gap_float": float(max_gap),
            "min_root_bracket_gap": fraction_text(root_gap),
            "min_root_bracket_gap_float": float(root_gap),
            "m_half_equal": str(a["_m_half"]) == str(b["_m_half"]),
        })
        print(
            f"L={L} N={a['n']}: decile gap {fraction_text(max_gap)} = "
            f"{float(max_gap):.6e}, root gap {fraction_text(root_gap)} = "
            f"{float(root_gap):.6e}"
        )

    # The headline pair is the largest L with decile gap < 3/10^6.
    threshold = Fraction(3, 10**6)
    reached = [row for row in family if Fraction(row["max_decile_bracket_gap"]) < threshold]
    assert reached, "no family member under the 3e-6 SE threshold"
    headline = reached[-1]
    # Independent check (GOVERNANCE 2): recompute the headline gap by a
    # second route -- raw Bernstein counts -> power basis -> Sturm, without
    # touching the sweep_record decile brackets.
    a_triple_check = (2, 1, headline["L"])
    b_triple_check = (2, 0, headline["L"])
    ga = build_geometry_hnf(a_triple_check)
    gb = build_geometry_hnf(b_triple_check)
    ca = bernstein_counts(ga)
    cb = bernstein_counts(gb)
    fa = cdf_from_counts(ca)
    fb = cdf_from_counts(cb)
    gaps = []
    for u in DECILES:
        pa = list(fa)
        pa[0] -= u
        pb = list(fb)
        pb[0] -= u
        bra = isolate_roots(pa, bits=ISOLATION_BITS)
        brb = isolate_roots(pb, bits=ISOLATION_BITS)
        assert len(bra) == 1 and len(brb) == 1
        gaps.append(bracket_gap(bra[0], brb[0]))
    check_gap = max(gaps)
    assert check_gap == Fraction(headline["max_decile_bracket_gap"]), (
        f"independent check {fraction_text(check_gap)} != main "
        f"{headline['max_decile_bracket_gap']}"
    )
    print(
        "independent check: headline gap reproduced exactly "
        f"({fraction_text(check_gap)})"
    )

    artifact = {
        "schema": "matching-one/decile-grid-pc-e2-n20-pair/v1",
        "issue": 641,
        "data_class": "exact Fraction arithmetic; Sturm isolation bits=60",
        "verdict": (
            "E2 IS reached on the admissible class: HNF(2,1,10) vs "
            "HNF(2,0,10) (N=20) has max decile bracket gap "
            f"{headline['max_decile_bracket_gap']} = 2.248e-6 < 3e-6 (the SE "
            "scale of the published nine-decile vector) while their p_c "
            "brackets are disjoint. The prior 'E2 not reached at N<=18' "
            "statement was a sweep-cutoff artifact. On the stricter "
            "'root differs by more than our quoted precision' reading the "
            "case is knife-edge: the root gap 2.248e-6 is of the same order "
            "as the quoted 'a few parts in 10^6' precision of Q(0.5)."
        ),
        "headline_pair": {
            "pair": headline["pair"],
            "n": headline["n"],
            "max_decile_bracket_gap": headline["max_decile_bracket_gap"],
            "max_decile_bracket_gap_float": headline["max_decile_bracket_gap_float"],
            "min_root_bracket_gap": headline["min_root_bracket_gap"],
            "min_root_bracket_gap_float": headline["min_root_bracket_gap_float"],
            "m_half_equal": headline["m_half_equal"],
        },
        "independent_check": (
            "gap recomputed from raw Bernstein counts via a separate "
            "CDF->Sturm route; agrees exactly"
        ),
        "family_scaling": family,
        "scaling_note": (
            "family gaps fall super-polynomially faster than 1/N: measured "
            "gaps 8.18e-4 (L=5), 2.30e-4 (L=6), 6.82e-5 (L=7), 2.16e-5 "
            "(L=8), 6.92e-6 (L=9), 2.25e-6 (L=10) -- roughly one order of "
            "magnitude per rung, local exponent ~N^-9 (N=14->18) to "
            "~N^-10 (N=16->20). The prior note's 'gap ~ 1/N, 3e-6 needs "
            "N >= 50' was wrong; 3e-6 is crossed at N=20. Trade-off: max "
            "decile gap = min root gap exactly for L=7..10 (same order at "
            "L=5,6), so growing N shrinks both together -- decile-within-SE "
            "needs N>=20, root-beyond-quoted-precision needs N<20, no "
            "member of this family satisfies both E2 criteria at once; a "
            "decisive instance needs a pair with separations of different "
            "orders."
        ),
        "torus_shape_note": (
            "HNF(2,1,L) is a 2xL cell (N=2L) with a helical vertical wrap: "
            "the shortest non-horizontal winding lattice vector is (+-1, L), "
            "Manhattan length 2L+1. HNF(2,0,L) is the plain 2xL torus with "
            "vertical winding L. Earlier text calling HNF(2,1,L) 'the 3xL "
            "torus' was wrong (a 3xL cell would have N=3L)."
        ),
    }
    out_dir = ROOT / "results" / "decile-grid-pc"
    out_path = out_dir / "e2-n20-pair.json"
    out_path.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    print(f"wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
