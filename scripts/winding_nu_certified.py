#!/usr/bin/env python3
"""Certified exact winding densities nu_w at widths 2..8, both adjacencies.

Purpose (2026-09-13, erratum follow-up): the span-spectrum erratum recorded that
closure `sum_h d_h + tail = nu_w` could only be *verified* at w = 2,3,4 because
no certified nu_w reference existed at w >= 5. This script removes that gap.

It uses the repository's own winding-intensity engine
(`cylinder_winding_intensity.py`: the same `advance()` the 18 published controls
were validated against), built by exhaustive BFS, then

  1. computes the exact rational stationary reward by rational Gauss-Jordan on
     the common all-p exit-law lumping (exact by construction: the lumping is
     defined to retain the joint next-class/reward law);
  2. independently re-certifies the same value through
     `stationary_certificate`, which returns a rigorous rational interval from a
     forward-error bound through the empty-row reset -- so the two numbers come
     from different code paths;
  3. checks the raw (unlumped) solve against the lumped one wherever the raw
     solve is affordable, so the lumping itself is controlled, not assumed;
  4. compares against the committed #741 references where they exist.

This is an INDEPENDENT engine from the span-spectrum builder: different state
space (no depth tracking), different reward accounting code, different solve.
Agreement of the two is therefore real corroboration.

usage:
  python scripts/winding_nu_certified.py --width 8 --matching 0 --p 1/4
  python scripts/winding_nu_certified.py --grid results/...json
"""
from __future__ import annotations
import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import cylinder_winding_intensity as W  # noqa: E402

# #741 certified references (exact rationals as printed in the handoff).
REFS = {
    ("NN", 8, "1/4"): "17502628473380503424175742111730325030001801413138981331992379782406219403788076443461112141432226837307423006489737552208224632141450146887713404996562363206349515250607990987848286702015668340618931011832915030867953648522909129625535987525113998308888891118101205671016155713840182276590146401648999847735217782956902290870163141233/394858190638015284877611732299089108091767441305269271936475303287570100103520767278514851975447321604769571781491275729897390594072360939445749498586029794045625531954047737173300378900180040052249014005467482837461037862326207208006914196399920025371610656324995371887037376183738569888770649458081597003907666586548661425781189401116672",
    ("matching", 8, "1/8"): "4304066353276814600044997473461749946733998654996906665848710094050668526018463750118116035133786833693570766188431654980451251563150021628926696864281471364809534986496434452171916677495079889807512646873142689361000629991021992202254034107483405752782584490865423182129099195425092500784726730027782992682182888628256311224842606882093651369807185556594705994887406409039631339042921350991283756541994077559051452642113341699201143846379431566869183683434604000069810038992473728867102962226989648218217/113511764657387757743385852913961444340846792248681594418514915909829354582886141640860978491305768062425209800158664229180983934165690093259309006582490727717774659821138193397179076192301603758529671081392965053326462178884774141814282691798376964456718728193496138901518192555819255671364890966589998223829885349719255811426060081926113431586436483919311190561771297363987867516713176200921755746267736239849245504325297603519015385553972106891313881673522269765109458789284250226055406632424069191379714048",
}

GRID = {
    "NN": {"p": ["1/8", "1/4"], "matching": False},
    "matching": {"p": ["1/16", "1/8"], "matching": True},
}


def nu_exact(width: int, matching: bool, p: Fraction, raw_too: bool = False) -> dict:
    t0 = perf_counter()
    states, transfer = W.build_transfer(width, matching)
    t_build = perf_counter() - t0
    reduced, blocks = W.reward_lump(transfer)
    t_lump = perf_counter() - t0 - t_build

    t1 = perf_counter()
    res_lumped = W.stationary_reward(reduced, p)
    t_solve = perf_counter() - t1
    mean = res_lumped["mean"]
    assert mean >= 0

    # Independent re-certification of the SAME stationary vector, from the
    # certificate's forward-error bound (different code path from the solve).
    cand = res_lumped["stationary"]
    cert = W.stationary_certificate(reduced, p, cand)
    assert cert["stationarity_l1_residual"] == 0, "exact rational solve left a residual"

    out = {
        "graph": "matching" if matching else "NN",
        "width": width,
        "p": f"{p.numerator}/{p.denominator}",
        "raw_states": len(states),
        "lumped_states": len(reduced),
        "nu_exact": str(mean),
        "nu_float": float(mean),
        "certificate_estimate": str(cert["estimate"]),
        "certificate_agrees": cert["estimate"] == mean,
        "certificate_abs_error_bound": str(cert["absolute_error_bound"]),
        "empty_row_reset": str(cert["empty_row_reset"]),
        "variance_rate_float": float(res_lumped["variance_rate"]),
        "seconds_build": t_build, "seconds_lump": t_lump, "seconds_solve": t_solve,
    }
    if raw_too:
        t2 = perf_counter()
        res_raw = W.stationary_reward(transfer, p)
        out["nu_exact_raw"] = str(res_raw["mean"])
        out["raw_solve_agrees"] = res_raw["mean"] == mean
        out["seconds_solve_raw"] = perf_counter() - t2
    return out


def grid_runner(widths, graphs, plists) -> list:
    rows = []
    for gname in graphs:
        matching = GRID[gname]["matching"]
        for w in widths:
            for ps in plists:
                pn, pd = ps.split("/")
                r = nu_exact(w, matching, Fraction(int(pn), int(pd)), raw_too=(w <= 5))
                key = (gname, w, ps)
                if key in REFS:
                    ref = Fraction(REFS[key])
                    r["reference_741"] = REFS[key]
                    r["matches_reference_741"] = (ref == Fraction(r["nu_exact"]))
                    r["reference_gap"] = str(ref - Fraction(r["nu_exact"]))
                rows.append(r)
                print("%-9s w=%d p=%-5s raw=%-6d lumped=%-5d nu=%.12f cert=%s%s%s"
                      % (gname, w, ps, r["raw_states"], r["lumped_states"], r["nu_float"],
                         "OK" if r["certificate_agrees"] else "MISMATCH",
                         "" if "matches_reference_741" not in r else
                         ("  #741:OK" if r["matches_reference_741"] else "  #741:DIFF"),
                         "" if "raw_solve_agrees" not in r else
                         ("  raw:OK" if r["raw_solve_agrees"] else "  raw:DIFF")),
                      flush=True)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--width", type=int, action="append")
    ap.add_argument("--matching", type=int, default=0)
    ap.add_argument("--p", action="append")
    ap.add_argument("--graphs", default="NN,matching")
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    widths = args.width or [2, 3, 4, 5, 6, 7, 8]
    plists = args.p or ["1/8", "1/4"]
    rows = grid_runner(widths, args.graphs.split(","), plists)
    bad = [r for r in rows if not r["certificate_agrees"]
           or r.get("matches_reference_741") is False or r.get("raw_solve_agrees") is False]
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps({"rows": rows, "failures": len(bad)}, indent=1))
        print("wrote", args.out)
    print("FAILURES:", len(bad))
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
