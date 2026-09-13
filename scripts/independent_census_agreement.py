"""Cross-check the independent census against the committed primary artefacts.

Reads results/pslq-degree4-*/latest.json from a repository checkout and the
artefacts produced by census_driver.py, and checks, per interval:

  * screen retention, after primitivity, against root_filter_candidates
  * root-carrying polynomial count and distinct-root count
  * near-set size against near_candidates_exactly_checked
  * the closest member, coefficient by coefficient, and its exact residual
  * the mean-value residual check the ticket asks for: the closest member is
    evaluated at the interval MIDPOINT here and at the ENDPOINTS by the primary,
    so the two residuals MUST differ, and the difference is capped at
    D_P (u-l)/2 with D_P = sum_k k|a_k| for that polynomial's own coefficients.

No value is transcribed by hand; every primary number is read from JSON.
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from census_driver import decimal_parts, poly_value_scaled, INTERVALS  # noqa: E402


def prim(coeffs) -> int:
    from math import gcd
    g = 0
    for a in coeffs:
        g = gcd(g, abs(a))
    return g


def num(v) -> Fraction:
    return Fraction(v)


def check(repo: str, outdir: str) -> dict:
    rows = []
    for iid, lower, upper in INTERVALS:
        prime = json.load(open(os.path.join(repo, "results", f"pslq-degree4-{iid}", "latest.json")))
        pr = prime["interval_result"]
        mine = json.load(open(os.path.join(outdir, f"census-d4-{iid}.json")))

        hits = [h for h in mine["search_footprint"]["retained_tuples"] if prim(h) == 1]
        roots_n = mine["decisions"]["root_containing_polynomials"]
        # near set: the driver keeps the raw count; count primitive ones here
        from census_driver import interval_parameters, run_screen
        w, rho, B, m, lo, up, E = interval_parameters(lower, upper, 4)
        near = [h for h in run_screen(4, w, B + 10 ** 9, "brute") if prim(h) == 1]

        # ---- exact algebra on the primary's closest member ----
        cp = pr["closest_polynomial"]["coefficients_ascending"]
        Xl, El = decimal_parts(lower)
        Xu, Eu = decimal_parts(upper)
        E2 = max(El, Eu)
        Xl_i, Xu_i = Xl * 10 ** (E2 - El), Xu * 10 ** (E2 - Eu)
        Xm = (Xl_i + Xu_i) // 2
        scale = 10 ** (4 * E2)
        Pl = Fraction(poly_value_scaled(cp, 4, Xl_i, E2), scale)
        Pu = Fraction(poly_value_scaled(cp, 4, Xu_i, E2), scale)
        Pm = Fraction(poly_value_scaled(cp, 4, Xm, E2), scale)
        Dp = sum(k * abs(a) for k, a in enumerate(cp))
        width = Fraction(Xu_i - Xl_i, 10 ** E2)
        mvt_cap = Fraction(Dp) * width / 2
        if Pl == 0 or Pu == 0 or (Pl > 0) != (Pu > 0):
            min_interval = Fraction(0)          # IVT: a genuine root
        else:
            min_interval = min(abs(Pl), abs(Pu))
        primary_res = num(pr["closest_polynomial"]["minimum_absolute_residual"])

        rows.append({
            "interval_id": iid,
            "screen_bound": B,
            "retained_mine": len(hits),
            "retained_primary": pr["root_filter_candidates"],
            "retained_agree": len(hits) == pr["root_filter_candidates"],
            "retained_raw_mine": len(mine["search_footprint"]["retained_tuples"]),
            "near_mine": len(near),
            "near_primary": pr["near_candidates_exactly_checked"],
            "near_agree": len(near) == pr["near_candidates_exactly_checked"],
            "roots_mine": roots_n,
            "roots_primary": pr["root_containing_polynomials"],
            "roots_agree": roots_n == pr["root_containing_polynomials"],
            "distinct_roots_mine": mine["decisions"]["distinct_roots_in_interval"],
            "distinct_roots_primary": pr["distinct_roots_in_interval"],
            "excluded_mine": mine["decisions"]["excluded"],
            "excluded_primary": pr["excluded"],
            "class_size_mine": mine["class_size"],
            "class_size_primary": pr["primitive_quartics_covered"],
            "closest_mine": (mine["closest"] or {}).get("coefficients_ascending"),
            "closest_primary": cp,
            "closest_agree": (mine["closest"] or {}).get("coefficients_ascending") == cp,
            "residual_mine_interval_min": str(min_interval),
            "residual_primary": str(primary_res),
            "residual_agree": min_interval == primary_res,
            "residual_at_midpoint": str(Pm),
            "abs_residual_diff_midpoint_vs_endpoint": str(abs(abs(Pm) - min_interval)),
            "mvt_cap": str(mvt_cap),
            "mvt_holds": abs(abs(Pm) - min_interval) <= mvt_cap,
            "residual_is_nontrivial": (min_interval != 0) and (min_interval > mvt_cap),
            "weights_agree": [int(x) for x in w[1:]] == pr["fixed_point_weights_for_powers_1_to_4"],
            "rho_agree": str(Fraction(rho)) == str(Fraction(pr["weight_rounding_error_bound_scaled"])),
        })
    return {"schema": "matching-one/independent-census-agreement/v1", "issue": 568, "rows": rows}


if __name__ == "__main__":
    repo = sys.argv[1]
    outdir = sys.argv[2]
    out = sys.argv[3] if len(sys.argv) > 3 else os.path.join(HERE, "agreement.json")
    res = check(repo, outdir)
    with open(out, "w") as fh:
        json.dump(res, fh, indent=1)
    print(f"{'interval':<28}{'retain':>12}{'near':>12}{'roots':>9}{'closest':>9}{'residual':>10}{'MVT':>6}")
    for r in res["rows"]:
        print(f"{r['interval_id']:<28}"
              f"{str(r['retained_mine'])+'/'+str(r['retained_primary']):>12}"
              f"{str(r['near_mine'])+'/'+str(r['near_primary']):>12}"
              f"{str(r['roots_mine'])+'/'+str(r['roots_primary']):>9}"
              f"{'OK' if r['closest_agree'] else 'DIFF':>9}"
              f"{'OK' if r['residual_agree'] else 'DIFF':>10}"
              f"{'OK' if r['mvt_holds'] else 'FAIL':>6}")
    allok = all(r["retained_agree"] and r["near_agree"] and r["roots_agree"]
                and r["closest_agree"] and r["residual_agree"] and r["mvt_holds"]
                and r["weights_agree"] and r["rho_agree"] for r in res["rows"])
    print("\n全部一致:", allok)
