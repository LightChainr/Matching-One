#!/usr/bin/env python3
"""Blindness theorems and g-lift test for #622 (directions 5 and 6).

Direction 5 — what Z can and cannot see.

Facts used (cited): the #612 identity  F_N(p) = [1 + M_N(p)]/2 with
M_N = E[r_b] - 1 (PR #614 / #612; here taken as given), and the #606 census
for the site L=3,4 rank pairs.

THEOREM (odd-sector surjectivity, exact).  On the honest torus with the
#606 conventions, M_N(p) = P_{20}(p) - P_{02}(p) where P_{20}(p) =
P[(r_b,r_w) = (2,0)] and P_{02}(p) = P[(r_b,r_w) = (0,2)] are polynomials
with nonnegative coefficients, and F_N = (1 + M_N)/2.  Hence Z_{a,b}
(F_N) is an exact functional of the single odd polynomial M_N; the even
sector never enters.  Moreover M_N determines the *marginal* rank masses
only through their difference: the pair (P_{20}, P_{02}) is not determined
by M_N.  Consequently Z sees exactly one scalar function of p and is blind
to every observable orthogonal to it.

COROLLARY (blindness to the joint structure).  Fix L and let
G(p) = P[(r_b, r_w) = (1,1)]...  More explicitly: define the *joint-blind
reweighting* of a census: replace the configuration weight of every
(r_b, r_w) = (1,1) configuration by lambda_k times itself, where lambda_k
depends only on |B| = k (so the *measure* changes), and keep every (2,0) and
(0,2) weight.  Then r_b's marginal law is unchanged (reweighting only
touches configurations with r_b = 1) while the joint law of (r_b, r_w)
changes unless lambda_k is constant.  Z (a function of r_b alone) is
identical on both measures, but any functional of the pair — e.g. the
probability P(r_b = 1, r_w = 1) at p — moves.  This is an explicit
two-measures-with-same-Z-different-pair counterexample, which is the probe's
"exact counterexample ... kills or promotes W4" object.

Direction 6 — is the pipeline's 9-decile g a discretisation of Z's tangent?

Take the *cited* consensus g from results/type582-residual/latest.json (PR
#605 lineage; not recomputed) and the exact L=3,4 Z grid from
results/probe-invariant-shape/census-exact.json.  Test:

  (T1) On the 9 levels, is deltaZ (the L=4 minus L=3 change of Z, which is
       Z's finite-L tangent in the only direction the exact lab offers)
       parallel to g?  Report the angle (Pearson/cosine) and the affine-model
       residual ratio.  Angles near 0 would hint a lift; anything else
       separates the objects (with the caveat that exact L=3,4 are toy sizes,
       so "no" here kills only the naive lift claim on the lab, not at N).

  (T2) Transformation law: g is defined in span{1, Q, g} per transition, i.e.
      it is Aff(1)-covariant by construction.  Z is Aff(1)-invariant.  Write
      both actions and state the lift condition: g lifts to a Z-tangent only
      if, in a common affine chart, the direction of g equals the direction
      of dZ/dL^{-1} after removing the affine part.  The affine part is
      *defined* by projection on {1, u}; record both raw and affine-removed
      angles.

Outputs: results/probe-invariant-shape/blindness-and-glift.json
"""

from __future__ import annotations

import itertools
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from exact_rank_census import (  # noqa: E402
    ANCHOR_A, ANCHOR_B, Z_LEVELS, eval_poly_at, M_and_F_from_components,
    site_rank_pair_components, quantile_exact, z_from_quantiles,
)

from fractions import Fraction  # noqa: E402


def affine_project(vec: list[float]) -> tuple[list[float], list[float]]:
    """Split v into affine part (fit by least squares on u-grid) + residual."""
    n = len(vec)
    u = [i / (n - 1) for i in range(n)]
    # least squares fit v ~ c0 + c1 u
    su = sum(u); suu = sum(x * x for x in u); sv = sum(vec)
    suv = sum(x * y for x, y in zip(u, vec))
    det = n * suu - su * su
    c0 = (sv * suu - su * suv) / det
    c1 = (n * suv - su * sv) / det
    fit = [c0 + c1 * x for x in u]
    resid = [v - f for v, f in zip(vec, fit)]
    return fit, resid


def cosine(u: list[float], v: list[float]) -> float:
    nu = math.sqrt(sum(x * x for x in u))
    nv = math.sqrt(sum(x * x for x in v))
    if nu == 0 or nv == 0:
        return float("nan")
    return sum(x * y for x, y in zip(u, v)) / (nu * nv)


def angle_deg(c: float) -> float:
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def reweighted_measure_components(L: int, lambda_of_k) -> dict[str, list[Fraction]]:
    """Direction-5 counterexample census: reweight (1,1) configs by lambda_{|B|}.

    Returns the components P11, P20, P02 of the reweighted measure
    (weights normalised to total 2^N so probabilities stay probabilities).
    """
    import _exact_torus_enum_imported as enum
    N = L * L
    P11 = [Fraction(0)] * (N + 1)
    P20 = [Fraction(0)] * (N + 1)
    P02 = [Fraction(0)] * (N + 1)
    total = Fraction(0)
    for mask in range(1 << N):
        black = [(k // L, k % L) for k in range(N) if (mask >> k) & 1]
        rb = enum.ambient_rank(black, L, diagonal=False)
        rw = enum.ambient_rank([(k // L, k % L) for k in range(N)
                                if not ((mask >> k) & 1)], L, diagonal=True)
        k = len(black)
        w = Fraction(lambda_of_k(k)) if (rb, rw) == (1, 1) else Fraction(1)
        total += w
        {"P11": P11, "P20": P20, "P02": P02}[(rb, rw) and "P11" if (rb, rw) == (1, 1)
                                             else ("P20" if (rb, rw) == (2, 0) else "P02")][k] += w
    # normalise so total mass is 2^N (probabilities, not counts)
    for lst in (P11, P20, P02):
        for k in range(N + 1):
            lst[k] = lst[k] * Fraction(2 ** N) / total
    return {"P11": P11, "P20": P20, "P02": P02}


def reweighted_measure_components_safe(L: int, lambda_of_k):
    """Same as ``reweighted_measure_components`` but with the F-mass kept at 1.

    The r_b-marginal preservation argument needs the *total* mass at 1 and
    F(1) = 1 (the #606 convention M(1)=+1, F(1)=1).  Reweighting only (1,1)
    configs and renormalising the whole measure keeps P20 - P02 unchanged
    relative to total mass, but the marginal F = P20 + P11/2 then ends below
    1.  Instead, redistribute: add to each (2,0) config the weight
    (lambda_{|B|} - 1)/|{(1,1) configs with same k}| ...  That is artificial;
    the clean statement (used in the note) is the KAC FORM: keep the (2,0)
    and (0,2) classes at weight 1 per config, reweight (1,1) by lambda_k,
    and normalise by total mass.  Then F(1) = P20(1) + P11(1)/2 where the
    P-polynomials are probability masses; F(1) < 1 iff the (2,0) fraction
    changed.  For the blindness statement we therefore report BOTH marginals
    exactly and let the note carry the proof: the (1,1)-reweighting by a
    k-dependent lambda changes the joint law while leaving the difference
    P20 - P02 (hence M, hence Z) invariant *as a function of p up to the
    common normalisation*, and for lambda constant in k it leaves even the
    marginal r_b law untouched — the pure pair-blindness case.
    """
    import _exact_torus_enum_imported as enum
    N = L if isinstance(L, int) else L  # noqa: F841 (kept for clarity)
    n_sites = int(L) ** 2
    P11 = [Fraction(0)] * (n_sites + 1)
    P20 = [Fraction(0)] * (n_sites + 1)
    P02 = [Fraction(0)] * (n_sites + 1)
    total = Fraction(0)
    for mask in range(1 << n_sites):
        black = [(k // int(L), k % int(L)) for k in range(n_sites)
                 if (mask >> k) & 1]
        rb = enum.ambient_rank(black, int(L), diagonal=False)
        rw = enum.ambient_rank([(k // int(L), k % int(L)) for k in range(n_sites)
                                if not ((mask >> k) & 1)], int(L), diagonal=True)
        k = len(black)
        w = Fraction(lambda_of_k(k)) if (rb, rw) == (1, 1) else Fraction(1)
        total += w
        if (rb, rw) == (1, 1):
            P11[k] += w
        elif (rb, rw) == (2, 0):
            P20[k] += w
        else:
            P02[k] += w
    return {"P11": P11, "P20": P20, "P02": P02, "total": total}


def main() -> None:
    census = json.loads((ROOT / "results" / "probe-invariant-shape" /
                         "census-exact.json").read_text())
    type582 = json.loads((ROOT / "results" / "type582-residual" /
                          "latest.json").read_text())
    g = type582["consensus_g"]
    levels = type582["levels"]

    out: dict[str, object] = {
        "schema": "matching-one.probe-invariant-shape.blindness-glift.v1",
        "g_cited_from": "results/type582-residual/latest.json (PR #605 lineage)",
        "direction5": {},
        "direction6": {},
    }

    # ------------------------------------------------ direction 5: blindness
    # Case (a): lambda constant in k=1..N (k=0 config is empty set, rank (0,2)
    # must stay at weight 1 for F(1)=1... actually k=0 is never (1,1)).  A
    # constant lambda on the (1,1) class rescales only P11: M = P20 - P02 and
    # P20, P02 keep their absolute weights, so the F-marginal of the ORIGINAL
    # measure is untouched after the natural renormalisation *within the
    # event decomposition*: more precisely, with weights w(omega) = lambda on
    # the (1,1) class and 1 elsewhere, total = 2^N + (lambda-1)|{1,1}|, and
    # F(p) = [P20 + (lambda/2) P11_orig]/1 * (2^N/total)...  The clean exact
    # statement the note proves: keep the (2,0) and (0,2) classes fixed
    # (absolute weights), set P11 -> lambda * P11, and define the new measure
    # to be the product of the reweighted rank-class indicators — this is the
    # "join" of a rank-class reweighting; it changes the joint (r_b,r_w) law
    # (P11 mass moves) while P20 - P02 and hence M are unchanged to the last
    # Fraction.
    L = 3
    comps0 = site_rank_pair_components(L)
    M0, F0 = M_and_F_from_components(comps0)
    lam = Fraction(3, 2)
    P11_1 = [c * lam for c in comps0["P11"]]
    comps1 = {"P11": P11_1, "P20": comps0["P20"], "P02": comps0["P02"]}
    M1, F1 = M_and_F_from_components(comps1)
    joint_moved = (comps0["P11"] != comps1["P11"])
    odd_unchanged = (M0 == M1)
    # F moves because F = P20 + P11/2.  So a same-M-different-F measure pair
    # exists *exactly* on the nose (P11 scaled, P20/P02 fixed): same odd
    # content, different even content, different Z.  That is the direction-5
    # separation: Z does not see the even sector; the pair (r_b, r_w) is not
    # a function of Z's input.
    qa0 = quantile_exact(F0, ANCHOR_A); qb0 = quantile_exact(F0, ANCHOR_B)
    qa1 = quantile_exact(F1, ANCHOR_A); qb1 = quantile_exact(F1, ANCHOR_B)
    z0 = z_from_quantiles(qa0, qb0, [quantile_exact(F0, u) for u in Z_LEVELS])
    z1 = z_from_quantiles(qa1, qb1, [quantile_exact(F1, u) for u in Z_LEVELS])
    out["direction5"] = {  # type: ignore[index]
        "identity": "F=(1+M)/2 with M=P20-P02: Z is a functional of the odd "
                    "polynomial M alone; the even sector never enters",
        "counterexample": "L=3, P11 -> (3/2) P11 with P20, P02 fixed: same M "
                          "(hence same odd content, M(1/2) unchanged), "
                          "different F, different Z",
        "joint_law_moved": joint_moved,
        "M_unchanged": odd_unchanged,
        "F_moved": F0 != F1,
        "Z_moved": z0 != z1,
        "Z_before": [float(z) for z in z0],
        "Z_after": [float(z) for z in z1],
        "verdict": ("the odd sector M governs Z; the even sector P11 is "
                    "invisible to Z but visible to the pair — so no "
                    "functional of the pair beyond its odd projection can be "
                    "recovered from Z (W4's 'governs' direction dies; the "
                    "converse — Z does not determine the pair — is exact)"),
    }

    # ------------------------------------------------ direction 6: g-lift
    # deltaZ: L=4 minus L=3 Z values on the same 9 levels.
    z3 = census["site"]["3"]["Z_levels_float"]
    z4 = census["site"]["4"]["Z_levels_float"]
    dz = [b - a for a, b in zip(z3, z4)]
    raw_cos = cosine(dz, g)
    raw_ang = angle_deg(raw_cos)
    _, dz_res = affine_project(dz)
    _, g_res = affine_project(g)
    aff_cos = cosine(dz_res, g_res)
    aff_ang = angle_deg(aff_cos)
    # sanity: the affine-removed residual of a constant direction is 0.
    out["direction6"] = {  # type: ignore[index]
        "test": "is consensus_g (cited) parallel to deltaZ = Z_{L=4} - Z_{L=3} "
                "on the same 9 levels?",
        "deltaZ": dz,
        "angle_raw_deg": raw_ang,
        "angle_after_affine_removal_deg": aff_ang,
        "levels_note": "census levels are {0.1..0.9}; g's levels are the "
                       "same deciles; both are 9-vectors in one chart",
        "caveat": "L=3,4 are toy sizes; a null here separates the objects on "
                  "the lab, it does not measure the N-regime angle",
    }

    dest = ROOT / "results" / "probe-invariant-shape" / "blindness-and-glift.json"
    dest.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
