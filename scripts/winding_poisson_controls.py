#!/usr/bin/env python3
"""Finite controls for one-anchor-per-winding-component Poisson reduction.

Standard library only. Physical lifted edges are retained, including parallel
edges at width two. No simulated samples or fitted percolation constants.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from itertools import product
import json
import math
from pathlib import Path
from typing import Iterable


def steps(matching: bool) -> tuple[tuple[int, int], ...]:
    return ((1, 0), (-1, 0), (0, 1), (0, -1)) + (
        ((1, 1), (1, -1), (-1, 1), (-1, -1)) if matching else ()
    )


def adjacency(w: int, m: int, matching: bool = False, periodic_y: bool = True):
    if w < 2 or m < 2:
        raise ValueError('width and height must both be at least two')
    out = [[] for _ in range(w * m)]
    for y in range(m):
        for x in range(w):
            for dx, dy in steps(matching):
                yy = y + dy
                if not periodic_y and not 0 <= yy < m:
                    continue
                out[y * w + x].append((((yy % m) * w + (x + dx) % w), dx, dy))
    return out


def components(mask: int, w: int, m: int, adj):
    """Independent graph-potential traversal; return vertices and lift cycles."""
    unseen = mask
    answer = []
    while unseen:
        bit = unseen & -unseen
        root = bit.bit_length() - 1
        unseen ^= bit
        lift = {root: (0, 0)}
        todo = [root]
        verts = []
        cycles = set()
        while todo:
            v = todo.pop()
            verts.append(v)
            vx, vy = lift[v]
            for z, dx, dy in adj[v]:
                if not mask & (1 << z):
                    continue
                q = (vx + dx, vy + dy)
                if z not in lift:
                    lift[z] = q
                    unseen &= ~(1 << z)
                    todo.append(z)
                else:
                    gain = (q[0] - lift[z][0], q[1] - lift[z][1])
                    if gain != (0, 0):
                        cycles.add(gain)
        answer.append((tuple(verts), tuple(sorted(cycles))))
    return answer


def ambient_rank(comps) -> int:
    vecs = [z for _, cyc in comps for z in cyc]
    if not vecs:
        return 0
    ax, ay = vecs[0]
    return 2 if any(ax * by != ay * bx for bx, by in vecs[1:]) else 1


def global_anchors(comps, w: int, m: int, height: int) -> tuple[int, ...]:
    """One lexicographic anchor for each full component using <=height rows."""
    found = []
    for verts, cyc in comps:
        if not any(dx for dx, _ in cyc):
            continue
        rows = {v // w for v in verts}
        if len(rows) > height or len(rows) == m:
            continue
        starts = [y for y in rows if (y - 1) % m not in rows]
        if len(starts) != 1:
            raise AssertionError('connected component has non-contiguous row projection')
        bottom = starts[0]
        if any((y - bottom) % m >= height for y in rows):
            continue
        if any(dy for _, dy in cyc):
            raise AssertionError('short component has vertical winding')
        x = min(v % w for v in verts if v // w == bottom)
        found.append(bottom * w + x)
    return tuple(sorted(found))


def local_anchor_bits(mask: int, w: int, height: int, local_adj) -> tuple[int, ...]:
    """Window rows 0..height+1, bottom row 1, two guard rows.

    Does not call global_anchors and does not use vertical periodic edges.
    """
    ans = []
    for verts, cyc in components(mask, w, height + 2, local_adj):
        ys = {v // w for v in verts}
        if 0 in ys or height + 1 in ys or 1 not in ys:
            continue
        if not any(dx for dx, _ in cyc):
            continue
        if any(dy for _, dy in cyc):
            raise AssertionError('open vertical window has vertical winding')
        ans.append(min(v % w for v in verts if v // w == 1))
    return tuple(sorted(ans))


def local_all_anchors(mask: int, w: int, m: int, height: int, local_adj):
    if not 1 <= height < m / 2 or height + 2 >= m:
        raise ValueError('need 1 <= cutoff < m/2 and cutoff+2 < m')
    rowmask = (1 << w) - 1
    ans = []
    for bottom in range(m):
        win = 0
        for k in range(height + 2):
            row = (bottom - 1 + k) % m
            win |= ((mask >> (w * row)) & rowmask) << (w * k)
        ans.extend(bottom * w + x for x in local_anchor_bits(win, w, height, local_adj))
    return tuple(sorted(ans))


def support_rows(anchor: int, w: int, m: int, height: int) -> frozenset[int]:
    y = anchor // w
    return frozenset((y - 1 + k) % m for k in range(height + 2))


def dependency(w: int, m: int, height: int):
    rr = [support_rows(i, w, m, height) for i in range(w * m)]
    return [tuple(j for j in range(w * m) if rr[i] & rr[j]) for i in range(w * m)]


def bernstein_value(coef: Iterable[int], p: Fraction) -> Fraction:
    v = list(coef)
    n = len(v) - 1
    return sum((Fraction(c) * p ** k * (1 - p) ** (n - k)
                for k, c in enumerate(v)), Fraction())


def local_intensity_coefficients(w: int, height: int, matching: bool):
    n = w * (height + 2)
    adj = adjacency(w, height + 2, matching, False)
    cs = [0] * (n + 1)
    for mask in range(1 << n):
        cs[mask.bit_count()] += len(local_anchor_bits(mask, w, height, adj))
    return cs


def fstr(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f'{x.numerator}/{x.denominator}'


def univariate_tv(law: dict[int, Fraction], lam: Fraction) -> float:
    z = float(lam)
    overlap = sum(min(float(prob), math.exp(-z) * z ** k / math.factorial(k))
                  for k, prob in law.items())
    return max(0.0, 1.0 - overlap)


def enumerate_one(w: int, m: int, height: int, matching: bool, ps=(Fraction(1, 10), Fraction(1, 4))):
    n = w * m
    adj = adjacency(w, m, matching)
    ladj = adjacency(w, height + 2, matching, False)
    dep = dependency(w, m, height)
    singles = [[0] * (n + 1) for _ in range(n)]
    pairs = {(i, j): [0] * (n + 1) for i in range(n) for j in dep[i] if j != i}
    laws = Counter()
    absent = [0] * (n + 1)
    mismatch = [0] * (n + 1)
    large = [0] * (n + 1)
    maxz = 0
    # Cache the independent local-window classification, not global labels.
    lcache = {b: local_anchor_bits(b, w, height, ladj)
              for b in range(1 << (w * (height + 2)))}
    rowmask = (1 << w) - 1
    for mask in range(1 << n):
        cc = components(mask, w, m, adj)
        aa = global_anchors(cc, w, m, height)
        loc = []
        for y in range(m):
            win = 0
            for k in range(height + 2):
                win |= (((mask >> (w * ((y - 1 + k) % m))) & rowmask)
                        << (w * k))
            loc.extend(y * w + x for x in lcache[win])
        if aa != tuple(sorted(loc)):
            raise AssertionError(('local/global mismatch', w, m, height, matching, mask))
        k = mask.bit_count()
        r = ambient_rank(cc)
        z = len(aa)
        maxz = max(maxz, z)
        laws[z, k] += 1
        absent[k] += (r == 0)
        is_bad = any(len({v // w for v in verts}) > height for verts, _ in cc)
        large[k] += is_bad
        mismatch[k] += ((r == 0) != (z == 0))
        if ((r == 0) != (z == 0)) and not is_bad:
            raise AssertionError('void mismatch not covered by localization failure')
        for i in aa:
            singles[i][k] += 1
            for j in aa:
                if (i, j) in pairs:
                    pairs[i, j][k] += 1
    local_cs = local_intensity_coefficients(w, height, matching)
    cases = []
    for p in ps:
        pp = [bernstein_value(c, p) for c in singles]
        lam = sum(pp)
        nu = bernstein_value(local_cs, p)
        assert lam == m * nu
        b1 = sum((pp[i] * pp[j] for i in range(n) for j in dep[i]), Fraction())
        b2 = sum((bernstein_value(c, p) for c in pairs.values()), Fraction())
        law = {z: sum((Fraction(laws[z, k]) * p ** k * (1-p) ** (n-k)
                       for k in range(n+1)), Fraction()) for z in range(maxz+1)}
        assert sum(law.values()) == 1
        assert sum(z * prob for z, prob in law.items()) == lam
        actual_void = bernstein_value(absent, p)
        miss = bernstein_value(mismatch, p)
        bad = bernstein_value(large, p)
        assert abs(actual_void - law[0]) <= miss <= bad
        tv = univariate_tv(law, lam)
        # AGG process bound in sup-event total-variation convention.
        bound = min(Fraction(1), 2 * (b1 + b2))
        if tv > float(bound) + 1e-13:
            raise AssertionError(('Poisson bound failed', tv, bound))
        cases.append({'p': fstr(p), 'lambda': fstr(lam), 'nu_per_row': fstr(nu),
                      'b1': fstr(b1), 'b2': fstr(b2),
                      'anchored_count_law': {str(z): fstr(prob) for z, prob in law.items()},
                      'rank_zero_probability': fstr(actual_void),
                      'void_mismatch_probability': fstr(miss),
                      'localization_failure_probability': fstr(bad),
                      'poisson_tv_float_diagnostic': tv,
                      'agg_process_tv_bound': fstr(bound),
                      'rank_void_poisson_error_float_diagnostic': abs(float(actual_void)-math.exp(-float(lam)))})
    return {'width': w, 'length': m, 'cutoff_rows': height, 'matching': matching,
            'configurations': 1 << n, 'local_configurations': 1 << (w*(height+2)),
            'local_global_anchor_failures': 0, 'max_anchored_count': maxz,
            'local_intensity_bernstein_counts': local_cs,
            'cases': cases}


def two_colour_control(w=2, m=4, height=1):
    """Exact categorical coupling: low=1/4, middle=1/2, high=1/4."""
    n = w * m
    adjs = [adjacency(w, m, False), adjacency(w, m, True)]
    deps = dependency(w, m, height)
    dset = [set(d) for d in deps]
    one = [0] * (2*n)
    pair = Counter()
    law = Counter()
    void = Counter()
    total = 4 ** n
    lowq = highq = Fraction(1, 4)
    tables = []
    for adj in adjs:
        table = []
        for mask in range(1 << n):
            cc = components(mask, w, m, adj)
            table.append((global_anchors(cc, w, m, height), ambient_rank(cc)))
        tables.append(table)
    for labels in product(range(3), repeat=n):
        low = high = 0
        mid = 0
        for i, label in enumerate(labels):
            if label == 0:
                low |= 1 << i
            elif label == 2:
                high |= 1 << i
            else:
                mid += 1
        weight = 1 << mid
        al, rl = tables[0][low]
        ah, rh = tables[1][high]
        aa = list(al) + [n+i for i in ah]
        law[len(al), len(ah)] += weight
        void[rl == 0, rh == 0] += weight
        for i in aa:
            one[i] += weight
            for j in aa:
                if j != i and j % n in dset[i % n]:
                    pair[i,j] += weight
    assert sum(law.values()) == total
    pp = [Fraction(v, total) for v in one]
    lams = [sum(pp[:n]), sum(pp[n:])]
    b1 = sum((pp[i] * pp[j] for i in range(2*n) for j in range(2*n)
              if j % n in dset[i % n]), Fraction())
    b2 = Fraction(sum(pair.values()), total)
    pmf = {key: Fraction(v, total) for key, v in law.items()}
    lp, lq = map(float, lams)
    tv = 1 - sum(min(float(prob), math.exp(-lp-lq) * lp**z * lq**v /
                        (math.factorial(z)*math.factorial(v))) for (z,v),prob in pmf.items())
    assert tv <= float(2*(b1+b2)) + 1e-13
    cov = sum((z*v*pr for (z,v),pr in pmf.items()), Fraction()) - lams[0]*lams[1]
    # Coarse positive-winding events, unlike component anchors, are oppositely
    # monotone functions of the shared uniform labels, hence negatively associated.
    p_lo = sum(Fraction(v,total) for (zl,zh),v in void.items() if not zl)
    p_hi = sum(Fraction(v,total) for (zl,zh),v in void.items() if not zh)
    p_joint = Fraction(void[False,False],total)
    assert p_joint <= p_lo*p_hi
    return {'width':w, 'length':m,'cutoff_rows':height,'categorical_configurations':3**n,
            'low_probability':fstr(lowq),'high_probability':fstr(highq),
            'lambda_low':fstr(lams[0]),'lambda_high':fstr(lams[1]),
            'b1':fstr(b1),'b2':fstr(b2),'count_covariance':fstr(cov),
            'count_law': {f'{z},{v}':fstr(pr) for (z,v),pr in sorted(pmf.items())},
            'opposite_monotone_joint':fstr(p_joint),
            'opposite_monotone_product':fstr(p_lo*p_hi),
            'joint_poisson_tv_float_diagnostic':tv,
            'agg_process_tv_bound':fstr(min(Fraction(1),2*(b1+b2)))}


def exact_row_contact_control(w=3, m=5, p=Fraction(1,4)):
    """Two isolated full-row clusters may be positively correlated.

    They share a CLOSED guard row; applying BK to the anchors would be wrong.
    """
    single = p**w * (1-p)**(2*w)
    joint = p**(2*w) * (1-p)**(3*w)
    assert joint > single*single
    return {'width': w, 'length_at_least':5,'p':fstr(p),
            'single_anchor_probability':fstr(single),
            'two_anchors_two_rows_apart':fstr(joint),
            'product_of_marginals':fstr(single*single),
            'ratio':fstr(joint/(single*single)),
            'bk_applies_to_disjoint_increasing_winding_witnesses_not_anchors':True}


def geometry_checks():
    nchecks = 0
    for matching in (False, True):
        w,m,h=3,12,2
        adj=adjacency(w,m,matching)
        ladj=adjacency(w,h+2,matching,False)
        # All rows are either entirely open or entirely closed: tests seam,
        # multiple components, and cutoff coverage without random sampling.
        for rowpattern in range(1<<m):
            mask=sum(((1<<w)-1)<<(w*y) for y in range(m) if rowpattern & (1<<y))
            a=global_anchors(components(mask,w,m,adj),w,m,h)
            b=local_all_anchors(mask,w,m,h,ladj)
            assert a==b
            nchecks+=1
        ds=dependency(w,m,h)
        for i in range(w*m):
            for j in range(w*m):
                if j not in ds[i]:
                    assert not support_rows(i,w,m,h)&support_rows(j,w,m,h)
    return {'deterministic_row_masks_both_graphs':nchecks,
            'nonoverlap_neighbourhoods_checked':True}



def component_activity_counts(w: int, height: int, matching: bool):
    """Enumerate component SHAPES, not configurations of their surroundings.

    Return counts by (occupied volume, distinct external boundary volume).
    The shapes live in the interior rows 1..height of the open cylinder window.
    This is independent of local_anchor_bits and of the window probability sum.
    """
    adj = adjacency(w, height+2, matching, False)
    answer = Counter()
    for interior in range(1 << (w*height)):
        mask = interior << w
        cc = components(mask, w, height+2, adj)
        if len(cc) != 1:
            continue
        verts, cycles = cc[0]
        if not any(v // w == 1 for v in verts) or not any(dx for dx,dy in cycles):
            continue
        vs = set(verts)
        boundary = {z for v in verts for z,dx,dy in adj[v] if z not in vs}
        assert len(boundary) <= len(steps(matching))*len(vs)
        answer[len(vs),len(boundary)] += 1
    return answer


def component_activity_control(w: int, height: int, matching: bool):
    shapes = component_activity_counts(w, height, matching)
    coef = local_intensity_coefficients(w, height, matching)
    nw = w*(height+2)
    cases = []
    for p in (Fraction(1,5), Fraction(1,3), Fraction(2,5)):
        q = 1-p
        act = {nb: Fraction(c)*p**nb[0]*q**nb[1] for nb,c in shapes.items()}
        nu = sum(act.values(), Fraction())
        assert nu == bernstein_value(coef,p)
        # First path: component activity differentiation with distinct boundary.
        mean_n = sum((nb[0]*v for nb,v in act.items()),Fraction())/nu
        mean_b = sum((nb[1]*v for nb,v in act.items()),Fraction())/nu
        score = {nb: Fraction(nb[0],1)/p-Fraction(nb[1],1)/q for nb in act}
        es = sum((score[nb]*v for nb,v in act.items()),Fraction())/nu
        es2 = sum((score[nb]**2*v for nb,v in act.items()),Fraction())/nu
        contact = mean_n/p**2+mean_b/q**2
        curvature = es2-es**2-contact
        # Second path: differentiate the entire surrounding-window polynomial.
        d1 = d2 = Fraction()
        for k,c in enumerate(coef):
            weight = Fraction(c)*p**k*q**(nw-k)
            score_all = Fraction(k,1)/p-Fraction(nw-k,1)/q
            contact_all = Fraction(k,1)/p**2+Fraction(nw-k,1)/q**2
            d1 += weight*score_all
            d2 += weight*(score_all**2-contact_all)
        assert es == d1/nu
        assert curvature == d2/nu-(d1/nu)**2
        # Independent logit-coordinate score and chain rule.
        zscore = {nb: q*nb[0]-p*nb[1] for nb in act}
        ez = sum((zscore[nb]*v for nb,v in act.items()),Fraction())/nu
        ez2 = sum((zscore[nb]**2*v for nb,v in act.items()),Fraction())/nu
        zcurv = ez2-ez**2-p*q*(mean_n+mean_b)
        assert ez == p*q*es
        assert zcurv == (p*q)**2*curvature+p*q*(1-2*p)*es
        assert curvature >= -contact
        cases.append({'p':fstr(p),'nu':fstr(nu),'mean_component_volume':fstr(mean_n),
                      'mean_distinct_boundary_volume':fstr(mean_b),
                      'log_nu_derivative':fstr(es),'log_nu_second_derivative':fstr(curvature),
                      'curvature_lower_bound':fstr(-contact),
                      'logit_log_nu_second_derivative':fstr(zcurv),
                      'two_differentiation_routes_agree':True})
    return {'width':w,'cutoff_rows':height,'matching':matching,
            'interior_subsets':1<<(w*height),'surrounding_window_configurations':1<<nw,
            'component_shapes':sum(shapes.values()),
            'activity_terms':[{ 'volume':n,'boundary_volume':b,'multiplicity':c}
                              for (n,b),c in sorted(shapes.items())], 'cases':cases}


def generate():
    cases = [enumerate_one(w,m,h,g) for (w,m,h) in ((2,5,1),(2,6,2),(3,5,1))
             for g in (False, True)]
    return {'schema':'matching-one.winding-poisson-controls.v1',
            'scope':'finite exact anchoring and dependency controls; not asymptotic simulations',
            'total_graph_configurations':sum(c['configurations'] for c in cases),
            'local_window_controls':cases,
            'two_colour':two_colour_control(),
            'closed_boundary_contact':exact_row_contact_control(),
            'component_activities':[component_activity_control(w,h,g)
                                    for w,h in ((2,1),(2,2),(3,2)) for g in (False,True)],
            'geometry':geometry_checks(),
            'probability_accuracy':'all probabilities, intensities and b1/b2 are exact fractions; TV against exp is a float diagnostic',
            'nonclaims':['no numerical mass kappa', 'no prefactor A or beta determination',
                         'no proof by finite enumeration', 'no full repository CI']}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        raise SystemExit(f'refusing to overwrite {args.output}')
    data=generate()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(f'wrote {args.output}; {data["total_graph_configurations"]} graph/configuration controls')


if __name__=='__main__':
    main()
