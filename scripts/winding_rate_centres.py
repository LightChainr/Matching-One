#!/usr/bin/env python3
"""Finite controls for the microscopic winding-cost / two-centre theorem.

Only deterministic finite product measures are computed. No number returned
here is an estimate of the infinite-lattice inverse correlation length.
Run from the repository root; output files are never overwritten.
"""
from __future__ import annotations
import argparse
from collections import deque
from fractions import Fraction
import json
from math import comb
from pathlib import Path
from typing import Sequence

STEPS4 = ((1, 0), (-1, 0), (0, 1), (0, -1))
STEPS8 = STEPS4 + ((1, 1), (1, -1), (-1, 1), (-1, -1))


def steps(kind: str) -> tuple[tuple[int, int], ...]:
    if kind not in ('NN', 'matching'):
        raise ValueError('kind must be NN or matching')
    return STEPS4 if kind == 'NN' else STEPS8


def fraction_record(value: Fraction) -> dict[str, str]:
    return {'numerator': str(value.numerator), 'denominator': str(value.denominator)}


def decimal(value: Fraction, places: int = 24) -> str:
    from decimal import Decimal, localcontext
    with localcontext() as ctx:
        ctx.prec = places + 12
        return format(Decimal(value.numerator) / Decimal(value.denominator), f'.{places}f')


def bernstein_count_value(counts: Sequence[int], p: Fraction) -> Fraction:
    """counts[k] is a count, not an already normalized Bernstein coefficient."""
    if not 0 <= p <= 1:
        raise ValueError('probability outside [0,1]')
    n = len(counts) - 1
    return sum((Fraction(c) * p**k * (1-p)**(n-k) for k, c in enumerate(counts)), Fraction())


def finite_adjacency(width: int, height: int, kind: str) -> tuple[tuple[int, ...], ...]:
    result = []
    for y in range(height):
        for x in range(width):
            result.append(tuple(yy * width + xx for dx, dy in steps(kind)
                                if 0 <= (xx := x+dx) < width
                                and 0 <= (yy := y+dy) < height))
    return tuple(result)


def reached(mask: int, adjacency: Sequence[Sequence[int]], source: int) -> set[int]:
    if not (mask >> source) & 1:
        return set()
    seen = {source}
    stack = [source]
    while stack:
        u = stack.pop()
        for v in adjacency[u]:
            if ((mask >> v) & 1) and v not in seen:
                seen.add(v)
                stack.append(v)
    return seen


def seed_counts(kind: str) -> list[int]:
    """P((0,1)<->(2,1) inside 3x3 | (0,1) occupied).

    Eight independent bits remain. This is an explicit finite witness and
    hence only provides an UPPER bound -log(q)/2 on the planar mass.
    """
    adjacency = finite_adjacency(3, 3, kind)
    source, target = 3, 5
    others = [i for i in range(9) if i != source]
    counts = [0] * 9
    for word in range(1 << 8):
        mask = 1 << source
        for j, i in enumerate(others):
            mask |= ((word >> j) & 1) << i
        if target in reached(mask, adjacency, source):
            counts[word.bit_count()] += 1
    return counts


def seed_root_interval(counts: Sequence[int], target: Fraction = Fraction(1, 16),
                       bits: int = 64) -> tuple[Fraction, Fraction]:
    """Enclose the finite seed root q(p)=target; NOT the torus birth centre."""
    lo, hi = Fraction(0), Fraction(1)
    for _ in range(bits):
        mid = (lo+hi)/2
        if bernstein_count_value(counts, mid) < target:
            lo = mid
        else:
            hi = mid
    assert bernstein_count_value(counts, lo) <= target <= bernstein_count_value(counts, hi)
    return lo, hi


def torus_adjacency(width: int, height: int, kind: str):
    """Keep parallel periodic edges and their actual lifted steps."""
    return tuple(tuple((((x+dx) % width) + width*((y+dy) % height), dx, dy)
                       for dx, dy in steps(kind))
                 for y in range(height) for x in range(width))


def winding_walk(mask: int, width: int, height: int, kind: str,
                 adjacency=None) -> list[tuple[int, int]] | None:
    """Independent graph-potential traversal returning a lifted closed walk."""
    adjacency = adjacency or torus_adjacency(width, height, kind)
    potentials: dict[int, tuple[int, int]] = {}
    parents: dict[int, int | None] = {}
    for root in range(width*height):
        if not ((mask >> root) & 1) or root in potentials:
            continue
        potentials[root] = (root % width, root // width)
        parents[root] = None
        queue = deque([root])
        while queue:
            u = queue.popleft()
            ux, uy = potentials[u]
            for v, dx, dy in adjacency[u]:
                if not ((mask >> v) & 1):
                    continue
                candidate = (ux+dx, uy+dy)
                if v not in potentials:
                    potentials[v] = candidate
                    parents[v] = u
                    queue.append(v)
                    continue
                vx, vy = potentials[v]
                gain = (candidate[0]-vx, candidate[1]-vy)
                if gain == (0, 0):
                    continue
                assert gain[0] % width == gain[1] % height == 0
                up = []
                node = u
                while node is not None:
                    up.append(node)
                    node = parents[node]
                walk = [potentials[node] for node in reversed(up)]
                walk.append(candidate)
                node = parents[v]
                while node is not None:
                    nx, ny = potentials[node]
                    walk.append((nx+gain[0], ny+gain[1]))
                    node = parents[node]
                assert walk[-1] != walk[0]
                assert ((walk[-1][0]-walk[0][0]) % width == 0 and
                        (walk[-1][1]-walk[0][1]) % height == 0)
                return walk
    return None


def cut_witness(walk: Sequence[tuple[int, int]], width: int):
    """First range reaching w-1, not a radius-w/2 estimate.

    The entire retained prefix lies in ONE injecting w-by-w vertex square.
    All used edges are its planar edges, never the extra torus seam edges.
    """
    if width < 2 or not walk:
        raise ValueError('nonempty walk and width >=2 required')
    x0, y0 = walk[0]
    xmin = xmax = x0
    ymin = ymax = y0
    for j, (x, y) in enumerate(walk):
        xmin, xmax = min(xmin, x), max(xmax, x)
        ymin, ymax = min(ymin, y), max(ymax, y)
        if max(xmax-xmin, ymax-ymin) == width-1:
            return list(walk[:j+1]), (xmin, ymin), ('x' if xmax-xmin == width-1 else 'y')
        assert max(xmax-xmin, ymax-ymin) < width-1
    raise AssertionError('nonzero winding did not span the injecting square')


def verify_cut(mask: int, width: int, height: int, kind: str, walk) -> None:
    prefix, (xmin, ymin), axis = cut_witness(walk, width)
    for x, y in prefix:
        assert xmin <= x < xmin+width and ymin <= y < ymin+width
        assert (mask >> ((x % width)+width*(y % height))) & 1
    for (x, y), (xx, yy) in zip(prefix, prefix[1:]):
        assert (xx-x, yy-y) in steps(kind)
    # Independent planar finite-box connectivity check on the projected mask.
    boxmask = 0
    for by in range(width):
        for bx in range(width):
            v = ((xmin+bx) % width)+width*((ymin+by) % height)
            boxmask |= ((mask >> v) & 1) << (by*width+bx)
    adjacency = finite_adjacency(width, width, kind)
    side0 = [width*j for j in range(width)] if axis == 'x' else list(range(width))
    side1 = {width*j+width-1 for j in range(width)} if axis == 'x' else set(range(width*(width-1), width*width))
    assert any(reached(boxmask, adjacency, source) & side1 for source in side0)


def rectangle_connection(mask: int, width: int, height: int, kind: str,
                         xleft: int, length: int, row: int = 1) -> bool:
    """A length+1 by 3 planar seed, with coordinates mapped periodically."""
    localmask = 0
    for y in range(3):
        for x in range(length+1):
            v = (xleft+x) % width + width*((row-1+y) % height)
            localmask |= ((mask >> v) & 1) << (y*(length+1)+x)
    adj = finite_adjacency(length+1, 3, kind)
    return 2*(length+1)-1 in reached(localmask, adj, length+1)


def seed_ring(mask: int, width: int, height: int, kind: str) -> bool:
    """Concatenate length-2 seeds, then force at most one leftover edge."""
    if width < 3 or height < 3:
        raise ValueError('individual 3x3 seeds must inject')
    k, remainder = divmod(width, 2)
    if not all(rectangle_connection(mask, width, height, kind, 2*j, 2) for j in range(k)):
        return False
    if remainder:
        return bool((mask >> (width+width-1)) & 1 and (mask >> width) & 1)
    return True


def seed_and_cluster_controls() -> dict:
    result = {}
    for kind in ('NN', 'matching'):
        counts = seed_counts(kind)
        lo, hi = seed_root_interval(counts)
        # A direct finite-cluster comparison, independent of the asymptotic proof.
        adj = finite_adjacency(3, 3, kind)
        p, q, cap = Fraction(1, 5), Fraction(1, 4), 5
        tau_p = Fraction()
        tau_q_small = Fraction()
        connected_count = 0
        for mask in range(1 << 9):
            cluster = reached(mask, adj, 3)
            if 5 not in cluster:
                continue
            connected_count += 1
            n = mask.bit_count()
            tau_p += p**n*(1-p)**(9-n)
            if len(cluster) <= cap:
                tau_q_small += q**n*(1-q)**(9-n)
        assert tau_q_small <= (q/p)**cap*tau_p
        qquarter = bernstein_count_value(counts, Fraction(1,4))
        assert qquarter > Fraction(1,16)
        result[kind] = {
            'conditional_origin': True,
            'seed_vertex_count': 9,
            'independent_remaining_sites': 8,
            'source': [0,1], 'target': [2,1],
            'conditional_connection_counts_by_occupied_other_sites': counts,
            'q_at_one_quarter': fraction_record(qquarter),
            'seed_root_target': '1/16',
            'seed_root_interval': [fraction_record(lo), fraction_record(hi)],
            'seed_root_midpoint_decimal': decimal((lo+hi)/2),
            'strict_upper_bound_on_first_centre_at_d_log4': decimal(hi),
            'cluster_comparison': {'cap': cap, 'p': str(p), 'q': str(q),
                                   'lhs': fraction_record(tau_q_small),
                                   'rhs': fraction_record((q/p)**cap*tau_p),
                                   'passed': True},
        }
    # NN centre a is <= its seed root; b is >= one minus matching seed root.
    result['certified_d_log4_centre_brackets'] = {
        'a_lower': '1/12',
        'a_upper': result['NN']['seed_root_interval'][1],
        'b_lower': fraction_record(1-Fraction(int(result['matching']['seed_root_interval'][1]['numerator']),
                                             int(result['matching']['seed_root_interval'][1]['denominator']))),
        'b_upper': '27/28',
        'warning': 'Seed roots bound the infinite centres; they are NOT centre estimates.'
    }
    return result


def census_control(width: int, height: int, kind: str, ring_check: bool = False) -> dict:
    n = width*height
    if height < width:
        raise ValueError('cut control assumes w<=m')
    adj = torus_adjacency(width, height, kind)
    winding_counts = [0]*(n+1)
    ring_counts = [0]*(n+1)
    windings = 0
    for mask in range(1 << n):
        walk = winding_walk(mask, width, height, kind, adj)
        if walk is not None:
            windings += 1
            winding_counts[mask.bit_count()] += 1
            verify_cut(mask, width, height, kind, walk)
        if ring_check and seed_ring(mask, width, height, kind):
            assert walk is not None
            ring_counts[mask.bit_count()] += 1
    finite_checks = []
    for p in (Fraction(1,100), Fraction(1,4), Fraction(1,2)):
        f = bernstein_count_value(winding_counts, p)
        item = {'p': str(p), 'winding_probability': fraction_record(f)}
        branch = len(steps(kind))-1
        if branch*p < 1:
            # kappa >= -log(branch*p) from nonbacktracking path counting.
            upper = 2*n*width**2*p*(branch*p)**(width-1)
            assert f <= upper
            item['cut_upper_using_path_mass_lower'] = fraction_record(upper)
        if ring_check:
            qseed = bernstein_count_value(seed_counts(kind), p)
            k, remainder = divmod(width,2)
            lower = p**(remainder+1)*qseed**k
            ring = bernstein_count_value(ring_counts, p)
            assert lower <= ring <= f
            item['conditional_harris_lower'] = fraction_record(lower)
            item['ring_probability'] = fraction_record(ring)
        finite_checks.append(item)
    return {'width': width, 'height': height, 'graph': kind,
            'configurations': 1 << n, 'nonzero_winding_configurations': windings,
            'first_span_cut_witnesses_checked': windings,
            'winding_counts': winding_counts,
            'ring_counts': ring_counts if ring_check else None,
            'rational_probability_checks': finite_checks}


def run_controls() -> dict:
    seeds = seed_and_cluster_controls()
    censuses = [census_control(w,h,kind,ring_check=(w==3 and h==3))
                for w,h in ((3,3),(3,4),(4,4)) for kind in ('NN','matching')]
    return {
        'schema': 'matching-one.winding-rate-centres.v1',
        'date': '2026-09-13',
        'claim_boundary': 'Finite deterministic controls only; no infinite correlation length numerically computed.',
        'seed_bounds': seeds,
        'censuses': censuses,
        'total_graph_configuration_checks': sum(c['configurations'] for c in censuses),
        'total_cut_witnesses': sum(c['first_span_cut_witnesses_checked'] for c in censuses),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    # Reserve before computing to ensure no result file is silently overwritten.
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(run_controls(), handle, ensure_ascii=False, indent=2)
        handle.write('\n')

if __name__ == '__main__':
    main()
