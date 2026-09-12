#!/usr/bin/env python3
"""Independent lifted-graph checks of rank, local site sources and root response.

No boundary-state update is used to obtain graph ranks. Integer polynomials are
unnormalized Bernstein counts. Enumeration is bounded to existing tiny tori.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
from itertools import combinations
import json
from math import comb, gcd
from pathlib import Path
import time
from width4_site_sources import DEFAULT_CERTIFICATE, EXPECTED_BLOB, follow, load_certificate


def lifted_graph(word, width=4):
    """Return ambient rank and component cycle gain generators from physical edges."""
    length = len(word)
    if width < 2 or length < 2 or any(not 0 <= b < 1 << width for b in word):
        raise ValueError('honest rectangular periods and valid masks required')
    occupied = sum(mask << (width*y) for y, mask in enumerate(word))
    positions = {}
    images = []
    all_gains = []
    for root in range(width*length):
        if root in positions or not ((occupied >> root) & 1):
            continue
        positions[root] = (0, 0)
        todo = [root]
        gains = set()
        while todo:
            u = todo.pop()
            x, y = u % width, u // width
            px, py = positions[u]
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                v = ((y+dy) % length)*width + (x+dx) % width
                if not (occupied >> v) & 1:
                    continue
                candidate = (px+dx, py+dy)
                if v not in positions:
                    positions[v] = candidate
                    todo.append(v)
                else:
                    wx, wy = candidate[0]-positions[v][0], candidate[1]-positions[v][1]
                    if wx % width or wy % length:
                        raise AssertionError('invalid lifted displacement')
                    gain = (wx//width, wy//length)
                    if gain != (0, 0):
                        gains.add(gain)
        gains = sorted(gains)
        index = 0
        for a, b in combinations(gains, 2):
            index = gcd(index, abs(a[0]*b[1]-a[1]*b[0]))
        if index:
            rank = 2
            saturated = index == 1
        elif gains:
            rank = 1
            divisor = 0
            for x, y in gains:
                divisor = gcd(divisor, gcd(abs(x), abs(y)))
            saturated = divisor == 1
        else:
            rank, saturated = 0, True
        images.append({'rank': rank, 'saturated': saturated, 'rank_two_index': index,
                       'cycle_gains': gains})
        all_gains.extend(gains)
    rank = 0 if not all_gains else 1
    if any(a[0]*b[1]-a[1]*b[0] for a, b in combinations(all_gains, 2)):
        rank = 2
    return rank, images


def unnormalized_bernstein(coefficients, p):
    p = Fraction(p)
    degree = len(coefficients)-1
    return sum((c*p**k*(1-p)**(degree-k) for k, c in enumerate(coefficients)), Fraction(0))


def bernstein_derivative(coefficients):
    n = len(coefficients)-1
    return [(k+1)*coefficients[k+1]-(n-k)*coefficients[k] for k in range(n)]


def bernstein_interval(coefficients, lo, hi):
    """Outward exact termwise interval, not a floating root approximation."""
    n = len(coefficients)-1
    lower = upper = Fraction(0)
    for k, c in enumerate(coefficients):
        a = lo**k*(1-hi)**(n-k)
        b = hi**k*(1-lo)**(n-k)
        lower += c*(a if c >= 0 else b)
        upper += c*(b if c >= 0 else a)
    return lower, upper


def physical_controls(source):
    start = time.perf_counter()
    tables = {}
    total = 0
    saturation_components = 0
    for length in (2, 3, 4):
        n = 4*length
        counts = [[0]*(n+1) for _ in range(3)]
        mixed = {distance: [0]*(n-3) for distance in range(1, min(3, length-1))}
        for mask in range(1 << n):
            word = [(mask >> (4*y)) & 15 for y in range(length)]
            rank, images = lifted_graph(word)
            if rank != source['quotient_rank_output'][follow(source, word)]:
                raise AssertionError('graph and automaton disagree')
            if any(not image['saturated'] for image in images):
                raise AssertionError('embedded component has nonsaturated image')
            if sum(image['rank'] == 2 for image in images) > 1:
                raise AssertionError('disjoint independent rank-two components')
            if rank == 2 and not any(image['rank'] == 2 for image in images):
                raise AssertionError('disjoint rank-one components generate rank two')
            saturation_components += len(images)
            k = mask.bit_count()
            counts[rank][k] += 1
            for distance, coefficients in mixed.items():
                a, b = word[1], word[1+distance]
                sign = ((a & 1)-((a >> 2) & 1))*((b & 1)-((b >> 2) & 1))
                if sign:
                    coefficients[k-2] += sign*(rank-1)
            total += 1
        if any(sum(counts[j][k] for j in range(3)) != comb(n, k) for k in range(n+1)):
            raise AssertionError('occupation-count conservation failed')
        matching = [b-a for a, b in zip(counts[0], counts[2])]
        tables[str(length)] = {'sites': n, 'configurations': 1 << n,
                              'rank_bernstein_counts': counts, 'M_bernstein_counts': matching,
                              'mixed_degree_N_minus_4_coefficients': mixed,
                              'M_half': str(unnormalized_bernstein(matching, Fraction(1, 2))),
                              'mixed_at_half': {d: str(unnormalized_bernstein(c, Fraction(1, 2)))
                                                for d, c in mixed.items()}}
    matching = tables['4']['M_bernstein_counts']
    derivative = bernstein_derivative(matching)
    lo, hi = Fraction(0), Fraction(1)
    for _ in range(90):
        mid = (lo+hi)/2
        if unnormalized_bernstein(matching, mid) < 0:
            lo = mid
        else:
            hi = mid
    if not unnormalized_bernstein(matching, lo) < 0 < unnormalized_bernstein(matching, hi):
        raise AssertionError('root bracket signs failed')
    dlo, dhi = bernstein_interval(derivative, lo, hi)
    if not 0 < dlo <= dhi:
        raise AssertionError('positive root derivative not certified')
    root_response = []
    for distance, coefficients in tables['4']['mixed_degree_N_minus_4_coefficients'].items():
        clo, chi = bernstein_interval(coefficients, lo, hi)
        if not 0 < clo <= chi:
            raise AssertionError('positive mixed source response not certified')
        bound = [-chi/dlo, -clo/dhi]
        root_response.append({'row_distance': distance,
                              'mixed_M_interval': [str(clo), str(chi)],
                              'mixed_root_interval': list(map(str, bound)),
                              'mixed_M_diagnostic': float(unnormalized_bernstein(coefficients, (lo+hi)/2)),
                              'mixed_root_diagnostic': float(sum(bound)/2)})
    # A rank-one spiral is the essential counterexample to projection-label bookkeeping.
    spiral = [3, 6, 12, 9]
    spiral_rank, spiral_images = lifted_graph(spiral)
    if spiral_rank != 1 or not any(x and y for im in spiral_images for x, y in im['cycle_gains']):
        raise AssertionError('spiral control failed')
    return {'schema': 'matching-one.site-source-physical-checks.v1', 'source_blob': EXPECTED_BLOB,
            'physical_configurations_checked': total,
            'component_saturation_checks': saturation_components,
            'all_component_images_saturated': True,
            'all_rank2_images_carried_by_one_component': True,
            'tables': tables,
            'finite_4x4_root_bracket': [str(lo), str(hi)],
            'finite_4x4_root_diagnostic': float((lo+hi)/2),
            'root_Mprime_interval': [str(dlo), str(dhi)],
            'root_source_responses': root_response,
            'spiral_control': {'row_masks': spiral, 'rank': spiral_rank,
                               'components': spiral_images, 'critical_polynomial_sector': '1D'},
            'limits': ['Exact finite responses, not continuum field assignments.',
                       'Root mixed derivative is -mixed_M/Mprime only at the baseline root.',
                       'The elementary Bernstein inequalities certify signs over the whole root bracket.'],
            'elapsed_seconds': time.perf_counter()-start}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--certificate', type=Path, default=DEFAULT_CERTIFICATE)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    result = physical_controls(load_certificate(args.certificate))
    text = json.dumps(result, indent=2, allow_nan=False)+'\n'
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open('x', encoding='utf-8') as f:
            f.write(text)
    else:
        print(text, end='')


if __name__ == '__main__':
    main()
