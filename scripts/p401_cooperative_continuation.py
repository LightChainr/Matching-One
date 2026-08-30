#!/usr/bin/env python3
"""Exact finite-horizon continuation certificates for Matching One issue #401.

Uses the unchanged Gaussian rank oracle from fee33287 (blob 62e06795...).
This module adds no topology implementation and performs no Monte Carlo.
Python 3.10+, standard library. Run from the repository root:
  python scripts/p401_cooperative_continuation.py --output results/p401-cooperative-continuation/exact.json
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, permutations
import json
from math import comb, factorial, gcd
from pathlib import Path
from typing import Any, Sequence

from p334_birth_age_collision_review_20260830 import enumerate_states

State = tuple[int, tuple[int, int] | None]
SCHEMA = 'matching-one/p401-cooperative-continuation/v1'


def validate_states(n: int, states: Sequence[State]) -> None:
    """Check the rank/line contracts needed by all continuation arguments."""
    if not isinstance(n, int) or not 2 <= n <= 20 or len(states) != 1 << n:
        raise ValueError('Require 2<=N<=20 and exactly 2**N ordered states.')
    for rank, line in states:
        if type(rank) is not int or rank not in (0, 1, 2):
            raise ValueError('Rank must be an integer in {0,1,2}.')
        if rank == 1:
            if (not isinstance(line, tuple) or len(line) != 2
                    or any(type(x) is not int for x in line)):
                raise ValueError('Rank one requires an integer primitive line.')
            a, b = line
            if gcd(abs(a), abs(b)) != 1 or a < 0 or (a == 0 and b <= 0):
                raise ValueError('Line must be primitive and sign-normalized.')
        elif line is not None:
            raise ValueError('Only rank-one states carry a line.')
    if states[0][0] != 0 or states[-1][0] != 2:
        raise ValueError('Require empty rank zero and full rank two.')
    for mask, (rank, line) in enumerate(states):
        for v in range(n):
            if mask >> v & 1:
                continue
            after, next_line = states[mask | (1 << v)]
            if after < rank:
                raise ValueError('Nonmonotone rank cache.')
            if rank == after == 1 and line != next_line:
                raise ValueError('A rank-one inclusion changed its projective line.')


def continuation(n: int, states: Sequence[State], mask: int) -> dict[str, Any]:
    """Count minimal completions through order three and direct horizon events.

    The caller must validate the full cache once before calling repeatedly.
    All edges/triples are unordered; normalization uses binomial coefficients.
    """
    if not 0 <= mask < (1 << n) or states[mask][0] != 1:
        raise ValueError('Continuation starts at a valid rank-one mask.')
    vacant = [v for v in range(n) if not mask >> v & 1]
    d = len(vacant)
    hits = {v for v in vacant if states[mask | (1 << v)][0] == 2}
    safe = [v for v in vacant if v not in hits]
    edges = {pair for pair in combinations(safe, 2)
             if states[mask | sum(1 << v for v in pair)][0] == 2}
    degree = Counter(v for pair in edges for v in pair)
    wedges = sum(comb(value, 2) for value in degree.values())
    triangles = minimal_triples = 0
    for triple in combinations(safe, 3):
        edge_count = sum(pair in edges for pair in combinations(triple, 2))
        triangles += edge_count == 3
        minimal_triples += (edge_count == 0 and
                            states[mask | sum(1 << v for v in triple)][0] == 2)
    successes = [0]
    for horizon in range(1, min(3, d) + 1):
        successes.append(sum(states[mask | sum(1 << v for v in added)][0] == 2
                             for added in combinations(vacant, horizon)))
    if successes[1] != len(hits):
        raise AssertionError('One-step count mismatch.')
    if d >= 2:
        predicted = comb(d, 2) - comb(len(safe), 2) + len(edges)
        if successes[2] != predicted:
            raise AssertionError('Two-step completion identity failed.')
    if d >= 3:
        predicted = (comb(d, 3) - comb(len(safe), 3)
                     + len(edges) * (len(safe) - 2)
                     - wedges + triangles + minimal_triples)
        if successes[3] != predicted:
            raise AssertionError('Three-step overlap identity failed.')
    return {'mask': mask, 'k': mask.bit_count(), 'line': list(states[mask][1]),
            'vacancies': d, 'x': len(hits), 'c2': len(edges),
            'pair_wedges': wedges, 'pair_triangles': triangles,
            'minimal_triples': minimal_triples,
            'trigger_sites': sorted(hits), 'cooperative_pairs': [list(e) for e in sorted(edges)],
            'successful_subsets_by_horizon': successes,
            'exit_probability': {str(m): str(Fraction(successes[m], comb(d, m)))
                                 for m in range(1, len(successes))},
            'bernoulli_second_derivative_at_zero': 2 * (len(edges) - comb(len(hits), 2))}


def survival_signature(n: int, states: Sequence[State], mask: int) -> list[int]:
    """Exact b_m, not a claim that this vector is Markov under its own updates."""
    if not 0 <= mask < (1 << n) or states[mask][0] != 1:
        raise ValueError('Signature starts at a valid rank-one mask.')
    available = ((1 << n) - 1) ^ mask
    counts = [0] * (available.bit_count() + 1)
    sub = available
    while True:
        if states[mask | sub][0] == 1:
            counts[sub.bit_count()] += 1
        if sub == 0:
            return counts
        sub = (sub - 1) & available


def prefix_weights(n: int, states: Sequence[State]) -> dict[int, Counter]:
    """Count actual uniformly ordered prefixes by first essential birth."""
    result: dict[int, Counter] = {}
    for mask, (rank, _) in enumerate(states):
        if rank != 1:
            continue
        k = mask.bit_count()
        counts: Counter = Counter()
        for v in range(n):
            if not mask >> v & 1:
                continue
            previous = mask ^ (1 << v)
            if states[previous][0] == 0:
                counts[k] += factorial(k - 1)
            else:
                counts.update(result[previous])
        if sum(counts.values()) != factorial(k):
            raise AssertionError('Prefix recurrence does not reproduce k!.')
        result[mask] = counts
    return result


def age_stratum(n: int, records: dict[int, dict], weights: dict[int, Counter],
                k: int, line: tuple[int, int], x: int, horizon: int,
                c2: int | None = None) -> dict[str, Any]:
    if not 1 <= horizon <= 3 or n-k < horizon:
        raise ValueError('Unsupported finite horizon.')
    totals: Counter = Counter()
    weighted: Counter = Counter()
    cooperation: Counter = Counter()
    masks = []
    for mask, row in records.items():
        if (row['k'], tuple(row['line']), row['x']) != (k, line, x):
            continue
        if c2 is not None and row['c2'] != c2:
            continue
        masks.append(mask)
        for birth, count in weights[mask].items():
            totals[birth] += count
            weighted[birth] += count * row['successful_subsets_by_horizon'][horizon]
            cooperation[birth] += count * row['c2']
    if not totals:
        raise ValueError('Empty requested stratum.')
    rows = [{'K1': j, 'prefixes': totals[j], 'successful_subset_weight': weighted[j],
             'cooperative_pair_weight': cooperation[j],
             'mean_c2': str(Fraction(cooperation[j], totals[j])),
             'probability': str(Fraction(weighted[j], comb(n-k, horizon)*totals[j]))}
            for j in sorted(totals)]
    probabilities = [Fraction(row['probability']) for row in rows]
    return {'N': n, 'k': k, 'line': list(line), 'x': x, 'c2_control': c2,
            'horizon': horizon, 'current_subsets': len(masks), 'rows': rows,
            'latest_minus_earliest': str(probabilities[-1] - probabilities[0])}


def brute_n10(states: Sequence[State]) -> dict[str, Any]:
    """Independent of prefix DP and the cooperative-completion count formula."""
    totals: Counter = Counter()
    pair_outcomes: Counter = Counter()
    checked = 0
    for prefix in permutations(range(10), 5):
        checked += 1
        mask, birth = 0, None
        for k, vertex in enumerate(prefix, 1):
            mask |= 1 << vertex
            if birth is None and states[mask][0] >= 1:
                birth = k
        if states[mask] != (1, (0, 1)):
            continue
        vacant = [v for v in range(10) if not mask >> v & 1]
        if sum(states[mask | (1 << v)][0] == 2 for v in vacant) != 1:
            continue
        totals[birth] += 1
        pair_outcomes[birth] += sum(states[mask | (1 << u) | (1 << v)][0] == 2
                                    for u, v in combinations(vacant, 2))
    if totals != {4: 960, 5: 2640} or pair_outcomes != {4: 6240, 5: 17760}:
        raise AssertionError('Independent N10 prefix certificate failed.')
    return {'ordered_prefixes_checked': checked,
            'prefix_counts': dict(sorted(totals.items())),
            'successful_pair_weight': dict(sorted(pair_outcomes.items()))}


def build_certificate() -> dict[str, Any]:
    gates = []
    data: dict[int, tuple] = {}
    for a, b in [(2, 1), (3, 0), (3, 1), (3, 2)]:
        n, states = enumerate_states(a, b)
        validate_states(n, states)
        records = {mask: continuation(n, states, mask)
                   for mask, (rank, _) in enumerate(states) if rank == 1}
        weights = prefix_weights(n, states)
        for mask, row in records.items():
            sig = survival_signature(n, states, mask)
            d = row['vacancies']
            for m in range(1, min(3, d) + 1):
                if sig[m] + row['successful_subsets_by_horizon'][m] != comb(d, m):
                    raise AssertionError('Survival signature count mismatch.')
        gates.append({'generator': [a, b], 'N': n, 'states': len(states),
                      'rank_one_states': len(records),
                      'two_step_identity_cases': sum(r['vacancies'] >= 2 for r in records.values()),
                      'three_step_identity_cases': sum(r['vacancies'] >= 3 for r in records.values()),
                      'full_survival_signatures_checked': len(records)})
        data[n] = (states, records, weights)
    n10, r10, w10 = data[10]
    n13, r13, w13 = data[13]
    witnesses = [age_stratum(10, r10, w10, 5, (0, 1), 1, 2),
                 age_stratum(13, r13, w13, 6, (0, 1), 0, 2),
                 age_stratum(13, r13, w13, 6, (0, 1), 0, 3, c2=5)]
    if [x['latest_minus_earliest'] for x in witnesses] != ['1/44', '1/154', '2/315']:
        raise AssertionError('History witness mismatch.')
    qs = [Fraction(row['exit_probability']['2']) for row in r10.values()
          if (row['k'], row['line'], row['x']) == (5, [0, 1], 1)]
    mean = sum(qs, Fraction()) / len(qs)
    second = sum((q*q for q in qs), Fraction()) / len(qs)
    variance = second - mean*mean
    if variance != Fraction(1, 450):
        raise AssertionError('Two-future residual variance mismatch.')
    fixtures = []
    for n, masks in [(10, [155, 157]), (13, [655, 693])]:
        states, records, _ = data[n]
        for mask in masks:
            fixtures.append({'N': n, **records[mask],
                             'survival_signature': survival_signature(n, states, mask)})
    return {'schema': SCHEMA, 'status': 'exact_finite_volume',
            'topology_source': {'commit': 'fee33287cf4830e07ccef6177f43034add02256e',
                                'path': 'scripts/p334_birth_age_collision_review_20260830.py',
                                'git_blob': '62e06795fdfa91a956aedd62b7344e84aa5efc5c',
                                'changed': False},
            'gates': gates, 'history_witnesses': witnesses, 'configuration_witnesses': fixtures,
            'independent_prefix_check': brute_n10(n10),
            'two_future_control': {'N': 10, 'k': 5, 'line': [0, 1], 'x': 1,
                                   'current_subsets': len(qs), 'mean_q2': str(mean),
                                   'mean_q2_squared': str(second),
                                   'variance_q2_given_k_line_x': str(variance),
                                   'variance_after_also_controlling_c2': '0'},
            'boundaries': ['No continuum-memory or asymptotic exponent is proved.',
                           'Fixed-checkpoint reliability sufficiency does not prove Markov closure.',
                           'No Q4/Jordan identification or BA/TM/ULC counterexample.',
                           'No Monte Carlo, new production field, or repository-wide CI was run.']}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    text = json.dumps(build_certificate(), indent=2) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding='utf-8')
    print(text, end='')


if __name__ == '__main__':
    main()
