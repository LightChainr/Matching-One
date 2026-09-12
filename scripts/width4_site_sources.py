#!/usr/bin/env python3
"""Source-compatible rank quotients for the exact width-four site automaton.

Standard library only. Upstream certificate: PR #708, immutable blob checked.
Column probability groups may vary independently at every spatial row.
Coarsest strong lumping is not an arbitrary positive-realization minimum.
"""
from __future__ import annotations
import argparse
from collections import Counter, deque
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CERTIFICATE = ROOT / 'results/research-control-20260912/width4-rank-closure-certificate.json'
EXPECTED_BLOB = '50b7297deefe7c50215aea2ed534ca5810461af3'
PROFILES = {
    'homogeneous': ((0, 1, 2, 3),),
    'one_marked': ((0,), (1, 2, 3)),
    'opposite_pairs': ((0, 2), (1, 3)),
    'adjacent_pairs': ((0, 1), (2, 3)),
    'two_marked_adjacent': ((0,), (1,), (2, 3)),
    'two_marked_opposite': ((0,), (2,), (1, 3)),
    'all_independent': ((0,), (1,), (2,), (3,)),
}
EXPECTED_COUNTS = (94, 303, 179, 262, 509, 303, 509)
D4 = tuple(tuple((sign*j+k) % 4 for j in range(4))
           for sign in (1, -1) for k in range(4))


def load_certificate(path=DEFAULT_CERTIFICATE):
    raw = Path(path).read_bytes()
    blob = hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
    if blob != EXPECTED_BLOB:
        raise ValueError(f'wrong certificate blob: {blob}; expected {EXPECTED_BLOB}')
    source = json.loads(raw)
    rows = source['quotient_transitions']
    if len(rows) != 509 or any(len(row) != 16 for row in rows):
        raise ValueError('unexpected automaton dimensions')
    return source


def canonical_labels(values):
    lookup = {}
    return [lookup.setdefault(value, len(lookup)) for value in values]


def validate_groups(groups):
    groups = tuple(tuple(g) for g in groups)
    if (any(not g for g in groups) or any(type(j) is not int for g in groups for j in g)
            or sorted(j for g in groups for j in g) != list(range(4))):
        raise ValueError('groups must partition columns 0,1,2,3')
    return groups


def mask_signature(mask, groups):
    return tuple(sum((mask >> j) & 1 for j in group) for group in groups)


def refine(source, groups):
    """Polynomial coefficient refinement for all independently varying group p's."""
    groups = validate_groups(groups)
    rows = source['quotient_transitions']
    keys = [mask_signature(mask, groups) for mask in range(16)]
    labels = canonical_labels(source['quotient_rank_output'])
    ladder = [len(set(labels))]
    while True:
        lookup = {}
        new = []
        for i, row in enumerate(rows):
            totals = Counter((keys[mask], labels[j]) for mask, j in enumerate(row))
            key = (labels[i], tuple(sorted(totals.items())))
            new.append(lookup.setdefault(key, len(lookup)))
        ladder.append(len(lookup))
        if new == labels:
            return {'groups': groups, 'labels': labels, 'ladder': ladder,
                    'classes': len(lookup)}
        labels = new


def transform_mask(mask, permutation):
    return sum(((mask >> j) & 1) << permutation[j] for j in range(4))


def follow(source, word):
    if not word or any(type(mask) is not int or not 0 <= mask < 16 for mask in word):
        raise ValueError('word must contain row masks in [0,15]')
    state = source['quotient_initial'][word[0]]
    for mask in word[1:]:
        state = source['quotient_transitions'][state][mask]
    return state


def histories(source):
    """A physical (at least two rows) word for every deterministic rank class."""
    result = {}
    todo = deque()
    for a in range(16):
        for b in range(16):
            state = follow(source, [a, b])
            if state not in result:
                result[state] = [a, b]
                todo.append(state)
    while todo:
        state = todo.popleft()
        for mask, target in enumerate(source['quotient_transitions'][state]):
            if target not in result:
                result[target] = result[state]+[mask]
                todo.append(target)
    if len(result) != 509:
        raise AssertionError('not all rank classes have physical histories')
    return [result[i] for i in range(509)]


def group_actions(source, words=None):
    words = histories(source) if words is None else words
    rows = source['quotient_transitions']
    output = source['quotient_rank_output']
    actions = []
    for permutation in D4:
        action = [follow(source, [transform_mask(b, permutation) for b in word])
                  for word in words]
        if sorted(action) != list(range(509)):
            raise AssertionError('group action is not a permutation')
        for i, row in enumerate(rows):
            if output[action[i]] != output[i]:
                raise AssertionError('rank not invariant')
            for mask, j in enumerate(row):
                if action[j] != rows[action[i]][transform_mask(mask, permutation)]:
                    raise AssertionError('transition action fails')
        actions.append(action)
    return actions


def row_integer_weights(probabilities):
    """Integer row weights with one common denominator; endpoints are supported."""
    if len(probabilities) != 4:
        raise ValueError('exactly four probabilities required')
    probabilities = [Fraction(x) for x in probabilities]
    if any(p < 0 or p > 1 for p in probabilities):
        raise ValueError('site probabilities must lie in [0,1]')
    denominator = 1
    for p in probabilities:
        denominator *= p.denominator
    weights = []
    for mask in range(16):
        value = 1
        for j, p in enumerate(probabilities):
            value *= p.numerator if (mask >> j) & 1 else p.denominator-p.numerator
        weights.append(value)
    if sum(weights) != denominator:
        raise AssertionError('row weights not normalized')
    return weights, denominator


def signed_output(source, weights_by_row, denominators):
    """Exact M or a signed-source coefficient; does not require positive weights."""
    if len(weights_by_row) < 2 or len(denominators) != len(weights_by_row):
        raise ValueError('at least two equally specified rows required')
    n = len(source['quotient_transitions'])
    vector = [0]*n
    denominator = 1
    for mask, j in enumerate(source['quotient_initial']):
        vector[j] += weights_by_row[0][mask]
    for d in denominators:
        denominator *= d
    for weights in weights_by_row[1:]:
        moved = [0]*n
        for i, value in enumerate(vector):
            if value:
                for mask, j in enumerate(source['quotient_transitions'][i]):
                    moved[j] += value*weights[mask]
        vector = moved
    return Fraction(sum(value*(r-1) for value, r in
                        zip(vector, source['quotient_rank_output'])), denominator)


def product_measure_M(source, probability_rows):
    entries = [row_integer_weights(row) for row in probability_rows]
    return signed_output(source, [a for a, _ in entries], [d for _, d in entries])


def dipole_response(source, length, source_rows, p=Fraction(1, 2)):
    """Mixed spatial-row dipole derivatives at any rational interior baseline p."""
    if type(length) is not int or length < 2:
        raise ValueError('length must be an integer at least two')
    if len(set(source_rows)) != len(source_rows) or any(type(y) is not int or y < 0 or y >= length for y in source_rows):
        raise ValueError('distinct valid source rows required')
    p = Fraction(p)
    if not 0 < p < 1:
        raise ValueError('interior probability required for signed local sources')
    a, d = p.numerator, p.denominator
    base, denominator = row_integer_weights([p]*4)
    odd = []
    for mask, weight in enumerate(base):
        numerator = weight*d*d*((mask & 1)-((mask >> 2) & 1))
        value, remainder = divmod(numerator, a*(d-a))
        if remainder:
            raise AssertionError('dipole derivative should have integer coefficients')
        odd.append(value)
    return signed_output(source, [odd if y in source_rows else base for y in range(length)],
                         [denominator]*length)


def dipole_response_half(source, length, source_rows):
    return dipole_response(source, length, source_rows, Fraction(1, 2))


def four_sign(source, length, row1, row2, eps1, eps2, p=Fraction(1, 2)):
    if row1 == row2 or not 0 <= row1 < length or not 0 <= row2 < length:
        raise ValueError('two distinct source rows required')
    eps1, eps2, p = Fraction(eps1), Fraction(eps2), Fraction(p)
    if not eps1 or not eps2:
        raise ValueError('nonzero amplitudes required')
    total = Fraction(0)
    corners = {}
    for s in (-1, 1):
        for t in (-1, 1):
            ps = [[p]*4 for _ in range(length)]
            for y, sign, eps in ((row1, s, eps1), (row2, t, eps2)):
                ps[y][0] += sign*eps
                ps[y][2] -= sign*eps
            value = product_measure_M(source, ps)
            corners[f'{s},{t}'] = str(value)
            total += s*t*value
    return total/(4*eps1*eps2), corners


def profile_report(source):
    words = histories(source)
    actions = group_actions(source, words)
    for a, pa in enumerate(D4):
        for b, pb in enumerate(D4):
            c = D4.index(tuple(pa[pb[j]] for j in range(4)))
            if any(actions[a][actions[b][i]] != actions[c][i] for i in range(509)):
                raise AssertionError('D4 multiplication law fails')
    profiles = {}
    for (name, groups), expected in zip(PROFILES.items(), EXPECTED_COUNTS):
        result = refine(source, groups)
        eligible = [k for k, permutation in enumerate(D4)
                    if all({permutation[j] for j in group} == set(group) for group in groups)]
        orbit_labels = canonical_labels(tuple(sorted({actions[k][i] for k in eligible}))
                                        for i in range(509))
        if result['labels'] != orbit_labels or result['classes'] != expected:
            raise AssertionError(f'orbit/refinement disagreement: {name}')
        result['label_preserving_D4_indices'] = eligible
        result['equals_orbit_partition_blockwise'] = True
        profiles[name] = result
    # Check the 94-class failure without using orbit names as proof.
    a, b = follow(source, [1, 1]), follow(source, [2, 2])
    hom = profiles['homogeneous']['labels']
    if hom[a] != hom[b]:
        raise AssertionError('witness not in one homogeneous block')
    p, q = Fraction(2, 5), Fraction(3, 5)
    weights, denominator = row_integer_weights([q, p, p, p])
    output = source['quotient_rank_output']
    rows = source['quotient_transitions']
    difference = Fraction(sum(weights[mask]*(output[rows[a][mask]]-output[rows[b][mask]])
                              for mask in range(16)), denominator)
    if difference != q-p:
        raise AssertionError('marked-column witness failed')
    response = []
    for m in (3, 4, 6, 8, 12):
        one = dipole_response_half(source, m, (1,))
        if one:
            raise AssertionError('odd linear response not zero')
        for distance in range(1, min(4, m-1)):
            mixed = dipole_response_half(source, m, (1, 1+distance))
            response.append({'length': m, 'row_distance': distance, 'linear': '0',
                             'mixed_exact': str(mixed), 'mixed_diagnostic': float(mixed)})
    finite = []
    for eps1, eps2 in ((Fraction(1, 20), Fraction(1, 7)),
                       (Fraction(1, 3), Fraction(1, 4))):
        value, corners = four_sign(source, 4, 1, 2, eps1, eps2)
        if value != Fraction(327, 1024):
            raise AssertionError('finite-amplitude identity failed')
        finite.append({'amplitudes': [str(eps1), str(eps2)], 'corners': corners,
                       'four_sign_divided_difference': str(value)})
    return {'schema': 'matching-one.width4-site-sources.v1', 'source_blob': EXPECTED_BLOB,
            'scope': 'coarsest common strong rank lumping on PR708; no arbitrary positive minimum',
            'profiles': profiles, 'D4_permutations': D4,
            'D4_actions': actions, 'physical_class_histories': words,
            'group_transition_checks': 8*509*16,
            'witness': {'histories': [[1, 1], [2, 2]], 'states': [a, b],
                        'same_homogeneous_class': hom[a], 'same_occupation': 2,
                        'same_immediate_rank': 1, 'next_row_probabilities': ['q', 'p', 'p', 'p'],
                        'rank_difference_polynomial': 'q-p', 'exact_check': str(difference)},
            'dipole_responses_at_p_half': response, 'finite_amplitude_checks': finite,
            'limits': ['Seven column set partitions modulo D4 exhaust the grouping cases.',
                       'Spatial occupation sources, not continuous-time P398 generators.',
                       'No continuum spin, original-U candidate, or new sampling verdict.']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--certificate', type=Path, default=DEFAULT_CERTIFICATE)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    start = time.perf_counter()
    result = profile_report(load_certificate(args.certificate))
    result['elapsed_seconds'] = time.perf_counter()-start
    text = json.dumps(result, indent=2, allow_nan=False)+'\n'
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open('x', encoding='utf-8') as f:
            f.write(text)
    else:
        print(text, end='')


if __name__ == '__main__':
    main()
