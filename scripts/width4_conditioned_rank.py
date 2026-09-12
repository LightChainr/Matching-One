#!/usr/bin/env python3
"""Exact finite-horizon conditioning on a final torus rank, using PR708 states.

Positive integer backward messages implement the elementary Doob transform.
This is not a new general sampling theorem and not a large-width cost claim.
No burn-in, topology-event rejection, or floating-point sampling weights.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from random import Random
import time
from width4_site_sources import (DEFAULT_CERTIFICATE, EXPECTED_BLOB, follow,
                                 load_certificate, row_integer_weights)
from site_source_graph_checks import lifted_graph


def integer_choice(weights, rng):
    if any(type(x) is not int or x < 0 for x in weights):
        raise ValueError('nonnegative integer weights required')
    total = sum(weights)
    if not total:
        raise ValueError('cannot sample zero total weight')
    draw = rng.randrange(total)
    for i, weight in enumerate(weights):
        if draw < weight:
            return i
        draw -= weight
    raise AssertionError('integer choice failed')


class RankBridge:
    """Condition independent site occupations on final rank 0, 1, or 2."""
    def __init__(self, source, probability_rows, target):
        if type(target) is not int or target not in (0, 1, 2):
            raise ValueError('rank target must be 0,1,2')
        if len(probability_rows) < 2:
            raise ValueError('at least two physical rows required')
        self.source = source
        self.target = target
        entries = [row_integer_weights(row) for row in probability_rows]
        self.weights = [w for w, _ in entries]
        self.denominator = 1
        for _, d in entries:
            self.denominator *= d
        self.length = len(entries)
        # h[y](state after row y-1) includes pending spatial rows y,...,m-1.
        self.h = [None]*(self.length+1)
        self.h[self.length] = [int(r == target) for r in source['quotient_rank_output']]
        for y in range(self.length-1, 0, -1):
            hnext = self.h[y+1]
            weights = self.weights[y]
            self.h[y] = [sum(weights[mask]*hnext[j] for mask, j in enumerate(row))
                         for row in source['quotient_transitions']]
        self.first_weights = [self.weights[0][mask]*self.h[1][j]
                              for mask, j in enumerate(source['quotient_initial'])]
        self.mass = sum(self.first_weights)
        if not self.mass:
            raise ValueError('conditioning event has zero probability')
        self.probability = Fraction(self.mass, self.denominator)

    def sample(self, rng):
        """Exact under a uniform integer random-bit source; seeded PRNG for controls."""
        mask = integer_choice(self.first_weights, rng)
        word = [mask]
        state = self.source['quotient_initial'][mask]
        for y in range(1, self.length):
            row = self.source['quotient_transitions'][state]
            weights = [self.weights[y][b]*self.h[y+1][j] for b, j in enumerate(row)]
            if sum(weights) != self.h[y][state]:
                raise AssertionError('backward normalizer inconsistent')
            mask = integer_choice(weights, rng)
            word.append(mask)
            state = row[mask]
        if self.source['quotient_rank_output'][state] != self.target:
            raise AssertionError('wrong final rank')
        return word

    def sequential_word_probability(self, word):
        if len(word) != self.length or any(type(b) is not int or not 0 <= b < 16 for b in word):
            raise ValueError('word has incorrect length or masks')
        mask = word[0]
        numerator = self.first_weights[mask]
        if not numerator:
            return Fraction(0)
        result = Fraction(numerator, self.mass)
        state = self.source['quotient_initial'][mask]
        for y in range(1, self.length):
            target = self.source['quotient_transitions'][state][word[y]]
            numerator = self.weights[y][word[y]]*self.h[y+1][target]
            if not numerator:
                return Fraction(0)
            result *= Fraction(numerator, self.h[y][state])
            state = target
        return result

    def direct_word_probability(self, word):
        if len(word) != self.length or any(type(b) is not int or not 0 <= b < 16 for b in word):
            raise ValueError('word has incorrect length or masks')
        state = follow(self.source, word)
        if self.source['quotient_rank_output'][state] != self.target:
            return Fraction(0)
        mass = 1
        for y, mask in enumerate(word):
            mass *= self.weights[y][mask]
        return Fraction(mass, self.mass)


def sector_moments(source, probability_rows):
    """Exact total occupation moments K,K^2; no numerical differentiation."""
    if len(probability_rows) < 2:
        raise ValueError('at least two rows required')
    entries = [row_integer_weights(row) for row in probability_rows]
    size = len(source['quotient_transitions'])
    mass, first, second = ([0]*size for _ in range(3))
    for b, state in enumerate(source['quotient_initial']):
        value = entries[0][0][b]
        k = b.bit_count()
        mass[state] += value
        first[state] += k*value
        second[state] += k*k*value
    denominator = entries[0][1]
    for weights, den in entries[1:]:
        denominator *= den
        new, new1, new2 = ([0]*size for _ in range(3))
        for i, value in enumerate(mass):
            if not value:
                continue
            for b, j in enumerate(source['quotient_transitions'][i]):
                w, k = weights[b], b.bit_count()
                new[j] += w*value
                new1[j] += w*(first[i]+k*value)
                new2[j] += w*(second[i]+2*k*first[i]+k*k*value)
        mass, first, second = new, new1, new2
    result = []
    for rank in range(3):
        indices = [i for i, r in enumerate(source['quotient_rank_output']) if r == rank]
        z = sum(mass[i] for i in indices)
        z1 = sum(first[i] for i in indices)
        z2 = sum(second[i] for i in indices)
        result.append({'mass': z, 'denominator': denominator, 'first_weighted': z1,
                       'second_weighted': z2, 'probability': Fraction(z, denominator),
                       'mean': Fraction(z1, z) if z else None,
                       'variance': Fraction(z2, z)-Fraction(z1, z)**2 if z else None})
    if sum(r['probability'] for r in result) != 1:
        raise AssertionError('total probability is not one')
    ps = [Fraction(p) for row in probability_rows for p in row]
    mean = sum(ps)
    total_first = sum(r['first_weighted'] for r in result)
    total_second = sum(r['second_weighted'] for r in result)
    if Fraction(total_first, denominator) != mean:
        raise AssertionError('unconditional mean mismatch')
    if Fraction(total_second, denominator)-mean**2 != sum(p*(1-p) for p in ps):
        raise AssertionError('unconditional variance mismatch')
    return result


def bridge_report(source):
    started = time.perf_counter()
    small_rows = [[Fraction(1, 3), Fraction(2, 5), Fraction(3, 7), Fraction(4, 9)],
                  [Fraction(3, 4), Fraction(2, 3), Fraction(3, 5), Fraction(1, 2)],
                  [Fraction(2, 7), Fraction(4, 5), Fraction(1, 4), Fraction(5, 6)]]
    bridges = [RankBridge(source, small_rows, r) for r in range(3)]
    moments = sector_moments(source, small_rows)
    sums = [[0, 0, 0] for _ in range(3)]
    checked = 0
    for mask in range(1 << 12):
        word = [(mask >> (4*y)) & 15 for y in range(3)]
        rank = lifted_graph(word)[0]
        mass = 1
        for y, b in enumerate(word):
            mass *= bridges[rank].weights[y][b]
        k = mask.bit_count()
        sums[rank][0] += mass
        sums[rank][1] += k*mass
        sums[rank][2] += k*k*mass
        for bridge in bridges:
            if bridge.sequential_word_probability(word) != bridge.direct_word_probability(word):
                raise AssertionError('conditional path law fails')
            checked += 1
    for rank, (z, z1, z2) in enumerate(sums):
        if (z, z1, z2) != (moments[rank]['mass'], moments[rank]['first_weighted'],
                           moments[rank]['second_weighted']):
            raise AssertionError('independent conditional occupation moments fail')
        if bridges[rank].mass != z:
            raise AssertionError('independent rank normalizer fails')
    # Fixed rational p close to q4, not the exact root and not a new p_c value.
    p = Fraction(591417, 1000000)
    length = 128
    ps = [[p]*4 for _ in range(length)]
    begin = time.perf_counter()
    rare_moments = sector_moments(source, ps)
    rare = []
    for rank in (0, 2):
        bridge = RankBridge(source, ps, rank)
        if bridge.mass != rare_moments[rank]['mass']:
            raise AssertionError('forward/backward rare normalizers disagree')
        controls = []
        for seed in (20260912, 20260913, 20260914, 20260915):
            word = bridge.sample(Random(seed))
            graph_rank = lifted_graph(word)[0]
            if graph_rank != rank:
                raise AssertionError('conditioned sample violates independent rank oracle')
            if bridge.sequential_word_probability(word) != bridge.direct_word_probability(word):
                raise AssertionError('long path telescoping fails')
            controls.append({'seed': seed, 'physical_rank': graph_rank,
                             'occupied_sites': sum(b.bit_count() for b in word),
                             'row_masks': word,
                             'row_mask_sha256': hashlib.sha256(bytes(word)).hexdigest()})
        rare.append({'rank': rank, 'probability_exact': str(bridge.probability),
                     'probability_diagnostic': float(bridge.probability),
                     'mean_occupation_exact': str(rare_moments[rank]['mean']),
                     'mean_occupation_diagnostic': float(rare_moments[rank]['mean']),
                     'variance_occupation_diagnostic': float(rare_moments[rank]['variance']),
                     'counting_bits_in_normalizer': bridge.mass.bit_length(),
                     'stored_backward_state_values': (length)*509,
                     'bounded_algorithm_control_draws': controls})
    score = (rare_moments[2]['mean']-rare_moments[0]['mean'])/(p*(1-p))
    union_probability = rare_moments[0]['probability']+rare_moments[2]['probability']
    return {'schema': 'matching-one.width4-exact-rank-bridge.v1', 'source_blob': EXPECTED_BLOB,
            'small_inhomogeneous_probability_rows': [[str(p) for p in row] for row in small_rows],
            'small_physical_configurations': 4096, 'conditional_path_equalities': checked,
            'small_sector_probabilities': [str(b.probability) for b in bridges],
            'small_conditional_mean_occupations': [str(r['mean']) for r in moments],
            'rare_control': {'length': length, 'sites': 4*length, 'p_exact': str(p),
                             'note': 'rational near q4; not an exact balance root',
                             'sectors': rare,
                             'probability_rank_not_1': float(union_probability),
                             'rank_odds_thermal_derivative': float(score),
                             'rank_odds_thermal_derivative_exact_hex': {
                                 'numerator': hex(score.numerator), 'denominator': hex(score.denominator)},
                             'elapsed_seconds': time.perf_counter()-begin},
            'limits': ['Eight fixed-seed draws are algorithm controls, not a new probability estimate.',
                       'Each completed configuration is checked by physical lifted-graph traversal.',
                       'Uniform integer bits give the exact conditional law; a seeded PRNG is used in controls.',
                       'No rejection of a whole percolation configuration, burn-in, or asymptotic conditioning.',
                       'Conditioning requires exact backward normalizers on the certified finite state space.',
                       'No cost guarantee as circumference increases, and no comparison with optimized transfer codes.'],
            'elapsed_seconds': time.perf_counter()-started}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--certificate', type=Path, default=DEFAULT_CERTIFICATE)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    result = bridge_report(load_certificate(args.certificate))
    text = json.dumps(result, indent=2, allow_nan=False)+'\n'
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open('x', encoding='utf-8') as f:
            f.write(text)
    else:
        print(text, end='')


if __name__ == '__main__':
    main()
