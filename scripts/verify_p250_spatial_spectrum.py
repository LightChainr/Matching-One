#!/usr/bin/env python3
"""Independent verifier for the N505 spatial-spectrum witness certificate.

Uses cyclic quotient labels and weighted union-find, not the generator's
adjugate-pair keys or lifted BFS. Does not import the discovery executable.
The written proof, not a floating DFT, justifies the cyclotomic rank bound.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from math import gcd
from pathlib import Path

SCHEMA = "matching-one/p250-spatial-spectrum-certificate/v1"
NN_FORWARD = ((1, 0), (0, 1))
MATCHING_FORWARD = NN_FORWARD + ((1, 1), (1, -1))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


class LiftUnionFind:
    def __init__(self, n: int) -> None:
        self.parent = list(range(n))
        self.size = [1] * n
        self.dx = [0] * n
        self.dy = [0] * n
        self.basis = [[] for _ in range(n)]

    def find(self, v: int):
        if self.parent[v] != v:
            old = self.parent[v]
            root, x, y = self.find(old)
            self.dx[v] += x
            self.dy[v] += y
            self.parent[v] = root
        return self.parent[v], self.dx[v], self.dy[v]

    def add_basis(self, root: int, vector) -> None:
        if vector == (0, 0) or len(self.basis[root]) == 2:
            return
        if not self.basis[root]:
            self.basis[root].append(vector)
        elif self.basis[root][0][0] * vector[1] != self.basis[root][0][1] * vector[0]:
            self.basis[root].append(vector)

    def edge(self, u: int, v: int, dx: int, dy: int) -> None:
        ru, ux, uy = self.find(u)
        rv, vx, vy = self.find(v)
        wx, wy = dx + ux - vx, dy + uy - vy
        if ru == rv:
            self.add_basis(ru, (wx, wy))
            return
        if self.size[ru] < self.size[rv]:
            ru, rv, wx, wy = rv, ru, -wx, -wy
        self.parent[rv] = ru
        self.dx[rv], self.dy[rv] = wx, wy
        self.size[ru] += self.size[rv]
        for vector in self.basis[rv]:
            self.add_basis(ru, vector)


def rank_map(n: int, y_label: int, enabled: set, steps):
    uf = LiftUnionFind(n)
    for vertex in sorted(enabled):
        for dx, dy in steps:
            target = (vertex + dx + y_label * dy) % n
            if target in enabled:
                uf.edge(vertex, target, dx, dy)
    roots = {v: uf.find(v)[0] for v in enabled}
    counts = Counter(roots.values())
    components = sorted(({'size': size, 'ambient_rank': len(uf.basis[root])}
                         for root, size in counts.items()),
                        key=lambda row: (row['size'], row['ambient_rank']))
    return components, {v: len(uf.basis[root]) for v, root in roots.items()}


def reduce_phi5(polynomial):
    """Integer long division by 1+t+t^2+t^3+t^4."""
    values = list(polynomial) + [0] * max(0, 4 - len(polynomial))
    for degree in range(len(values) - 1, 3, -1):
        top = values[degree]
        for offset in range(5):
            values[degree - offset] -= top
    return values[:4]


def verify_witness(row: dict, hand: str, a: int, b: int) -> None:
    n = a * a + b * b
    require(row['hand'] == hand and row['gaussian_period'] == [a, b], 'wrong hand geometry')
    require(row['order'] == n == 505, 'wrong child order')
    require(row['period_matrix_columns'] == [[a, b], [-b, a]], 'wrong period columns')
    require(gcd(b, n) == 1, 'axial verifier requires invertible b')
    y_label = (-a * pow(b, -1, n)) % n
    label = lambda x, y: (x + y_label * y) % n
    require(label(a, b) == label(-b, a) == 0, 'period quotient mismatch')
    sign = 1 if b > 0 else -1
    expected_lifts = [[x, 0] for x in range(a)] + [[a, sign * y] for y in range(abs(b))]
    require(row['occupied_lifts'] == expected_lifts, 'occupied staircase changed')
    black = {label(*point) for point in row['occupied_lifts']}
    white = set(range(n)) - black
    require(row['occupied_count'] == len(black) == a + abs(b), 'wrong occupied count')
    require(row['vacant_count'] == len(white), 'wrong vacant count')
    bc, br = rank_map(n, y_label, black, NN_FORWARD)
    wc, wr = rank_map(n, y_label, white, MATCHING_FORWARD)
    require(row['black_NN_components'] == bc, 'black topology mismatch')
    require(row['white_matching_components'] == wc, 'white topology mismatch')
    require(bc == [{'size': len(black), 'ambient_rank': 1}], 'black witness not rank one')
    require(wc == [{'size': len(white), 'ambient_rank': 1}], 'white witness not rank one')
    points = [label(j + 10 * f, f) for j in range(101) for f in range(5)]
    require(len(set(points)) == row['parent_fiber_label_count'] == 505, 'fiber labels collapse')
    require(row['charged_field_basis'] == '5*F_r=sum_{j=0}^3 coefficients[j]*zeta5^j',
            'field normalization changed')
    scalar = {v: (int(br[v] == 1) if v in black else -int(wr[v] == 1)) for v in range(n)}
    require(set(row['charges']) == {'1', '2', '3', '4'}, 'charge set changed')
    for r in range(1, 5):
        field = []
        for j in range(101):
            polynomial = [0] * 5
            for f in range(5):
                polynomial[(-r * f) % 5] += scalar[label(j + 10 * f, f)]
            field.append(reduce_phi5(polynomial))
        record = row['charges'][str(r)]
        zeros = [i for i, value in enumerate(field) if value == [0] * 4]
        nonzeros = [i for i, value in enumerate(field) if value != [0] * 4]
        require(record['zero_parents'] == len(zeros) == 72, 'zero support mismatch')
        require(record['nonzero_parents'] == len(nonzeros) == 29, 'nonzero support mismatch')
        require(record['zero_parent'] == zeros[0], 'zero witness index changed')
        require(record['nonzero_parent'] == nonzeros[0], 'nonzero witness index changed')
        require(record['zero_coefficients'] == field[zeros[0]], 'zero witness changed')
        require(record['nonzero_coefficients'] == field[nonzeros[0]], 'nonzero witness changed')
        serial = json.dumps(field, sort_keys=True, separators=(',', ':')).encode('utf-8')
        require(record['all_coefficients_sha256'] == hashlib.sha256(serial).hexdigest(),
                'charged coefficient hash mismatch')
        require(record['all_fifth_root_gauges_preserve_zero_support'] is True, 'missing gauge check')
        for value in field:
            for t in range(5):
                require(any(value) == any(reduce_phi5([0] * t + value)), 'gauge changed support')


def verify_certificate(data: dict) -> dict:
    require(data['schema'] == SCHEMA and data['issue'] == 406, 'wrong certificate type')
    require(data['data_class'] == 'exact finite witness; no Monte Carlo', 'wrong data class')
    require(data['source_semantics_commit'] == '33c557b9aebed1bc9c07019b9cd5cee6c04be947',
            'source contract changed')
    require(data['parent'] == {'gaussian_period': [10, 1], 'order': 101,
                               'label': '(x-10*y) mod 101', 'fiber_step': [10, 1]},
            'parent geometry changed')
    require(data['cyclotomic_degree'] == {'phi_5': 4, 'phi_101': 100,
                                         'phi_505': 400, 'extension_degree': 100},
            'cyclotomic degree data changed')
    require(all(101 % q for q in range(2, 11)), '101 must be prime')
    require(sum(gcd(k, 505) == 1 for k in range(505)) == 400, 'totient identity failed')
    require(len(data['witnesses']) == 2, 'two witness configurations required')
    for row, specification in zip(data['witnesses'], [('plus', 19, 12), ('minus', 21, -8)]):
        verify_witness(row, *specification)
    require(len(data['endpoint_aliases']) == 3, 'three radius audits required')
    for row, radius, distinct in zip(data['endpoint_aliases'], (4, 5, 6), (41, 61, 77)):
        fibers = {}
        for x in range(-radius, radius + 1):
            for y in range(-radius, radius + 1):
                if abs(x) + abs(y) <= radius:
                    fibers.setdefault((x - 10 * y) % 101, []).append([x, y])
        expected = [{'parent_label': key, 'displacements': values}
                    for key, values in sorted(fibers.items()) if len(values) > 1]
        require(row == {'radius': radius, 'displacement_labels': 1 + 2 * radius * (radius + 1),
                        'distinct_parent_vertices': distinct, 'repeated_classes': expected},
                'endpoint alias table changed')
    control = data['projection_leakage_control']
    require(set(control['checks']) == {'microscopic_commutation', 'orthogonal_projector',
                                      'compression_commutator_nonzero', 'leakage_identity'},
            'projection control checks changed')
    require(all(value is True for value in control['checks'].values()), 'projection check failed')
    require(control['commutator'] == [['0', '-4/9', '4/9'], ['4/9', '0', '-4/9'],
                                     ['-4/9', '4/9', '0']], 'wrong compression commutator')
    conclusion = data['conclusion']
    require(conclusion['complete_spatial_rank_lower_bound'] == 100 and
            conclusion['positive_nonzero_spatial_frequencies_per_hand_charge'] == 100,
            'wrong exact spatial rank bound')
    require(conclusion['ensemble'] == '0<p<1 independent Bernoulli and uniform parent-anchor target',
            'ensemble assumption changed')
    require(conclusion['boundary'] == 'Not a statistical effective-rank or physical-field count; no production spectrum estimated.',
            'claim boundary changed')
    return {'verified': True, 'independent_algorithm': 'cyclic quotient plus weighted union-find',
            'exact_spatial_rank_lower_bound': 100, 'witnesses': 2,
            'production_data_read': False}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate', type=Path)
    args = parser.parse_args()
    try:
        result = verify_certificate(json.loads(args.certificate.read_text(encoding='utf-8')))
    except (ValueError, KeyError, TypeError, IndexError) as error:
        parser.exit(1, 'certificate rejected: ' + str(error) + '\n')
    print(json.dumps(result, sort_keys=True, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
