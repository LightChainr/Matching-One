#!/usr/bin/env python3
"""Exact trace-equivalence versus branching/Markov-closure certificates.

Reuses the unchanged fee33287 Gaussian rank oracle. This is a bounded exact
analysis, not a production engine. Independent N16 verification is supplied
in verify_p429_n16.cpp. Python 3.10+; standard library only.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from math import comb, factorial, gcd
from pathlib import Path
import json
from typing import Any, Sequence

try:
    from .p334_birth_age_collision_review_20260830 import enumerate_states
except ImportError:
    from p334_birth_age_collision_review_20260830 import enumerate_states

State = tuple[int, tuple[int, int] | None]
SCHEMA = 'matching-one/p429-branching-continuation/v1'
GEOMETRIES = ((2, 1), (3, 0), (3, 1), (3, 2), (4, 0), (4, 1))
TARGET_COUNTS = (1, 7, 18, 20, 8, 0, 0, 0, 0)
N16_COORDINATES = (
    ((0, 0), (1, 0), (2, 0), (3, 0), (1, 1), (3, 1), (0, 3), (1, 3)),
    ((0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (3, 1), (0, 3)),
)


class ContinuationKernel:
    """Killed uniform-insertion kernel over an explicit, validated rank cache."""

    def __init__(self, n: int, states: Sequence[State]) -> None:
        if type(n) is not int or not 1 <= n <= 20 or len(states) != 1 << n:
            raise ValueError('require a complete Boolean rank cache, 1 <= N <= 20')
        self.n, self.states = n, tuple(states)
        self.full = (1 << n) - 1
        for mask, (rank, line) in enumerate(self.states):
            if type(rank) is not int or rank not in (0, 1, 2):
                raise ValueError('rank must be 0, 1 or 2')
            if rank == 1:
                if (not isinstance(line, tuple) or len(line) != 2
                        or any(type(x) is not int for x in line)
                        or gcd(*map(abs, line)) != 1
                        or line[0] < 0 or (line[0] == 0 and line[1] <= 0)):
                    raise ValueError('rank-one line must be canonical and primitive')
            elif line is not None:
                raise ValueError('only rank-one states have a plateau line')
            for child in self.children(mask):
                after, child_line = self.states[child]
                if after < rank or (rank == after == 1 and child_line != line):
                    raise ValueError('rank monotonicity or plateau-line nesting fails')
        if self.states[0][0] != 0 or self.states[-1][0] != 2:
            raise ValueError('empty/full endpoint ranks must be zero/two')

        self.counts: dict[int, tuple[int, ...]] = {}
        self.classes: dict[int, int] = {}
        self.groups: dict[tuple, list[int]] = defaultdict(list)
        class_ids: dict[tuple, int] = {}
        ordered = sorted((s for s, st in enumerate(self.states) if st[0] == 1),
                         key=lambda s: (-s.bit_count(), s))
        for s in ordered:
            k, d = s.bit_count(), n - s.bit_count()
            children = self.children(s)
            safe = [t for t in children if self.states[t][0] == 1]
            counts = [1]
            for m in range(1, d + 1):
                total = sum(self.counts[t][m - 1] for t in safe)
                if total % m:
                    raise ArithmeticError('safe-subset double count is not integral')
                counts.append(total // m)
            self.counts[s] = tuple(counts)
            # -1 denotes one absorbing cemetery state, excluded from counts.
            histogram = tuple(sorted(Counter(self.classes.get(t, -1)
                                             for t in children).items()))
            key = (k, self.states[s][1], histogram)
            if key not in class_ids:
                class_ids[key] = len(class_ids)
            self.classes[s] = class_ids[key]
            self.groups[self.signature(s)].append(s)
        self.class_count = len(class_ids)

    def children(self, mask: int) -> tuple[int, ...]:
        if type(mask) is not int or not 0 <= mask <= self.full:
            raise ValueError('mask outside the Boolean cube')
        return tuple(mask | (1 << v) for v in range(self.n) if not mask & (1 << v))

    def signature(self, mask: int) -> tuple:
        if mask not in self.counts:
            raise ValueError('a survival signature requires rank one')
        return mask.bit_count(), self.states[mask][1], self.counts[mask]

    def survival(self, mask: int, horizon: int) -> Fraction:
        if mask not in self.counts:
            raise ValueError('survival requires a rank-one checkpoint')
        d = self.n - mask.bit_count()
        if type(horizon) is not int or not 0 <= horizon <= d:
            raise ValueError('invalid survival horizon')
        return Fraction(self.counts[mask][horizon], comb(d, horizon))

    def direct_counts(self, mask: int) -> tuple[int, ...]:
        """Independent safe-subset enumeration, not the recurrence."""
        if mask not in self.counts:
            raise ValueError('direct counts require rank one')
        d = self.n - mask.bit_count()
        out = [0] * (d + 1)
        rest = self.full ^ mask
        added = rest
        while True:
            if self.states[mask | added][0] == 1:
                out[added.bit_count()] += 1
            if not added:
                return tuple(out)
            added = (added - 1) & rest

    def fork_after_one(self, mask: int, horizon: int = 1) -> tuple[Fraction, Fraction]:
        """Mean single-branch survival and joint two-branch survival."""
        self.signature(mask)
        d = self.n - mask.bit_count()
        if type(horizon) is not int or not 0 <= horizon <= d - 1:
            raise ValueError('branch horizon exceeds post-update vacancies')
        values = [self.survival(t, horizon) if t in self.counts else Fraction(0)
                  for t in self.children(mask)]
        return sum(values, Fraction(0)) / d, sum((v*v for v in values), Fraction(0)) / d

    def next_exit_histogram(self, mask: int) -> dict[str, int]:
        d = self.n - mask.bit_count()
        if d < 2:
            raise ValueError('requires at least two remaining sites')
        out: Counter[str] = Counter()
        for t in self.children(mask):
            out['absorbed' if t not in self.counts else str(d-1-self.counts[t][1])] += 1
        return dict(sorted(out.items()))

    def summary(self, a: int, b: int) -> dict[str, Any]:
        split = sum(len({self.classes[s] for s in group}) > 1
                    for group in self.groups.values())
        return dict(gaussian=[a, b], N=self.n, rank_one_states=len(self.counts),
                    survival_classes=len(self.groups), strong_markov_classes=self.class_count,
                    split_survival_classes=split, cemetery_included=False)


def inherited_coordinates(a: int, b: int) -> tuple[tuple[int, int], ...]:
    """Site registry of the unchanged oracle; no topology calculation here."""
    n = a*a + b*b
    def key(x: int, y: int) -> tuple[int, int]:
        return (a*x+b*y) % n, (-b*x+a*y) % n
    reps = [(0, 0)]
    seen = {key(0, 0)}
    for x, y in reps:
        for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            point = key(x+dx, y+dy)
            if point not in seen:
                seen.add(point)
                reps.append((x+dx, y+dy))
    if len(reps) != n:
        raise ArithmeticError('wrong quotient cardinality')
    return tuple(reps)


def n16_mask(coordinates: Sequence[tuple[int, int]]) -> int:
    lookup = {(x % 4, y % 4): j for j, (x, y) in enumerate(inherited_coordinates(4, 0))}
    canonical = [(x % 4, y % 4) for x, y in coordinates]
    if len(set(canonical)) != len(canonical):
        raise ValueError('duplicate physical site')
    return sum(1 << lookup[p] for p in canonical)


def witness_record(kernel: ContinuationKernel, mask: int) -> dict[str, Any]:
    k, line, counts = kernel.signature(mask)
    if kernel.direct_counts(mask) != counts:
        raise ArithmeticError('direct subsets disagree with recurrence')
    mean, branch = kernel.fork_after_one(mask)
    return dict(k=k, line=list(line), counts=list(counts),
                survival=[str(kernel.survival(mask, j)) for j in range(len(counts))],
                next_exit_choices=kernel.next_exit_histogram(mask),
                unbranched_two_step=str(mean), shared_step_fork=str(branch),
                successor_prediction_variance=str(branch-mean*mean))


def n16_history_certificate(kernel: ContinuationKernel) -> dict[str, Any]:
    key = (8, (1, 0), TARGET_COUNTS)
    group = set(kernel.groups[key])
    prefixes: dict[int, Counter[int]] = {}
    totals: Counter[int] = Counter()
    weighted: dict[int, Counter[str]] = defaultdict(Counter)
    for s, (rank, _) in enumerate(kernel.states):
        k = s.bit_count()
        if rank != 1 or k > 8:
            continue
        counts: Counter[int] = Counter()
        for v in range(kernel.n):
            if not s & (1 << v):
                continue
            before = s ^ (1 << v)
            previous_rank = kernel.states[before][0]
            if previous_rank == 0:
                counts[k] += factorial(k-1)
            elif previous_rank == 1:
                counts.update(prefixes[before])
        if sum(counts.values()) != factorial(k):
            raise ArithmeticError('birth-prefix weights fail factorial normalization')
        prefixes[s] = counts
        if s in group:
            histogram = kernel.next_exit_histogram(s)
            for birth, count in counts.items():
                totals[birth] += count
                for next_x, choices in histogram.items():
                    weighted[birth][next_x] += count*choices
    rows = []
    for birth, total in sorted(totals.items()):
        weights = dict(sorted(weighted[birth].items()))
        rows.append(dict(K1=birth, prefixes=total, next_exit_weights=weights,
                         next_exit_probabilities={x: str(Fraction(w, 8*total))
                                                  for x, w in weights.items()}))
    return dict(stratum_states=len(group), selected_prefixes=sum(totals.values()),
                line=[1, 0], k=8, full_counts=list(TARGET_COUNTS), rows=rows,
                event='rank(S9)=1 and x(S9)=3', gap_K1_8_minus_4='1/66')


def build_artifact() -> dict[str, Any]:
    summaries = []
    witness = history = None
    for a, b in GEOMETRIES:
        n, states = enumerate_states(a, b)
        kernel = ContinuationKernel(n, states)
        summaries.append(kernel.summary(a, b))
        if (a, b) == (4, 0):
            masks = [n16_mask(p) for p in N16_COORDINATES]
            witness = [dict(name=name, coordinates=[list(point) for point in coords],
                            **witness_record(kernel, mask))
                       for name, coords, mask in zip(('A', 'B'), N16_COORDINATES, masks)]
            history = n16_history_certificate(kernel)
    artifact = dict(schema=SCHEMA, status='exact_finite_volume', issue=429,
                    inherited_oracle_commit='fee33287cf4830e07ccef6177f43034add02256e',
                    inherited_oracle_blob='62e06795fdfa91a956aedd62b7344e84aa5efc5c',
                    census=summaries, n16_witness=witness, n16_history=history,
                    claim_boundary='No all-HNF minimality, asymptotic complexity, or CFT field claim.')
    verify_expected_artifact(artifact)
    return artifact


def verify_expected_artifact(artifact: dict[str, Any]) -> None:
    expected = [(5, 10, 2, 2, 0), (9, 162, 10, 10, 0), (10, 310, 16, 16, 0),
                (13, 2340, 62, 62, 0), (16, 19932, 210, 214, 4), (17, 38896, 346, 390, 42)]
    keys = ('N', 'rank_one_states', 'survival_classes', 'strong_markov_classes', 'split_survival_classes')
    if artifact.get('schema') != SCHEMA or [tuple(r[k] for k in keys) for r in artifact['census']] != expected:
        raise ValueError('census certificate mismatch')
    a, b = artifact['n16_witness']
    for record in (a, b):
        if tuple(record['counts']) != TARGET_COUNTS or record['unbranched_two_step'] != '9/14':
            raise ValueError('full-survival witness mismatch')
    if (Fraction(a['shared_step_fork']) != Fraction(95, 196)
            or Fraction(b['shared_step_fork']) != Fraction(93, 196)):
        raise ValueError('fork certificate mismatch')
    history = artifact['n16_history']
    expected_totals = {4: 110592, 5: 442368, 6: 1198080, 7: 2442240, 8: 3548160}
    expected_event = {4: '1/6', 5: '1/6', 6: '2/13', 7: '8/53', 8: '2/11'}
    if history['stratum_states'] != 192 or history['selected_prefixes'] != 7741440:
        raise ValueError('history stratum mismatch')
    if {r['K1']: r['prefixes'] for r in history['rows']} != expected_totals:
        raise ValueError('prefix counts mismatch')
    for r in history['rows']:
        if r['next_exit_probabilities']['3'] != expected_event[r['K1']]:
            raise ValueError('history-dependent transition mismatch')
        for x, w in r['next_exit_weights'].items():
            if str(Fraction(w, 8*r['prefixes'])) != r['next_exit_probabilities'][x]:
                raise ValueError('inconsistent transition normalization')
        if sum(r['next_exit_weights'].values()) != 8*r['prefixes']:
            raise ValueError('incomplete transition weights')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    text = json.dumps(build_artifact(), indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding='utf-8')
    else:
        print(text, end='')


if __name__ == '__main__':
    main()
