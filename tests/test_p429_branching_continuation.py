from __future__ import annotations
import copy
import hashlib
import itertools
import json
from pathlib import Path
from collections import Counter
from fractions import Fraction
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import p429_branching_continuation as p


class BranchingCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.artifact = p.build_artifact()
        n, st = p.enumerate_states(4, 0)
        cls.kernel = p.ContinuationKernel(n, st)
        cls.a, cls.b = [p.n16_mask(c) for c in p.N16_COORDINATES]

    def test_inherited_oracle_is_original_blob(self):
        content = (ROOT / 'scripts/p334_birth_age_collision_review_20260830.py').read_bytes()
        git_blob = hashlib.sha1(b'blob ' + str(len(content)).encode() + b'\0' + content).hexdigest()
        self.assertEqual(git_blob, '62e06795fdfa91a956aedd62b7344e84aa5efc5c')

    def test_census_matches_frozen_certificate(self):
        p.verify_expected_artifact(self.artifact)
        stored = json.loads((ROOT / 'results/p429-branching-continuation/exact.json').read_text())
        self.assertEqual(self.artifact, stored)

    def test_complete_survival_and_direct_subsets(self):
        for s in (self.a, self.b):
            self.assertEqual(self.kernel.counts[s], p.TARGET_COUNTS)
            self.assertEqual(self.kernel.direct_counts(s), p.TARGET_COUNTS)
        self.assertEqual(self.kernel.signature(self.a), self.kernel.signature(self.b))

    def test_all_future_permutations_agree(self):
        expected = Counter({1: 5040, 2: 9360, 3: 11520, 4: 9792, 5: 4608})
        for s in (self.a, self.b):
            vacant = [v for v in range(16) if not s & (1 << v)]
            exits = Counter()
            for order in itertools.permutations(vacant):
                t = s
                for m, v in enumerate(order, 1):
                    t |= 1 << v
                    if self.kernel.states[t][0] == 2:
                        exits[m] += 1
                        break
            self.assertEqual(exits, expected)

    def test_shared_update_fork_by_392_choices(self):
        for s, expected in ((self.a, 190), (self.b, 186)):
            success = 0
            for v in range(16):
                if s & (1 << v):
                    continue
                t = s | (1 << v)
                for u in range(16):
                    if t & (1 << u):
                        continue
                    for w in range(16):
                        if t & (1 << w):
                            continue
                        success += (self.kernel.states[t | (1 << u)][0] == 1
                                    and self.kernel.states[t | (1 << w)][0] == 1)
            self.assertEqual(success, expected)
            mean, fork = self.kernel.fork_after_one(s)
            self.assertEqual(mean, Fraction(9, 14))
            self.assertEqual(fork, Fraction(expected, 392))
        self.assertEqual(self.kernel.fork_after_one(self.a)[1]
                         - self.kernel.fork_after_one(self.b)[1], Fraction(1, 98))

    def test_root_forks_cannot_distinguish(self):
        for m in range(9):
            for n in range(9):
                self.assertEqual(self.kernel.survival(self.a, m)*self.kernel.survival(self.a, n),
                                 self.kernel.survival(self.b, m)*self.kernel.survival(self.b, n))

    def test_successor_covariance_is_exact(self):
        self.assertEqual(self.kernel.fork_after_one(self.a)[1]-Fraction(9, 14)**2, Fraction(1, 14))
        self.assertEqual(self.kernel.fork_after_one(self.b)[1]-Fraction(9, 14)**2, Fraction(3, 49))
        self.assertNotEqual(self.kernel.classes[self.a], self.kernel.classes[self.b])

    def test_markov_partition_has_constant_successor_laws(self):
        seen = {}
        for s, label in self.kernel.classes.items():
            histogram = tuple(sorted(Counter(self.kernel.classes.get(t, -1)
                                             for t in self.kernel.children(s)).items()))
            value = (s.bit_count(), self.kernel.states[s][1], histogram)
            if label in seen:
                self.assertEqual(value, seen[label])
            seen[label] = value
        self.assertEqual(len(seen), 214)

    def test_actual_history_gap(self):
        hist = self.artifact['n16_history']
        rows = {r['K1']: r for r in hist['rows']}
        self.assertEqual(sum(r['prefixes'] for r in rows.values()), 192*40320)
        self.assertEqual(Fraction(rows[8]['next_exit_probabilities']['3'])
                         - Fraction(rows[4]['next_exit_probabilities']['3']), Fraction(1, 66))
        for r in rows.values():
            self.assertEqual(r['next_exit_probabilities']['absorbed'], '1/8')

    def test_physical_translation_invariance(self):
        for points, mask in zip(p.N16_COORDINATES, (self.a, self.b)):
            for dx, dy in itertools.product(range(4), repeat=2):
                moved = p.n16_mask([(x+dx, y+dy) for x, y in points])
                self.assertEqual(self.kernel.signature(mask), self.kernel.signature(moved))
                self.assertEqual(self.kernel.classes[mask], self.kernel.classes[moved])

    def test_corrupted_artifact_rejected(self):
        broken = copy.deepcopy(self.artifact)
        broken['n16_witness'][0]['shared_step_fork'] = '93/196'
        with self.assertRaises(ValueError):
            p.verify_expected_artifact(broken)
        broken = copy.deepcopy(self.artifact)
        broken['n16_history']['rows'][0]['next_exit_weights']['3'] += 1
        with self.assertRaises(ValueError):
            p.verify_expected_artifact(broken)

    def test_invalid_cache_and_horizon_rejected(self):
        n, states = p.enumerate_states(2, 1)
        with self.assertRaises(ValueError):
            p.ContinuationKernel(n, states[:-1])
        bad = list(states)
        bad[0] = (1, (1, 0))
        with self.assertRaises(ValueError):
            p.ContinuationKernel(n, bad)
        bad = list(states)
        first = next(i for i, st in enumerate(states) if st[0] == 1)
        bad[first] = (1, (2, 0))
        with self.assertRaises(ValueError):
            p.ContinuationKernel(n, bad)
        with self.assertRaises(ValueError):
            self.kernel.survival(self.a, 9)
        with self.assertRaises(ValueError):
            self.kernel.fork_after_one(self.a, 8)
        with self.assertRaises(ValueError):
            self.kernel.signature(0)
        with self.assertRaises(ValueError):
            p.n16_mask([(0, 0), (4, 0)])

    @unittest.skipUnless(shutil.which('g++'), 'g++ required for the independent N16 verifier')
    def test_independent_cpp_topology_and_all_selected_prefixes(self):
        with tempfile.TemporaryDirectory() as tmp:
            binary = Path(tmp) / 'verify'
            table = Path(tmp) / 'rank.txt'
            compile_result = subprocess.run(
                ['g++', '-std=c++17', '-O2', '-Wall', '-Wextra', '-Wpedantic',
                 str(ROOT / 'scripts/verify_p429_n16.cpp'), '-o', str(binary)],
                text=True, capture_output=True, check=True, timeout=30)
            self.assertEqual(compile_result.stderr, '')
            completed = subprocess.run([str(binary), str(table)], text=True,
                                       capture_output=True, check=True, timeout=30)
            report = json.loads(completed.stdout)
            registry = p.inherited_coordinates(4, 0)
            permutation = [(x % 4)+4*(y % 4) for x, y in registry]
            reference = {}
            for row in table.read_text().splitlines():
                mask, rank, x, y = map(int, row.split())
                reference[mask] = rank, (x, y) if rank == 1 else None
            self.assertEqual(len(reference), 65536)
            for s, state in enumerate(self.kernel.states):
                rowmajor = sum(1 << permutation[v] for v in range(16) if s & (1 << v))
                self.assertEqual(state, reference[rowmajor])
            self.assertEqual(report['rank_one_states'], 19932)
            self.assertEqual(report['stratum_states'], 192)
            self.assertEqual([r['fork_success_of_392'] for r in report['witnesses']], [190, 186])
            by_birth = {r['K1']: r for r in report['cohorts']}
            for r in self.artifact['n16_history']['rows']:
                independent = by_birth[r['K1']]
                self.assertEqual(independent['prefixes'], r['prefixes'])
                self.assertEqual(independent['next_three_weight'], r['next_exit_weights']['3'])


if __name__ == '__main__':
    unittest.main()
