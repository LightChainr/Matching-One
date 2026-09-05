"""Focused scientific checks only; no production runner or random stream."""
import copy
import json
import sys
import unittest
from fractions import Fraction
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from p487_rank_one_cut_network import HnfSquareTorus, cut_rank_one
from p487_check_cut_network import INPUT, checkpoint, exhaustive_control


class CutNetworkTests(unittest.TestCase):
    def setUp(self):
        self.t = HnfSquareTorus(4, 0, 4)
        self.a = {i for i in range(16) if 12463 >> i & 1}
        self.b = {i for i in range(16) if 4343 >> i & 1}

    def test_exact_hnf_period_coordinates(self):
        t = HnfSquareTorus(425, 268, 1)
        self.assertEqual(t.vertex(268, 1), 0)
        self.assertEqual(t.vertex(425, 0), 0)
        self.assertEqual(t.winding(425*12-268*19, -19), (12, -19))
        with self.assertRaises(ValueError):
            t.winding(1, 0)

    def test_reject_bad_geometry(self):
        for args in ((0, 0, 4), (4, 4, 4), (4, -1, 4), (2, 0, 2)):
            with self.assertRaises(ValueError):
                HnfSquareTorus(*args)
        with self.assertRaises(TypeError):
            HnfSquareTorus(True, 0, 4)

    def test_reject_bad_sites(self):
        for sites in ({-1}, {16}, {True}, {1.5}):
            with self.assertRaises(ValueError):
                self.t.rank_bfs(sites)

    def test_rank_one_scope(self):
        for sites in (set(), set(range(16))):
            with self.assertRaises(ValueError):
                cut_rank_one(self.t, sites)

    def test_cycle_corruption_rejected(self):
        vertices, directions, w = self.t.essential_cycle(self.a)
        with self.assertRaises(ValueError):
            cut_rank_one(self.t, self.a, (vertices, directions, (0, 0)))
        with self.assertRaises(ValueError):
            cut_rank_one(self.t, self.a, (vertices+vertices[:1], directions, w))
        with self.assertRaises(ValueError):
            cut_rank_one(self.t, self.a, (vertices, tuple([9]*len(directions)), w))

    def test_reversing_cut_exchanges_sides(self):
        r = cut_rank_one(self.t, self.a)
        v, d, w = r.cycle
        reverse = (v[:1] + v[:0:-1], tuple((z+2) % 4 for z in d[::-1]), (-w[0], -w[1]))
        s = cut_rank_one(self.t, self.a, reverse)
        self.assertEqual(r.left_sites, s.right_sites)
        self.assertEqual(r.right_sites, s.left_sites)
        self.assertEqual(r.neutral_sites, s.neutral_sites)
        self.assertEqual(r.minimal_pairs(), s.minimal_pairs())
        self.assertEqual(r.minimal_triples(), s.minimal_triples())

    def test_d4_rotation_covariance_of_events(self):
        def rotate(v):
            x, y = self.t.coordinates[v]
            return self.t.vertex(-y, x)
        r = cut_rank_one(self.t, self.a)
        s = cut_rank_one(self.t, {rotate(v) for v in self.a})
        self.assertEqual({tuple(sorted(map(rotate, p))) for p in r.minimal_pairs()}, s.minimal_pairs())
        self.assertEqual({tuple(sorted(map(rotate, p))) for p in r.minimal_triples()}, s.minimal_triples())

    def test_two_terminal_full_survival_on_n16_witnesses(self):
        for A in (self.a, self.b):
            r = cut_rank_one(self.t, A)
            q = sorted(r.network.vacancies)
            counts = []
            for m in range(len(q)+1):
                count = 0
                for U in combinations(q, m):
                    exit_cut = r.network.connects(U)
                    self.assertEqual(exit_cut, self.t.rank_union_find(A | set(U)) == 2)
                    count += not exit_cut
                counts.append(count)
            self.assertEqual(counts, [1, 7, 18, 20, 8, 0, 0, 0, 0])

    def test_network_preserves_existing_branching_counterexample(self):
        values = []
        for A in (self.a, self.b):
            network = cut_rank_one(self.t, A).network
            total = Fraction(0)
            for v in network.vacancies:
                child = network.activate(v)
                if child is None:
                    continue
                alive = sum(not child.connects({w}) for w in child.vacancies)
                total += Fraction(alive, len(child.vacancies))**2
            values.append(total/len(network.vacancies))
        self.assertEqual(values, [Fraction(95, 196), Fraction(93, 196)])
        self.assertEqual(values[0]-values[1], Fraction(1, 98))

    def test_contraction_has_one_variable_per_remaining_site(self):
        network = cut_rank_one(self.t, self.a).network
        for v in network.vacancies:
            child = network.activate(v)
            if child is not None:
                self.assertEqual(child.vacancies, network.vacancies-{v})
                self.assertTrue(all(x < 0 for x in set(child.adjacency)-child.vacancies))
                self.assertFalse(child.connects())
                self.assertTrue(child.connects(child.vacancies))

    def test_reject_bad_update(self):
        network = cut_rank_one(self.t, self.a).network
        with self.assertRaises(ValueError):
            network.activate(next(iter(self.a)))
        with self.assertRaises(ValueError):
            network.connects({next(iter(self.a))})

    def test_n10_complete_continuation_and_update(self):
        result = exhaustive_control("N10", (10, 3, 1))
        self.assertEqual(result["rank_one_checkpoints"], 310)
        self.assertEqual(result["full_future_checks"], 7200)
        self.assertEqual(result["safe_update_future_checks"], 11200)
        self.assertEqual(result["failures"], 0)

    def test_saved_checkpoints_without_large_replay(self):
        rows = json.loads(INPUT.read_text())["rows"]
        results = [checkpoint(row, False) for row in rows]
        self.assertEqual([r["pair_statistics"]["wedges"] for r in results], [926, 1466])
        self.assertEqual([r["minimal_triples"] for r in results], [583, 509])
        self.assertTrue(all(r["exact_archived_edge_equality"] for r in results))
        self.assertTrue(all(r["alternative_cut_pair_and_triple_sets_equal"] for r in results))

    def test_corrupted_checkpoint_certificate_rejected(self):
        row = copy.deepcopy(json.loads(INPUT.read_text())["rows"][0])
        row["graph_source_git_blob"] = "0"*40
        with self.assertRaises(AssertionError):
            checkpoint(row, False)
        row = copy.deepcopy(json.loads(INPUT.read_text())["rows"][0])
        row["expected_minimal_triples"] += 1
        with self.assertRaises(AssertionError):
            checkpoint(row, False)


if __name__ == "__main__":
    unittest.main()
