from collections import Counter
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from noncrossing_connectivity_codec import canonical_rgs
from p398_positive_cylinder import build_model, result


class PositiveCylinderTests(unittest.TestCase):
    def test_row_is_actual_eight_bond_sum(self):
        states, transfer, _, _ = build_model()
        index = {state: i for i, state in enumerate(states)}
        for col, state in enumerate(states):
            counts = Counter()
            for mask in range(256):
                parent = list(range(8))

                def find(x):
                    while parent[x] != x:
                        x = parent[x]
                    return x

                def union(a, b):
                    parent[find(a)] = find(b)

                for a in range(4):
                    for b in range(a):
                        if state[a] == state[b]:
                            union(a, b)
                    if mask & (1 << a):
                        union(a, a+4)
                    if mask & (1 << (a+4)):
                        union(a+4, (a+1) % 4+4)
                following = canonical_rgs(find(a+4) for a in range(4))
                counts[index[following]] += 1
            self.assertEqual([F(counts[row], 256) for row in range(14)],
                             [transfer[row][col] for row in range(14)])

    def test_existing_marks_are_physical_frontier_functions(self):
        states, _, _, readouts = build_model()
        for state, row in zip(states, readouts):
            adjacent = [int(state[i] == state[(i+1) % 4]) for i in range(4)]
            singleton = [int(state.count(state[i]) == 1) for i in range(4)]
            self.assertEqual(row, [adjacent[0]-adjacent[2], adjacent[1]-adjacent[3],
                                   singleton[0]-singleton[2], singleton[1]-singleton[3]])

    def test_two_propagating_modes_and_locked_artifact(self):
        value = result()
        self.assertEqual(value["characteristic_polynomial_high_to_low"], [1, "-3/32", "1/1024"])
        self.assertEqual(value["d1_cross_asymmetry_C_AL_minus_conjugate_C_LA"], ["1/1393", "-1/1393"])
        self.assertTrue(all(row["rank_over_C"] == 2 for row in value["correlations"][1:]))
        self.assertEqual(value, json.loads((ROOT / "results/p398-positive-cylinder/latest.json").read_text()))


if __name__ == "__main__":
    unittest.main()
