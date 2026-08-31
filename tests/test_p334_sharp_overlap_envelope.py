"""Focused checks for the new sharp extremal formula, not old graph tests."""
from collections import defaultdict
from itertools import product
from math import comb
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from p334_sharp_overlap_envelope import build, sharp_wedges


def test_sharp_bounds_against_all_tiny_bipartite_graphs():
    for L, R in ((1, 4), (2, 3), (3, 3), (3, 4)):
        by_edges = defaultdict(list)
        for bits in product((0, 1), repeat=L * R):
            rows = [sum(bits[i * R:(i + 1) * R]) for i in range(L)]
            columns = [sum(bits[i * R + j] for i in range(L)) for j in range(R)]
            by_edges[sum(bits)].append(sum(comb(d, 2) for d in rows + columns))
        for m, values in by_edges.items():
            answer = sharp_wedges(L, R, m)
            assert (answer["minimum"], answer["maximum"]) == (min(values), max(values))


def test_integer_maximizing_partitions_and_transposition():
    for L, R, m in ((14, 12, 108), (5, 29, 108), (8, 11, 55)):
        result = sharp_wedges(L, R, m)
        degrees = result["maximum_ferrers_partition"]
        assert sum(degrees) == m
        assert degrees == sorted(degrees, reverse=True)
        maximum = sum(comb(d, 2) + i * d for i, d in enumerate(degrees))
        assert maximum == result["maximum"]
        assert result["maximum"] == sharp_wedges(R, L, m)["maximum"]


def test_real_pair_has_disjoint_sharp_capacity_envelopes():
    result = build()
    A, B = result["records"][:2]
    assert (A["W2"]["minimum"], A["W2"]["maximum"]) == (796, 1046)
    assert (B["W2"]["minimum"], B["W2"]["maximum"]) == (1263, 1578)
    assert B["W2"]["minimum"] - A["W2"]["maximum"] == 217
    assert result["selected_summary"]["graphs"] == 22
    assert result["selected_summary"]["zero_width_envelopes"] == 8
