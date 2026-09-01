#!/usr/bin/env python3
"""Canonical Kreg: Bell8 closure via equality diagrams, not colour coarsenings.

Uses only Python's standard library. No colour assignments, coarsening sums,
falling factorials or repository computation modules occur in this algorithm.
"""
from __future__ import annotations

import csv
from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import time

HERE = Path(__file__).resolve().parent
SOURCE_COMMIT = "2ba8863f75e0ced211b7b5442e8cddbe2fbd3deb"
SOURCE_NOTE = "notes/local-pair-two-insertion-algebra.md"


def partitions(n):
    """Restricted growth strings, lexicographic and canonical."""
    word = [0] * n

    def visit(i, largest):
        if i == n:
            yield tuple(word)
            return
        for label in range(largest + 2):
            word[i] = label
            yield from visit(i + 1, max(largest, label))

    yield from visit(1, 0)


def canonical(labels):
    seen = {}
    return tuple(seen.setdefault(x, len(seen)) for x in labels)


def edges_for_partition(partition, offset=0):
    first = {}
    edges = []
    for vertex, label in enumerate(partition):
        vertex += offset
        if label in first:
            edges.append((first[label], vertex))
        else:
            first[label] = vertex
    return tuple(edges)


def join_partition(n, edges):
    parent = list(range(n))

    def find(x):
        while x != parent[x]:
            x = parent[x]
        return x

    for x, y in edges:
        x, y = find(x), find(y)
        if x != y:
            parent[y] = x
    return canonical(tuple(find(x) for x in range(n)))


def closed_loops(exterior_parent, first_edges, second_edges):
    """Sew equality diagrams to the exterior: one factor Q per component."""
    parent = list(exterior_parent)

    def find(x):
        while parent[x] != x:
            x = parent[x]
        return x

    for x, y in first_edges + second_edges:
        x, y = find(x), find(y)
        if x != y:
            parent[y] = x
    return sum(parent[x] == x for x in range(len(parent)))


def equality_diagrams():
    # Areg = (1/2)(1-d_ab)(1-d_cd) *
    #  [d_ac*d_bd+d_ad*d_bc -(d_ac+d_ad+d_bc+d_bd)/(Q-2)
    #                            +4/(Q*(Q-2))].
    # Kreg = (Areg + R Areg)/2. Each coefficient is stored as
    #  (u + v/(Q-2) + w/(Q*(Q-2)))/4.
    # This expands products of delta tensors directly; equalities impose
    # connectivity only and never demand distinct colours across blocks.
    bracket = [
        (((0, 2), (1, 3)), (1, 0, 0)),
        (((0, 3), (1, 2)), (1, 0, 0)),
        (((0, 2),), (0, -1, 0)),
        (((0, 3),), (0, -1, 0)),
        (((1, 2),), (0, -1, 0)),
        (((1, 3),), (0, -1, 0)),
        ((), (0, 0, 4)),
    ]
    terms = defaultdict(lambda: [0, 0, 0])
    raw_terms = 0
    for rotation in (0, 1):
        for edge_list, basis in bracket:
            for ab in (0, 1):
                for cd in (0, 1):
                    raw_terms += 1
                    edges = list(edge_list)
                    if ab:
                        edges.append((0, 1))
                    if cd:
                        edges.append((2, 3))
                    rotated = tuple(((x + rotation) % 4, (y + rotation) % 4)
                                    for x, y in edges)
                    pi = join_partition(4, rotated)
                    sign = -1 if (ab + cd) % 2 else 1
                    for i, value in enumerate(basis):
                        terms[pi][i] += sign * value
    result = []
    for pi, (u, v, w) in sorted(terms.items()):
        if not (u or v or w):
            continue
        # Exact Taylor jets at Q=1. d[1/(Q*(Q-2))]/dQ vanishes there.
        c4, dc4 = u - v - w, -v
        result.append(dict(partition=pi, numerator_basis=(u, v, w),
                           c4=c4, dc4=dc4,
                           edges=edges_for_partition(pi)))
    assert sum(x["c4"] for x in result) == 0
    assert sum(x["dc4"] for x in result) == 0
    return result, raw_terms


def exterior_parent(pi):
    first = {}
    return tuple(first.setdefault(label, i) for i, label in enumerate(pi))


def two_point(pi, diagrams):
    b = max(pi) + 1
    parent = exterior_parent(pi)
    b0_16 = 0
    da_16 = 0
    # For each ordered pair of equality diagrams,
    # B(Q) += c_i(Q)*c_j(Q)*Q**(loops - |pi|).
    # Differentiate *all three factors*, including the normalization.
    for x in diagrams:
        for y in diagrams:
            loops = closed_loops(parent, x["edges"],
                                 tuple((i + 4, j + 4) for i, j in y["edges"]))
            product = x["c4"] * y["c4"]
            b0_16 += product
            da_16 += (x["dc4"] * y["c4"] + x["c4"] * y["dc4"]
                      + (loops - b) * product)
    return b0_16, da_16


def single_point(pi, diagrams):
    b = max(pi) + 1
    parent = exterior_parent(pi)
    value4 = deriv4 = 0
    for term in diagrams:
        loops = closed_loops(parent, term["edges"], ())
        value4 += term["c4"]
        deriv4 += term["dc4"] + (loops - b) * term["c4"]
    return value4, deriv4


def expected_single(pi):
    blocks = Counter(pi)
    shape = sorted(blocks.values())
    if shape == [1, 1, 1, 1]:
        return 4
    if shape == [1, 1, 2]:
        pair = [i for i, label in enumerate(pi) if blocks[label] == 2]
        return 4 if abs(pair[1] - pair[0]) == 2 else 2
    if shape == [2, 2]:
        return -2 if pi[0] == pi[2] else -1
    return 0


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    started = time.process_time()
    diagrams, raw_terms = equality_diagrams()
    singles = []
    for pi in partitions(4):
        value4, deriv4 = single_point(pi, diagrams)
        assert value4 == 0
        assert deriv4 == expected_single(pi), (pi, deriv4)
        singles.append(dict(partition="".join(map(str, pi)), value4=value4,
                            derivative4=deriv4))
    by_shared = defaultdict(Counter)
    rows = []
    table = {}
    for pi in partitions(8):
        value16, deriv16 = two_point(pi, diagrams)
        shared = len(set(pi[:4]) & set(pi[4:]))
        assert value16 == 0
        if shared <= 1:
            assert deriv16 == 0, (pi, deriv16)
        by_shared[shared][deriv16] += 1
        table[pi] = deriv16
        rows.append(dict(partition="".join(map(str, pi)),
                         exterior_blocks=max(pi) + 1,
                         shared_components=shared,
                         value_numerator16=value16,
                         derivative_numerator16=deriv16,
                         derivative=str(Fraction(deriv16, 16))))
    assert len(rows) == 4140
    witness = (0, 1, 2, 3, 0, 3, 2, 1)
    assert table[witness] == 26
    # C4 invariance at each site and site exchange are direct symmetries.
    for pi, derivative in table.items():
        transformed = (pi[1], pi[2], pi[3], pi[0]) + pi[4:]
        assert table[canonical(transformed)] == derivative
        transformed = pi[:4] + (pi[5], pi[6], pi[7], pi[4])
        assert table[canonical(transformed)] == derivative
        assert table[canonical(pi[4:] + pi[:4])] == derivative
    csv_path = HERE / "BELL8_DIAGRAM_RESULTS.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    coefficients = [{k: value for k, value in item.items() if k != "edges"}
                    for item in diagrams]
    (HERE / "EQUALITY_DIAGRAMS.json").write_text(json.dumps({
        "source_commit": SOURCE_COMMIT, "source_path": SOURCE_NOTE,
        "source_formula": "Section 2 equation (5), then average with quarter turn",
        "coefficient_basis": "c(Q)=(u+v/(Q-2)+w/(Q*(Q-2)))/4",
        "c4": "4*c(1)=u-v-w", "dc4": "4*c'(1)=-v",
        "raw_delta_product_terms": raw_terms, "terms": coefficients,
        "single_closure_controls": singles,
    }, indent=2) + "\n")
    spectrum = Counter(table.values())
    summary = {
        "source_commit": SOURCE_COMMIT,
        "source_path": SOURCE_NOTE,
        "algorithm": "equality-diagram product expansion; DSU sewing; exact first jets",
        "colour_coarsening_used": False,
        "partition_count": len(rows),
        "raw_delta_product_terms": raw_terms,
        "distinct_equality_diagrams": len(diagrams),
        "all_B_at_Q1_zero": True,
        "shared_at_most_one_all_zero": True,
        "positive": sum(count for value, count in spectrum.items() if value > 0),
        "zero": spectrum[0],
        "negative": sum(count for value, count in spectrum.items() if value < 0),
        "minimum": str(Fraction(min(spectrum), 16)),
        "maximum": str(Fraction(max(spectrum), 16)),
        "minimum_examples": [row for row in rows if row["derivative_numerator16"] == min(spectrum)][:8],
        "maximum_examples": [row for row in rows if row["derivative_numerator16"] == max(spectrum)][:8],
        "spectrum_numerator_over_16": dict(sorted(spectrum.items())),
        "by_shared_components": {key: dict(sorted(value.items()))
                                 for key, value in sorted(by_shared.items())},
        "four_path_witness": {"partition": "".join(map(str, witness)),
                              "derivative_numerator16": table[witness],
                              "derivative": "13/8"},
        "single_closure_controls_all_pass": True,
        "independent_C4_and_site_exchange_controls_all_pass": True,
        "script_sha256": sha(Path(__file__)),
        "results_sha256": sha(csv_path),
        "cpu_seconds": time.process_time() - started,
        "new_lattice_populations": 0,
    }
    (HERE / "SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
