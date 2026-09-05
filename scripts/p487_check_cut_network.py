#!/usr/bin/env python3
"""Exact bounded and saved-checkpoint checks for Issue #487; no random samples."""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from itertools import combinations
from pathlib import Path

from p487_rank_one_cut_network import HnfSquareTorus, cut_rank_one, pair_statistics

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results/p487-cut-network/inputs.json"


def git_blob(data: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()


def exhaustive_control(name: str, parameters: tuple) -> dict:
    t = HnfSquareTorus(*parameters)
    states = [{v for v in range(t.n) if mask >> v & 1} for mask in range(1 << t.n)]
    ranks = [t.rank_bfs(A) for A in states]
    if any(r != t.rank_union_find(A) for A, r in zip(states, ranks)):
        raise AssertionError("independent rank implementations disagree")
    count = checkpoints = pairs_checked = triples_checked = updates = 0
    full = (1 << t.n)-1
    for mask, (A, rank) in enumerate(zip(states, ranks)):
        if rank != 1:
            continue
        result = cut_rank_one(t, A)
        checkpoints += 1
        free = full ^ mask
        subset = free
        while True:
            if result.network.connects(states[subset]) != (ranks[mask | subset] == 2):
                raise AssertionError((name, mask, subset, "full continuation"))
            count += 1
            if subset == 0:
                break
            subset = (subset-1) & free
        singletons = {v for v in result.network.vacancies if ranks[mask | (1 << v)] == 2}
        if singletons != result.singleton_triggers:
            raise AssertionError("singleton mismatch")
        safe = sorted(result.network.vacancies-singletons)
        expected_pairs = set()
        for pair in combinations(safe, 2):
            pairs_checked += 1
            if ranks[mask | (1 << pair[0]) | (1 << pair[1])] == 2:
                expected_pairs.add(pair)
        if expected_pairs != result.minimal_pairs():
            raise AssertionError("pair predicate mismatch")
        expected_triples = set()
        for triple in combinations(safe, 3):
            if any(pair in expected_pairs for pair in combinations(triple, 2)):
                continue
            triples_checked += 1
            m = mask | sum(1 << v for v in triple)
            if ranks[m] == 2:
                expected_triples.add(triple)
        if expected_triples != result.minimal_triples():
            raise AssertionError("minimal-triple predicate mismatch")
        for v in sorted(result.network.vacancies):
            child = result.network.activate(v)
            if (child is None) != (v in singletons):
                raise AssertionError("absorption mismatch")
            if child is not None:
                rest = free ^ (1 << v)
                subset = rest
                while True:
                    expected = ranks[mask | (1 << v) | subset] == 2
                    if child.connects(states[subset]) != expected:
                        raise AssertionError("contraction-update mismatch")
                    updates += 1
                    if subset == 0:
                        break
                    subset = (subset-1) & rest
    return dict(name=name, period_matrix=[[t.h, t.shear], [0, t.height]],
                all_rank_checks=len(states), rank_one_checkpoints=checkpoints,
                full_future_checks=count, pair_checks=pairs_checked,
                pair_free_triple_checks=triples_checked,
                safe_update_future_checks=updates, failures=0)


def checkpoint(row: dict, independent_large_rank_check: bool = True) -> dict:
    t = HnfSquareTorus(row["N"], row["h12"], 1)
    mask = int(row["occupied_mask_hex"], 16)
    if mask < 0 or mask >= 1 << t.n:
        raise ValueError("invalid occupied mask")
    A = {v for v in range(t.n) if mask >> v & 1}
    if len(A) != row["k0"]:
        raise ValueError("occupied count mismatch")
    result = cut_rank_one(t, A)
    pairs = result.minimal_pairs()
    raw = {key: row[key] for key in ("N", "h12", "k0", "seed", "replica_counter", "ell")}
    raw["safe_sites"] = sorted(result.network.vacancies-result.singleton_triggers)
    raw["minimal_trigger_pairs"] = [list(edge) for edge in sorted(pairs)]
    raw_bytes = (json.dumps(raw, separators=(",", ":")) + "\n").encode()
    reproduced = git_blob(raw_bytes)
    if reproduced != row["graph_source_git_blob"]:
        raise AssertionError("cut-generated graph is not the archived Git blob")
    triples = result.minimal_triples()
    stats = pair_statistics(pairs)
    if stats["wedges"] != row["expected_wedges"] or len(triples) != row["expected_minimal_triples"]:
        raise AssertionError("saved checkpoint count mismatch")
    nonisolated = {v for pair in pairs for v in pair}
    by_middle = Counter(next(v for v in triple if v in result.neutral_sites) for triple in triples)
    bicliques = result.bicliques()
    for block in bicliques:
        if any(v >= t.n for v in block["old_sites"]):
            raise AssertionError("a neutral block contains a split cycle site")
        block["old_ambient_rank"] = t.rank_union_find(block["old_sites"])
        block["side_sizes"] = [len(block["left"]), len(block["right"])]
    direct = sorted(tuple(sorted((u, v))) for u in result.left_sites for v in result.right_sites
                    if v in result.network.adjacency[u])
    alternative = cut_rank_one(t, A, t.essential_cycle(A, reverse_search=True))
    if alternative.minimal_pairs() != pairs or alternative.minimal_triples() != triples:
        raise AssertionError("changing the occupied cut changes an event")
    checked = 0
    if independent_large_rank_check:
        for v in result.network.vacancies:
            if (t.rank_union_find(A | {v}) == 2) != (v in result.singleton_triggers):
                raise AssertionError("large singleton mismatch")
            checked += 1
        for pair in combinations(raw["safe_sites"], 2):
            if (t.rank_union_find(A | set(pair)) == 2) != (pair in pairs):
                raise AssertionError("large pair mismatch")
            checked += 1
        for triple in triples:
            if t.rank_union_find(A | set(triple)) != 2:
                raise AssertionError("large predicted triple does not trigger")
            checked += 1
    return dict(name=row["name"], counter=row["replica_counter"],
                graph_git_blob=reproduced, exact_archived_edge_equality=True,
                original_cut_cycle=[list(result.cycle[0]), list(result.cycle[1]), list(result.cycle[2])],
                alternative_cycle_length=len(alternative.cycle[0]),
                alternative_cut_pair_and_triple_sets_equal=True,
                network_vertices=len(result.network.adjacency),
                network_edges=sum(map(len, result.network.adjacency.values())) // 2,
                permanent_components=len(result.old_components),
                left_sites=sorted(result.left_sites), right_sites=sorted(result.right_sites),
                neutral_sites=sorted(result.neutral_sites),
                nonisolated_left=sorted(result.left_sites & nonisolated),
                nonisolated_right=sorted(result.right_sites & nonisolated),
                singleton_triggers=sorted(result.singleton_triggers),
                pair_statistics=stats, direct_pair_edges=direct, bicliques=bicliques,
                minimal_triples=len(triples), triples_by_neutral_middle={str(k): v for k, v in sorted(by_middle.items())},
                independent_large_rank_checks=checked)


def run(include_large_checks: bool = True) -> dict:
    specs = [("N5", (5, 2, 1)), ("N8", (4, 2, 2)), ("N9", (3, 0, 3)),
             ("N10", (10, 3, 1)), ("N13", (13, 8, 1))]
    inputs = json.loads(INPUT.read_text())
    return dict(issue=487, new_random_samples=0,
                scope="Embedded square-NN rank-one checkpoints; fixed-cut vertex updates. Not a CFT, H4, or continuum-dimension claim.",
                source_commit=inputs["source_commit"], input_sha256=hashlib.sha256(INPUT.read_bytes()).hexdigest(),
                tiny_controls=[exhaustive_control(name, params) for name, params in specs],
                real_checkpoints=[checkpoint(row, include_large_checks) for row in inputs["rows"]])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--verify", type=Path)
    args = parser.parse_args()
    report = run()
    if args.verify:
        if report != json.loads(args.verify.read_text()):
            raise SystemExit("stored report differs from exact regeneration")
        print("exact regeneration verified")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if not args.output and not args.verify:
        print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
