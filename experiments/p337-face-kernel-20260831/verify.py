#!/usr/bin/env python3
"""Deterministic N50 gain-graph checks; no enumeration or random sampling."""
from __future__ import annotations

from collections import deque
from datetime import datetime, timezone
from fractions import Fraction as F
from itertools import combinations
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

ROOT = Path(__file__).resolve().parent
STEPS = ((1, 0), (0, 1), (-1, 0), (0, -1))
NN_HALF = ((1, 0), (0, 1))
MATCH_HALF = NN_HALF + ((1, 1), (1, -1))


def key(x, y):
    # Quotient periods (5,5), (-5,5); 50 classes with equal u/v parity.
    return ((x+y) % 10, (y-x) % 10)


def coordinate(vertex):
    u, v = vertex
    return ((u-v)//2, (u+v)//2)


ALL = {(u, v) for u in range(10) for v in range(10) if (u-v) % 2 == 0}
A = {v for v in ALL if v[0] % 2 == 0}
B = ALL-A
ROW = {key(1+i, i) for i in range(5)}


def edge_set(vertices, steps):
    edges = []
    for v in sorted(vertices):
        x, y = coordinate(v)
        for dx, dy in steps:
            w = key(x+dx, y+dy)
            if w in vertices:
                edges.append((v, w, (dx, dy)))
    return edges


def graph_stats(vertices, edges):
    adjacency = {v: [] for v in vertices}
    for v, w, (dx, dy) in edges:
        adjacency[v].append((w, (dx, dy)))
        adjacency[w].append((v, (-dx, -dy)))
    lifts = {}
    components = 0
    windings = []
    for start in sorted(vertices):
        if start in lifts:
            continue
        components += 1
        lifts[start] = (0, 0)
        queue = deque([start])
        while queue:
            v = queue.popleft()
            x, y = lifts[v]
            for w, (dx, dy) in adjacency[v]:
                proposed = (x+dx, y+dy)
                if w not in lifts:
                    lifts[w] = proposed
                    queue.append(w)
                else:
                    gx, gy = proposed[0]-lifts[w][0], proposed[1]-lifts[w][1]
                    assert (gx+gy) % 10 == (gy-gx) % 10 == 0
                    h = ((gx+gy)//10, (gy-gx)//10)
                    if h != (0, 0):
                        windings.append(h)
    rank = int(bool(windings))
    if windings:
        x, y = windings[0]
        if any(x*v-y*u != 0 for u, v in windings[1:]):
            rank = 2
    return {"V": len(vertices), "edges": len(edges), "components": components,
            "beta": len(edges)-len(vertices)+components, "rank": rank}


def clique_graph(occupied_A, occupied_B):
    edges = []
    counts = [0]*5
    correction = 0
    for a in sorted(occupied_A):
        x, y = coordinate(a)
        ports = [(key(x+dx, y+dy), (dx, dy)) for dx, dy in STEPS
                 if key(x+dx, y+dy) in occupied_B]
        m = len(ports)
        counts[m] += 1
        correction += (m-1)*(m-2)//2 if m >= 1 else 0
        # A separate edge per face and port pair; never simplify parallel edges.
        for (v, dv), (w, dw) in combinations(ports, 2):
            edges.append((v, w, (dw[0]-dv[0], dw[1]-dv[1])))
    return graph_stats(occupied_B, edges), counts, correction


def check(name, occupied_A, occupied_B):
    occupied = occupied_A | occupied_B
    vacant = ALL-occupied
    parent = graph_stats(occupied, edge_set(occupied, NN_HALF))
    white = graph_stats(vacant, edge_set(vacant, MATCH_HALF))
    full_faces = 0
    for x, y in map(coordinate, ALL):
        full_faces += all(key(x+dx, y+dy) in occupied
                          for dx, dy in ((0, 0), (1, 0), (0, 1), (1, 1)))
    vacant_edges = len(edge_set(vacant, NN_HALF))
    direct_S = parent["components"] + white["components"] + full_faces + vacant_edges
    q = parent["rank"]-1
    identity_S = 2*parent["beta"]-3*len(occupied)-q+100
    clique, counts, correction = clique_graph(occupied_A, occupied_B)
    clique_S = 2*(clique["beta"]-correction)-3*len(occupied)-q+100
    assert parent["beta"] == clique["beta"]-correction
    assert parent["rank"] == clique["rank"]
    assert direct_S == identity_S == clique_S
    return {"name": name, "K": len(occupied), "q": q, "E": q*q,
            "S_direct": direct_S, "S_cycle_identity": identity_S, "S_clique": clique_S,
            "parent": parent, "clique": clique, "occupied_faces_by_m": counts,
            "cycle_correction": correction, "full_faces": full_faces, "vacant_edges": vacant_edges}


def main():
    started = time.perf_counter()
    start_utc = datetime.now(timezone.utc).isoformat()
    a, b = key(2, 0), key(1, 1)
    witnesses = [("row_no_hole", A, ROW), ("row_a", A-{a}, ROW),
                 ("row_b", A-{b}, ROW), ("row_ab", A-{a, b}, ROW),
                 ("row_plus_one_B", A, ROW | {key(3, 0)}),
                 ("row_plus_two_B", A, ROW | {key(3, 0), key(2, -1)}),
                 ("A_only", A, set()), ("one_B", A, {key(1, 0)}),
                 ("all_occupied", A, B), ("empty", set(), set())]
    checked = [check(*witness) for witness in witnesses]
    assert [(w["K"], w["parent"]["rank"], w["S_direct"]) for w in checked[:4]] == [
        (30, 1, 22), (29, 1, 23), (29, 1, 23), (28, 0, 25)]
    mixed = {field: checked[3][field]-checked[1][field]-checked[2][field]+checked[0][field]
             for field in ("q", "E", "S_direct")}
    # Three face bits only, to verify a derived local probability at r=1/2.
    # This is 8 states of the declared motif, not graph/defect-pattern enumeration.
    r = F(1, 2)
    joint = F(0)
    marginal = F(0)
    for h0 in (0, 1):
        for h1 in (0, 1):
            for h2 in (0, 1):
                probability = (r if h0 else 1-r)*(r if h1 else 1-r)*(r if h2 else 1-r)
                X, Y = 1-h0*h1, 1-h0*h2
                marginal += probability*X
                joint += probability*X*Y
    assert joint == 1-2*r*r+r**3 and joint-marginal**2 == r**3*(1-r)
    result = {"schema": "matching-one.p337-face-kernel.deterministic.v1",
        "periods": [[5, 5], [-5, 5]], "parent_N": 50, "witnesses": checked,
        "review_two_hole_mixed_differences": mixed,
        "local_side_edge_check": {"hole_probability": str(r), "marginal_open": str(marginal),
            "joint_open": str(joint), "independent_same_marginal_joint": str(marginal**2),
            "covariance": str(joint-marginal**2)},
        "new_samples": 0, "new_exhaustive_graph_enumerations": 0,
        "declared_graph_configurations_checked": len(witnesses), "local_3bit_states_checked": 8}
    (ROOT/"latest.json").write_text(json.dumps(result, indent=2)+"\n")
    receipt = {"started_utc": start_utc, "finished_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.perf_counter()-started, "command": sys.argv, "python": sys.version,
        "machine": platform.machine(), "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "result_sha256": hashlib.sha256((ROOT/"latest.json").read_bytes()).hexdigest(),
        "run_count": 1, "cloud_jobs": 0, "new_samples": 0, "exit_code": 0}
    (ROOT/"run.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print(json.dumps({"checked_graphs": len(witnesses), "mixed_differences": mixed,
                      "local_edge_covariance": str(joint-marginal**2)}))


if __name__ == "__main__":
    main()
