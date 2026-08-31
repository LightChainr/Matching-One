#!/usr/bin/env python3
"""Two safe-disjoint essential white cycles certify trigger bipartiteness.

No new random paths: replay the 22 already selected checkpoint counters, then
optionally exhaust a bounded set of tiny quotients for this new certificate.
The theorem is in notes/p334-trigger-bipartite-theorem.md.  The min-cost
circulation only constructs independently inspectable cycles; it is not the
proof of the annular packing statement.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

import networkx as nx

from integer_period_torus import classify_configuration, integer_torus_geometry
from p334_checkpoint_scalar_collision import archived_permutation
from rank_one_survival_certificate import RankCache, trigger_layers

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "results/local-20260831/P334-cooperative-closure/trigger_graph_raw"
OUTPUT = ROOT / "results/p334-trigger-bipartite"


def bezout(a, b):
    """Return integral dual covector to a primitive period-basis line."""
    old_r, r, old_s, s, old_t, t = abs(a), abs(b), 1, 0, 0, 1
    while r:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t
    if old_r != 1:
        raise ValueError("line must be primitive")
    return old_s * (1 if a >= 0 else -1), old_t * (1 if b >= 0 else -1)


def white_spine(geometry, white, safe, line):
    """Fixed planar site/face incidence spine, with per-site deletion semantics.

    Face f is the unit square rooted at geometry.coordinates[f].  A spoke
    joins its centre to each white corner; parallel spokes on small quotients
    are retained with exact deck addresses.
    """
    graph = nx.MultiDiGraph()
    n = geometry.n
    capacities = {v: (1 if v in safe else 2) for v in white}
    capacities.update({n + f: 2 for f in range(n)})
    for vertex, capacity in capacities.items():
        graph.add_edge((vertex, 0), (vertex, 1), key="capacity",
                       capacity=capacity, weight=0, winding=(0, 0))
    covector = bezout(*line)
    for face, (x, y) in enumerate(geometry.coordinates):
        for dx, dy in ((0, 0), (1, 0), (1, 1), (0, 1)):
            site = geometry.vertex((x + dx, y + dy))
            if site not in white:
                continue
            sx, sy = geometry.coordinates[site]
            winding = geometry.periods.winding((x + dx - sx, y + dy - sy))
            charge = sum(a * b for a, b in zip(covector, winding))
            key = f"{face}:{dx}:{dy}"
            graph.add_edge((n + face, 1), (site, 0), key=key,
                           capacity=2, weight=-charge, winding=winding,
                           spoke=[face, dx, dy], direction="face-to-site")
            graph.add_edge((site, 1), (n + face, 0), key=key,
                           capacity=2, weight=charge,
                           winding=tuple(-v for v in winding),
                           spoke=[face, dx, dy], direction="site-to-face")
    return graph


def decompose_circulation(graph, flow):
    remaining = {(u, v, key): value
                 for u, nbrs in flow.items() for v, keys in nbrs.items()
                 for key, value in keys.items() if value > 0}
    cycles = []
    while remaining:
        start = next(iter(remaining))[0]
        position, path, vertex = {}, [], start
        while vertex not in position:
            position[vertex] = len(path)
            edge = next(edge for edge in remaining if edge[0] == vertex)
            path.append(edge)
            vertex = edge[1]
        cycle = path[position[vertex]:]
        count = min(remaining[edge] for edge in cycle)
        for edge in cycle:
            remaining[edge] -= count
            if remaining[edge] == 0:
                del remaining[edge]
        winding = tuple(sum(graph[u][v][key]["winding"][axis]
                            for u, v, key in cycle) for axis in (0, 1))
        if winding != (0, 0):
            for _ in range(count):
                cycles.append((cycle, winding))
    return cycles


def certificate(matrix, occupied_labels, safe_labels, trigger_pairs):
    geometry = integer_torus_geometry(tuple(map(tuple, matrix)))
    n = geometry.n
    # Inputs here use the backend's vertex ids. Archive labels are converted
    # explicitly by archived_certificates; tiny mode already uses these ids.
    occupied, safe = set(occupied_labels), set(safe_labels)
    white = set(range(n)) - occupied
    black_channels, _ = classify_configuration(geometry, [v in occupied for v in range(n)])
    white_channels, components = classify_configuration(
        geometry, [v in white for v in range(n)], matching=True)
    if black_channels.max_rank != 1 or white_channels.max_rank != 1:
        raise AssertionError("not a rank-one matching pair")
    lines = {c.basis[0] for c in components if c.rank == 1}
    if len(lines) != 1:
        raise AssertionError("white essential components do not share one line")
    line = next(iter(lines))
    graph = white_spine(geometry, white, safe, line)
    flow = nx.min_cost_flow(graph)
    cycles = decompose_circulation(graph, flow)
    if len(cycles) < 2:
        raise AssertionError("no two-cycle packing")
    records = []
    for edges, winding in cycles[:2]:
        if winding != line:
            raise AssertionError(f"nonprimitive or wrongly oriented cycle: {winding}, {line}")
        vertices = [u[0] for u, v, key in edges if key == "capacity"]
        if len(vertices) != len(set(vertices)):
            raise AssertionError("cycle not simple in incidence spine")
        sites = sorted(v for v in vertices if v < n)
        records.append({"winding": list(winding), "white_sites": sites,
                        "safe_sites": sorted(set(sites) & safe),
                        "spoke_walk": [{"spoke": graph[u][v][key]["spoke"],
                                        "direction": graph[u][v][key]["direction"],
                                        "deck": list(graph[u][v][key]["winding"])}
                                       for u, v, key in edges if key != "capacity"]})
    left, right = (set(item["safe_sites"]) for item in records)
    if left & right:
        raise AssertionError("cycle packing shares a safe site")
    uncovered = [(u, v) for u, v in trigger_pairs
                 if not ((u in left and v in right) or (v in left and u in right))]
    if uncovered:
        raise AssertionError(f"trigger pair not cross-cycle: {uncovered[0]}")
    return {"N": n, "period_matrix": matrix, "white_line": list(line),
            "occupied_vertices": sorted(occupied), "safe_vertices": sorted(safe),
            "minimal_trigger_pairs": [list(edge) for edge in trigger_pairs],
            "cycles": records, "trigger_edges": len(trigger_pairs),
            "packing_winding": len(cycles),
            "all_trigger_pairs_cross_cycles": True,
            "safe_cycle_intersection": [],
            "outside_cycles_incident_trigger_vertices": []}


def archived_certificates():
    rows = []
    for path in sorted(ARCHIVE.glob("*.json")):
        record = json.loads(path.read_text())
        n = record["N"]
        matrix = [[n, record["h12"]], [0, 1]]
        geometry = integer_torus_geometry(tuple(map(tuple, matrix)))
        vertex = lambda label: geometry.vertex((label, 0))
        prefix = archived_permutation(n, record["seed"], record["replica_counter"])[:record["k0"]]
        item = certificate(matrix, [vertex(v) for v in prefix],
                           [vertex(v) for v in record["safe_sites"]],
                           [(vertex(u), vertex(v)) for u, v in record["minimal_trigger_pairs"]])
        # Convert independent backend ids back to the production HNF site labels.
        inverse = {vertex(v): v for v in range(n)}
        item["site_label_by_backend_vertex"] = [inverse[v] for v in range(n)]
        item["source"] = str(path.relative_to(ROOT))
        item["seed"] = record["seed"]
        item["replica_counter"] = record["replica_counter"]
        item["k0"] = record["k0"]
        rows.append(item)
        print(f"archive {path.stem}: two cycles, {item['trigger_edges']} cross-edges", flush=True)
    return rows


def tiny_census():
    matrices = []
    for n in range(1, 7):
        for a in range(1, n + 1):
            if n % a == 0:
                for b in range(a):
                    matrices.append(((a, b), (0, n // a)))
    matrices.extend([((3, 0), (0, 3)), ((3, -1), (1, 3)), ((3, -2), (2, 3))])
    results = []
    for matrix in matrices:
        geometry = integer_torus_geometry(matrix)
        cache = RankCache(geometry)
        counts = Counter()
        first = None
        for mask in range(1 << geometry.n):
            if cache.rank(mask) != 1:
                continue
            counts["rank_one_states"] += 1
            singleton, pairs = trigger_layers(cache, mask)
            if not pairs:
                continue
            safe = [v for v in range(geometry.n) if not (mask >> v) & 1 and v not in singleton]
            # New property: exact pair graph bipartiteness, not old rank/Ferrers tests.
            graph = nx.Graph()
            graph.add_nodes_from(safe)
            graph.add_edges_from(pairs)
            if not nx.is_bipartite(graph):
                return {"counterexample": {"period_matrix": matrix, "mask": mask,
                                            "trigger_pairs": pairs}}
            counts["nonempty_trigger_graphs"] += 1
            if first is None:
                first = certificate(matrix, [v for v in range(geometry.n) if (mask >> v) & 1],
                                    safe, pairs)
                first["mask"] = mask
        results.append({"period_matrix": matrix, "N": geometry.n, **dict(counts),
                        "first_two_cycle_certificate": first})
        print(f"tiny {matrix}: {dict(counts)}", flush=True)
    return {"quotients": results, "counterexample": None}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scope", choices=("archive", "tiny", "both"), default="both")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.scope in ("archive", "both"):
        rows = archived_certificates()
        (args.output / "archived_two_cycles.json").write_text(json.dumps({
            "parent_commit": "1b5a9dea07e1c62f69798fddbf4899ff986c0b72",
            "new_random_paths": 0, "checkpoints": rows,
            "certificate_count": len(rows)}, indent=2, sort_keys=True) + "\n")
    if args.scope in ("tiny", "both"):
        (args.output / "tiny_census.json").write_text(json.dumps(tiny_census(), indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
