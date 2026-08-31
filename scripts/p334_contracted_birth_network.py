#!/usr/bin/env python3
"""Exact whole-event gain-graph and two-port mapping for two fixed N425 states."""
from collections import defaultdict, deque
import json
from pathlib import Path

import networkx as nx

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/p334-contracted-full-clock"


def build(row):
    n = 425
    occupied = set(row["occupied_prefix_labels"])
    vacant = set(range(n)) - occupied
    edges = [(v, (v + 1) % n, 1, 0) for v in range(n)]
    edges += [(v, (v - 268) % n, 0, 1) for v in range(n)]
    adjacency = defaultdict(list)
    for u, v, dx, dy in edges:
        adjacency[u].append((v, dx, dy))
        adjacency[v].append((u, -dx, -dy))
    root, potential, occupied_components, essential = {}, {}, {}, []
    for start in sorted(occupied):
        if start in root:
            continue
        root[start], potential[start] = start, (0, 0)
        queue, vertices, cycles = deque([start]), [], set()
        while queue:
            u = queue.popleft()
            vertices.append(u)
            for v, dx, dy in adjacency[u]:
                if v not in occupied:
                    continue
                proposed = (potential[u][0] + dx, potential[u][1] + dy)
                if v not in root:
                    root[v], potential[v] = start, proposed
                    queue.append(v)
                else:
                    cycle = (proposed[0] - potential[v][0], proposed[1] - potential[v][1])
                    if cycle != (0, 0):
                        if 19 * cycle[0] + 8 * cycle[1] != 0:
                            raise ValueError("prefix not on the archived ambient line")
                        cycles.add(cycle)
        occupied_components[start] = sorted(vertices)
        if cycles:
            essential.append(start)
    for v in vacant:
        root[v], potential[v] = v, (0, 0)
    fixed = set(occupied_components)
    contracted = set()
    for u, v, dx, dy in edges:
        a, b = root[u], root[v]
        gx = potential[u][0] + dx - potential[v][0]
        gy = potential[u][1] + dy - potential[v][1]
        gain = 19 * gx + 8 * gy
        if a == b:
            if gain:
                raise ValueError("nonzero fixed loop")
            continue
        if a > b:
            a, b, gain = b, a, -gain
        contracted.add((a, b, gain))
    # Remove one actual essential occupied component. Certify that the entire
    # remaining graph is balanced, then record the port gains to that root.
    terminal_root = min(essential)
    rest = nx.Graph()
    rest.add_nodes_from((fixed | vacant) - {terminal_root})
    incident = defaultdict(list)
    ports = []
    for u, v, gain in sorted(contracted):
        if u == terminal_root:
            ports.append((v, gain))
        elif v == terminal_root:
            ports.append((u, -gain))
        else:
            rest.add_edge(u, v)
            incident[u].append((v, gain))
            incident[v].append((u, -gain))
    components = []
    for nodes in sorted(nx.connected_components(rest), key=lambda c: min(c)):
        start = min(nodes)
        gauge = {start: 0}
        queue = deque([start])
        while queue:
            u = queue.popleft()
            for v, gain in incident[u]:
                proposed = gauge[u] + gain
                if v in gauge:
                    if gauge[v] != proposed:
                        raise ValueError("full graph minus essential component is unbalanced")
                else:
                    gauge[v] = proposed
                    queue.append(v)
        component_ports = [(v, gain - gauge[v]) for v, gain in ports if v in nodes]
        addresses = sorted({g for _, g in component_ports})
        record = {"vertices": sorted(nodes), "vacant_sites": sorted(nodes & vacant),
                  "fixed_components": sorted(nodes & fixed), "gauge": gauge,
                  "ports": component_ports, "addresses": addresses,
                  "ordinary_edges": sorted([sorted(e) for e in rest.subgraph(nodes).edges])}
        if len(addresses) == 2:
            graph = rest.subgraph(nodes).copy()
            terminals = (-2, -1)
            for v, address in component_ports:
                graph.add_edge(terminals[addresses.index(address)], v)
            # Only blocks along the unique block-cut-tree terminal path can
            # participate in an s-t path. All off-path variables are free.
            blocks = list(nx.biconnected_components(graph))
            articulations = set(nx.articulation_points(graph))
            tree = nx.Graph()
            home = {}
            for j, block in enumerate(blocks):
                bnode = ("block", j)
                tree.add_node(bnode)
                for v in block:
                    if v in articulations:
                        vnode = ("cut", v)
                        tree.add_edge(bnode, vnode)
                        home[v] = vnode
                    else:
                        home[v] = bnode
            route = nx.shortest_path(tree, home[terminals[0]], home[terminals[1]])
            core = set().union(*(blocks[node[1]] for node in route if node[0] == "block"))
            core_graph = graph.subgraph(core)
            record["two_terminal_network"] = {
                "terminals": list(terminals), "vertices": sorted(core),
                "edges": sorted([sorted(e) for e in core_graph.edges]),
                "vacant_sites": sorted(core & vacant), "fixed_components": sorted(core & fixed),
                "irrelevant_vacant_sites": sorted((nodes - core) & vacant),
                "block_sizes_on_terminal_route": [len(blocks[node[1]]) for node in route if node[0] == "block"],
                "min_degree_treewidth_upper_bound": nx.approximation.treewidth_min_degree(core_graph)[0],
                "maximal_exact_twin_group_sizes": sorted((len(vs) for vs in twin_groups(core_graph, core & vacant)), reverse=True),
            }
        components.append(record)
    return {"counter": row["replica_counter"], "seed": row["seed"], "N": n, "k0": row["k0"],
            "period_matrix": row["period_matrix"], "ell": row["ell"],
            "physical_line": [8, -19], "physical_transverse_covector": [19, 8],
            "occupied_components": occupied_components, "essential_component_roots": essential,
            "chosen_essential_root": terminal_root, "vacant_sites": sorted(vacant),
            "contracted_edges_with_transverse_gain": sorted(contracted),
            "balanced_after_removing_essential_root": True, "port_components": components}


def twin_groups(graph, sites):
    groups = defaultdict(list)
    for v in sites:
        groups[tuple(sorted(graph.neighbors(v)))].append(v)
    return list(groups.values())


def main():
    source = json.loads((ROOT / "results/p334-quartic-clock/full_quartics.json").read_text())["checkpoints"]
    records = [build(row) for row in source]
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT / "whole_event_networks.json").write_text(json.dumps({"parent_commit": "1614a17e10997656fdf2d5520846fff2a228a5cd", "new_samples": 0, "records": records}, indent=2, sort_keys=True) + "\n")
    for row in records:
        print(row["counter"], "fixed components", len(row["occupied_components"]), "essential", row["essential_component_roots"], "contracted edges", len(row["contracted_edges_with_transverse_gain"]))
        for c in row["port_components"]:
            net = c.get("two_terminal_network")
            print(" component", len(c["vertices"]), "addresses", c["addresses"],
                  "network", None if net is None else {k: v for k, v in net.items() if k not in ("vertices", "edges", "vacant_sites", "irrelevant_vacant_sites")},
                  "core vacant", None if net is None else len(net["vacant_sites"]))


if __name__ == "__main__":
    main()
