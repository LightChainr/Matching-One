#!/usr/bin/env python3
"""Locate minimal-triple middle sites in two already-solved physical networks.

No subset census or reliability solve: eliminate fixed nodes, then count
missing cross-edges between each interior site's L/R neighbors.
"""
from itertools import combinations, product
import json
from pathlib import Path

import networkx as nx

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/p334-all147-prefix-clocks/prefixes"
OUTPUT = ROOT / "results/p334-marked-triple-bridges"
COUNTERS = (43042508631, 43042514803)


def eliminate_fixed(network):
    graph = nx.Graph()
    graph.add_nodes_from(network["vertices"])
    graph.add_edges_from(network["edges"])
    witness = {}
    for u, v in graph.edges:
        witness[u, v], witness[v, u] = [u, v], [v, u]
    eliminations = []
    for fixed in sorted(network["fixed_components"]):
        neighbors = sorted(graph.neighbors(fixed))
        added = []
        for u, v in combinations(neighbors, 2):
            if not graph.has_edge(u, v):
                path = witness[u, fixed] + witness[fixed, v][1:]
                graph.add_edge(u, v)
                witness[u, v], witness[v, u] = path, list(reversed(path))
                added.append([u, v])
        graph.remove_node(fixed)
        eliminations.append({"fixed_component_root": fixed,
                             "neighbors_at_elimination": neighbors,
                             "new_clique_edges": added})
    return graph, witness, eliminations


def main():
    records = []
    for counter in COUNTERS:
        source = json.loads((SOURCE / f"{counter}.json").read_text())
        mapping = source["mapping"]
        records_by_factor = []
        for component in mapping["port_components"]:
            if "two_terminal_network" not in component:
                continue
            network = component["two_terminal_network"]
            graph, edge_paths, eliminations = eliminate_fixed(network)
            s, t = network["terminals"]
            # Both endpoints are retained throughout fixed-node elimination.
            direct = set(graph.neighbors(s)) & set(graph.neighbors(t))
            graph.remove_nodes_from(direct)
            L = set(graph.neighbors(s))
            R = set(graph.neighbors(t))
            I = set(graph) - L - R - {s, t}
            pair_edges = sorted([x, z] for x, z in product(L, R) if graph.has_edge(x, z))
            middles = []
            for y in sorted(I):
                left = sorted(set(graph.neighbors(y)) & L)
                right = sorted(set(graph.neighbors(y)) & R)
                exclusions = [[x, z] for x, z in product(left, right) if graph.has_edge(x, z)]
                bridges = [[x, y, z] for x, z in product(left, right) if not graph.has_edge(x, z)]
                if not bridges:
                    continue
                certificates = []
                for x, _, z in bridges:
                    certificates.append({"ordered_sites_L_I_R": [x, y, z],
                                         "fixed_expanded_path": edge_paths[s, x]
                                         + edge_paths[x, y][1:]
                                         + edge_paths[y, z][1:]
                                         + edge_paths[z, t][1:]})
                middles.append({"middle_site_id": y, "left_neighbors": left,
                                "right_neighbors": right, "potential_cross_pairs": len(left)*len(right),
                                "excluded_existing_pair_edges": exclusions,
                                "minimal_triple_contribution": len(bridges),
                                "ordered_triples": bridges,
                                "contracted_path_certificates": certificates})
            records_by_factor.append({"factor":len(records_by_factor),
                                      "port_addresses":component["addresses"],
                                      "fixed_component_roots":network["fixed_components"],
                                      "fixed_eliminations":eliminations,
                                      "reduced_edges_after_removing_direct": sorted(sorted(e) for e in graph.edges),
                                      "direct_sites":sorted(direct),
                                      "L":sorted(L),"R":sorted(R),"I":sorted(I),
                                      "minimal_pair_edges":pair_edges,
                                      "positive_middle_sites":middles,
                                      "minimal_triples":sum(m["minimal_triple_contribution"] for m in middles)})
        positive = [{"factor":f["factor"],"port_addresses":f["port_addresses"],**m}
                    for f in records_by_factor for m in f["positive_middle_sites"]]
        record = {"counter":counter,"source_artifact":str((SOURCE/f"{counter}.json").relative_to(ROOT)),
                  "original_row":source["clock"]["original_row"],
                  "chosen_essential_root":mapping["chosen_essential_root"],
                  "site_label_semantics":"Original quotient site IDs 0..424; HNF representative (id,0) modulo columns of [[425,268],[0,1]]. Fixed path nodes are occupied-component root IDs; -2/-1 denote the two address ports.",
                  "factors":records_by_factor,
                  "direct_sites":sorted(v for f in records_by_factor for v in f["direct_sites"]),
                  "minimal_pair_count":sum(len(f["minimal_pair_edges"]) for f in records_by_factor),
                  "genuine_minimal_triples":sum(f["minimal_triples"] for f in records_by_factor),
                  "positive_middle_sites":positive,
                  "positive_middle_site_count":len(positive)}
        records.append(record)
        print(counter, "triples",record["genuine_minimal_triples"],"middle sites",record["positive_middle_site_count"])
        for m in positive:
            print("  factor",m["factor"],"y",m["middle_site_id"],"L",m["left_neighbors"],
                  "R",m["right_neighbors"],"excluded",m["excluded_existing_pair_edges"],
                  "contribution",m["minimal_triple_contribution"])
    OUTPUT.mkdir(parents=True,exist_ok=True)
    result = {"parent_commit":"87b6ca5b39084c06143f31cafdaba53f90012e27",
              "selection":"Only the unique all147 same-(H2,b2,degree-square) different-clock witness pair; no expansion",
              "new_samples":0,"new_reliability_solves":0,"global_triple_subset_enumerations":0,
              "method":"Fixed-node neighbor-clique elimination; local marked middle-site L/R nonedge count",
              "records":records}
    (OUTPUT/"middle_site_bridges.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")


if __name__ == "__main__":
    main()
