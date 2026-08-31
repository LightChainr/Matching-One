"""Read exact new cycle witnesses, without repeating old pair/Ferrers scans."""
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from integer_period_torus import integer_torus_geometry


def check_certificate(record):
    geometry = integer_torus_geometry(tuple(map(tuple, record["period_matrix"])))
    occupied, safe = set(record["occupied_vertices"]), set(record["safe_vertices"])
    safe_sets = []
    for cycle in record["cycles"]:
        endpoints, winding, nodes = [], [0, 0], []
        for step in cycle["spoke_walk"]:
            face, dx, dy = step["spoke"]
            x, y = geometry.coordinates[face]
            site = geometry.vertex((x + dx, y + dy))
            sx, sy = geometry.coordinates[site]
            deck = geometry.periods.winding((x + dx - sx, y + dy - sy))
            source, target = ("face", face), ("site", site)
            if step["direction"] == "site-to-face":
                source, target, deck = target, source, tuple(-v for v in deck)
            assert list(deck) == step["deck"]
            assert site not in occupied
            endpoints.append((source, target))
            nodes.append(source)
            winding = [winding[i] + deck[i] for i in (0, 1)]
        assert all(endpoints[i][1] == endpoints[(i + 1) % len(endpoints)][0]
                   for i in range(len(endpoints)))
        assert len(nodes) == len(set(nodes))
        assert winding == record["white_line"] == cycle["winding"]
        sites = {v for kind, v in nodes if kind == "site"}
        assert sites == set(cycle["white_sites"])
        assert sites & safe == set(cycle["safe_sites"])
        safe_sets.append(sites & safe)
    left, right = safe_sets
    assert not left & right
    assert all((u in left and v in right) or (v in left and u in right)
               for u, v in record["minimal_trigger_pairs"])


def test_all_22_archived_two_cycle_certificates():
    archive = json.loads((ROOT / "results/p334-trigger-bipartite/archived_two_cycles.json").read_text())
    assert archive["new_random_paths"] == 0
    assert len(archive["checkpoints"]) == archive["certificate_count"] == 22
    for record in archive["checkpoints"]:
        check_certificate(record)


def test_tiny_quotient_cycle_certificates():
    tiny = json.loads((ROOT / "results/p334-trigger-bipartite/tiny_census.json").read_text())
    assert tiny["counterexample"] is None
    assert {row["N"] for row in tiny["quotients"]} == {1, 2, 3, 4, 5, 6, 9, 10, 13}
    for row in tiny["quotients"]:
        certificate = row["first_two_cycle_certificate"]
        if certificate is not None:
            check_certificate(certificate)
