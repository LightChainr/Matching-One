import csv
import sys


sys.path.insert(0, "scripts")
from score_p334_local_active_boundary_motif import load_fresh  # noqa: E402


def test_complement_typed_local_motif_coordinates(tmp_path):
    path = tmp_path / "pilot.csv"
    fields = [
        "n", "orientation", "batch", "ell_u", "ell_v", "age_steps",
        "next_exit", "essential_size", "essential_carriers",
        "occupied_frontier", "vacant_frontier", "boundary_cut_edges",
        "boundary_multicontact_sites", "boundary_contact_pairs", "core_vertices",
        "core_edges", "articulation_vertices", "bridges",
        "birth_r1_occupied", "birth_r1_essential",
        "birth_r1_same_side_pairs", "birth_r1_opposite_side_pairs",
        "exit_r1_occupied", "exit_r1_essential",
        "exit_r1_same_side_pairs", "exit_r1_opposite_side_pairs",
    ]
    row = {name: "0" for name in fields}
    row.update({"n": "325", "orientation": "first", "ell_u": "1",
                "birth_r1_occupied": "2", "exit_r1_occupied": "6",
                "birth_r1_essential": "1", "exit_r1_essential": "3"})
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader(); writer.writerow(row)
    loaded = load_fresh(path)[0]
    assert loaded["S_occupied"] == 4.0
    assert loaded["D_occupied"] == 2.0
    assert loaded["S_essential"] == 2.0
    assert loaded["D_essential"] == 1.0
