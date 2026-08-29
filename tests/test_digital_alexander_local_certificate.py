import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "local_certificate", ROOT / "scripts" / "digital_alexander_local_certificate.py"
)
assert SPEC and SPEC.loader
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


def test_all_six_pattern_classes_have_expected_multiplicity():
    payload = MOD.build_certificate()
    assert payload["category_counts"] == {
        "all_black": 1,
        "one_white": 4,
        "two_adjacent_white": 4,
        "two_diagonal_white": 2,
        "three_white": 4,
        "all_white": 1,
    }


def test_only_isolated_diagonal_pair_requires_diagonal_spine_edge():
    rows = MOD.build_certificate()["patterns"]
    diagonal_rows = [row for row in rows if row["category"] == "two_diagonal_white"]
    assert len(diagonal_rows) == 2
    assert all(len(row["spine_edges"]) == 1 for row in diagonal_rows)
    for row in rows:
        if row["category"] in {"three_white", "all_white"}:
            assert not any(
                set(edge) in ({"SW", "NE"}, {"SE", "NW"})
                for edge in row["spine_edges"]
            )


def test_every_redundant_diagonal_has_local_nn_replacement():
    rows = MOD.build_certificate()["patterns"]
    assert sum(len(row["removed_diagonals_with_NN_replacement"]) for row in rows) == 6
    assert sum(
        len(row["removed_diagonals_with_NN_replacement"])
        for row in rows
        if row["category"] == "three_white"
    ) == 4
    assert sum(
        len(row["removed_diagonals_with_NN_replacement"])
        for row in rows
        if row["category"] == "all_white"
    ) == 2


def test_two_diagonal_masks_are_the_expected_bit_patterns():
    rows = MOD.build_certificate()["patterns"]
    masks = {row["mask"] for row in rows if row["category"] == "two_diagonal_white"}
    assert masks == {"0101", "1010"}


def test_face_center_incidence_star_has_one_spoke_per_white_corner():
    rows = MOD.build_certificate()["patterns"]
    for row in rows:
        assert len(row["incidence_spokes"]) == len(row["white_corners"])
        assert {edge[0] for edge in row["incidence_spokes"]} == set(row["white_corners"])
