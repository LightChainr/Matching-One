#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_conditional_replica_filtration import render  # noqa: E402

PAYLOAD = render()


def geometry(n: int) -> dict:
    return next(item for item in PAYLOAD["geometries"] if item["N"] == n)


def test_every_level_passes_replica_four_replica_psd_and_telescope() -> None:
    for item in PAYLOAD["geometries"]:
        for vector in item["vectors"].values():
            assert vector["telescoping_exact"] is True
            assert vector["covariance_Y"] == vector["sum_Gamma"]
            assert vector["krawtchouk_F0_control"][
                "matches_PR245_p_biased_subset_moments"
            ] is True
            assert vector["krawtchouk_F0_control"]["C0_from_krawtchouk"] == vector[
                "levels"
            ][0]["C_j"]
            for level in vector["levels"]:
                assert level["C_j"] == level["C_j_conditional_replica"]
                assert level["Gamma_j"] == level["Gamma_j_four_replica"]
                assert level["rank_psd"]["psd_exact"] is True
            for cross in vector["martingale_orthogonality"].values():
                assert cross == [["0", "0"], ["0", "0"]]


def test_n10_topology_plane_is_rank_two_and_rotates_at_both_refinements() -> None:
    topology = geometry(10)["vectors"]["topology_pair"]
    assert [level["rank_psd"]["rank"] for level in topology["levels"]] == [2, 2, 2]
    assert [level["rank_psd"]["offdiagonal_sign"] for level in topology["levels"]] == [
        "negative",
        "positive",
        "negative",
    ]
    assert topology["levels"][1]["state_count"] == 3**5
    rotations = [
        item["principal_axis_rotation_degrees"]
        for item in topology["consecutive_principal_axis_rotations"]
    ]
    assert rotations[0] > 60.0
    assert rotations[1] > 60.0


def test_cross_irrep_pair_separates_radial_and_orientation_information() -> None:
    n5 = geometry(5)["vectors"]["symmetry_pair"]
    n10 = geometry(10)["vectors"]["symmetry_pair"]
    assert [level["rank_psd"]["rank"] for level in n5["levels"]] == [1, 0, 1]
    assert [level["rank_psd"]["rank"] for level in n10["levels"]] == [1, 1, 2]
    for vector in (n5, n10):
        for level in vector["levels"]:
            gamma = level["Gamma_j"]
            assert gamma[0][1] == gamma[1][0] == "0"


def test_n5_topology_is_completely_radial_but_n10_is_not() -> None:
    n5 = geometry(5)["vectors"]["topology_pair"]
    n10 = geometry(10)["vectors"]["topology_pair"]
    assert [level["rank_psd"]["rank"] for level in n5["levels"]] == [2, 0, 0]
    assert n5["levels"][0]["cumulative_predictable_trace_fraction"] == "1"
    assert n10["levels"][0]["cumulative_predictable_trace_fraction"] != "1"
    assert n10["levels"][-1]["cumulative_predictable_trace_fraction"] == "1"
