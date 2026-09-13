from scripts.p14_w5_joint_relative_dual import (
    build_report,
    checkerboard_periodic_embedding,
    relative_state,
    spherical_partition_from_state,
)
from scripts.p14_w5_terminal_duality import (
    partition_key,
    primal_partition,
    spherical_transform,
)


def test_relative_state_reconstructs_spherical_partition_for_all_configurations():
    for mask in range(256):
        bits = tuple((mask >> index) & 1 for index in range(8))
        expected = partition_key(primal_partition(spherical_transform(bits)))
        observed = partition_key(spherical_partition_from_state(relative_state(bits)))
        assert observed == expected


def test_enriched_state_is_a_192_state_connectivity_quotient():
    manifest = {
        "schema": "matching-one/p14-w5-joint-relative-dual/v1",
        "parent_certificate": "results/terminal-reliability/p14-w5-terminal-duality.json",
        "torus_length": 4,
        "claim_boundary": "test",
    }
    finite = build_report(manifest)["finite_cell"]
    assert finite["labelled_enriched_state_count"] == 192
    assert finite["d4_orbit_count"] == 41
    assert finite["configuration_multiplicity_histogram"] == {"1": 176, "5": 16}
    assert finite["terminal_pair_ambiguous_for_spherical_output"] > 0
    assert finite["enriched_state_determines_primal_partition"] is True
    assert finite["enriched_state_determines_spherical_partition_after_outer_gluing"] is True


def test_checkerboard_gluing_is_periodic_dual_but_not_graph_self_dual():
    periodic = checkerboard_periodic_embedding(4)
    assert periodic["black_w5_cells"] == 8
    assert periodic["primal"] == {
        "vertices": 24,
        "edges": 64,
        "degree_histogram": {"4": 8, "6": 16},
    }
    assert periodic["disk_relative_dual"] == {
        "vertices": 40,
        "edges": 64,
        "degree_histogram": {"3": 32, "4": 8},
    }
    assert periodic["edge_bijection_count"] == 64
    assert periodic["every_grid_rim_owned_once"] is True
    assert periodic["torus_euler_characteristic"] == 0
    assert periodic["graph_self_isomorphism_obstructed"] is True


def test_checkerboard_counts_scale_to_next_even_torus():
    periodic = checkerboard_periodic_embedding(6)
    assert periodic["primal"]["vertices"] == 54
    assert periodic["primal"]["edges"] == 144
    assert periodic["disk_relative_dual"]["vertices"] == 90
    assert periodic["edge_bijection_count"] == 144
