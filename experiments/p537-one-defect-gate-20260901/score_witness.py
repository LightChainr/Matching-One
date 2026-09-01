#!/usr/bin/env python3
"""Score one physical #537 edge with the globally frozen pooled-root Schur data."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction as F
from pathlib import Path

from vendor_siteflip_schur import (
    GEOMETRIES,
    Interval,
    baseline_packet,
    interval_record,
    offsite_weight,
    parse_fraction,
    read_aggregates,
    read_baseline,
    read_root,
    sha256,
    source_global_packets,
)


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASELINE = ROOT / "experiments" / "p537-landing-matrix-preflight-20260901"


def bell_from_joint(key: int) -> int:
    """Restrict a canonical x4+y4+z4 key to x4+y4 and recanonicalize it."""

    canonical: dict[int, int] = {}
    bell = 0
    for port in range(8):
        label = (key >> (4 * port)) & 15
        if label not in canonical:
            canonical[label] = len(canonical)
        bell |= canonical[label] << (3 * port)
    return bell


def excludes_zero(value: Interval) -> bool:
    return value.lo > 0 or value.hi < 0


def display_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def score(args: argparse.Namespace) -> dict[str, object]:
    witness = json.loads(args.witness.read_text())
    if witness.get("schema") != "matching-one/p537-one-defect-witness/v1":
        raise ValueError("unexpected witness schema")
    if witness.get("geometry", {}).get("id") != "axis":
        raise ValueError("the frozen witness must be in the axis geometry")

    nodes = witness.get("nodes")
    if not isinstance(nodes, list) or len(nodes) != 2:
        raise ValueError("witness needs exactly two endpoint nodes")
    before, after = nodes
    for node in nodes:
        if bell_from_joint(int(node["joint_C"])) != int(node["bell"]):
            raise ValueError("Bell key is not the x/y restriction of joint C")
    if int(before["rank_index_q_plus_1"]) == int(after["rank_index_q_plus_1"]):
        raise ValueError("edge does not change landing/rank")
    if int(before["bell"]) == int(after["bell"]):
        raise ValueError("edge does not change source Bell")
    if int(before["g16"]) == int(after["g16"]):
        raise ValueError("edge does not change the actual source value")

    n = int(witness["N"])
    p = read_root(args.baseline_root)
    coefficients = {
        geometry: read_baseline(path, n)
        for geometry, path in {
            "axis": args.baseline_axis,
            "tilted": args.baseline_tilted,
        }.items()
    }
    baseline = {
        geometry: baseline_packet(coefficients[geometry], n, p)
        for geometry in GEOMETRIES
    }
    rows = read_aggregates(args.global_source, n)
    source_packets, beta = source_global_packets(
        rows, baseline, n, p, args.a_raw_denominator
    )

    component = str(witness["source_ports"]["source_component"])
    if set(beta) != {component}:
        raise ValueError(f"global source file must contain only {component!r}")
    delta = parse_fraction(args.delta)
    c = {"axis": F(1, 1) / delta, "tilted": -F(1, 1) / delta}
    mt = (baseline["axis"]["q_t"] + baseline["tilted"]["q_t"]) / 2
    yt = (baseline["axis"]["e_t"] - baseline["tilted"]["e_t"]) / delta
    root_response = yt / mt
    mu_h = {
        geometry: 2 * c[geometry] * baseline[geometry]["e"]
        - root_response * baseline[geometry]["q"]
        for geometry in GEOMETRIES
    }

    geometry = "axis"
    k = int(witness["k_minus"])
    off_weight = offsite_weight(k, n, p)
    mu_a = source_packets[(geometry, component)]["mu_a"]
    beta_component = beta[component]
    stats = witness["sufficient_statistics"]
    source_total = Interval.of(0)
    counterterm_total = Interval.of(0)
    state_records: list[dict[str, object]] = []
    for state in (0, 1):
        wi = (1 - p) if state == 0 else p
        ui = state - p
        s_minus = k - (n - 1) * p
        si = s_minus + ui
        bi = ui * si - p * (1 - p)
        count = F(stats["count"])
        sum_q = F(stats[f"sum_q{state}"])
        sum_e = F(stats[f"sum_E{state}"])
        divisor = n * args.a_raw_denominator
        sum_a = F(stats[f"sum_a16_{state}"]) / divisor
        sum_qa = F(stats[f"sum_q{state}_a16_{state}"]) / divisor
        sum_ea = F(stats[f"sum_E{state}_a16_{state}"]) / divisor
        sum_h = 2 * c[geometry] * sum_e - root_response * sum_q - mu_h[geometry] * count
        sum_ha = (
            2 * c[geometry] * (sum_ea - mu_a * sum_e)
            - root_response * (sum_qa - mu_a * sum_q)
            - mu_h[geometry] * (sum_a - mu_a * count)
        )

        # The factor 1/2 is the geometry pool.  The individual representative
        # uses one fixed z direction; the orbit record restores all four.
        positive_mass = off_weight * wi / 2
        source_rep = positive_mass * ui * sum_ha
        counterterm_rep = -positive_mass * beta_component * bi * sum_h
        source_orbit = args.z_orbit_multiplicity * source_rep
        counterterm_orbit = args.z_orbit_multiplicity * counterterm_rep
        source_total += source_orbit
        counterterm_total += counterterm_orbit
        state_records.append(
            {
                "state": state,
                "positive_mass_Pi_fixed_frame_pooled": interval_record(positive_mass),
                "source_midpoint_Si_fixed_frame_pooled": interval_record(source_rep),
                "root_counterterm_Si_fixed_frame_pooled": interval_record(counterterm_rep),
                "full_Si_fixed_frame_pooled": interval_record(source_rep + counterterm_rep),
                "source_midpoint_Si_C4_orbit_pooled": interval_record(source_orbit),
                "root_counterterm_Si_C4_orbit_pooled": interval_record(counterterm_orbit),
                "full_Si_C4_orbit_pooled": interval_record(source_orbit + counterterm_orbit),
            }
        )

    full_total = source_total + counterterm_total
    delta_a = F(int(after["g16"]) - int(before["g16"]), n * args.a_raw_denominator)
    robust = excludes_zero(source_total) and delta_a != 0
    if not robust or not excludes_zero(full_total):
        raise AssertionError("frozen edge failed the allocation-robust nonzero stop rule")

    payload = {
        "schema": "matching-one/p537-one-defect-score/v1",
        "status": "allocation_robust_physical_diagonal_edge_nonzero",
        "graph_contract": {
            "vertex": "(geometry,eta,z,state,rank_index,Bell,C_joint,B_outer,W_outer)",
            "physical_edge": "the Bernoulli pair obtained by changing only X_z from 0 to 1",
            "landing_rank_slow_variable": "rank_index=q+1 together with B/W landing labels",
            "source_Bell_slow_variable": "the x4+y4 restriction of the common C_joint map",
            "positive_mass": "P_i=(1/2) nu_{g,-z}(eta) w_i for one fixed-frame physical edge",
            "signed_mass": "S_i=P_i Htilde_i[(a_i-mu_a)u_i-beta_lambda b_i]",
            "edge_weight": "W(e)=S_0+S_1; reported both for one fixed direction and its frozen C4 orbit",
            "global_freeze": "p,R,mu_H,mu_a,beta_lambda are computed before cell disintegration from both complete geometries",
        },
        "transition_id": witness["transition_id"],
        "topology": {
            "rank_before_after": [before["rank_index_q_plus_1"], after["rank_index_q_plus_1"]],
            "bell_before_after": [before["bell"], after["bell"]],
            "joint_C_before_after": [before["joint_C"], after["joint_C"]],
            "g16_before_after": [before["g16"], after["g16"]],
            "delta_a_exact": str(delta_a),
            "off_port": witness["off_port"],
            "collar": witness["collar"],
        },
        "global": {
            "root_p": interval_record(p),
            "M_t": interval_record(mt),
            "Y_t": interval_record(yt),
            "R": interval_record(root_response),
            "mu_H_axis": interval_record(mu_h["axis"]),
            "mu_a_axis_component": interval_record(mu_a),
            "beta_component": interval_record(beta_component),
        },
        "states": state_records,
        "edge_weight_C4_orbit_pooled": {
            "source_midpoint_part": interval_record(source_total),
            "root_counterterm_part": interval_record(counterterm_total),
            "full": interval_record(full_total),
            "allocation_robust": robust,
            "reason": "Delta a is nonzero and the beta-free source-midpoint sum already excludes zero",
        },
        "stop_decision": {
            "rule": "first physical edge changing both slow variables with nonzero full Schur weight",
            "triggered": True,
            "blanket_full_graph_two_independent_defect_route": "falsified",
            "contact_collision_sector": "must remain in the leading signed functional",
            "separated_sector": "open after the fixed distance-at-most-one contact split; a six-arm gain still requires disjoint-annulus localization",
            "next_object": "the contact/collision contribution to the surviving leading four-arm signed functional",
            "no_further_graph_enumeration_required": True,
        },
        "margin_boundary": {
            "full_graph_zero_margins_checked": False,
            "reason": "Stop rule A is existential; S*1 and 1^T*S are full-graph checks needed only if no diagonal edge survives",
        },
        "scope": [
            "This exact finite N25 edge disproves an automatic two-spatial-defect or six-arm gain for the full physical decomposition.",
            "It is a local contact/collision witness (arm_mask=3), not an ordinary separated four-arm event or an asymptotic lower bound.",
            "The beta-free source part is nonzero, so the decision does not rely on counterterm allocation across Bell cells.",
        ],
        "inputs": {
            "witness": {"path": display_path(args.witness), "sha256": sha256(args.witness)},
            "global_source": {"path": display_path(args.global_source), "sha256": sha256(args.global_source)},
            "baseline_axis": {"path": display_path(args.baseline_axis), "sha256": sha256(args.baseline_axis)},
            "baseline_tilted": {"path": display_path(args.baseline_tilted), "sha256": sha256(args.baseline_tilted)},
            "baseline_root": {"path": display_path(args.baseline_root), "sha256": sha256(args.baseline_root)},
        },
    }
    carrier_scope = witness.get("carrier_scope")
    if carrier_scope is not None:
        if carrier_scope.get("classification") != "joint_incidence_typed_carrier":
            raise ValueError("unexpected nonadjacent carrier classification")
        if carrier_scope.get("annular_separation_certified") is not False:
            raise ValueError("fixed N5 witness must not claim annular separation")
        payload["topology"]["pairwise_NN_distances"] = witness["pairwise_NN_distances"]
        payload["topology"]["carrier_scope"] = carrier_scope
        payload["stop_decision"] = {
            "rule": "fixed physical edge changing both slow variables with nonzero full Schur weight",
            "triggered": True,
            "blanket_full_graph_two_independent_defect_route": "falsified",
            "distance_at_most_one_contact_split": "insufficient: the distance-two edge still changes joint terminal incidence",
            "carrier_classification": "joint-incidence/typed carrier",
            "separated_sector": "open: pairwise distance two on N5 is not a disjoint-annulus certificate",
            "next_object": "the joint-incidence/typed-carrier contribution to the surviving leading signed functional",
            "no_further_graph_enumeration_required": True,
        }
        payload["scope"] = [
            "This exact finite N25 edge disproves an automatic two-spatial-defect or six-arm gain for the full physical decomposition.",
            "All three marked centers have NN distance two, yet joint terminal incidence drops from two to one; metric distance at most one therefore does not define the full contact carrier.",
            "Distance two on the N5 quotient is not annular separation, so this witness proves neither a separated-sector statement nor an asymptotic lower bound.",
            "The beta-free source part is nonzero, so the decision does not rely on counterterm allocation across Bell cells.",
        ]
    return payload


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--witness", type=Path, default=ROOT / "results" / "p537-one-defect-gate-20260901" / "witness.json")
    result.add_argument("--global-source", type=Path, default=HERE / "global-diag1.csv")
    result.add_argument("--baseline-axis", type=Path, default=BASELINE / "baseline-axis.csv")
    result.add_argument("--baseline-tilted", type=Path, default=BASELINE / "baseline-tilted.csv")
    result.add_argument("--baseline-root", type=Path, default=BASELINE / "baseline-root.json")
    result.add_argument("--output", type=Path, default=ROOT / "results" / "p537-one-defect-gate-20260901" / "result.json")
    result.add_argument("--delta", default="1152/625")
    result.add_argument("--a-raw-denominator", type=int, default=16)
    result.add_argument("--z-orbit-multiplicity", type=int, default=4)
    return result


def main() -> None:
    args = parser().parse_args()
    payload = score(args)
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
