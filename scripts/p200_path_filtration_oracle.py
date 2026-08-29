#!/usr/bin/env python3
"""Exact path-filtration oracle for Issue #200 Phase E.

The final R2/R5 join is order independent.  This oracle instead retains the
intermediate ambient-H1 rank and asks when a final rank-two event first becomes
visible along the two filtrations.  It also records the corresponding exact
Doob and conditional-variance diagnostics under the uniform measure on the
smallest honest Gaussian N10 lift.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
from typing import Callable

sys.path.insert(0, str(Path(__file__).resolve().parent))

from p200_typed_mixed_homology_oracle import (  # noqa: E402
    FINAL_MATRIX,
    GEOMETRY,
    PERIOD_2,
    PERIOD_5,
    column_hnf,
    configuration_record,
)


LAYERS = ("black_NN", "white_matching")


def fraction_record(value: Fraction) -> dict[str, str | float]:
    return {"exact": str(value), "decimal": float(value)}


def fraction_histogram(values: list[Fraction]) -> dict[str, int]:
    return {
        str(value): count
        for value, count in sorted(Counter(values).items())
    }


def digest_fractions(values: list[Fraction]) -> str:
    encoded = json.dumps([str(value) for value in values], separators=(",", ":"))
    return hashlib.sha256(encoded.encode("ascii")).hexdigest()


def first_rank_two_time(h0: int, h_first: int, h25: int) -> int | None:
    """First stage in {0,1,2} at which ambient H1 has rank two."""

    if h25 != 2:
        return None
    if h0 == 2:
        return 0
    if h_first == 2:
        return 1
    return 2


def activation_row(record: dict, layer: str) -> dict[str, int | None]:
    ranks = record[layer]["ambient_H1_ranks"]
    tau_2_then_5 = first_rank_two_time(ranks["0"], ranks["2"], ranks["25"])
    tau_5_then_2 = first_rank_two_time(ranks["0"], ranks["5"], ranks["25"])
    contrast = 0
    if tau_2_then_5 is not None:
        contrast = tau_5_then_2 - tau_2_then_5
    return {
        "final_rank_two_event": int(ranks["25"] == 2),
        "tau_R2_then_R5": tau_2_then_5,
        "tau_R5_then_R2": tau_5_then_2,
        "C_activation": contrast,
    }


def conditional_expectations(
    rows: list[dict], key: Callable[[dict], tuple], target: Callable[[dict], int]
) -> dict[tuple, Fraction]:
    totals: dict[tuple, list[int]] = defaultdict(lambda: [0, 0])
    for row in rows:
        bucket = totals[key(row)]
        bucket[0] += target(row)
        bucket[1] += 1
    return {
        state: Fraction(total, count)
        for state, (total, count) in totals.items()
    }


def ambient_rank_doob(rows: list[dict], layer: str) -> dict:
    def ranks(row: dict) -> dict[str, int]:
        return row[layer]["ambient_H1_ranks"]

    target = lambda row: int(ranks(row)["25"] == 2)
    key_0 = lambda row: (ranks(row)["0"],)
    key_2 = lambda row: (ranks(row)["0"], ranks(row)["2"])
    key_5 = lambda row: (ranks(row)["0"], ranks(row)["5"])
    key_25 = lambda row: (
        ranks(row)["0"],
        ranks(row)["2"],
        ranks(row)["5"],
        ranks(row)["25"],
    )
    means = {
        "0": conditional_expectations(rows, key_0, target),
        "2": conditional_expectations(rows, key_2, target),
        "5": conditional_expectations(rows, key_5, target),
        "25": conditional_expectations(rows, key_25, target),
    }
    output = []
    for row in rows:
        m0 = means["0"][key_0(row)]
        m2 = means["2"][key_2(row)]
        m5 = means["5"][key_5(row)]
        m25 = means["25"][key_25(row)]
        if m25 != target(row):
            raise AssertionError("final rank sigma-algebra did not reveal target")
        q_2_then_5 = (m2 - m0) ** 2 + (m25 - m2) ** 2
        q_5_then_2 = (m5 - m0) ** 2 + (m25 - m5) ** 2
        loss_2 = m2 * (1 - m2)
        loss_5 = m5 * (1 - m5)
        output.append(
            {
                "M0": m0,
                "M2": m2,
                "M5": m5,
                "M25": m25,
                "Q_R2_then_R5": q_2_then_5,
                "Q_R5_then_R2": q_5_then_2,
                "C_Doob_Q": q_2_then_5 - q_5_then_2,
                "conditional_variance_after_R2": loss_2,
                "conditional_variance_after_R5": loss_5,
                "R2_information_advantage": loss_5 - loss_2,
            }
        )
    return {"rows": output, "conditional_means": means}


def partition_residual_doob(rows: list[dict], layer: str) -> list[Fraction]:
    """Natural scalar-rank filtration for the final Rc != 0 event."""

    def ranks(row: dict) -> dict[str, int]:
        return row[layer]["partition_ranks"]

    target = lambda row: int(row[layer]["R_nonlocal"] != 0)
    key_0 = lambda row: (ranks(row)["0"],)
    key_2 = lambda row: (ranks(row)["0"], ranks(row)["2"])
    key_5 = lambda row: (ranks(row)["0"], ranks(row)["5"])
    key_25 = lambda row: (
        ranks(row)["0"],
        ranks(row)["2"],
        ranks(row)["5"],
        ranks(row)["25"],
        row[layer]["J_local"],
    )
    means = [
        conditional_expectations(rows, key_0, target),
        conditional_expectations(rows, key_2, target),
        conditional_expectations(rows, key_5, target),
        conditional_expectations(rows, key_25, target),
    ]
    contrasts = []
    for row in rows:
        m0 = means[0][key_0(row)]
        m2 = means[1][key_2(row)]
        m5 = means[2][key_5(row)]
        mf = means[3][key_25(row)]
        if mf != target(row):
            raise AssertionError("final partition sigma-algebra did not reveal Rc event")
        q2 = (m2 - m0) ** 2 + (mf - m2) ** 2
        q5 = (m5 - m0) ** 2 + (mf - m5) ** 2
        contrasts.append(q2 - q5)
    return contrasts


def exact_summary(values: list[Fraction]) -> dict:
    mean = sum(values, Fraction()) / len(values)
    return {
        "nonzero_configurations": sum(value != 0 for value in values),
        "mean_at_p_half": fraction_record(mean),
        "minimum": str(min(values)),
        "maximum": str(max(values)),
        "histogram": fraction_histogram(values),
        "configuration_order_sha256": digest_fractions(values),
    }


def compact_layer_record(record: dict, layer: str) -> dict:
    return {
        "partition_ranks": record[layer]["partition_ranks"],
        "ambient_H1_ranks": record[layer]["ambient_H1_ranks"],
        "J_local": record[layer]["J_local"],
        "R_nonlocal": record[layer]["R_nonlocal"],
    }


def render() -> dict:
    records = [configuration_record(mask) for mask in range(1 << 10)]
    activation = {
        layer: [activation_row(record, layer) for record in records]
        for layer in LAYERS
    }
    activation_values = {
        layer: [Fraction(row["C_activation"]) for row in activation[layer]]
        for layer in LAYERS
    }
    activation_even = [
        (black + white) / 2
        for black, white in zip(
            activation_values["black_NN"],
            activation_values["white_matching"],
        )
    ]
    activation_odd = [
        (black - white) / 2
        for black, white in zip(
            activation_values["black_NN"],
            activation_values["white_matching"],
        )
    ]

    doob = {layer: ambient_rank_doob(records, layer) for layer in LAYERS}
    doob_values = {
        layer: [row["C_Doob_Q"] for row in doob[layer]["rows"]]
        for layer in LAYERS
    }
    information_values = {
        layer: [row["R2_information_advantage"] for row in doob[layer]["rows"]]
        for layer in LAYERS
    }
    if any(sum(doob_values[layer], Fraction()) for layer in LAYERS):
        raise AssertionError("Doob quadratic-variation isometry failed")

    residual_no_go = {
        layer: partition_residual_doob(records, layer) for layer in LAYERS
    }
    if any(value for layer in LAYERS for value in residual_no_go[layer]):
        raise AssertionError("scalar Rc rank filtration unexpectedly became ordered")

    activation_mask = next(
        mask
        for mask in range(1 << 10)
        if mask.bit_count() == 5
        and (activation_even[mask] != 0 or activation_odd[mask] != 0)
    )
    doob_mask = next(
        mask
        for mask in range(1 << 10)
        if mask.bit_count() == 5
        and any(doob_values[layer][mask] != 0 for layer in LAYERS)
    )
    endpoint_collision = None
    seen: dict[tuple[int, int], tuple[int, Fraction]] = {}
    for mask, record in enumerate(records):
        ranks = record["black_NN"]["ambient_H1_ranks"]
        if ranks["25"] != 2:
            continue
        key = ranks["0"], ranks["25"]
        value = doob_values["black_NN"][mask]
        if key in seen and seen[key][1] != value:
            endpoint_collision = {
                "shared_endpoint_ranks_h0_h25": list(key),
                "first_mask": seen[key][0],
                "first_C_Doob_Q": str(seen[key][1]),
                "second_mask": mask,
                "second_C_Doob_Q": str(value),
            }
            break
        seen[key] = mask, value
    if endpoint_collision is None:
        raise AssertionError("Doob contrast collapsed to an endpoint function")

    final_hnf, _ = column_hnf(FINAL_MATRIX)
    hnf_2, _ = column_hnf(PERIOD_2.matrix)
    hnf_5, _ = column_hnf(PERIOD_5.matrix)

    activation_witness = records[activation_mask]
    doob_witness = records[doob_mask]
    return {
        "schema": "matching-one.p200-path-filtration-oracle.v1",
        "issue": 200,
        "status": "exact_tiny_path_ordered_witness_production_stopped",
        "geometry": {
            "lineage": "N1 -> N2/N5 -> N10 over Z[i]",
            "final_period_matrix": [list(row) for row in FINAL_MATRIX],
            "final_column_HNF": [list(row) for row in final_hnf],
            "N2_column_HNF": [list(row) for row in hnf_2],
            "N5_column_HNF": [list(row) for row in hnf_5],
            "coordinates_in_reference_engine_order": [
                list(value) for value in GEOMETRY.coordinates
            ],
            "lift_convention": "raw displacement between exact reference representatives; ambient-H1 path rows are convention-labelled",
        },
        "filtrations": {
            "target": "A_c=1{ambient H1 rank after R2 join R5 is two}",
            "R2_then_R5": "sigma(h0) subset sigma(h0,h2) subset sigma(h0,h2,h5,h25)",
            "R5_then_R2": "sigma(h0) subset sigma(h0,h5) subset sigma(h0,h2,h5,h25)",
            "common_endpoint": True,
            "activation_time": "first stage k in {0,1,2} at which the target rank two is visible; undefined off A_c",
            "C_activation": "tau(R5 then R2)-tau(R2 then R5)=A_c*(1{h2=2}-1{h5=2})",
            "C_Doob_Q": "sum_k (Delta E[A_c|F_k^(2,5)])^2 - sum_k (Delta E[A_c|F_k^(5,2)])^2",
            "R2_information_advantage": "Var(A_c|sigma(h0,h5))-Var(A_c|sigma(h0,h2)); positive means R2 reveals more at stage one",
        },
        "exhaustive_checks": {
            "configurations": len(records),
            "join_endpoint_equal_every_configuration_and_colour": True,
            "activation": {
                "black_NN": exact_summary(activation_values["black_NN"]),
                "white_matching": exact_summary(activation_values["white_matching"]),
                "typed_complement_even": exact_summary(activation_even),
                "typed_complement_odd": exact_summary(activation_odd),
            },
            "Doob_quadratic_variation": {
                "black_NN": exact_summary(doob_values["black_NN"]),
                "white_matching": exact_summary(doob_values["white_matching"]),
                "mean_zero_reason": "exact martingale isometry; the configurationwise contrast is nevertheless nonzero",
            },
            "conditional_variance_information": {
                layer: exact_summary(information_values[layer])
                for layer in LAYERS
            },
            "partition_Rc_scalar_rank_no_go": {
                "target": "1{R_nonlocal != 0}",
                "same_rank_summary_filtrations": True,
                "C_Doob_Q_zero_all_black": True,
                "C_Doob_Q_zero_all_white": True,
                "interpretation": "the natural scalar ranks do not chronologize Rc on N10; a nonzero Rc path mark needs richer intermediate state",
            },
        },
        "symmetry_typing": {
            "path_orientation_R2_swap_R5": "C_activation, C_Doob_Q, and R2_information_advantage are odd by definition",
            "typed_layer_swap": {
                "even": "(C_black_NN+C_white_matching)/2 is invariant",
                "odd": "(C_black_NN-C_white_matching)/2 changes sign",
            },
            "geometric_N650_orientation": "not inferred from the one-fiber N10 oracle; S/D needs both real N650 geometries or a second exact finite geometry",
        },
        "witnesses": {
            "balanced_first_H1_activation": {
                "mask": activation_mask,
                "occupied_sites": activation_witness["occupied_sites"],
                "black_NN": compact_layer_record(activation_witness, "black_NN"),
                "white_matching": compact_layer_record(
                    activation_witness, "white_matching"
                ),
                "black_path": activation["black_NN"][activation_mask],
                "white_path": activation["white_matching"][activation_mask],
                "C_typed_even": str(activation_even[activation_mask]),
                "C_typed_odd": str(activation_odd[activation_mask]),
            },
            "balanced_Doob_quadratic_variation": {
                "mask": doob_mask,
                "occupied_sites": doob_witness["occupied_sites"],
                "black_NN": compact_layer_record(doob_witness, "black_NN"),
                "white_matching": compact_layer_record(
                    doob_witness, "white_matching"
                ),
                "black_Doob": {
                    key: str(value)
                    for key, value in doob["black_NN"]["rows"][doob_mask].items()
                },
                "white_Doob": {
                    key: str(value)
                    for key, value in doob["white_matching"]["rows"][doob_mask].items()
                },
            },
            "same_endpoint_different_Doob_contrast": endpoint_collision,
        },
        "minimal_extra_data": {
            "for_H1_activation": "retain per-colour h0,h2,h5,h25 before batch aggregation; do not retain only the final join or mixed Delta25",
            "for_partition_Rc_chronology": "dynamic edge order or marked cluster lineage at the intermediate join; scalar r0,r2,r5,r25,J_local is exactly insufficient in this oracle",
        },
        "decision": {
            "exact_toy_result": "a genuine path-filtration witness exists even though final joins commute",
            "mechanism": "norm-five and norm-two intermediate ranks reveal the same final rank-two event at different stopping times",
            "production": "stopped; do not run the frozen 100M mixed-join job and do not reinterpret the Phase D Rc rejection as path order",
            "not_claimed": "continuum memory, Jordan structure, or a resolved N650 S/D channel",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = render()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
