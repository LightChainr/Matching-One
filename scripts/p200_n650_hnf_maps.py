#!/usr/bin/env python3
"""Exact HNF quotient maps for the two Issue #200 N650 lineages."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gaussian_crt_commutator import (  # noqa: E402
    ALPHA,
    BETA,
    SOURCES,
    gaussian_multiply,
    multiplication_matrix,
)
from integer_period_torus import IntegerPeriods  # noqa: E402


Gaussian = tuple[int, int]
Matrix = tuple[tuple[int, int], tuple[int, int]]


def extended_gcd(a: int, b: int) -> tuple[int, int, int]:
    old_r, r = abs(a), abs(b)
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t
    return old_r, (-old_s if a < 0 else old_s), (-old_t if b < 0 else old_t)


def column_hnf(matrix: Matrix) -> tuple[Matrix, dict[str, int]]:
    """Match threshold_rank_integer_period_mc.cpp's upper column-HNF."""

    (a, b), (c, d) = matrix
    order = abs(a * d - b * c)
    h22, bezout_c, bezout_d = extended_gcd(c, d)
    h11 = order // h22
    h12 = (a * bezout_c + b * bezout_d) % h11
    hnf = ((h11, h12), (0, h22))
    original = IntegerPeriods(matrix)
    reduced = IntegerPeriods(hnf)
    for column in ((hnf[0][0], hnf[1][0]), (hnf[0][1], hnf[1][1])):
        original.winding(column)
    for column in ((a, c), (b, d)):
        reduced.winding(column)
    return hnf, {"h11": h11, "h12": h12, "h22": h22, "order": order}


def hnf_label(point: tuple[int, int], hnf: dict[str, int]) -> int:
    x, y = point
    quotient_y = (y - (y % hnf["h22"])) // hnf["h22"]
    ry = y - quotient_y * hnf["h22"]
    rx = (x - quotient_y * hnf["h12"]) % hnf["h11"]
    return rx + hnf["h11"] * ry


def hnf_representative(label: int, hnf: dict[str, int]) -> tuple[int, int]:
    return label % hnf["h11"], label // hnf["h11"]


def digest(values: list) -> str:
    encoded = json.dumps(values, separators=(",", ":")).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def lattice_record(value: Gaussian) -> dict:
    matrix = multiplication_matrix(value)
    hnf_matrix, hnf = column_hnf(matrix)
    return {
        "gaussian": list(value),
        "period_matrix": [list(row) for row in matrix],
        "column_HNF": [list(row) for row in hnf_matrix],
        "hnf_parameters": hnf,
    }


def projection_map(source: Gaussian) -> dict:
    after_2 = gaussian_multiply(source, ALPHA)
    after_5 = gaussian_multiply(source, BETA)
    final = gaussian_multiply(after_2, BETA)
    records = {
        "source_N65": lattice_record(source),
        "intermediate_N130": lattice_record(after_2),
        "intermediate_N325": lattice_record(after_5),
        "final_N650": lattice_record(final),
    }
    final_hnf = records["final_N650"]["hnf_parameters"]
    hnf_130 = records["intermediate_N130"]["hnf_parameters"]
    hnf_325 = records["intermediate_N325"]["hnf_parameters"]
    source_hnf = records["source_N65"]["hnf_parameters"]
    map_650_to_130, map_650_to_325, map_650_to_65 = [], [], []
    for label in range(650):
        representative = hnf_representative(label, final_hnf)
        map_650_to_130.append(hnf_label(representative, hnf_130))
        map_650_to_325.append(hnf_label(representative, hnf_325))
        map_650_to_65.append(hnf_label(representative, source_hnf))
    if set(Counter(map_650_to_130).values()) != {5}:
        raise AssertionError("N650 -> N130 fibers are not all size five")
    if set(Counter(map_650_to_325).values()) != {2}:
        raise AssertionError("N650 -> N325 fibers are not all size two")
    if set(Counter(map_650_to_65).values()) != {10}:
        raise AssertionError("N650 -> N65 fibers are not all size ten")

    map_130_to_65 = [hnf_label(hnf_representative(i, hnf_130), source_hnf) for i in range(130)]
    map_325_to_65 = [hnf_label(hnf_representative(i, hnf_325), source_hnf) for i in range(325)]
    if [map_130_to_65[i] for i in map_650_to_130] != map_650_to_65:
        raise AssertionError("N130 composition failed")
    if [map_325_to_65[i] for i in map_650_to_325] != map_650_to_65:
        raise AssertionError("N325 composition failed")
    for source_label in range(65):
        pairs = {
            (map_650_to_130[i], map_650_to_325[i])
            for i in range(650)
            if map_650_to_65[i] == source_label
        }
        if len(pairs) != 10:
            raise AssertionError("actual HNF labels failed the CRT pair gate")
    return {
        "lattices": records,
        "map_cardinalities": {
            "N650_to_N130": {"targets": 130, "fiber_size": 5},
            "N650_to_N325": {"targets": 325, "fiber_size": 2},
            "N650_to_N65": {"targets": 65, "fiber_size": 10},
        },
        "map_sha256_in_final_HNF_label_order": {
            "N650_to_N130": digest(map_650_to_130),
            "N650_to_N325": digest(map_650_to_325),
            "N650_to_N65": digest(map_650_to_65),
            "paired_N130_N325": digest(list(map(list, zip(map_650_to_130, map_650_to_325)))),
        },
        "first_16_map_rows": [
            {"final_label": i, "N130_label": map_650_to_130[i], "N325_label": map_650_to_325[i], "N65_label": map_650_to_65[i]}
            for i in range(16)
        ],
        "composition_to_N65_equal_by_both_paths": True,
        "N130_N325_label_pair_is_bijective_on_each_N65_fiber": True,
    }


def render() -> dict:
    return {
        "schema": "matching-one.p200-n650-hnf-map.v1",
        "issue": 200,
        "base_exact_result": "results/exact-cover-character-oracles/n650_gaussian_crt_commutator.json",
        "status": "exact_configuration_label_operationalization",
        "lineages": [projection_map(source) for source in SOURCES],
        "configuration_semantics": {
            "well_typed_per_final_site": ["N650 HNF label", "N130 quotient label", "N325 quotient label", "N65 quotient label"],
            "not_well_typed": "subtracting N130 and N325 labels or treating them as duplicate endpoint samples",
            "linear_or_full_join_order_mark": "exact implementation null",
        },
        "mixed_join_acquisition_gate": {
            "canonical_algebra": "Delta25 h = h(Pi join R2 join R5)-h(Pi join R2)-h(Pi join R5)+h(Pi)",
            "required_same_configuration_rows": ["h0", "h2", "h5", "h25"],
            "required_batch_output": "four means plus full 4x4 covariance before forming Delta25",
            "current_threshold_histogram_sufficient": False,
            "current_runner_obstruction": "rank histograms discard the typed connectivity partition and lifted homology state needed to apply both joins",
            "q_specific_obstruction": "ordinary binary pushdown does not preserve black/white complementarity, so matching charge q needs an explicit typed two-colour transport",
            "implementation_decision": "no C++ production extension until h and its typed lift transport are frozen",
        },
        "high_information_candidate": {
            "name": "typed two-colour mixed homology defect",
            "definition": "h_R=(rank_H1(black_NN join R)-rank_H1(white_matching join R))/2 in the common parent-period basis; Delta25=h_R2R5-h_R2-h_R5+h_empty",
            "configuration_level_construction": [
                "build black NN and white matching connectivity separately on the final lift",
                "for each R in {empty,R2,R5,R2 join R5}, add only same-colour fiber-identification edges with exact lifted deck displacement",
                "measure both ambient winding ranks before discarding either DSU",
                "emit four correlated rows and their full covariance",
            ],
            "exact_symmetry_target": "complement swaps the two typed layers and negates every h_R and Delta25",
            "status": "exploratory mechanism; requires a tiny exact complement/winding oracle before C++ production",
            "why_not_ordinary_q": "mixed-colour fibers can occupy both quotient layers, so this is a typed defect observable rather than the standard matching charge of one binary quotient mask",
        },
        "evidence_boundary": {
            "exact": "both real N650 HNF quotient squares commute and have fiber sizes 2,5,10",
            "exact_toy": "the symmetric partition-rank Delta25 equals -4 in the base artifact",
            "open": "whether typed ambient-H1 or matching charge has a nonzero Delta25 on percolation configurations",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = render()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
