#!/usr/bin/env python3
"""Exact tiny oracle for the production marked-birth permutation path.

At pre-insertion size k, the next permutation site is uniform among the N-k
absent sites.  Multiplication by N-k is therefore an exact Horvitz estimator
of the absent-site insertion sum conditional on the occupied set.  This file
exhausts both sides of that identity at every k.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import itertools
import json
from math import comb, factorial
from pathlib import Path
from typing import Any, Iterable, Sequence

from digital_alexander_filtration_oracle import rank_mark
from homology_rank_birth_insertion import rank_birth_insertion
from integer_period_torus import axis_integer_torus, gaussian_integer_torus
from rank_birth_parity_channels import complement_pair, spin4_character


METRICS = (
    "q",
    "q2",
    "active_S",
    "active_D",
    "inactive_S",
    "inactive_D",
    "site_S",
    "site_D",
    "J_S_re",
    "J_S_im",
    "J_D_re",
    "J_D_im",
    "q_J_D_re",
    "q_J_D_im",
)


def _fraction(value: str | int | None) -> Fraction:
    return Fraction(0) if value is None else Fraction(value)


def insertion_values(geometry, active: Sequence[bool], vertex: int) -> dict[str, Fraction]:
    primal, matching, _ = complement_pair(geometry, active, vertex)
    rank_before = int(primal["rank_before"])
    active_s = int(primal["even"])
    active_d = int(primal["odd"])
    inactive_s = int(matching["even"])
    inactive_d = int(matching["odd"])
    if active_s != inactive_s or active_d != -inactive_d:
        raise AssertionError("active/inactive insertion source did not complement-pair")
    s_value = Fraction(active_s + inactive_s, 2)
    d_value = Fraction(active_d - inactive_d, 2)
    line = primal["ell"]
    if line is None:
        chi_re = chi_im = Fraction(0)
    else:
        chi_re, chi_im = spin4_character(geometry.periods.period_vector(line))
    q_value = rank_before - 1
    if primal["gate_01"] and primal["gate_12"]:
        if line is not None or d_value != 0 or s_value != 2:
            raise AssertionError("direct 0->2 did not map to null-line S=2,D=0")
    return {
        "q": Fraction(q_value),
        "q2": Fraction(q_value * q_value),
        "active_S": Fraction(active_s),
        "active_D": Fraction(active_d),
        "inactive_S": Fraction(inactive_s),
        "inactive_D": Fraction(inactive_d),
        "site_S": s_value,
        "site_D": d_value,
        "J_S_re": s_value * chi_re,
        "J_S_im": s_value * chi_im,
        "J_D_re": d_value * chi_re,
        "J_D_im": d_value * chi_im,
        "q_J_D_re": q_value * d_value * chi_re,
        "q_J_D_im": q_value * d_value * chi_im,
    }


def _zero() -> dict[str, Fraction]:
    return {name: Fraction(0) for name in METRICS}


def _add(target: dict[str, Fraction], source: dict[str, Fraction], factor: int = 1) -> None:
    for name in METRICS:
        target[name] += factor * source[name]


def _texts(values: dict[str, Fraction]) -> dict[str, str]:
    return {name: str(value) for name, value in values.items()}


def exact_site_sum(geometry, k: int) -> dict[str, Fraction]:
    total = _zero()
    vertices = tuple(range(geometry.n))
    for occupied in itertools.combinations(vertices, k):
        active = [False] * geometry.n
        for vertex in occupied:
            active[vertex] = True
        rank = rank_mark(geometry, active, matching=False)[0]
        q_value = rank - 1
        total["q"] += q_value
        total["q2"] += q_value * q_value
        for vertex in vertices:
            if not active[vertex]:
                values = insertion_values(geometry, active, vertex)
                for name in METRICS[2:]:
                    total[name] += values[name]
    denominator = comb(geometry.n, k)
    return {name: value / denominator for name, value in total.items()}


def exact_path_horvitz(geometry, k: int) -> dict[str, Fraction]:
    total = _zero()
    permutations = factorial(geometry.n)
    absent = geometry.n - k
    for order in itertools.permutations(range(geometry.n)):
        active = [False] * geometry.n
        for vertex in order[:k]:
            active[vertex] = True
        values = insertion_values(geometry, active, order[k])
        total["q"] += values["q"]
        total["q2"] += values["q2"]
        for name in METRICS[2:]:
            total[name] += absent * values[name]
    return {name: value / permutations for name, value in total.items()}


def summarize_geometry(name: str, geometry) -> dict[str, Any]:
    rows = []
    direct_mass = Fraction(0)
    lifted_frame_examples: set[tuple[str, str, str]] = set()
    for k in range(geometry.n):
        site = exact_site_sum(geometry, k)
        path = exact_path_horvitz(geometry, k)
        residual = {metric: path[metric] - site[metric] for metric in METRICS}
        if any(residual.values()):
            raise AssertionError(f"Horvitz identity failed for {name}, k={k}: {residual}")
        rows.append(
            {
                "k": k,
                "exact_absent_site_sum": _texts(site),
                "path_Horvitz_average": _texts(path),
                "residual": _texts(residual),
            }
        )

        # Record the direct-birth mass and examples proving that chi4 is made
        # from the physical lifted vector, not the period-coordinate line.
        for occupied in itertools.combinations(range(geometry.n), k):
            active = [False] * geometry.n
            for vertex in occupied:
                active[vertex] = True
            for vertex in range(geometry.n):
                if active[vertex]:
                    continue
                record = rank_birth_insertion(geometry, active, vertex)
                if record["gate_0_to_1"] and record["gate_1_to_2"]:
                    direct_mass += Fraction(1, comb(geometry.n, k))
                for birth in record["births"]:
                    if birth["ell"] is not None:
                        mark = birth["homology_h4"]
                        lifted_frame_examples.add(
                            (mark["physical_vector"], mark["cos4"], mark["sin4"])
                        )
    return {
        "id": name,
        "N": geometry.n,
        "period_matrix": [list(row) for row in geometry.periods.matrix],
        "microcanonical_rows": rows,
        "direct_0_to_2_absent_site_mass_summed_over_k": str(direct_mass),
        "lifted_Euclidean_chi4_examples": [
            {"physical_vector": vector, "cos4": real, "sin4": imaginary}
            for vector, real, imaginary in sorted(lifted_frame_examples)
        ],
    }


def build_artifact() -> dict[str, Any]:
    geometries = [
        summarize_geometry("axis-L2", axis_integer_torus(2)),
        summarize_geometry("gaussian-2-1", gaussian_integer_torus(2, 1)),
    ]
    if Fraction(geometries[0]["direct_0_to_2_absent_site_mass_summed_over_k"]) <= 0:
        raise AssertionError("axis direct-birth control has zero direct mass")
    if not any(
        row["physical_vector"] in {"2,1", "-2,-1"}
        and row["cos4"] == "-7/25"
        and row["sin4"] == "24/25"
        for row in geometries[1]["lifted_Euclidean_chi4_examples"]
    ):
        raise AssertionError("Gaussian lifted-frame chi4 certificate is missing")
    return {
        "schema": "matching-one/marked-birth-path-oracle/v1",
        "issues": [215, 269, 275, 276],
        "status": "tiny_exact_microcanonical_Horvitz_oracle",
        "identity": "E_path[(N-k) I(next site) | |A|=k] = E_A[sum_(v absent) I_v(A)]",
        "canonical_scorer_bridge": "multiply stored site sums by N/(N-k), then convolve against Bin(N-1,k)",
        "geometries": geometries,
        "exact_boundaries": {
            "full_source": "S=(S_active+S_inactive)/2 and D=(D_active-D_inactive)/2; every tiny insertion has S_active=S_inactive and D_active=-D_inactive",
            "direct_birth": "0->2 has ell=null, S=2 and D=0",
            "first_birth": "strict 0->1 uses the post-insertion rank-one line/index",
            "second_birth": "strict 1->2 uses the pre-insertion plateau line/index",
            "chi4_frame": "primitive lifted Euclidean P*ell",
            "production_iota": "C++ runner preserves raw winding coefficients before primitive reduction; tiny geometric examples happen to have iota=1",
        },
    }


def render_markdown(artifact: dict[str, Any]) -> str:
    lines = [
        "# Exact marked-birth permutation-path oracle",
        "",
        "For every microcanonical size `k`, exhaustive permutation paths agree",
        "coefficient-for-coefficient with exhaustive absent-site insertion sums:",
        "",
        "```text",
        "E[(N-k) I(next site)] = E[sum_(v absent) I_v].",
        "```",
        "",
        "The checked vector is `(q,q2,S_active/D_active,S_inactive/D_inactive,",
        "S_full,D_full,chi4*S,chi4*D,q*chi4*D)` with both",
        "real and imaginary spin-four components. All residuals in the JSON are",
        "exactly zero. Axis `L=2` supplies nonzero direct `0->2` mass and verifies",
        "`ell=null,D=0,S=2`. Gaussian `(2,1)` verifies that `chi4` uses the lifted",
        "Euclidean direction `P*ell`; the `(2,1)` line gives `(-7+24i)/25`.",
        "",
        "The production path stores Horvitz absent-site sums. A canonical Russo",
        "score multiplies them by `N/(N-k)` and convolves with `Bin(N-1,k)`.",
        "For the common-field product, the stored root-deleted `q_before*J_D` is",
        "completed to full-configuration `q*J_D` by adding `p*J_D`, because",
        "`q_after=q_before+S` and `S*J_D=J_D`.",
        "No random-root replay and no per-step all-site scan is required.",
        "",
        "## Frozen production schema",
        "",
        "The sparse row is one aligned-batch histogram of",
        "`(K1,K2,site01,site12,ell,iota01,iota12,P*ell,local-H4 marks)`.",
        "The per-`k` path retains `q,q2`, both active-primal and inactive-matching",
        "reverse gates, their canonical full source, `chi4(ell)S/D`, `q*chi4(ell)D`,",
        "and local-H4 `S/D`. The full matching-function source is",
        "`S_full=(S_active+S_inactive)/2` and",
        "`D_full=(D_active-D_inactive)/2`.",
        "",
        "The raw sides stay in the CSV even though exact complement pairing makes",
        "`S_active=S_inactive` and `D_active=-D_inactive`. The angular mark uses the",
        "primitive physical lift `P*ell`, never the period-coordinate line as an",
        "angle. The saturation index is the gcd of raw same-line winding",
        "coefficients before primitive reduction; rank two stores zero.",
        "",
        "The prerevealed discovery pilot is `N=65` and its exact q2 child `N=130`",
        "at 20,000 samples each, plus max-leverage `P50 N=145` at 10,000 samples.",
        "",
        "## Scientific card",
        "",
        "1. MECHANISM SPACE: separates active/inactive topology, even/odd birth sources, line polarization, and local landing geometry.",
        "2. NOT PROVED: a finite pilot cannot identify `Q4 epsilon` or establish `x=21/4` scaling.",
        "3. OBSERVER-SECTOR-SOURCE-GEOMETRY: `A_top` | Alexander odd/even | `chi4(ell)D/S` and local H4 | Gaussian orientation pairs.",
        "4. DEPENDENCY GROUP: all scores and covariances reuse the same counter-coupled permutation batches.",
        "5. UPWEIGHT OBSERVATION: complement-clean q2 sign/phase transfer of connected `q-J_D4` and `gamma_D4`.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    artifact = build_artifact()
    text = (
        json.dumps(artifact, indent=2, sort_keys=True) + "\n"
        if args.format == "json"
        else render_markdown(artifact)
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
