#!/usr/bin/env python3
"""Exact Q->1 relative-cluster-fugacity formula for the matching observable."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

from euler_motif_controls import cluster_difference, configuration_identity
from integer_period_torus import IntegerTorusGeometry, gaussian_integer_torus


def fraction_record(value: Fraction) -> dict[str, int | str]:
    return {"numerator": value.numerator, "denominator": value.denominator, "text": str(value)}


def active_mask(mask: int, n: int) -> list[bool]:
    return [bool((mask >> vertex) & 1) for vertex in range(n)]


def plaquette_pattern_histogram(geometry: IntegerTorusGeometry, active: list[bool]) -> tuple[int, ...]:
    """Counts of all 16 oriented 2x2 binary patterns in covering coordinates."""
    histogram = [0] * 16
    for x, y in geometry.coordinates:
        vertices = (
            geometry.vertex((x, y)),
            geometry.vertex((x + 1, y)),
            geometry.vertex((x, y + 1)),
            geometry.vertex((x + 1, y + 1)),
        )
        pattern = sum(int(active[vertex]) << bit for bit, vertex in enumerate(vertices))
        histogram[pattern] += 1
    return tuple(histogram)


def enumerate_source(geometry: IntegerTorusGeometry) -> dict:
    counts: dict[int, list[int]] = {-1: [0] * (geometry.n + 1), 0: [0] * (geometry.n + 1), 1: [0] * (geometry.n + 1)}
    matching_coefficients = [0] * (geometry.n + 1)
    factorization_failures = 0
    channel_failures = 0
    for mask in range(1 << geometry.n):
        active = active_mask(mask, geometry.n)
        record = configuration_identity(geometry, active, mask)
        black_clusters, white_clusters = cluster_difference(geometry, active)
        local_euler = record.motifs["V"] - record.motifs["E"] + record.motifs["F0"]
        source_exponent = black_clusters - white_clusters - local_euler
        factorization_failures += source_exponent != record.q
        channel_failures += len(set(record.wrapping.values())) != 1
        occupied = record.motifs["V"]
        counts[record.q][occupied] += 1
        matching_coefficients[occupied] += record.q

    half_sector = {
        q: Fraction(sum(values), 1 << geometry.n) for q, values in counts.items()
    }
    half_matching = sum(Fraction(q) * value for q, value in half_sector.items())
    half_second = sum(Fraction(q * q) * value for q, value in half_sector.items()) - half_matching**2
    return {
        "geometry": geometry.name,
        "N": geometry.n,
        "configuration_count": 1 << geometry.n,
        "sector_Bernstein_coefficients_by_q": {str(q): values for q, values in counts.items()},
        "matching_Bernstein_coefficients": matching_coefficients,
        "factorization_failures": factorization_failures,
        "topology_channel_failures": channel_failures,
        "p_half": {
            "Z_minus": fraction_record(half_sector[-1]),
            "Z_zero": fraction_record(half_sector[0]),
            "Z_plus": fraction_record(half_sector[1]),
            "first_logQ_derivative": fraction_record(half_matching),
            "second_logQ_derivative": fraction_record(half_second),
        },
    }


def obstruction_witness() -> dict:
    geometry = gaussian_integer_torus(4, 3)
    masks = (0x1D24768, 0x0F6ACA0)
    rows = []
    for mask in masks:
        active = active_mask(mask, geometry.n)
        record = configuration_identity(geometry, active, mask)
        black_clusters, white_clusters = cluster_difference(geometry, active)
        rows.append({
            "mask_hex": hex(mask),
            "black_clusters": black_clusters,
            "white_matching_clusters": white_clusters,
            "q_matching": record.q,
            "plaquette_pattern_histogram": list(plaquette_pattern_histogram(geometry, active)),
            "Euler_motifs": {name: record.motifs[name] for name in ("V", "E", "F0")},
        })
    return {
        "geometry": geometry.name,
        "N": geometry.n,
        "declared_score_class": "a*k_black + sum_x f(oriented 2x2 plaquette pattern at x), for arbitrary a and f",
        "configurations": rows,
        "same_black_cluster_count": rows[0]["black_clusters"] == rows[1]["black_clusters"],
        "same_complete_plaquette_histogram": rows[0]["plaquette_pattern_histogram"] == rows[1]["plaquette_pattern_histogram"],
        "different_matching_charge": rows[0]["q_matching"] != rows[1]["q_matching"],
        "cause": "the complementary white matching-cluster count differs",
        "no_go": (
            "no ordinary one-colour FK cluster derivative plus any translation-invariant plaquette-local "
            "counterterm can equal q on every configuration"
        ),
    }


def build_oracle() -> dict:
    source = enumerate_source(gaussian_integer_torus(2, 1))
    witness = obstruction_witness()
    return {
        "schema": "matching-one.p114-relative-cluster-fugacity.v1",
        "issue": 114,
        "finite_partition_function": {
            "definition": (
                "Z_rel(Q,p)=sum_omega P_p(omega) "
                "Q^k_black Q^(-k_white) Q^[-(V-E+F0)]"
            ),
            "normalization": "Z_rel(1,p)=1",
            "local_counterterm_factorization": "Q^[-(V-E+F0)]=product_sites Q^-n product_NN_edges Q^(n_i n_j) product_faces Q^(-prod_corner n)",
            "topological_resolution": "Z_rel(Q,p)=Q Z_+(p)+Z_0(p)+Q^-1 Z_-(p)",
            "sector_definition": "Z_r(p)=sum_{omega:q(omega)=r} P_p(omega), r in {-1,0,+1}",
            "exact_derivatives": {
                "first": "(Q d_Q) log Z_rel | Q=1 = E[q] = M_N(p)",
                "second": "(Q d_Q)^2 log Z_rel | Q=1 = Var(q)",
            },
            "configuration_identity": "q=k_black-k_white-(V-E+F0)=I_black(H)-I_white_matching(H)",
            "typed_topology": "H may be any repository channel whose primal-minus-matching difference equals q on the declared quotient",
            "continuum_modular_gate": (
                "use cross/either for a basis-independent modular-scalar continuum source; finite direction/both "
                "identities do not promote those labels to modular scalars"
            ),
        },
        "tiny_exact_oracle": source,
        "ordinary_FK_local_obstruction": witness,
        "derivative_ledger": {
            "fixed_p_relative_Q_direction": [
                "black NN cluster-fugacity derivative +k_black",
                "white matching cluster-fugacity derivative -k_white",
                "explicit local Euler derivative -(V-E+F0)",
            ],
            "issue_258_relation": (
                "at fixed p and fixed topology semantics these are the complete source terms; moving along a "
                "Q-dependent critical manifold adds the separate measure-score covariance and must not be folded into this identity"
            ),
            "projector_derivative": (
                "not needed for the finite equality when q is used as a fixed sector charge; it becomes a separate "
                "term if the generic-Q homology projector itself is continued"
            ),
        },
        "structural_connections": {
            "issue_233": (
                "the insertion D=k_black-C k_white C-local is the first relative-Q pull-through obstruction; "
                "complement/matching exchange sends Q to Q^-1"
            ),
            "issue_123": (
                "the paired black/white connectivity state is not optional: local composition generates partial "
                "boundary partitions, so the source cannot close on the all-or-none Bernoulli tensor alone"
            ),
            "issue_120_257": (
                "the unmarked source is Potts-colour singlet; a charged [2] endpoint requires an additional character insertion"
            ),
        },
        "locality_status": {
            "positive": (
                "for every fixed transfer width the source has a finite paired-connectivity representation tracking "
                "black NN and white NN+NNN partitions plus local site/edge/face weights"
            ),
            "negative": witness["no_go"],
            "interpretation": (
                "matching is a local coupled/two-connectivity cluster source, but not an ordinary one-colour FK Q derivative"
            ),
        },
        "claim_boundary": {
            "proved": [
                "the finite partition-function derivative formula",
                "its topological-sector resolution",
                "the N5 exhaustive factorization",
                "the N25 plaquette-local one-colour FK obstruction",
            ],
            "not_proved": [
                "a width-independent finite-bond-dimension tensor representation",
                "a closed continuum CFT matrix element for the relative source",
                "that generic-Q homology projector derivatives vanish",
                "a no-go for local terms of arbitrarily large range",
            ],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    payload = json.dumps(build_oracle(), indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
