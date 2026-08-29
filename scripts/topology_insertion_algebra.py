#!/usr/bin/env python3
"""Exact finite insertion algebra and all-order topology-observer no-go.

For q=r-1 in {-1,0,1}, insertion has S=Delta q and D=Delta(q^2).
Every function of q reduces to a quadratic.  On the support of a nonzero
insertion source, q_before=-(S-D)/2 and q_after=(S+D)/2.  The resulting
finite algebra closes every connected coupling of f(q) to a marked gate
source in terms of marked gate means.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import itertools
import json
from pathlib import Path
from typing import Any, Callable, Optional

from digital_alexander_filtration_oracle import rank_mark
from integer_period_torus import axis_integer_torus, gaussian_integer_torus
from rank_birth_parity_channels import _gate_record, spin4_character


Matrix = list[list[str]]


TRANSITIONS = (
    {"name": "inactive_minus", "q0": -1, "q1": -1, "S": 0, "D": 0},
    {"name": "inactive_zero", "q0": 0, "q1": 0, "S": 0, "D": 0},
    {"name": "inactive_plus", "q0": 1, "q1": 1, "S": 0, "D": 0},
    {"name": "strict_01", "q0": -1, "q1": 0, "S": 1, "D": -1},
    {"name": "strict_12", "q0": 0, "q1": 1, "S": 1, "D": 1},
    {"name": "direct_02", "q0": -1, "q1": 1, "S": 2, "D": 0},
)


def _fraction(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _mask(mask: int, n: int) -> list[bool]:
    return [bool(mask & (1 << vertex)) for vertex in range(n)]


def _gate_monomials(s: int, d: int) -> dict[str, int]:
    return {"S": s, "D": d, "S2": s * s, "SD": s * d, "D2": d * d}


def truth_table() -> list[dict[str, Any]]:
    rows = []
    for raw in TRANSITIONS:
        row = dict(raw)
        s, d = row["S"], row["D"]
        row.update(
            {
                "I01": Fraction(s - d, 2),
                "I12": Fraction(s + d, 2),
                "strict01_projector": Fraction(d * d - d, 2),
                "strict12_projector": Fraction(d * d + d, 2),
                "direct02_projector": Fraction(s * s - d * d, 4),
                "Delta_q_check": row["q1"] - row["q0"],
                "Delta_q2_check": row["q1"] ** 2 - row["q0"] ** 2,
            }
        )
        rows.append({key: _fraction(value) if isinstance(value, Fraction) else value for key, value in row.items()})
    return rows


def algebra_matrices() -> dict[str, Any]:
    # Matrices use column convention: column j is multiplication applied to
    # basis element j. Gate basis is (1,S,D,D^2).
    q_matrix: Matrix = [
        ["0", "0", "0"],
        ["1", "0", "1"],
        ["0", "1", "0"],
    ]
    s_matrix: Matrix = [
        ["0", "0", "0", "0"],
        ["1", "2", "0", "0"],
        ["0", "0", "1", "0"],
        ["0", "-1", "0", "1"],
    ]
    d_matrix: Matrix = [
        ["0", "0", "0", "0"],
        ["0", "0", "0", "0"],
        ["1", "1", "0", "1"],
        ["0", "0", "1", "0"],
    ]
    delta_matrix: Matrix = [
        ["0", "0", "0"],
        ["0", "1", "0"],
        ["0", "0", "1"],
        ["0", "0", "0"],
    ]
    bernoulli_restriction: Matrix = [
        ["1", "0", "0"],
        ["0", "p-1/2", "1/2"],
        ["0", "1/2", "p-1/2"],
        ["0", "0", "0"],
    ]
    return {
        "matrix_convention": "columns are images of basis vectors",
        "observer_basis": ["1", "q", "q^2"],
        "Q_multiplication": q_matrix,
        "Q_minimal_polynomial": "Q*(Q-1)*(Q+1)=0",
        "gate_basis": ["1", "S", "D", "D^2"],
        "S_multiplication": s_matrix,
        "D_multiplication": d_matrix,
        "S_minimal_polynomial": "S*(S-1)*(S-2)=0",
        "D_minimal_polynomial": "D*(D-1)*(D+1)=0",
        "joint_relations": ["S*D=D", "D^2=2*S-S^2", "S*D^2=D^2", "D^3=D"],
        "insertion_difference_observer_to_gate": delta_matrix,
        "insertion_difference_statement": "Delta(a+b*q+c*q^2)=b*S+c*D",
        "bernoulli_root_restriction_observer_to_gate": bernoulli_restriction,
        "restriction_statement": "E_x[f(q_x)|environment]=a+(b*(p-1/2)+c/2)S+(b/2+c*(p-1/2))D on insertion support",
        "endpoint_multiplication_on_gate_ideal": {
            "Q_before": "(-S+D)/2",
            "Q_after": "(S+D)/2",
            "Q_before_minimal_polynomial": "Q_before*(Q_before+1)=0",
            "Q_after_minimal_polynomial": "Q_after*(Q_after-1)=0",
        },
    }


def _line_mark(geometry, record: dict[str, Any]) -> tuple[Fraction, Fraction]:
    if record["ell"] is None:
        return Fraction(0), Fraction(0)
    return spin4_character(geometry.periods.period_vector(record["ell"]))


def tiny_geometry_oracle(name: str, geometry) -> dict[str, Any]:
    states = [rank_mark(geometry, _mask(mask, geometry.n), matching=False) for mask in range(1 << geometry.n)]
    checked = 0
    transition_counts: dict[str, int] = {}
    failures: list[dict[str, Any]] = []
    for vertex in range(geometry.n):
        bit = 1 << vertex
        for environment in range(1 << geometry.n):
            if environment & bit:
                continue
            record = _gate_record(states[environment], states[environment | bit])
            q0, q1 = int(record["rank_before"]) - 1, int(record["rank_after"]) - 1
            s, d = int(record["even"]), int(record["odd"])
            transition_counts[f"{q0}->{q1}"] = transition_counts.get(f"{q0}->{q1}", 0) + 1
            marks = {"unmarked": (Fraction(1), Fraction(0)), "line_chi4": _line_mark(geometry, record)}
            for p in (Fraction(2, 5), Fraction(3, 5)):
                t = p - Fraction(1, 2)
                for mark_name, mark_pair in marks.items():
                    for component, mark in zip(("re", "im"), mark_pair):
                        for gate_name, gate in _gate_monomials(s, d).items():
                            for observer, f0, f1, rhs_factor in (
                                ("q", q0, q1, t * s + Fraction(d, 2)),
                                ("q2", q0 * q0, q1 * q1, Fraction(s, 2) + t * d),
                            ):
                                lhs = ((1 - p) * f0 + p * f1) * mark * gate
                                rhs = rhs_factor * mark * gate
                                checked += 1
                                if lhs != rhs:
                                    failures.append(
                                        {
                                            "vertex": vertex,
                                            "environment": environment,
                                            "p": _fraction(p),
                                            "mark": mark_name + "_" + component,
                                            "gate": gate_name,
                                            "observer": observer,
                                            "lhs": _fraction(lhs),
                                            "rhs": _fraction(rhs),
                                        }
                                    )
    if failures:
        raise AssertionError(f"finite insertion algebra failed for {name}: {failures[:2]}")
    return {
        "id": name,
        "N": geometry.n,
        "period_matrix": [list(row) for row in geometry.periods.matrix],
        "transition_counts": dict(sorted(transition_counts.items())),
        "exact_pointwise_checks": checked,
        "failures": failures,
    }


def build_artifact() -> dict[str, Any]:
    geometries = [
        tiny_geometry_oracle("axis-L2", axis_integer_torus(2)),
        tiny_geometry_oracle("gaussian-2-1", gaussian_integer_torus(2, 1)),
    ]
    return {
        "schema": "matching-one/topology-insertion-algebra/v1",
        "issues": [215, 258, 275],
        "status": "exact_all_order_topology_observer_no_go",
        "truth_table": truth_table(),
        "finite_algebra": algebra_matrices(),
        "theorem": {
            "observer": "every f:{-1,0,1}->C has unique f(q)=a+b*q+c*q^2",
            "marked_source": "X=H(environment)*G(S,D), with H independent of root occupation and G(0,0)=0",
            "raw_closure": "E[f(q)X]=a*mu_G+(b*(p-1/2)+c/2)*mu_SG+(b/2+c*(p-1/2))*mu_DG",
            "connected_closure": "Cov(f(q),X)=b*((p-1/2)mu_SG+mu_DG/2-<q>mu_G)+c*(mu_SG/2+(p-1/2)mu_DG-<q^2>mu_G)",
            "all_order_reason": "Q*(Q-1)*(Q+1)=0 reduces arbitrary polynomial, analytic, or lookup-table f to degree at most two",
            "no_go": "topology-only f(q) cannot supply an independent field matrix element against any root-occupation-independent marked insertion source",
        },
        "tiny_exact_oracle": geometries,
        "observer_frontier": {
            "current_fields_sufficient_for": [
                "all f(q) couplings to unmarked S/D",
                "all f(q) couplings to lifted-line chi4 S/D",
                "all f(q) couplings to same-root local-landing H4 S/D",
                "contact subtraction using q,q2,J_S,J_D and q_before*J_D",
            ],
            "external_classes": {
                "local_landing_observer": "landing or arm/motif value at a separately sampled site or separation, not merely a mark multiplying the same root gate",
                "bulk_Betti_observer": "component/cycle counts or Euler/Betti density that varies within a fixed ambient-rank q sector",
                "charged_seam_observer": "signed/character-valued seam crossing, winding charge, or defect partition observable resolving states within the same q sector",
                "modular_or_shape_observer": "independent stress/KdV response or wrapping channel not determined by ambient rank alone",
            },
            "minimal_missing_statistics": [
                "O_ext and O_ext^2 per microcanonical k/batch",
                "O_ext*J_D4 and O_ext*J_S4 on the same configuration",
                "for field Gram: two-root J_D4*conj(J_S4) and |J_S4|^2",
                "root/separation metadata when O_ext is local",
                "seam charge and orientation phase when O_ext is charged",
            ],
        },
        "claim_boundary": {
            "exact": "finite state/gate algebra, all-order closure, and tiny marked geometry checks",
            "not_ruled_out": "mean insertion sources themselves, couplings to observers varying inside a fixed q sector, or continuum identification from independent modular data",
        },
        "scientific_card": [
            "MECHANISM SPACE: every q-only connected insertion response collapses to a four-dimensional gate algebra.",
            "NOT PROVED: the theorem says nothing about observers that resolve configurations within one q sector.",
            "OBSERVER-SECTOR-SOURCE-GEOMETRY: f(q) | ambient-rank topology | arbitrary marked S/D gate polynomial | any finite torus.",
            "DEPENDENCY GROUP: previous A_top/J_D and q2 orthogonal scores are corollaries, not independent evidence.",
            "UPWEIGHT OBSERVATION: an independent bulk, separated-local, or charged-seam observer with same-field source products.",
        ],
    }


def render_markdown(artifact: dict[str, Any]) -> str:
    checked = sum(item["exact_pointwise_checks"] for item in artifact["tiny_exact_oracle"])
    return "\n".join(
        [
            "# All-order topology-insertion no-go",
            "",
            "For `q in {-1,0,1}`, every topology-only observer is exactly",
            "`f(q)=a+bq+cq^2`. Insertion defines `S=Delta q`, `D=Delta(q^2)`,",
            "and the gate algebra has basis `(1,S,D,D^2)` with",
            "",
            "```text",
            "S D = D,       D^2 = 2S-S^2,",
            "S(S-1)(S-2)=0, D(D-1)(D+1)=0.",
            "```",
            "",
            "For any root-occupation-independent mark `H` and insertion polynomial",
            "`G` with `G(0,0)=0`,",
            "write `mu_G=<HG>`, `mu_SG=<HSG>`, `mu_DG=<HDG>`. Then",
            "",
            "```text",
            "Cov(f(q),HG) =",
            " b[(p-1/2)mu_SG + mu_DG/2 - <q>mu_G]",
            "+c[mu_SG/2 + (p-1/2)mu_DG - <q^2>mu_G].",
            "```",
            "",
            "This is an all-order no-go: no lookup table, nonlinear function, or",
            "higher polynomial of ambient rank can create an independent insertion",
            "matrix element. The minimal polynomial `Q(Q-1)(Q+1)=0` reduces them all.",
            "",
            f"The tiny axis/Gaussian oracle performs {checked} exact marked pointwise checks with zero failures.",
            "",
            "The escape is precise: use an observer varying within a fixed q sector—",
            "a separated local landing/arm observable, bulk Betti/Euler data, a charged",
            "seam or winding character, or an independent modular/stress response.",
            "",
            "## Scientific card",
            "",
            *[f"{i}. {line}" for i, line in enumerate(artifact["scientific_card"], 1)],
            "",
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args()
    artifact = build_artifact()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_markdown(artifact), encoding="utf-8")


if __name__ == "__main__":
    main()
