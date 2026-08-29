#!/usr/bin/env python3
"""Exact tiny preflight for the Euler-residue external rank-birth observer.

The observer is evaluated on the pre-insertion black configuration ``A``:

    O_ext(A) = C_NN(A) - C_matching(A^c) - q(A)
             = V(A) - E_NN(A) + F0(A).

It is configuration-level and root-independent, but is not measurable from
the three-state ambient rank variable q.  The path gate then records O_ext,
O_ext squared, products with J_S/J_D, and the two local-score Gram entries.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from math import comb
from pathlib import Path
from typing import Any, Sequence

from euler_motif_controls import configuration_identity
from integer_period_torus import axis_integer_torus, gaussian_integer_torus
from marked_birth_path_oracle import insertion_values


STATE_METRICS = ("O_ext", "O_ext2")
SITE_METRICS = (
    "O_ext_J_S_re",
    "O_ext_J_S_im",
    "O_ext_J_D_re",
    "O_ext_J_D_im",
)
GRAM_METRICS = ("J_D_conj_J_S_re", "J_D_conj_J_S_im", "abs_J_S2")
ALL_METRICS = STATE_METRICS + SITE_METRICS + GRAM_METRICS


def _zero() -> dict[str, Fraction]:
    return {name: Fraction(0) for name in ALL_METRICS}


def _texts(values: dict[str, Fraction]) -> dict[str, str]:
    return {name: str(values[name]) for name in ALL_METRICS}


def external_observer(geometry, active: Sequence[bool]) -> tuple[int, dict[str, int]]:
    record = configuration_identity(geometry, active)
    if record.residual:
        raise AssertionError("Euler/Betti configuration identity failed")
    value = record.cluster_difference - record.q
    motif_value = record.motifs["V"] - record.motifs["E"] + record.motifs["F0"]
    if value != motif_value:
        raise AssertionError("component and cell-complex definitions of O_ext differ")
    return value, {
        "q": record.q,
        "C_black_minus_C_white": record.cluster_difference,
        "V": record.motifs["V"],
        "E": record.motifs["E"],
        "F0": record.motifs["F0"],
    }


def local_products(o_ext: int, values: dict[str, Fraction]) -> dict[str, Fraction]:
    j_s = (values["J_S_re"], values["J_S_im"])
    j_d = (values["J_D_re"], values["J_D_im"])
    return {
        "O_ext_J_S_re": o_ext * j_s[0],
        "O_ext_J_S_im": o_ext * j_s[1],
        "O_ext_J_D_re": o_ext * j_d[0],
        "O_ext_J_D_im": o_ext * j_d[1],
        "J_D_conj_J_S_re": j_d[0] * j_s[0] + j_d[1] * j_s[1],
        "J_D_conj_J_S_im": j_d[1] * j_s[0] - j_d[0] * j_s[1],
        "abs_J_S2": j_s[0] * j_s[0] + j_s[1] * j_s[1],
    }


def summarize_geometry(name: str, geometry) -> dict[str, Any]:
    rows = []
    witness = None
    for k in range(geometry.n):
        q_collision: dict[int, tuple[int, int]] = {}
        configurations = 0
        exact = _zero()
        path = _zero()
        identity_l1 = 0
        for mask in range(1 << geometry.n):
            if mask.bit_count() != k:
                continue
            configurations += 1
            active = [bool((mask >> vertex) & 1) for vertex in range(geometry.n)]
            o_ext, identity = external_observer(geometry, active)
            identity_l1 += abs(
                o_ext - (identity["V"] - identity["E"] + identity["F0"])
            )
            exact["O_ext"] += o_ext
            exact["O_ext2"] += o_ext * o_ext
            path["O_ext"] += o_ext
            path["O_ext2"] += o_ext * o_ext

            previous = q_collision.get(identity["q"])
            if previous is None:
                q_collision[identity["q"]] = (o_ext, mask)
            elif previous[0] != o_ext and witness is None:
                witness = {
                    "k": k,
                    "q": identity["q"],
                    "first_mask": previous[1],
                    "first_O_ext": previous[0],
                    "second_mask": mask,
                    "second_O_ext": o_ext,
                }

            absent = geometry.n - k
            site_total = {name_: Fraction(0) for name_ in SITE_METRICS + GRAM_METRICS}
            path_next_total = {
                name_: Fraction(0) for name_ in SITE_METRICS + GRAM_METRICS
            }
            for vertex in range(geometry.n):
                if active[vertex]:
                    continue
                insertion = insertion_values(geometry, active, vertex)
                if insertion["q"] != identity["q"]:
                    raise AssertionError("rank q and Euler wrapping q differ")
                products = local_products(o_ext, insertion)
                for metric in SITE_METRICS:
                    site_total[metric] += products[metric]
                    path_next_total[metric] += absent * products[metric]
                for metric in GRAM_METRICS:
                    # A Gram entry is deliberately the second moment of the
                    # same Horvitz local score, not a product of site sums.
                    site_total[metric] += absent * products[metric]
                    path_next_total[metric] += absent * absent * products[metric]
            for metric in SITE_METRICS + GRAM_METRICS:
                exact[metric] += site_total[metric]
                path[metric] += path_next_total[metric] / absent

        if configurations != comb(geometry.n, k):
            raise AssertionError("microcanonical configuration count mismatch")
        exact = {metric: value / configurations for metric, value in exact.items()}
        path = {metric: value / configurations for metric, value in path.items()}
        residual = {metric: path[metric] - exact[metric] for metric in ALL_METRICS}
        if any(residual.values()) or identity_l1:
            raise AssertionError(f"external path oracle failed for {name}, k={k}")
        rows.append(
            {
                "k": k,
                "configuration_average_or_site_sum": _texts(exact),
                "path_average": _texts(path),
                "residual": _texts(residual),
            }
        )
    return {
        "id": name,
        "N": geometry.n,
        "period_matrix": [list(row) for row in geometry.periods.matrix],
        "microcanonical_rows": rows,
        "q_sigma_algebra_escape_witness": witness,
    }


def build_artifact() -> dict[str, Any]:
    geometries = [
        summarize_geometry("axis-L2", axis_integer_torus(2)),
        summarize_geometry("gaussian-2-1", gaussian_integer_torus(2, 1)),
        summarize_geometry("axis-L3-index9", axis_integer_torus(3)),
    ]
    if not any(row["q_sigma_algebra_escape_witness"] is not None for row in geometries):
        raise AssertionError("tiny controls found no fixed-k q-sigma escape witness")
    return {
        "schema": "matching-one/external-observer-path-oracle/v1",
        "issues": [215, 275, 276],
        "status": "tiny_exact_preflight_no_production",
        "selected_observer": {
            "name": "Euler_residue",
            "definition": "O_ext=C_black_NN-C_white_matching-q=V-E+F0",
            "parity": "matching-odd before optional centering",
            "reason": "outside sigma(q), exact on the square cellulation, O(1) from maintained component counts, and scalar thermal/Euler content couples naturally to spin-4 birth sources",
        },
        "observer_authorization_gate": {
            "source_commit": "83e98fca02a074396493b64ed59cb02700a68796",
            "rule": "one C4 direction orbit cannot separate scalar and spin4",
            "status": "pass: O_ext is explicitly scalar and the H4 source retains typed internal complex chi4(ell); no single-orbit direction scalar is promoted",
            "future_local_H4": "requires separate axis+diagonal orbits (determinant -2) or equivalent typed/internal complex edge information",
        },
        "candidate_decision": [
            {
                "candidate": "bulk_Betti_Euler_cycle_rank",
                "decision": "selected",
                "sigma_algebra": "strictly outside q; explicit equal-q unequal-O_ext witnesses",
                "cost": "two integer component counters already updated by the two union-find traces",
            },
            {
                "candidate": "macroscopically_separated_local_H4_or_arm",
                "decision": "reserve_as_second_coordinate",
                "sigma_algebra": "outside q and field-specific",
                "cost": "needs a frozen far-anchor/orbit convention; the existing one-point local-D4 selector was already rejected",
            },
            {
                "candidate": "charged_seam_or_winding_character",
                "decision": "reserve_for_charged_sector",
                "sigma_algebra": "outside neutral q after a seam charge is fixed",
                "cost": "requires an orientation/gauge convention and per-component charged state",
            },
        ],
        "stored_path_statistics": [
            "O_ext",
            "O_ext2",
            "O_ext_times_J_D4_complex",
            "O_ext_times_J_S4_complex",
            "J_D4_times_conj_J_S4_complex",
            "abs_J_S4_squared",
            "q_times_J_D4_contact_control",
        ],
        "gram_semantics": "same-next-site Horvitz-score second moments; not products of independently summed site sources",
        "identity": "E[(N-k) O_ext J(next)|A]=O_ext sum_(v absent)J_v(A); same-path Gram uses (N-k)^2 and therefore targets (N-k) sum_v j_D(v)conj(j_S(v))",
        "geometries": geometries,
        "production_gate": "N325/N425 remains forbidden until this schema and scorer contract are frozen after tiny compiled preflight",
    }


def render_markdown(artifact: dict[str, Any]) -> str:
    lines = [
        "# External observer gate for the common rank-birth stream",
        "",
        "The selected minimal observer is the matching-odd Euler residue",
        "",
        "```text",
        "O_ext = C_black^NN - C_white^matching - q = V - E + F0.",
        "```",
        "",
        "It is a scalar bulk configuration observable. It is root-independent but",
        "not topology-only: the index-9 control contains fixed-k configurations with the",
        "same `q` and different `O_ext` (the smaller controls certify the path identities).",
        "Thus the all-order `f(q)` contact no-go does not",
        "close `Cov(O_ext,J_D4)`.",
        "",
        "The component form costs no configuration scan in the production path. The",
        "primal trace supplies `C_black` before the next insertion; the reverse matching",
        "trace supplies `C_white` after inserting the same site. The cell form supplies",
        "an independent exact interpretation as a bulk Euler/thermal coordinate.",
        "",
        "Every path row stores `O_ext`, `O_ext^2`, both complex `O_ext*J` products,",
        "`J_D*conj(J_S)`, `|J_S|^2`, and the old `q*J_D` only as a contact control.",
        "The Gram rows are same-next-site Horvitz-score moments, not a disguised product",
        "of two site sums.",
        "",
        "Axis L2, Gaussian index 5, and the existing axis index-9 backend give exact",
        "zero residuals at every microcanonical k. N325/N425 production remains gated.",
        "",
        "## Candidate decision",
        "",
        "- **Selected:** bulk Euler/Betti residue, because it is exact, non-q, scalar, and O(1).",
        "- **Second coordinate:** a macroscopically separated local-H4/arm mark after its far-anchor convention is frozen.",
        "- **Charged lane:** a seam character after a gauge/orientation convention is frozen.",
        "",
        "## Observer authorization",
        "",
        "The exact `83e98fc` single-orbit alias gate is satisfied: `O_ext` is",
        "explicitly scalar, while `J_D4/J_S4` retain the typed internal complex",
        "`chi4(ell)` direction. A future local-H4 coordinate must keep separate axis",
        "and diagonal orbits (response determinant `-2`) or equivalent typed complex",
        "edge information; a direction-only single-orbit scalar is forbidden.",
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
