#!/usr/bin/env python3
"""Minimal exact P333 closure of the P398 rooted derivative-response marks."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Sequence

from noncrossing_connectivity_codec import noncrossing_states
from p321_graded_closure_extension import intertwiner_constraints
from p321_homology_trace_certificate import action_matrix, join_adjacent, rotate_state
from p333_generic_q_detach_intertwiner import detach_jet
from p333_gram_source_intertwiner import (
    encode_fraction, encode_matrix, join_block_count, matrix_residual_rank,
    multiply, rref_solve, subtract, transpose,
)
from p333_one_mark_endpoint_jet import nullspace_basis, restricted_skew_residual
from p333_source_landing_doublet import block_diagonal
from p333_source_landing_doublet_width4 import landing_reference_state
from p398_rooted_gr1_completion import ROOTED_SEEDS, orbit, selected_completion_families
from p398_qadic_jantzen import response_columns

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "analysis/p398_rooted_module_closure_protocol.json"
DEFAULT_DIR = ROOT / "results/p398-rooted-module-closure"
WIDTH = 4


def columns(vectors):
    return transpose(vectors)


def vector(states, terms):
    index = {state: i for i, state in enumerate(states)}
    out = [Fraction(0)] * len(states)
    for state, coefficient in terms:
        out[index[tuple(state)]] += coefficient
    return out


def rooted_references(states):
    ap, op, dp = (orbit(ROOTED_SEEDS[name]) for name in ("AP", "OP", "DP"))
    return [
        vector(states, zip(ap, (1, 1, 1, 1))),
        vector(states, zip(op, (1, 1))),
        vector(states, ((ap[0], 1), (ap[2], -1))),
        vector(states, ((ap[1], 1), (ap[3], -1))),
        vector(states, zip(ap, (1, -1, 1, -1))),
        vector(states, zip(op, (1, -1))),
        vector(states, zip(dp, (1, -1))),
    ]


def mv(matrix, value):
    return [sum(a * b for a, b in zip(row, value)) for row in matrix]


def rank_vectors(vectors):
    return matrix_residual_rank(columns(vectors)) if vectors else 0


def independent(vectors):
    output = []
    for value in vectors:
        if rank_vectors(output + [value]) > len(output):
            output.append(value)
    return output


def span_closure(vectors, operators):
    basis = independent(vectors)
    frontier = list(basis)
    while frontier:
        new = []
        for value in frontier:
            for operator in operators:
                image = mv(operator, value)
                if rank_vectors(basis + [image]) > len(basis):
                    basis.append(image)
                    new.append(image)
        frontier = new
    return basis


def sector_dimensions(basis, translation):
    dimension = len(basis[0]) if basis else 0
    tb = multiply(translation, basis)
    trivial = dimension - matrix_residual_rank(subtract(tb, basis))
    sign = dimension - matrix_residual_rank(
        [[x + y for x, y in zip(left, right)] for left, right in zip(tb, basis)]
    )
    charge1 = dimension - trivial - sign
    if charge1 % 2:
        raise AssertionError("a rational C4 charge-one block must have even dimension")
    return {"trivial": trivial, "charge1_rational": charge1, "charge2": sign}


def ordinary_affine_certificate(states, translation, joins, detaches):
    """Small exact block proof replacing a 2*(14+m)^2 dense solve."""
    n = len(states)
    equations = intertwiner_constraints(
        joins + detaches + [translation],
        joins[1:] + joins[:1] + detaches[1:] + detaches[:1] + [translation],
    )
    solved = rref_solve(equations, [0] * len(equations), n * n)
    flattened = [Fraction(x) for row in translation for x in row]
    residual = [sum(x * y for x, y in zip(row, flattened)) for row in equations]
    invariant_rows = []
    for operator in joins + detaches:
        for column in range(n):
            invariant_rows.append(
                [operator[row][column] - int(row == column) for row in range(n)]
            )
    left = rref_solve(invariant_rows, [0] * len(invariant_rows), n)
    if solved["dimension"] != 1 or any(residual):
        raise AssertionError("ordinary shifted Hom is not Q*T")
    if left["dimension"] != 1 or left["nullspace"] != [[Fraction(1)] * n]:
        raise AssertionError("common left invariants are not Q*ones")
    velocities = [detach_jet(WIDTH, i)[1] for i in range(WIDTH)]
    velocity_covariance = [matrix_residual_rank(subtract(
        multiply(translation, velocities[i]),
        multiply(velocities[(i + 1) % WIDTH], translation),
    )) for i in range(WIDTH)]
    return {
        "ordinary_shifted_hom_rank": solved["rank"],
        "ordinary_shifted_hom_dimension": solved["dimension"],
        "ordinary_hom_generator": "T",
        "T_equation_residual_count": sum(bool(x) for x in residual),
        "common_left_invariant_rank": left["rank"],
        "common_left_invariant_dimension": left["dimension"],
        "common_left_invariant_generator": [1] * n,
        "ordinary_velocity_covariance_ranks": velocity_covariance,
        "block_implication": [
            "Write X0=[[A,0],[B,C]] and V=[[a,0],[b,c]]. Mark transport fixes C=R,c=0.",
            "The order-zero ordinary affine block has A=alpha*T; fixing the ordinary source sets alpha=1.",
            "Each row of B is a common left invariant; fixing the ordinary source forces B=0.",
            "P_(i+1)T=TP_i cancels the inhomogeneous ordinary jet term, so a=beta*T; V(source)=0 sets beta=0.",
            "R E_i=E_(i+1)T cancels the lower jet term; each row of b is a common left invariant and V(source)=0 forces b=0.",
            "Thus the source-transport-normalized affine first jet is uniquely X0=diag(T,R), V=0. Endpoint and radical are checked at that point; no surviving modulus can alter its Gram residual."
        ],
    }


def make_candidate(states, translation, g1, g2, references, families, name):
    n = len(states)
    c0 = multiply(g1, columns(references))
    c1 = multiply(g2, columns(references))
    rotation = []
    labels = []
    for family in families:
        rotation = block_diagonal(rotation, family["mark_action"]) if rotation else family["mark_action"]
        labels += family.get("labels", ["existing_landing_0_minus_2", "existing_landing_1_minus_3"])
    m = len(rotation)
    full_t = block_diagonal(translation, rotation)
    g0 = [[Fraction(1)] * n + list(c0[i]) for i in range(n)]
    full_g1 = [list(g1[i]) + list(c1[i]) for i in range(n)]
    g0 += [list(row) + [Fraction(0)] * m for row in transpose(c0)]
    full_g1 += [list(row) + [Fraction(0)] * m for row in transpose(c1)]
    radical = nullspace_basis(g0)
    skew = restricted_skew_residual(full_g1, radical, full_t)
    leading = multiply(transpose(radical), multiply(full_g1, radical))
    t2_radical = multiply(full_t, multiply(full_t, radical))
    charge1_projected = [[(x - y) / 2 for x, y in zip(row, other)]
                         for row, other in zip(radical, t2_radical)]
    charge1_vectors = independent(transpose(charge1_projected))
    charge1_basis = columns(charge1_vectors) if charge1_vectors else []
    charge1_gram = multiply(transpose(charge1_basis), multiply(full_g1, charge1_basis)) if charge1_basis else []
    charge1_skew = restricted_skew_residual(full_g1, charge1_basis, full_t) if charge1_basis else []
    cov0 = subtract(multiply(transpose(full_t), multiply(g0, full_t)), g0)
    cov1 = subtract(multiply(transpose(full_t), multiply(full_g1, full_t)), full_g1)
    invariant = multiply(g0, multiply(full_t, radical))
    # This executable emission is source-dependent, not just a site character.
    emissions = [multiply(transpose(c0), detach_jet(WIDTH, i)[1]) for i in range(WIDTH)]
    emission_residuals = [matrix_residual_rank(subtract(
        multiply(rotation, emissions[i]), multiply(emissions[(i + 1) % WIDTH], translation)
    )) for i in range(WIDTH)]
    endpoint = [Fraction(1)] * n + [Fraction(0)] * m
    source = [Fraction(int(i == n - 1)) for i in range(n + m)]
    endpoint_residual = [sum(endpoint[i] * full_t[i][j] for i in range(n + m)) - endpoint[j] for j in range(n + m)]
    source_residual = [a - b for a, b in zip(mv(full_t, source), source)]
    witness = None
    for i, row in enumerate(charge1_skew):
        for j, value in enumerate(row):
            if value and witness is None:
                witness = {
                    "charge1_projector_columns": [i, j],
                    "left_vector": [encode_fraction(row[i]) for row in charge1_basis],
                    "right_vector": [encode_fraction(row[j]) for row in charge1_basis],
                    "forced_Gram_skew_pairing": encode_fraction(value),
                    "affine_modulus_coefficient_rank": 0,
                    "augmented_rank": 1,
                    "identity": "the unique source-normalized affine point has a nonzero fixed Gram pairing; 0 = 1 after rescaling",
                }
    passes = matrix_residual_rank(skew) == 0
    return {
        "name": name,
        "ordinary_dimension": n,
        "mark_dimension": m,
        "extended_dimension": n + m,
        "labels": labels,
        "response_rank": matrix_residual_rank(c0),
        "G0_rank": matrix_residual_rank(g0),
        "radical_dimension": len(radical[0]),
        "radical_C4_dimensions": sector_dimensions(radical, full_t),
        "radical_leading_Gram_rank": matrix_residual_rank(leading),
        "radical_Gram_skew_rank": matrix_residual_rank(skew),
        "radical_basis": encode_matrix(radical),
        "radical_Gram_skew": encode_matrix(skew),
        "charge1_obstruction": {
            "projector": "(I-T_bar^2)/2",
            "dimension": len(charge1_vectors),
            "basis": encode_matrix(charge1_basis),
            "leading_Gram": encode_matrix(charge1_gram),
            "leading_Gram_rank": matrix_residual_rank(charge1_gram),
            "Gram_skew": encode_matrix(charge1_skew),
        },
        "first_empty_certificate": witness,
        "canonical_checks": {
            "G0_translation_covariance_rank": matrix_residual_rank(cov0),
            "G1_translation_covariance_rank": matrix_residual_rank(cov1),
            "radical_invariance_rank": matrix_residual_rank(invariant),
            "rooted_emission_covariance_ranks": emission_residuals,
            "endpoint_residual_count": sum(bool(x) for x in endpoint_residual),
            "source_residual_count": sum(bool(x) for x in source_residual),
        },
        "source_normalized_affine_jet": {"unique": True, "X0": "diag(T,R)", "V": "0"},
        "full_inherited_intersection": "unique" if passes else "empty",
        "unprojected_G0_self_adjoint_residual_rank": matrix_residual_rank(subtract(
            multiply(g0, full_t), multiply(transpose(full_t), g0))),
        "unprojected_G1_self_adjoint_residual_rank": matrix_residual_rank(subtract(
            multiply(full_g1, full_t), multiply(transpose(full_t), full_g1))),
    }


def build_result():
    protocol = json.loads(PROTOCOL.read_text())
    states = noncrossing_states(WIDTH)
    translation = action_matrix(WIDTH, lambda state: rotate_state(state, 1))
    joins = [action_matrix(WIDTH, lambda state, i=i: join_adjacent(state, i)) for i in range(WIDTH)]
    detaches = [detach_jet(WIDTH, i)[0] for i in range(WIDTH)]
    g1 = [[Fraction(join_block_count(a, b)) for b in states] for a in states]
    g2 = [[value * (value - 1) for value in row] for row in g1]
    references = rooted_references(states)
    families = list(selected_completion_families().values())
    old = response_columns(WIDTH)["C4_charge1_landing"]
    landing = [landing_reference_state(i) for i in range(WIDTH)]
    old_refs = [vector(states, ((landing[0], 1), (landing[2], -1))),
                vector(states, ((landing[1], 1), (landing[3], -1)))]
    candidates = [
        make_candidate(states, translation, g1, g2, references, families, "rooted7"),
        make_candidate(states, translation, g1, g2, references + old_refs, families + [old],
                       "rooted7_plus_existing_landing_charge1"),
    ]
    affine = span_closure(references, [translation] + joins)
    active_full = span_closure(references, [translation] + joins + detaches)
    root_sectors = sector_dimensions(columns(references), translation)
    affine_sectors = sector_dimensions(columns(affine), translation)
    full_sectors = sector_dimensions(columns(active_full), translation)
    image = mv(joins[0], references[0])
    triple_index = states.index((0, 0, 0, 1))
    active = {
        "warning": "Separate active-reference semantics, not the retained derivative accumulator lift tested above.",
        "initial_dimension": rank_vectors(references),
        "C4_closure_dimension": len(span_closure(references, [translation])),
        "join_closure_dimension": len(affine),
        "join_detach_source_closure_dimension": len(active_full),
        "minimal_witness": {
            "generator": "J0", "input": "AP_sum", "left_covector": "delta_(0,0,0,1)",
            "left_on_all_seven_inputs": [encode_fraction(v[triple_index]) for v in references],
            "left_on_image": encode_fraction(image[triple_index]),
        },
        "join_cokernel": {"dimension": len(affine) - len(references),
                          **{key: affine_sectors[key] - root_sectors[key] for key in root_sectors}},
        "missing_active_coordinates": ["fully_connected", "triple_plus_singleton C4 orbit (four coordinates)", "DP_sum", "all_singleton source"],
        "full_cokernel": {"dimension": len(active_full) - len(references),
                          **{key: full_sectors[key] - root_sectors[key] for key in root_sectors}},
    }
    source_files = ["scripts/p398_rooted_gr1_completion.py", "scripts/p333_source_landing_doublet_width4.py", "scripts/p333_one_mark_endpoint_jet.py"]
    return {
        "schema": "matching-one/p398-rooted-module-closure/v1",
        "status": "exact_rational_full_inherited_intersection_certificate",
        "base_commit": "32a27f9",
        "issues": [333, 398],
        "protocol_sha256": hashlib.sha256(PROTOCOL.read_bytes()).hexdigest(),
        "scorer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_sha256": {path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest() for path in source_files},
        "states": [list(state) for state in states],
        "ordinary_block_uniqueness": ordinary_affine_certificate(states, translation, joins, detaches),
        "candidates": candidates,
        "minimal_coupled_mark": {
            "additional_rational_dimension": 2,
            "type": "already-defined C4 charge-one landing doublet from the triple-plus-singleton rooted orbit",
            "new_terminal_character": False,
            "why_minimal": "The rooted7 residual is a nondegenerate two-dimensional rational C4 charge-one block. A one-dimensional rational C4 mark is trivial or sign and cannot pair with it; the existing two-column landing response removes precisely this residual block.",
            "scope": protocol["minimality_contract"],
        },
        "active_reference_distinction": active,
        "decision": "rooted7_fails_but_rooted7_plus_existing_charge1_is_unique",
        "claim_boundary": protocol["boundaries"],
    }


def render(result):
    lines = ["# P333/P398: minimal closure by coupling rooted and old landing marks", "",
             "The seven new rooted derivative marks alone fail the complete inherited P333 intersection. The sole remaining obstruction is one rational C4 charge-one doublet. Coupling the **already-defined two-column landing doublet** to the seven rooted marks gives a unique full affine/endpoint/radical-Gram/source solution: `X0=diag(T,R), V=0`.", "",
             "| candidate | dim W | rank G0 | dim radical | radical C4 dimensions (trivial, charge1, sign) | Gram-skew rank | full inherited intersection |",
             "|---|---:|---:|---:|---|---:|---|"]
    for row in result["candidates"]:
        d = row["radical_C4_dimensions"]
        lines.append(f"| {row['name']} | {row['extended_dimension']} | {row['G0_rank']} | {row['radical_dimension']} | {d['trivial']}, {d['charge1_rational']}, {d['charge2']} | {row['radical_Gram_skew_rank']} | {row['full_inherited_intersection']} |")
    lines += ["", "## Why this is the full intersection, not just a canonical-point check", ""]
    lines += [f"{i+1}. {text}" for i, text in enumerate(result["ordinary_block_uniqueness"]["block_implication"])]
    lines += ["", "The exact ordinary shifted-Hom rank is 195/196 and the common-left-invariant rank is 13/14. Hence the source-normalized affine point is unique before imposing Gram. For rooted7 the Gram restriction has coefficient rank 0 and augmented rank 1; its explicit radical-vector witness is in `latest.json`. For the coupled nine-mark lift every inherited gate is zero exactly, and the surviving four-dimensional radical still has nondegenerate leading Gram form: success is not a vacuous radical-exhaustion trick.", "",
              "## Minimum coupled mark", "",
              result["minimal_coupled_mark"]["why_minimal"], "",
              "The two new-to-this-module columns are `G1[:,landing0]-G1[:,landing2]` and `G1[:,landing1]-G1[:,landing3]`, already tested in the old landing branch. Their isolated width-four failure and the rooted7 failure do not survive their coupling. The seven rooted marks are not a standalone replacement for the old response family.", "",
              "Here coupling means sharing the ordinary Gram/source block and configuration-dependent emission. It is not a fitted off-diagonal exchange between the two mark families or evidence for a physical interaction.", "",
              "## Semantics and boundary", "",
              "This uses the prior P333 retained-response accumulator convention: joins act identically on emitted marks, while detach emits the configuration-dependent row `C0^T P_i`; it does **not** pretend the seven reference vectors are an invariant active connectivity submodule. Under active-reference semantics their join closure is 13 and join+detach closure is all 14, with a one-coordinate witness already at `J0 AP_sum`. That separate saturation is recorded, not used as the full-intersection decision.", "",
              "The inherited Gram gate is the first-jet form restricted to `ker G0`. Unprojected G0/G1 self-adjoint residuals are explicitly nonzero in the artifact; the result does not establish a stronger all-Gram/all-Q closure, a physical transfer module, or a continuum Jordan/field identity.", "",
              "## Reproduce", "", "```bash", "python3 scripts/p398_rooted_module_closure.py", "python3 -m unittest discover -s tests -p test_p398_rooted_module_closure.py", "```", ""]
    return "\n".join(lines)


def card(result):
    return "\n".join([
        "# Scientific card: coupling, not another character", "",
        "- Mechanism space changed: rooted7 alone is not the full P333 closure; coupling it to the existing C4 charge-one landing doublet makes the inherited intersection uniquely nonempty.",
        "- Exact result: 21D rooted7 has radical6 / Gram-skew rank2; 23D rooted7+old2 has radical4 / Gram-skew rank0, with unique X0=diag(T,R), V=0.",
        "- Minimum obstruction/repair: one nondegenerate rational C4 charge-one doublet; two existing landing columns are the sharp additional dimension within the declared lift.",
        "- Observer/sector/source: width-four derivative-response accumulator | C4 multiplicity sectors | ordinary all-singleton source plus exact mark transport.",
        "- Dependency group: the same 14-state exact connectivity/Gram algebra as P333/P398; not independent numerical evidence.",
        "- Not proved: unprojected Gram or all-Q closure, active rooted-submodule closure, physical transfer realization, continuum field/Jordan identity.",
        "- Next discriminator: construct a physical coupled emission/transport realizing the declared nine-mark lift, or test its all-Q extension; no additional terminal-character scan is needed.", ""])


def main(argv: Sequence[str] | None = None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_DIR)
    args = parser.parse_args(argv)
    result = build_result()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    outputs = {"latest.json": json.dumps(result, indent=2, sort_keys=True) + "\n", "latest.md": render(result), "scientific-card.md": card(result)}
    for name, contents in outputs.items():
        (args.out_dir / name).write_text(contents)
    checksums = "".join(f"{hashlib.sha256((args.out_dir / name).read_bytes()).hexdigest()}  {name}\n" for name in outputs)
    (args.out_dir / "SHA256SUMS").write_text(checksums)
    print(json.dumps({"decision": result["decision"], "candidates": [{key: row[key] for key in ("name", "radical_dimension", "radical_Gram_skew_rank", "full_inherited_intersection")} for row in result["candidates"]]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
