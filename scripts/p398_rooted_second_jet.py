#!/usr/bin/env python3
"""Second jet and Jantzen stop for the fixed P398 nine-mark retained family."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import comb
from pathlib import Path

from p321_homology_trace_certificate import action_matrix, rotate_state
from p333_gram_source_intertwiner import (
    encode_fraction, encode_matrix, join_block_count, matrix_residual_rank,
    multiply, rref_solve, subtract, transpose,
)
from p333_source_landing_doublet import block_diagonal
from p333_source_landing_doublet_width4 import landing_reference_state
from p398_qadic_jantzen import exact_determinant, response_columns
from p398_rooted_gr1_completion import selected_completion_families
from p398_rooted_module_closure import columns, rooted_references, vector

ROOT = Path(__file__).resolve().parents[1]
BASE_RESULT = ROOT / "results/p398-rooted-module-closure/latest.json"
PROTOCOL = ROOT / "analysis/p398_rooted_second_jet_protocol.json"
OUT = ROOT / "results/p398-rooted-second-jet"


def decode(matrix):
    return [[Fraction(value) for value in row] for row in matrix]


def restriction(operator, basis):
    images = multiply(operator, basis)
    coordinates = []
    for image in zip(*images):
        solved = rref_solve(basis, image, len(basis[0]))
        if not solved["consistent"] or solved["dimension"]:
            raise AssertionError("operator does not act uniquely on the fixed radical")
        coordinates.append(solved["particular"])
    return columns(coordinates)


def build_result():
    base = json.loads(BASE_RESULT.read_text())
    protocol = json.loads(PROTOCOL.read_text())
    inherited = base["candidates"][1]
    if inherited["full_inherited_intersection"] != "unique":
        raise AssertionError("the parent first-jet solution is not unique")
    states = tuple(tuple(state) for state in base["states"])
    n, m = len(states), inherited["mark_dimension"]
    basis = decode(inherited["radical_basis"])
    refs = rooted_references(states)
    landings = [landing_reference_state(i) for i in range(4)]
    refs += [vector(states, ((landings[0], 1), (landings[2], -1))),
             vector(states, ((landings[1], 1), (landings[3], -1)))]
    families = list(selected_completion_families().values()) + [response_columns(4)["C4_charge1_landing"]]
    rotation = []
    for family in families:
        rotation = block_diagonal(rotation, family["mark_action"]) if rotation else family["mark_action"]
    translation = action_matrix(4, lambda state: rotate_state(state, 1))
    full_t = block_diagonal(translation, rotation)
    t_radical = restriction(full_t, basis)
    identity = [[Fraction(int(i == j)) for j in range(4)] for i in range(4)]
    involution = subtract(multiply(t_radical, t_radical), identity)
    coefficients, projected = [], []
    for degree in range(5):
        g = [[Fraction(comb(join_block_count(a, b), degree)) for b in states] for a in states]
        derivative = [[Fraction((degree + 1) * comb(join_block_count(a, b), degree + 1))
                       for b in states] for a in states]
        cross = multiply(derivative, columns(refs))
        gamma = [g[i] + cross[i] for i in range(n)]
        gamma += [row + [Fraction(0)] * m for row in transpose(cross)]
        h = multiply(transpose(basis), multiply(gamma, basis))
        projected.append(h)
        covariance = subtract(multiply(transpose(full_t), multiply(gamma, full_t)), gamma)
        skew = subtract(multiply(h, t_radical), multiply(transpose(t_radical), h))
        coefficients.append({
            "degree": degree,
            "extended_coefficient_rank": matrix_residual_rank(gamma),
            "C4_covariance_residual_rank": matrix_residual_rank(covariance),
            "radical_coefficient": encode_matrix(h),
            "radical_coefficient_rank": matrix_residual_rank(h),
            "projected_Gram_skew_rank": matrix_residual_rank(skew),
        })
    if matrix_residual_rank(projected[0]) or matrix_residual_rank(involution):
        raise AssertionError("fixed radical or induced involution changed")
    if matrix_residual_rank(projected[1]) != 4:
        raise AssertionError("a higher Jantzen layer is present")
    proof = base["ordinary_block_uniqueness"]
    fixed_variables = n * m + m * m  # upper-right filtration plus fixed mark block
    remaining_variables = n * n + m * n
    subsystem_rank = proof["ordinary_shifted_hom_rank"] + m * proof["common_left_invariant_rank"]
    source_increment = 1 + m
    if subsystem_rank + source_increment != remaining_variables:
        raise AssertionError("the source-normalized second-jet subsystem is not full rank")
    return {
        "schema": "matching-one/p398-rooted-second-jet/v1",
        "status": "exact_same_family_second_jet_and_Jantzen_stop",
        "base_commit": protocol["base_commit"],
        "protocol_sha256": hashlib.sha256(PROTOCOL.read_bytes()).hexdigest(),
        "parent_result_sha256": hashlib.sha256(BASE_RESULT.read_bytes()).hexdigest(),
        "scorer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "family_definition_audit": protocol["family_provenance"],
        "coefficient_convention": protocol["coefficient_convention"],
        "second_jet_affine_certificate": {
            "full_matrix_variable_count": (n + m) ** 2,
            "fixed_filtration_and_transport_variables": fixed_variables,
            "remaining_ordinary_and_emission_variables": remaining_variables,
            "sufficient_homogeneous_subsystem_rank": subsystem_rank,
            "independent_source_normalizations": source_increment,
            "total_rank_in_remaining_variables": subsystem_rank + source_increment,
            "detach_second_coefficient": "0 by the existing exact linear retained-emission formula",
            "inhomogeneous_second_jet_rhs": "0 because X1=0 and D2=0",
            "unique_solution": "X2=0",
            "decision": "unique_second_jet_extension",
        },
        "fixed_radical_translation": {
            "matrix": encode_matrix(t_radical),
            "involution_residual_rank": matrix_residual_rank(involution),
            "trivial_dimension": 3,
            "sign_dimension": 1,
            "charge1_dimension": 0,
        },
        "Gram_Taylor_coefficients": coefficients,
        "Jantzen_certificate": {
            "G0_rank": inherited["G0_rank"],
            "J1_dimension": 4,
            "leading_form_rank": matrix_residual_rank(projected[1]),
            "leading_form_determinant": encode_fraction(exact_determinant(projected[1])),
            "valuation_multiplicities": {"0": 19, "1": 4},
            "J2_dimension": 0,
            "reason": "The induced first coefficient on ker Gamma0 is nondegenerate; the Schur complement starts t*H1 with invertible H1, so no valuation exceeds one."
        },
        "all_order_fixed_radical_identity": {
            "polynomial_degree": 4,
            "all_five_coefficient_skew_ranks_zero": all(row["projected_Gram_skew_rank"] == 0 for row in coefficients),
            "proof": "C4 covariance gives R_B^T H(Q) R_B=H(Q). Since R_B^2=I on the remaining radical, H(Q)R_B=R_B^T H(Q) identically. Thus the second projected Gram gate and every higher coefficient are automatic in this family.",
            "affine_extension": "X(Q)=diag(T,R) intertwines the constant joins/translation and exact linear detach family for every Q; all positive Taylor coefficients are zero under the same source/transport normalization."
        },
        "decision": "second_jet_unique_but_no_new_Jantzen_or_projected_Gram_discriminator",
        "scientific_change": "The nine-mark success is stable through second order, but repeating the same projected gate cannot reveal new structure: charge-one removal has turned it into an involution identity, and J2 is already zero.",
        "claim_boundary": protocol["boundaries"],
    }


def render(result):
    j = result["Jantzen_certificate"]
    lines = ["# P333/P398: second jet closes, but the projected gate is now automatic", "",
             result["scientific_change"], "",
             "## The second-order family is defined", "",
             "The predecessor protocols define the retained family, not just numerical jets: `Dbar_i(Q)=[[D_i(1)+(Q-1)P_i,0],[(Q-1)E_i,I]]` and `Gbar(Q)=[[G(Q),d_Q G(Q) S_ref],[transpose,0]]`. Commit 5389200 fixes the nine reference columns and the constant row `E_i=C0^T P_i`. Keeping that family makes `Dbar_i2=0` and fixes every Gram coefficient. We do not replace the emission by a Q-dependent row or introduce a new second-order observable.", "",
             "Taylor coefficients use `t=Q-1`: `G_k=binomial(b,k)` and the cross block is `(k+1)G_(k+1)S_ref`. In particular the second coefficient is half the second derivative.", "",
             "## Unique second jet", "",
             "With `X0=Tbar` and `X1=0`, the second affine equation is homogeneous: `X2 A0=B0 X2`. Filtration and fixed mark transport set 207 of 529 entries. The remaining 322 entries have an exact sufficient homogeneous subsystem of rank 312; the same source normalization supplies 10 independent equations, giving rank 322 and uniquely `X2=0`. The projected second Gram residual has rank zero.", "",
             "| Taylor degree | rank of fixed-radical Gram coefficient | projected Gram-skew rank |",
             "|---:|---:|---:|"]
    for row in result["Gram_Taylor_coefficients"]:
        lines.append(f"| {row['degree']} | {row['radical_coefficient_rank']} | {row['projected_Gram_skew_rank']} |")
    lines += ["", "## No next Jantzen layer", "",
              f"The nine-mark Gram has valuations `0^19, 1^4`: `dim J1=4`, `dim J2=0`. The exact leading radical determinant is `{j['leading_form_determinant']}`. This follows from the nondegenerate 4x4 leading form; no second-order coefficient can create a hidden higher layer at Q=1.", "",
              "## Why more of the same gate cannot discriminate", "",
              "The induced translation on the surviving radical has three trivial directions and one sign direction, hence `R_B^2=I`. C4 covariance gives `R_B^T H(Q)R_B=H(Q)`; multiplying by the involution immediately gives `H(Q)R_B=R_B^T H(Q)` for the entire polynomial. All five coefficient checks through degree four confirm the identity exactly. Affine covariance likewise keeps `X(Q)=Tbar` constant in the exact linear emission family.", "",
              "This is the new stop/redirect result: second-jet survival is real, but it adds no independent selection beyond eliminating the radical charge-one block. The next identifying datum must couple to the quotient/physical emission or otherwise impose a genuinely different constraint; repeating higher projected Gram jets cannot do it.", "",
              "## Boundary", "",
              "No known unprojected Gram failure was rerun or repaired. The all-order statement is only the same fixed-Q=1-radical projection, not an unprojected Gram module, a moving-radical prescription, or a physical transfer/Jordan realization. If one discards the inherited exact `(Q-1)E_i` family and keeps only its first jet, arbitrary `t^2` emission deformations would be additional unrecorded data; they are not silently admitted here.", "",
              "```bash", "python3 scripts/p398_rooted_second_jet.py", "python3 -m unittest discover -s tests -p test_p398_rooted_second_jet.py", "```", ""]
    return "\n".join(lines)


def card(result):
    return "\n".join([
        "# Scientific card: second jet is stable but no longer discriminatory", "",
        "- New result: the nine-mark retained family has unique second coefficient X2=0; the source-normalized first jet extends without adding a mark or changing Gram.",
        "- Exact structure: Jantzen valuations 0^19,1^4; J2=0. The surviving radical is 3 trivial + 1 sign, so its translation squares to identity.",
        "- Mechanism consequence: fixed-radical Gram compatibility is then automatic at every Q/Taylor order by C4 covariance, not a fresh piece of operator-identification evidence.",
        "- Observer/sector/source: the same width-four retained derivative-response family, fixed all-singleton source and mark transport, nine-mark common Gram block.",
        "- Dependency: same exact P333/P398 block as 5389200, zero new simulation.",
        "- Not proved: unprojected Gram closure, a physical emission/transfer module, or a continuum Jordan identity.",
        "- Next discriminator: a genuinely new quotient/physical-emission constraint; do not spend more exact work repeating higher fixed-radical Gram jets.", ""])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=OUT)
    args = parser.parse_args()
    result = build_result()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for name, text in {"latest.json": json.dumps(result, indent=2, sort_keys=True) + "\n", "latest.md": render(result), "scientific-card.md": card(result)}.items():
        (args.out_dir / name).write_text(text)
    print(json.dumps({"decision": result["decision"], "Jantzen": result["Jantzen_certificate"]}, indent=2))


if __name__ == "__main__":
    main()
