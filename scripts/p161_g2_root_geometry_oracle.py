#!/usr/bin/env python3
"""Exact geometry/no-go oracle for a pure Weierstrass-g2 root bias."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp

from derive_hexagonal_degree2_hecke import normalized_eisenstein


def modular_basis(weight: int) -> list[tuple[int, int]]:
    return [
        (a, b)
        for a in range(weight // 4 + 1)
        for b in range(weight // 6 + 1)
        if 4 * a + 6 * b == weight
    ]


def c4_allowed_ring(max_weight: int = 36) -> list[dict]:
    rows = []
    for weight in range(4, max_weight + 1, 2):
        for a, b in modular_basis(weight):
            if weight % 4:
                continue
            rows.append({
                "weight": weight,
                "monomial": f"E4^{a}*E6^{b}",
                "E4_power": a,
                "E6_power": b,
                "root_N_power": weight // 2,
                "hex_zero_order": a,
                "survives_exact_hex": a == 0,
                "hex_child_character": f"zeta^{a % 3}",
            })
    return rows


def numerical_oracle(dps: int = 90) -> dict:
    with mp.workdps(dps):
        i = mp.j
        rectangle = 2 * mp.j
        rho = (1 + mp.sqrt(3) * mp.j) / 2
        children = [2 * rho, rho / 2, (rho + 1) / 2]
        square = normalized_eisenstein(4, i, dps=dps)[0]
        rectangular = normalized_eisenstein(4, rectangle, dps=dps)[0]
        hexagonal = normalized_eisenstein(4, rho, dps=dps)[0]
        child_values = [normalized_eisenstein(4, tau, dps=dps)[0] for tau in children]
        zeta = mp.exp(2 * mp.pi * mp.j / 3)
        return {
            "dps": dps,
            "E4hat": {
                "i": {"real": mp.nstr(mp.re(square), 70), "imag": mp.nstr(mp.im(square), 70)},
                "2i": {"real": mp.nstr(mp.re(rectangular), 70), "imag": mp.nstr(mp.im(rectangular), 70)},
                "rho": {"real": mp.nstr(mp.re(hexagonal), 70), "imag": mp.nstr(mp.im(hexagonal), 70)},
                "hex_children": [
                    {"real": mp.nstr(mp.re(value), 70), "imag": mp.nstr(mp.im(value), 70)}
                    for value in child_values
                ],
            },
            "errors": {
                "rectangle_11_over_4": mp.nstr(abs(rectangular / square - mp.mpf(11) / 4), 12),
                "hex_zero": mp.nstr(abs(hexagonal), 12),
                "child_zeta": mp.nstr(abs(child_values[1] / child_values[0] - zeta), 12),
                "child_zeta2": mp.nstr(abs(child_values[2] / child_values[0] - zeta**2), 12),
                "child_sum": mp.nstr(abs(sum(child_values)), 12),
            },
        }


def analyze(dps: int = 90) -> dict:
    allowed = c4_allowed_ring(36)
    survivors = [row for row in allowed if row["survives_exact_hex"]]
    return {
        "schema": "matching-one/g2-root-geometry-no-go/v1",
        "issue": 161,
        "status": "exact_three_modulus_fingerprint_and_C4_hex_selection",
        "pure_g2_law": {
            "complex": "delta p_chiral(Lambda)=lambda4*g2(Lambda)",
            "real": "delta p_root=lambda4*Re[g2(Lambda)] for a reflection-even aligned square lattice",
            "area_normalized": "A^2 delta p_chiral=C*E4hat(tau), E4hat=Im(tau)^2 E4(tau)",
            "required_observable": "root or residual-to-slope response in which the common thermal-primary block cancels",
        },
        "three_modulus_fingerprint": {
            "moduli": ["i", "2i", "rho=exp(i*pi/3)"],
            "amplitude_vector": "(W(i),W(2i),W(rho))=A*(4,11,0)",
            "nulls": ["4 W(2i)-11 W(i)=0", "W(rho)=0"],
            "meaning": "one common amplitude, one enhanced rectangle, one exact elliptic zero",
        },
        "hexagonal_degree2_closure": {
            "child_order": ["2rho", "rho/2", "(rho+1)/2"],
            "complex_vector": "B*(1,zeta,zeta^2)",
            "complex_nulls": ["Y1-zeta*Y0=0", "Y2-zeta^2*Y0=0", "Y0+Y1+Y2=0"],
            "reflection_even_vector": "B*(1,-1/2,-1/2)",
            "real_nulls": ["2Y1+Y0=0", "2Y2+Y0=0"],
            "grading_use": "E6 contamination lies in DFT r=0 and E4^2 contamination in r=2, while pure g2 lies in r=1",
        },
        "C4_selection_at_exact_hex": {
            "microscopic_rule": "a C4-invariant square-lattice coupling admits only spins/weights divisible by 4",
            "ordinary_ring_rule": "among E4^a E6^b with weight divisible by 4, exact rho kills every a>0 term",
            "first_survivor": survivors[0],
            "conclusion": "after the weight-4 g2 zero, the first C4-allowed ordinary holomorphic term nonzero at exact rho is E6^2 of weight 12, giving root bias L^-12=N^-6",
            "pell_approximant": "if tau_L-rho=O(L^-2), the simple E4 zero leaks the leading g2 term as L^-6=N^-3 and dominates the intrinsic E6^2 N^-6 term",
            "allowed_spectrum": allowed,
        },
        "no_go_conditions": [
            "a nonzero leading N^-2 coefficient at exact rho",
            "failure of W(2i)/W(i)=11/4 after the same root/slope normalization",
            "hex-child DFT support outside r=1 at leading N^-2 order",
            "failure of the complex cube-root phase closure even when the real cos(4theta) projection appears correct",
            "a half-integer or non-ring correction exponent claimed as the same ordinary local g2 field",
        ],
        "identification_boundary": {
            "passing_does_not_select_module": "M4 is one-dimensional, so every ordinary scalar/chiral weight-4 modular response is proportional to E4; the test supports a local weight-4 law but cannot by itself distinguish thermal Q4 from another weight-4 module",
            "sector_valued_exception": "a vector-valued KdV/defect response may be nonzero at rho by occupying a nontrivial sector character; that is not a scalar g2 root law",
            "additional_terms": "Jordan/quasimodular derivatives, [2] four-leg fields, topological sectors, and quotient arithmetic can violate the nulls without contradicting the mathematical E4 identities",
        },
        "frozen_score": {
            "primary_vector_residual": ["4W_2i-11W_i", "W_rho"],
            "secondary_child_residual": ["Y1-zetaY0", "Y2-zeta^2Y0"],
            "covariance": "recompute centers, H4 projection, slope division, area normalization, and all residuals inside each synchronized replicate",
            "amplitude_parameters": 1,
        },
        "numerical_oracle": numerical_oracle(dps),
        "claim_boundary": {
            "exact": "E4/Hecke/elliptic identities, C3 child grading, and C4 modular-ring selection",
            "conditional_bridge": "one typed matching-root observable must share a common local weight-4 coupling and frame across all geometries",
            "not_claimed": "that the current lattice root bias already passes these nulls or that passing identifies a unique CFT module",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=90)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = json.dumps(analyze(args.dps), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
