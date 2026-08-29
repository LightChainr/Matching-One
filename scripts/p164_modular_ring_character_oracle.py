#!/usr/bin/env python3
"""Exact modular-ring/Hecke-character oracle for low irrelevant corrections."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp

from derive_hexagonal_degree2_hecke import normalized_eisenstein


def modular_basis(weight: int) -> list[tuple[int, int]]:
    """Exponent pairs (a,b) for E4^a E6^b of total weight."""
    return [
        (a, b)
        for a in range(weight // 4 + 1)
        for b in range(weight // 6 + 1)
        if 4 * a + 6 * b == weight
    ]


def ring_spectrum(max_weight: int = 24) -> list[dict]:
    rows = []
    for weight in range(2, max_weight + 1, 2):
        basis = modular_basis(weight)
        rows.append({
            "weight": weight,
            "basis": [f"E4^{a}*E6^{b}" for a, b in basis],
            "dimension": len(basis),
            "N_correction_power_if_x_equals_weight": f"N^-{(weight - 2) // 2}" if basis else None,
            "hex_child_characters": [f"zeta^{a % 3}" for a, _ in basis],
            "square_zero": [b > 0 for _, b in basis],
            "hex_parent_zero_order": [a for a, _ in basis],
        })
    return rows


def dft(vector: list[mp.mpc], *, dps: int) -> list[mp.mpc]:
    with mp.workdps(dps):
        zeta = mp.exp(2 * mp.pi * mp.j / 3)
        return [
            mp.fsum(zeta ** (-r * j) * vector[j] for j in range(3)) / 3
            for r in range(3)
        ]


def numerical_oracle(dps: int = 90) -> dict:
    with mp.workdps(dps):
        omega = (1 + mp.sqrt(3) * mp.j) / 2
        children = [2 * omega, omega / 2, (omega + 1) / 2]
        e4 = [normalized_eisenstein(4, tau, dps=dps)[0] for tau in children]
        e6 = [normalized_eisenstein(6, tau, dps=dps)[0] for tau in children]
        e4sq = [value * value for value in e4]
        vectors = {"E4": e4, "E6": e6, "E4^2": e4sq}
        payload = {}
        for name, vector in vectors.items():
            coeff = dft(vector, dps=dps)
            total = mp.sqrt(mp.fsum(abs(value) ** 2 for value in coeff))
            payload[name] = {
                "children": [
                    {"real": mp.nstr(mp.re(value), 70), "imag": mp.nstr(mp.im(value), 70)}
                    for value in vector
                ],
                "DFT": [
                    {"real": mp.nstr(mp.re(value), 70), "imag": mp.nstr(mp.im(value), 70)}
                    for value in coeff
                ],
                "normalized_DFT_power": [mp.nstr(abs(value) ** 2 / total**2, 20) for value in coeff],
            }
        zeta = mp.exp(2 * mp.pi * mp.j / 3)
        return {
            "dps": dps,
            "vectors": payload,
            "errors": {
                "E4_character": mp.nstr(max(abs(e4[1] / e4[0] - zeta), abs(e4[2] / e4[0] - zeta**2)), 12),
                "E6_constant": mp.nstr(max(abs(e6[1] / e6[0] - 1), abs(e6[2] / e6[0] - 1)), 12),
                "E4sq_conjugate_character": mp.nstr(max(abs(e4sq[1] / e4sq[0] - zeta**2), abs(e4sq[2] / e4sq[0] - zeta)), 12),
            },
        }


def analyze(dps: int = 90) -> dict:
    return {
        "schema": "matching-one/modular-ring-character-syzygy/v1",
        "issue": 164,
        "status": "exact_ring_grading_syzygy_and_frozen_three_child_nulls",
        "serre_derivative_syzygies": {
            "D4_E4": "D_4 E4=-E6/3",
            "D6_E6": "D_6 E6=-E4^2/2",
            "composed": "6 D_6 D_4 E4=E4^2",
            "product": "D_8(E4^2)=-2 E4 E6/3",
            "interpretation": "modular-covariant derivative corrections through weight 8 add no new shape directions beyond E6 and E4^2",
        },
        "hexagonal_degree2_ring_character": {
            "child_order": ["2omega", "omega/2", "(omega+1)/2"],
            "generator_images": {
                "E4": "A*(1,zeta,zeta^2)",
                "E6": "B*(1,1,1)",
            },
            "monomial_rule": "E4^a E6^b maps to one C3 child character zeta^(a mod 3)",
            "ring_homomorphism": "chi_child: C[E4,E6] -> Z/3Z with deg(E4)=1 and deg(E6)=0",
            "parent_gate": "at omega, every monomial with a>0 vanishes to order a; pure E6^b survives",
        },
        "three_child_discriminator": {
            "E4": {
                "ratios": ["1", "zeta", "zeta^2"],
                "DFT_support": 1,
                "nulls": ["y1-zeta*y0=0", "y2-zeta^2*y0=0"],
            },
            "E6_or_D4E4": {
                "ratios": ["1", "1", "1"],
                "DFT_support": 0,
                "nulls": ["y1-y0=0", "y2-y0=0"],
            },
            "E4^2_or_D6E6": {
                "ratios": ["1", "zeta^2", "zeta"],
                "DFT_support": 2,
                "nulls": ["y1-zeta^2*y0=0", "y2-zeta*y0=0"],
            },
            "real_projection_warning": "E4 and E4^2 both project to 1:-1/2:-1/2; their imaginary/chiral signs are opposite and must be retained.",
        },
        "missing_spectrum_rule": {
            "statement": "A pure holomorphic vacuum modular correction has even weight in {0,4,6,8,10,...}; weight 2 and every odd weight are absent.",
            "finite_size_consequence": "for a dimensionless torus response with x=weight, allowed N powers are integer N^-1,N^-2,N^-3,...; no half-integer power belongs to the ordinary C[E4,E6] vacuum ring",
            "first_shape_degeneracy": "weights 4,6,8,10 are one-dimensional; the first ordinary modular ambiguity is weight 12 (N^-5), span{E4^3,E6^2}",
            "weight12_elliptic_resolution": "E4^3 is visible at square and zero at hex; E6^2 is zero at square and visible at hex",
        },
        "spectrum_through_weight24": ring_spectrum(24),
        "frozen_score": {
            "input": "one complex observable over the three degree-2 hexagonal children in the registered order",
            "primary": "DFT support must be exactly r=1 for E4, r=0 for E6/D4E4, or r=2 for E4^2/D6E6",
            "mixture_test": "more than one nonzero DFT component is direct evidence for a modular mixture, quasimodular/Jordan tangent, or nonlocal sector",
            "amplitude_free": True,
        },
        "numerical_oracle": numerical_oracle(dps),
        "claim_boundary": {
            "exact": "Serre/Ramanujan identities, modular-ring weights, elliptic zeros, and degree-2 child characters",
            "conditional_bridge": "the lattice correction must be a typed complex/chiral torus response with the same frame across children",
            "not_claimed": "that every irrelevant correction is holomorphic or that a real-only child vector distinguishes conjugate characters",
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
