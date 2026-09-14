#!/usr/bin/env python3
"""Subleading spectrum and logit derivatives of the fixed-width charge transfer.

Imports the transparent safe frontier kernel from fixed_width_charge_transfer.py.
The outputs are finite-width numerical controls for the charge-sector scaling
notes; they are not new pc estimates.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.sparse.linalg import eigs

from fixed_width_charge_transfer import SafeTransfer


def logistic(h: float) -> float:
    if h >= 0:
        z = math.exp(-h)
        return 1.0 / (1.0 + z)
    z = math.exp(h)
    return z / (1.0 + z)


def leading_moduli(transfer: SafeTransfer, p: float, count: int = 2) -> list[float]:
    matrix = transfer.matrix(p)
    n = matrix.shape[0]
    if n <= count + 1:
        vals = np.linalg.eigvals(matrix.toarray())
    else:
        vals = eigs(
            matrix,
            k=count,
            which="LM",
            return_eigenvectors=False,
            tol=1e-12,
            maxiter=200000,
        )
    moduli = sorted((abs(complex(v)) for v in vals), reverse=True)
    return [float(x) for x in moduli[:count]]


def I_of_h(transfer: SafeTransfer, h: float) -> float:
    return transfer.perron(logistic(h))["I0"]


def charge_of_h(g4: SafeTransfer, g8: SafeTransfer, h: float) -> float:
    # Complement white logit is -h exactly.
    return I_of_h(g4, h) - I_of_h(g8, -h)


def central_derivatives(g4: SafeTransfer, g8: SafeTransfer, h: float, eps: float):
    f0 = charge_of_h(g4, g8, h)
    fm1 = charge_of_h(g4, g8, h - eps)
    fp1 = charge_of_h(g4, g8, h + eps)
    fm2 = charge_of_h(g4, g8, h - 2 * eps)
    fp2 = charge_of_h(g4, g8, h + 2 * eps)
    d1 = (fp1 - fm1) / (2 * eps)
    d2 = (fp1 - 2 * f0 + fm1) / eps**2
    d3 = (fp2 - 2 * fp1 + 2 * fm1 - fm2) / (2 * eps**3)
    return d1, d2, d3


def individual_second(transfer: SafeTransfer, h: float, eps: float) -> float:
    return (
        I_of_h(transfer, h + eps)
        - 2 * I_of_h(transfer, h)
        + I_of_h(transfer, h - eps)
    ) / eps**2


def record(width: int, *, eps: float) -> dict[str, object]:
    g4 = SafeTransfer(width, matching=False)
    g8 = SafeTransfer(width, matching=True)

    def theta_p(p: float) -> float:
        return g4.perron(p)["I0"] - g8.perron(1.0 - p)["I0"]

    root = brentq(theta_p, 0.5000000001, 0.999999999, xtol=2e-14)
    hroot = math.log(root / (1.0 - root))
    d1, d2, d3 = central_derivatives(g4, g8, hroot, eps)

    # Repeat at half the step as a finite-difference stability control.
    d1b, d2b, d3b = central_derivatives(g4, g8, hroot, eps / 2)

    c4 = individual_second(g4, hroot, eps)
    c8 = individual_second(g8, -hroot, eps)

    spec4 = leading_moduli(g4, root, 2)
    spec8 = leading_moduli(g8, 1.0 - root, 2)
    gap4 = math.log(spec4[0] / spec4[1])
    gap8 = math.log(spec8[0] / spec8[1])

    return {
        "width": width,
        "charge_root": root,
        "logit_root": hroot,
        "finite_difference_step": eps,
        "charge_logit_derivatives": {
            "d1": d1,
            "d2": d2,
            "d3": d3,
            "half_step": {"d1": d1b, "d2": d2b, "d3": d3b},
            "scaled": {
                "d1_times_w_quarter": d1 * width ** 0.25,
                "d2_times_sqrt_w": d2 * math.sqrt(width),
                "d3_over_w_5_over_4": d3 / width ** 1.25,
            },
        },
        "individual_logit_second_derivatives": {
            "G4": c4,
            "G8_complement": c8,
            "G4_over_sqrt_w": c4 / math.sqrt(width),
            "G8_over_sqrt_w": c8 / math.sqrt(width),
        },
        "safe_spectral_relaxation": {
            "G4_lambda1": spec4[0],
            "G4_lambda2_modulus": spec4[1],
            "G4_gap": gap4,
            "G4_gap_times_w": gap4 * width,
            "G8_lambda1": spec8[0],
            "G8_lambda2_modulus": spec8[1],
            "G8_gap": gap8,
            "G8_gap_times_w": gap8 * width,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-width", type=int, default=4)
    parser.add_argument("--max-width", type=int, default=8)
    parser.add_argument("--eps", type=float, default=7e-4)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    result = {
        "schema": "fixed-width-charge-spectrum-derivatives-v1",
        "claim_boundary": [
            "Finite differences are controls; exact first-derivative identities are in the accompanying notes.",
            "Subleading safe-kernel eigenvalues are not yet formally identified with every periodic TL trace eigenvalue.",
            "No continuum exponent is inferred solely from this file."
        ],
        "records": [record(w, eps=args.eps) for w in range(args.min_width, args.max_width + 1)],
    }
    text = json.dumps(result, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
