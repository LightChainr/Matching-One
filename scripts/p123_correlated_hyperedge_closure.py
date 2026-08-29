#!/usr/bin/env python3
"""Exact four-terminal duality closure and a two-block hyperedge no-go."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path


Poly = tuple[Fraction, ...]
MULTIPLICITIES = (1, 4, 4, 1, 2, 2)


def poly(*coefficients: int | Fraction) -> Poly:
    values = [Fraction(value) for value in coefficients]
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values)


def add(*values: Poly) -> Poly:
    degree = max(map(len, values))
    return poly(*(sum(value[i] if i < len(value) else 0 for value in values) for i in range(degree)))


def scale(value: Poly, scalar: int | Fraction) -> Poly:
    return poly(*(Fraction(scalar) * item for item in value))


def multiply(left: Poly, right: Poly) -> Poly:
    output = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            output[i + j] += a * b
    return poly(*output)


def evaluate(value: Poly, point: Fraction) -> Fraction:
    result = Fraction(0)
    for coefficient in reversed(value):
        result = result * point + coefficient
    return result


def ptext(value: Poly) -> str:
    terms = []
    for degree, coefficient in enumerate(value):
        if not coefficient:
            continue
        factor = "" if degree == 0 else ("t" if degree == 1 else f"t^{degree}")
        if degree == 0:
            terms.append(str(coefficient))
        elif coefficient == 1:
            terms.append(factor)
        elif coefficient == -1:
            terms.append(f"-{factor}")
        else:
            terms.append(f"{coefficient}*{factor}")
    return " + ".join(terms).replace("+ -", "- ") if terms else "0"


def dual(probabilities: tuple[Poly, ...]) -> tuple[Poly, ...]:
    p1, p2, p3, p4, p5, p6 = probabilities
    return p4, p3, p2, p1, p6, p5


def normalization(probabilities: tuple[Poly, ...]) -> Poly:
    return add(*(scale(value, multiplicity) for value, multiplicity in zip(probabilities, MULTIPLICITIES)))


def duality_residual(probabilities: tuple[Poly, ...]) -> tuple[Poly, Poly, Poly]:
    p1, p2, p3, p4, p5, p6 = probabilities
    return add(p4, scale(p1, -1)), add(p3, scale(p2, -1)), add(p6, scale(p5, -1))


def independent_two_block_family() -> tuple[Poly, ...]:
    """C4 average of two all-or-none four-hyperedges sharing two internal vertices."""
    t = poly(0, 1)
    one_minus_t = poly(1, -1)
    p1 = multiply(one_minus_t, one_minus_t)
    p2 = scale(multiply(t, one_minus_t), Fraction(1, 2))
    p4 = multiply(t, t)
    zero = poly(0)
    return p1, p2, zero, p4, zero, zero


def rational_record(value: Fraction) -> dict[str, int | str]:
    return {"numerator": value.numerator, "denominator": value.denominator, "text": str(value)}


def selfdual_simplex_point(u: Fraction, v: Fraction) -> tuple[Fraction, ...]:
    w = Fraction(1, 2) - 4 * u - 2 * v
    if min(w, u, v) < 0:
        raise ValueError("point lies outside the probability simplex")
    return w, u, u, w, v, v


def build_oracle() -> dict:
    all_none = (poly(1, -1), poly(0), poly(0), poly(0, 1), poly(0), poly(0))
    composed = independent_two_block_family()
    residual = duality_residual(composed)
    sample_points = []
    for u, v in ((Fraction(0), Fraction(0)), (Fraction(1, 32), Fraction(1, 32)), (Fraction(1, 16), Fraction(0))):
        point = selfdual_simplex_point(u, v)
        sample_points.append({
            "u": rational_record(u),
            "v": rational_record(v),
            "P": [rational_record(value) for value in point],
            "normalization": rational_record(sum(m * x for m, x in zip(MULTIPLICITIES, point))),
            "dual_fixed": point == (point[3], point[2], point[1], point[0], point[5], point[4]),
        })

    return {
        "schema": "matching-one.p123-correlated-hyperedge-closure.v1",
        "issue": 123,
        "canonical_coordinates": {
            "orbit_representatives": [
                "P1=A|B|C|D",
                "P2=one nearest-neighbour pair",
                "P3=three connected",
                "P4=ABCD",
                "P5=one diagonal pair",
                "P6=two noncrossing pairs",
            ],
            "multiplicities": list(MULTIPLICITIES),
            "normalization": "P1+4P2+4P3+P4+2P5+2P6=1",
            "duality": "(P1,P2,P3,P4,P5,P6)->(P4,P3,P2,P1,P6,P5)",
            "selfdual_equations": ["P4=P1", "P3=P2", "P6=P5"],
            "source": "Damavandi-Ziff arXiv:1506.06125v2, equations (2)-(6)",
        },
        "exact_square_site_embedding": {
            "tensor": [ptext(value) for value in all_none],
            "meaning": "a Bernoulli site is exactly the all-or-none correlated four-terminal hyperedge",
            "normalization": ptext(normalization(all_none)),
            "duality_residual": [ptext(value) for value in duality_residual(all_none)],
            "selfdual_intersection": "t=1/2 only",
            "transversality": "d_t(P4-P1)=2",
        },
        "minimal_composition_obstruction": {
            "construction": (
                "two all-or-none four-hyperedges share two internal vertices; average the AB|CD and BC|DA "
                "orientations to restore C4"
            ),
            "boundary_tensor": [ptext(value) for value in composed],
            "normalization": ptext(normalization(composed)),
            "generated_partial_state": "P2=t(1-t)/2 per nearest-neighbour orbit member",
            "duality_residual": [ptext(value) for value in residual],
            "no_common_selfdual_point": True,
            "proof": (
                "P3-P2=-t(1-t)/2 forces t=0 or 1, while P4-P1=2t-1 forces t=1/2"
            ),
        },
        "finite_local_correlated_no_go": {
            "family": (
                "for either C4 orientation let Prob(00)=a, Prob(10)=Prob(01)=b, Prob(11)=c, "
                "with a+2b+c=1"
            ),
            "boundary_coordinates": "P1=a, P2=b/2, P3=0, P4=c, P5=P6=0",
            "selfduality_solution": "b=0 and a=c=1/2",
            "meaning": (
                "any nonzero single-block state generates an unpaired P2 direction; selfduality forces perfect "
                "00/11 correlation and collapses the composite back to the all-or-none t=1/2 tensor"
            ),
        },
        "maximal_positive_correlated_family": {
            "coordinates": "P=(w,u,u,w,v,v), w=1/2-4u-2v",
            "domain": "u>=0, v>=0, 4u+2v<=1/2",
            "selfdual_configurationwise_in_coordinate_space": True,
            "sample_exact_points": sample_points,
            "all_connected_bound": "P4=w<=1/2",
            "consequence": (
                "adding dual-paired partial states yields a genuine two-dimensional selfdual correlated-cell "
                "manifold, but moves P4 downward from 1/2 and cannot preserve an all-connected marginal above 1/2"
            ),
        },
        "alternating_dual_embedding": {
            "family": "place any tensor P on one hyperedge colour and dual(P) on the other",
            "global_duality": "duality exchanges the two cell colours, so the alternating model is selfdual",
            "all_none_specialization": "A(t) alternates with A(1-t)",
            "homogeneous_intersection": "A(t)=A(1-t) iff t=1/2",
            "structural_result": (
                "a selfdual embedding of an arbitrary site probability exists only as an inhomogeneous dual pair; "
                "it does not impose an exact equation on the homogeneous square-site threshold"
            ),
        },
        "decision": {
            "result": "finite-local obstruction for the homogeneous Bernoulli slice, plus an exact correlated positive family",
            "square_site_implication": (
                "local four-edge duality either fixes the exact Bernoulli tensor at 1/2, generates partial "
                "connectivity couplings, or alternates p with 1-p; none yields a new exact constraint on homogeneous square-site pc"
            ),
            "next_live_route": (
                "a successful exact embedding must use a nonlocal projection or a larger cell whose terminal "
                "semantics no longer identify P4 directly with the site occupation probability"
            ),
        },
        "claim_boundary": {
            "proved": [
                "the all-or-none embedding and its unique local selfdual point",
                "the exact two-block composition tensor and no common selfdual root",
                "the no-go for the full correlated (a,b,b,c) two-block family",
                "the dual-fixed two-dimensional probability simplex",
            ],
            "not_proved": [
                "a no-go for arbitrarily large decorated cells or nonlocal projections",
                "a Yang-Baxter obstruction",
                "any new numerical or closed-form value of square-site pc",
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
