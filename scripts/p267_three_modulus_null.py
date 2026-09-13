#!/usr/bin/env python3
"""N100 three-modulus affine-shape nulls in Q(sqrt(2)); no simulation."""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
import json
from math import sqrt
from pathlib import Path


@dataclass(frozen=True)
class Q2:
    """Only the quadratic field needed for the three fixed shapes."""
    a: F = F(0)
    b: F = F(0)

    @staticmethod
    def cast(value):
        return value if isinstance(value, Q2) else Q2(F(value))

    def __add__(self, other):
        other = self.cast(other)
        return Q2(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self):
        return Q2(-self.a, -self.b)

    def __sub__(self, other):
        return self + -self.cast(other)

    def __rsub__(self, other):
        return self.cast(other) + -self

    def __mul__(self, other):
        other = self.cast(other)
        return Q2(self.a*other.a + 2*self.b*other.b,
                  self.a*other.b + self.b*other.a)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.cast(other)
        norm = other.a**2 - 2*other.b**2
        if not norm:
            raise ZeroDivisionError
        return self * Q2(other.a/norm, -other.b/norm)

    def as_float(self):
        return float(self.a) + float(self.b)*sqrt(2)

    def payload(self):
        return {"rational": str(self.a), "sqrt2": str(self.b),
                "decimal": self.as_float()}


S = Q2(F(0), F(1))
G2 = Q2(F(11, 4))
G4 = (91+60*S)/16
GS = (91-60*S)/16
SHAPES = {
    "affine_E4": (G2, G4, GS),
    "affine_height_only_E4": (G2, G4, Q2(F(1))),
    "affine_y_squared": tuple(map(Q2.cast, (4, 16, 1))),
}
WEIGHTS = {
    "affine_E4": (120*S, 47-60*S, -47-60*S),
    "affine_height_only_E4": (75+60*S, Q2(F(-28)), -47-60*S),
    "affine_y_squared": tuple(map(Q2.cast, (5, -1, -4))),
}
FIELDS = ("A_top", "E_top", "C", "W")


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), Q2())


def slope(model):
    a, b, c = SHAPES[model]
    return (c-a)/(b-a)


def normalized_weights(model):
    r = slope(model)
    return (r-1, -r, Q2(F(1)))


def projection_matrix(model):
    w = normalized_weights(model)
    return [[w[k//4].as_float() if k % 4 == j else 0.0
             for k in range(12)] for j in range(4)]


def project_joint(mean, covariance, model):
    """Apply the fixed four-vector null and full 12x12 covariance, no fitting."""
    if len(mean) != 12 or len(covariance) != 12 or any(len(r) != 12 for r in covariance):
        raise ValueError("need shape-major 12-vector and full 12x12 covariance")
    L = projection_matrix(model)
    residual = [sum(L[i][a]*mean[a] for a in range(12)) for i in range(4)]
    out_cov = [[sum(L[i][a]*covariance[a][b]*L[j][b]
                    for a in range(12) for b in range(12))
                for j in range(4)] for i in range(4)]
    return residual, out_cov


def certificate():
    # Normalize theta3(i)^4=1, theta4(i)^4=1/2.  Duplication gives these
    # theta3(2i)^4 and theta4(2i)^4.  The exact E4 formula then yields G4.
    u, v = (3+2*S)/8, S/2
    duplicated_G4 = (u*u + 14*u*v + v*v)*F(4, 3)
    assert duplicated_G4 == G4
    assert G4+GS == Q2(F(91, 8))
    output = {}
    for name, g in SHAPES.items():
        w = WEIGHTS[name]
        assert sum(w, Q2()) == Q2()
        assert dot(w, g) == Q2()
        nw = normalized_weights(name)
        assert sum(nw, Q2()) == dot(nw, g) == Q2()
        output[name] = {
            "shape_coordinates": [x.payload() for x in g],
            "unnormalized_exact_weights": [x.payload() for x in w],
            "normalized_weights_shear_coefficient_one": [x.payload() for x in nw],
            "signed_secant_ratio_r": slope(name).payload(),
            "offset_residual_exact": sum(w, Q2()).payload(),
            "shape_residual_exact": dot(w, g).payload(),
            "joint_projection_matrix_L": projection_matrix(name),
        }
    comparisons = []
    names = list(SHAPES)
    for i, first in enumerate(names):
        for second in names[i+1:]:
            gap = slope(first)-slope(second)
            assert gap != Q2()
            comparisons.append({
                "true_shape": first, "tested_shape": second,
                "normalized_wrong_model_residual_over_Y4_minus_Y2": gap.payload(),
                "intersection": "constant three-shape response only, independently in each field",
                "illustrative_3sigma_SE_over_absolute_span": abs(gap.as_float())/3,
            })
    return {
        "schema": "matching-one/p267-three-modulus-null/v1",
        "geometry_source_commit": "b9e4ea19bc585cbed18ec6ba1d13e85f2b5accc7",
        "N": 100, "shape_order": ["2i", "4i", "1/2+i"],
        "field_order": list(FIELDS),
        "stacking": "shape-major: [Y2_A,Y2_E,Y2_C,Y2_W,Y4_A,...,Ys_W]",
        "measure": "square-site fixed p_ref=0.59274605079; same signed chi4-normalized quotient contrasts as ce01e4d",
        "smith_pair": [[1, 100], [5, 20]],
        "coupling_assumption": "Y_j(tau)=a_j+b_j*f(tau), at this N and p, with arbitrary coordinate-specific a_j,b_j constant across only the three declared shapes",
        "models": output, "pairwise_discriminators": comparisons,
        "theta_duplication_G4": duplicated_G4.payload(),
        "area_offset_free_observable": "R_f=Ys-Y2-r_f*(Y4-Y2); use linear null, not a noisy empirical ratio",
        "covariance_contract": "R=L mean; Cov(R)=L Sigma L^T with full same-stream cross-shape/cross-field covariance; no independent-evidence multiplication",
        "degeneracy": "Y4-Y2=0 in every field makes the three affine shape hypotheses coincide; absence of a resolved span is underpowered, not acceptance",
        "boundary": "Conditional finite-lattice coupling test, not automatic transfer of E4 to A/E/C/W, not an exponent, field identity, or approval for acquisition. No N50 estimates or area scale are used.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    payload = json.dumps(certificate(), indent=2)+"\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(payload)
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
