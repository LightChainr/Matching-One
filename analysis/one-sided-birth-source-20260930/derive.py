#!/usr/bin/env python3
"""Transport existing exact responses into an actual moving-root E_top map.

Reads existing finite tables; no topology census or new sampling. Exact
positive coefficients prove the sign. Root decimals are illustrations only.
"""
from fractions import Fraction as F
import hashlib
import json
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = HERE.parent / "completion-pair-synergy-20260929/result.json"
KICK = HERE.parent / "birth-selection-mechanism-20260930/result.json"


def main():
    model = next(x for x in json.loads(BASE.read_text())["models"] if x["lattice"] == "square")
    start = next(x for x in json.loads(KICK.read_text())["rows"] if x["name"] == "A")
    points = [(v % 4, v//4) for v in start["occupied"]]
    orbit = set()
    for reflected in [False, True]:
        for rotation in range(4):
            transformed = []
            for x, y in points:
                if reflected:
                    x = -x
                for _ in range(rotation):
                    x, y = -y, x
                transformed.append((x, y))
            for tx in range(4):
                for ty in range(4):
                    orbit.add(tuple(sorted((x+tx) % 4+4*((y+ty) % 4) for x, y in transformed)))
    activation = F(len(orbit), comb(16, 6))
    h_to_derivative = {x["lag_insertions"]: F(x["first_kick_derivative_at_zero"])
                       for x in start["all_lags_single_kick"]}
    bernstein = [F(0) if k <= 6 else -activation*h_to_derivative[k-6] for k in range(17)]
    powers = [comb(16, k)*v for k, v in enumerate(bernstein)]
    assert all(v >= 0 for v in bernstein) and any(v > 0 for v in bernstein)
    counts = model["rank_counts"]
    def probabilities(p):
        prob = [0., 0., 0.]
        derivative = [0., 0., 0.]
        for k, rank, number in counts:
            prob[rank] += number*p**k*(1-p)**(16-k)
            derivative[rank] += number*((k*p**(k-1)*(1-p)**(16-k) if k else 0)
                                       - ((16-k)*p**k*(1-p)**(15-k) if k < 16 else 0))
        return prob, derivative
    lo, hi = .01, .99
    for _ in range(40):
        p = (lo+hi)/2
        prob, _ = probabilities(p)
        if prob[2]-prob[0] > 0:
            hi = p
        else:
            lo = p
    p = (lo+hi)/2
    prob, derivative = probabilities(p)
    g = sum(float(v)*p**k*(1-p)**(16-k) for k, v in enumerate(powers))
    mp = derivative[2]-derivative[0]
    normal_factor = -2*derivative[0]/mp
    output = {
        "schema": "matching-one.one-sided-birth-source.v1", "date": "2026-09-30",
        "source": "At count6 only, kick when occupied set is in the translation/dihedral orbit of A={0,1,4,5,8,12}; otherwise uniform. Uniform continuation after the kick.",
        "source_model": "square NN occupied-site L4; theta0 is the original uniform-permutation/iid-label model",
        "orbit_size": len(orbit), "prefix_activation_probability": str(activation),
        "delta_P2_Bernstein_degree16_coefficients": [str(v) for v in bernstein],
        "delta_P2_coefficients_of_p_k_times_1_minus_p_16_minus_k": [str(v) for v in powers],
        "exact_sign": "All coefficients nonnegative and some positive: delta_P2(p)>0 for every0<p<1. Delta_P0=0 exactly.",
        "root_illustration": {"p_balance": p, "baseline_Etop": prob[0]+prob[2],
                              "M_prime": mp, "P0_prime": derivative[0], "delta_M": g, "delta_Etop": g,
                              "root_derivative": -g/mp, "moving_root_Etop_derivative": g*normal_factor,
                              "normal_factor": normal_factor},
        "source_artifact_sha256": {str(f.relative_to(HERE.parent.parent)): hashlib.sha256(f.read_bytes()).hexdigest()
                                   for f in [BASE, KICK]},
        "scope": "Exact finite, symmetry-invariant constructed dynamic source. Not natural-archive evidence, original norm4 source, U six-coordinate normalization, large-size scaling or CFT identification."
    }
    (HERE / "result.json").write_text(json.dumps(output, indent=2)+"\n")
    lines = ["# One-sided birth source: moving-root E_top response", "",
             f"Square L4, one count6 kick on {len(orbit)} symmetry-related geometries; activation probability {activation}.",
             "First birth unchanged; deltaP0=0 and deltaM=deltaEtop=deltaP2=g(p)>0 at every interior p.", "",
             "Exact g(p), in positive p^k(1-p)^(16-k) terms:", ""]
    lines += [f"- k={k}: coefficient {v}" for k, v in enumerate(powers) if v]
    lines += ["", "Illustrative floating-point evaluation at the baseline balance root:", ""]
    lines += [f"- {key}: {value:.10g}" for key, value in output["root_illustration"].items()]
    lines += ["", "Sign proof uses exact coefficients and strict monotonicity, not decimal root precision.",
              "No old census was rerun. This is a constructed dynamic source, not the original norm4/U identification experiment.", ""]
    (HERE / "RESULT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
