#!/usr/bin/env python3
"""Exact two-parameter P398 charge-one propagation from the eight row bonds."""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import json
from math import acosh, comb, log, sqrt
from pathlib import Path

from noncrossing_connectivity_codec import noncrossing_states
from p321_homology_trace_certificate import action_matrix, identity, join_adjacent
from p333_generic_q_detach_intertwiner import detach_jet, detach_state
from p333_gram_source_intertwiner import multiply, rref_solve
from p333_source_landing_doublet_width4 import landing_gram_jet
from p398_rooted_gr1_completion import selected_completion_families
from p398_positive_cylinder import correlation, zdet, zmul

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/p398-anisotropic-cylinder/latest.json"


def add(a, b, scale=1):
    out = defaultdict(F, a)
    for key, value in b.items():
        out[key] += scale*value
    return {key: value for key, value in out.items() if value}


def mul(a, b):
    out = defaultdict(F)
    for (ah, av), x in a.items():
        for (bh, bv), y in b.items():
            out[ah+bh, av+bv] += x*y
    return {key: value for key, value in out.items() if value}


def encode(poly):
    return [{"h": h, "v": v, "coefficient": str(value)}
            for (h, v), value in sorted(poly.items())]


def polynomial_text(poly):
    return " + ".join(f"({value})*h^{h}*v^{v}" for (h,v), value in sorted(poly.items())) or "0"


def evaluate(poly, h, v):
    return sum(value*h**ph*v**pv for (ph,pv), value in poly.items())


def factored_certificate(polynomials):
    one, h, v = {(0,0):F(1)}, {(1,0):F(1)}, {(0,1):F(1)}
    w, z = add(one,h,-1), add(one,v,-1)
    q = mul(mul(w,w),mul(v,v))
    m00 = add(add(one,h),mul(h,v),-2)
    expected = {
        "a":mul(q,m00), "b":mul(q,add({},z,-1)),
        "c":mul(q,add({},mul(h,v),-1)), "d":mul(q,v)
    }
    difference = add(m00,v,-1)
    disc = mul(mul(q,q),add(mul(difference,difference),mul(mul(h,v),z),8))
    determinant = one
    for _ in range(5):
        determinant = mul(determinant,mul(w,v))
    expected.update({"trace":add(expected["a"],expected["d"]),
                     "determinant":determinant,"discriminant":disc})
    assert expected == polynomials
    return {"q":"(1-h)^2*v^2", "a":"q*(1+h-2*h*v)",
            "b":"-q*(1-v)", "c":"-q*h*v", "d":"q*v",
            "trace":"q*(1+h+v-2*h*v)", "determinant":"(1-h)^5*v^5",
            "discriminant":"q^2*((1+h-v-2*h*v)^2+8*h*v*(1-v))"}


def stationary_covariance(h, v):
    states = noncrossing_states(4)
    ap = selected_completion_families()["rooted_charge1"]["columns"]
    gram,_ = landing_gram_jet()
    functions = [ap[row]+gram[row][-2:] for row in range(14)]
    horizontal, vertical = identity(14), identity(14)
    for site in range(4):
        join = action_matrix(4,lambda state, site=site:join_adjacent(state,site))
        detach,_ = detach_jet(4,site)
        horizontal = multiply([[(1-h)*int(r==c)+h*join[r][c] for c in range(14)]
                               for r in range(14)],horizontal)
        vertical = multiply([[v*int(r==c)+(1-v)*detach[r][c] for c in range(14)]
                             for r in range(14)],vertical)
    transfer = multiply(horizontal,vertical)
    coefficients = [[transfer[r][c]-int(r==c) for c in range(14)] for r in range(14)]
    solution = rref_solve(coefficients+[[1]*14],[0]*14+[1],14)
    assert solution["consistent"] and solution["dimension"] == 0
    probability = solution["particular"]
    assert min(probability) > 0
    covariance = correlation(probability,functions,functions)
    assert zdet(covariance)[0] > 0
    return covariance


def critical_row(h, polynomials):
    v = 1-h
    a,b,c,d = [evaluate(polynomials[k],h,v) for k in ("a","b","c","d")]
    gap = sqrt(float((a-d)**2+8*b*c))
    lp,lm = (float(a+d)+gap)/2, (float(a+d)-gap)/2
    covariance = stationary_covariance(h,v)
    gamma = acosh(1+float(h*h/v))
    ratio = lm/lp
    modes = []
    for sign, eigenvalue in ((1,lp),(-1,lm)):
        # Eigenfunction A+kappa L, eigenvalue q*mu. Real coefficient times (1-i).
        k = (float(h)-sign*sqrt(1+float(v*v)))/2
        kl = complex(k,-k)
        caa = float(covariance[0][0][0])
        cll = float(covariance[1][1][0])
        cal = complex(*map(float,covariance[0][1]))
        variance = caa+abs(kl)**2*cll+2*(cal*kl.conjugate()).real
        modes.append({"eigenvalue":eigenvalue,"decay_length_rows":-1/log(eigenvalue),
                      "kappa_re_im":[k,-k],"stationary_variance_eigenfunction":variance,
                      "pure_mode_autocorrelation_at_d8":variance*eigenvalue**8})
    return {"h":str(h),"v":str(v),"gamma":gamma,"eigenmodes_plus_minus":modes,
            "lambda_minus_over_plus":ratio,"relative_decay_at_d8":ratio**8,
            "max_d_before_fast_relative_decay_below_10pct":int(log(.1)/log(ratio)),
            "C0_complex_re_im":[[[str(x) for x in z] for z in row] for row in covariance]}


def signed_jordan_witness(polynomials):
    h,v = F(-1,8),F(1,2)
    a,b,c,d = [evaluate(polynomials[k],h,v) for k in ("a","b","c","d")]
    eigenvalue = (a+d)/2
    n = [[(a-eigenvalue,F(0)),(b,b)],[(c,-c),(d-eigenvalue,F(0))]]
    for i in range(2):
        for j in range(2):
            x,y = zmul(n[i][0],n[0][j]),zmul(n[i][1],n[1][j])
            assert (x[0]+y[0],x[1]+y[1]) == (0,0)
    assert any(any(z) for row in n for z in row)
    return {"h":str(h),"v":str(v),"eigenvalue":str(eigenvalue),
            "B_minus_lambda_I_complex_re_im":[[[str(x) for x in z] for z in row] for row in n],
            "nilpotent_square_zero":True,"nonzero_nilpotent":True,
            "physical":False,"reason":"negative horizontal bond weight"}


def derive_sector():
    states = noncrossing_states(4)
    ap = selected_completion_families()["rooted_charge1"]["columns"]
    gram, _ = landing_gram_jet()
    functions = [ap[row]+gram[row][-2:] for row in range(14)]
    lookup = dict(zip(states, functions))
    bernstein = {}
    for nh in range(5):
        for nv in range(5):
            bernstein[nh,nv] = {
                (nh+a,nv+b): F((-1)**(a+b)*comb(4-nh,a)*comb(4-nv,b))
                for a in range(5-nh) for b in range(5-nv)
            }
    expected = [[{} for _ in range(4)] for _ in states]
    for row, state in enumerate(states):
        counts = defaultdict(lambda: [F(0)]*4)
        for mask in range(256):
            following = state
            for site in range(4):
                if not mask & (1 << site):
                    following = detach_state(following, site)
            for site in range(4):
                if mask & (1 << (site+4)):
                    following = join_adjacent(following, site)
            key = ((mask >> 4).bit_count(), (mask & 15).bit_count())
            counts[key] = [a+b for a,b in zip(counts[key], lookup[following])]
        for key, values in counts.items():
            for col, value in enumerate(values):
                expected[row][col] = add(expected[row][col], bernstein[key], value)
    real = [[{} for _ in range(4)] for _ in range(4)]
    for ph in range(5):
        for pv in range(5):
            for col in range(4):
                rhs = [expected[row][col].get((ph,pv), F(0)) for row in range(14)]
                solution = rref_solve(functions, rhs, 4)
                assert solution["consistent"] and solution["dimension"] == 0
                for row, value in enumerate(solution["particular"]):
                    if value:
                        real[row][col][ph,pv] = value
    # f=(A,L), E[f(next)|pi] = f(pi) B. Reflection fixes the two phases.
    a, b, c, d = real[0][0], real[0][2], real[2][0], real[2][2]
    assert real[0][1] == real[2][3] == {}
    assert real[0][3] == b and real[2][1] == add({}, c, -1)
    trace = add(a,d)
    determinant = add(mul(a,d), mul(b,c), -2)
    discriminant = add(mul(trace,trace), determinant, -4)
    return {"a":a,"b":b,"c":c,"d":d,"trace":trace,
            "determinant":determinant,"discriminant":discriminant}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, default=OUTPUT)
    args = parser.parse_args()
    polynomials = derive_sector()
    factored = factored_certificate(polynomials)
    value = {"schema":"matching-one/p398-anisotropic-cylinder/v1",
             "convention":"B=[[a,(1+i)b],[(1-i)c,d]], columns are output observables",
             "parameters":"horizontal h, vertical v, both in (0,1)",
             "polynomials": {k: encode(v) for k,v in polynomials.items()},
             "factorizations":factored,
             "positive_domain":"discriminant strictly positive and determinant positive; two distinct positive ordinary modes, no Jordan collision",
             "symmetrizer":"S=diag(1,sqrt(h*v/(1-v))); S^-1*(B/q)*S is Hermitian",
             "collision_set_closed_probability_square":"h=1 or v=0 gives B=0; (h,v)=(0,1) gives B=I; all collisions are scalar/semisimple",
             "signed_collision_curves":"h=-(1-v)/(1+sqrt(2*v))^2 or h=-(1-v)/(1-sqrt(2*v))^2 for 0<v<1; at v=1/2 only finite root h=-1/8",
             "signed_Jordan_witness":signed_jordan_witness(polynomials),
             "critical_line":{"equation":"h+v=1",
                 "eigenvalues":"(1-h)^4*(1-h+h^2 +/- h*sqrt(2-2*h+h^2))",
                 "rapidity":"gamma=acosh(1+h^2/(1-h)); lambda_pm=(1-h)^5*exp(+/-gamma)",
                 "generator_h_to_0":"(B-I)/h -> [[-5,-1-i],[-1+i,-5]], eigenvalues -5+/-sqrt(2)",
                 "eigenfunctions":"Psi_pm=A+(h -/+ sqrt(1+(1-h)^2))*(1-i)*L/2",
                 "samples":[critical_row(F(n,d),polynomials) for n,d in ((1,20),(1,10),(1,5),(1,2),(4,5))]},
             "input_sha256":{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
                 for path in [Path(__file__),*(ROOT/"scripts"/name for name in (
                     "p398_positive_cylinder.py", "noncrossing_connectivity_codec.py",
                     "p321_homology_trace_certificate.py", "p333_generic_q_detach_intertwiner.py",
                     "p333_gram_source_intertwiner.py", "p333_source_landing_doublet_width4.py",
                     "p398_rooted_gr1_completion.py"))]},
             "interpretation":"Anisotropy resolves row-distance mixtures, not a Jordan approach in rescaled physical time. At h->0 with t=h*d, two masses 5-/+sqrt(2) remain distinct; no cost-normalized Monte Carlo gain was measured.",
             "scope":"Same finite square-bond width-four positive cylinder; no continuum/Jordan identity in a different generator or source."}
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(value,indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"factorizations":factored,"positive_domain":value["positive_domain"]},indent=2))
    print("h v lambda+ lambda- ratio fast_C(d=8)")
    for row in value["critical_line"]["samples"]:
        plus,minus=row["eigenmodes_plus_minus"]
        print(row["h"],row["v"],plus["eigenvalue"],minus["eigenvalue"],
              row["lambda_minus_over_plus"],minus["pure_mode_autocorrelation_at_d8"])
    print("Signed Jordan witness: h=-1/8, v=1/2, lambda=243/1024 (nonphysical)")


if __name__ == "__main__":
    main()
