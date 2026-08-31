#!/usr/bin/env python3
"""P398 continuous-distance kernel and amplitude-free two-channel fingerprints."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from math import acosh, atan2, cosh, exp, log, sinh, sqrt
from pathlib import Path

from noncrossing_connectivity_codec import noncrossing_states
from p321_homology_trace_certificate import action_matrix, join_adjacent
from p333_generic_q_detach_intertwiner import detach_jet
from p333_gram_source_intertwiner import rref_solve
from p333_source_landing_doublet_width4 import landing_gram_jet
from p398_rooted_gr1_completion import selected_completion_families
from p398_positive_cylinder import correlation, encode, zadd, zconj, zdet, zmatmul, zmul
from p398_anisotropic_cylinder import derive_sector

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/"results/p398-continuous-kernel/latest.json"


def corner_jet(poly):
    """Constant and h coefficient after v=1-kappa*h, as exact kappa polynomials."""
    constant=F(0)
    slope=[F(0),F(0)]
    for (ph,pv),coefficient in poly.items():
        if ph==0:
            constant+=coefficient
            slope[1]-=pv*coefficient
        elif ph==1:
            slope[0]+=coefficient
    return constant,slope


def derive_generator():
    sector=derive_sector()
    jet={key:corner_jet(sector[key]) for key in ("a","b","c","d")}
    assert jet=={"a":(F(1),[F(-3),F(-2)]),"b":(F(0),[F(0),F(-1)]),
                 "c":(F(0),[F(-1),F(0)]),"d":(F(1),[F(-2),F(-3)])}
    return jet


def stationary_limit(kappa):
    """Exact stationary law of the existing 14-state continuous row generator."""
    states=noncrossing_states(4)
    generator=[[F(0) for _ in range(14)] for _ in range(14)]
    for site in range(4):
        join=action_matrix(4,lambda state,site=site:join_adjacent(state,site))
        detach,_=detach_jet(4,site)
        for r in range(14):
            for c in range(14):
                generator[r][c]+=join[r][c]+kappa*detach[r][c]-(1+kappa)*int(r==c)
    solution=rref_solve(generator+[[1]*14],[0]*14+[1],14)
    assert solution["consistent"] and solution["dimension"]==0
    probability=solution["particular"]
    assert min(probability)>0
    ap=selected_completion_families()["rooted_charge1"]["columns"]
    gram,_=landing_gram_jet()
    functions=[ap[r]+gram[r][-2:] for r in range(14)]
    covariance=correlation(probability,functions,functions)
    assert zdet(covariance)[0]>0
    # C'(0)=-C0 conjugate(M); record reversibility only in these readouts.
    mass=mass_matrix(kappa)
    product=zmatmul(covariance,[[zconj(z) for z in row] for row in mass])
    lag_skew=zadd(product[0][1],tuple(-x for x in zconj(product[1][0])))
    return {"stationary_probability":probability,"C0_complex_re_im":covariance,
            "C0_conjugate_M_offdiagonal_skew":lag_skew}


def mass_matrix(kappa):
    return [[(3+2*kappa,F(0)),(kappa,kappa)],[(F(1),F(-1)),(2+3*kappa,F(0))]]


def kernel(kappa,distance):
    """U=C0^-1 C(s)=exp(-conjugate(M)*s); decimal display, not derivation."""
    k=float(kappa)
    mean=5*(1+k)/2
    gap=sqrt(k*k+6*k+1)
    half=gap/2
    mass=[[3+2*k,k*(1-1j)],[1+1j,2+3*k]]
    scale=exp(-mean*distance)
    c,s=cosh(half*distance),sinh(half*distance)/half
    return [[scale*((c if r==col else 0)-s*(mass[r][col]-(mean if r==col else 0)))
             for col in range(2)] for r in range(2)]


def display_complex(matrix):
    return [[[complex(z).real,complex(z).imag] for z in row] for row in matrix]


def sample(kappa):
    k=float(kappa)
    gap=sqrt(k*k+6*k+1)
    slow,fast=(5*(1+k)-gap)/2,(5*(1+k)+gap)/2
    points=[]
    for distance in (.1,.5,1.):
        u=kernel(kappa,distance)
        trace=complex(u[0][0]+u[1][1]).real
        determinant=complex(u[0][0]*u[1][1]-u[0][1]*u[1][0]).real
        off_product=complex(u[0][1]*u[1][0]).real
        fingerprint=2*acosh(trace/(2*sqrt(determinant)))/(-log(determinant))
        mixing=4*off_product/(trace*trace-4*determinant)
        signed_axis=complex(u[0][0]-u[1][1]).real/sqrt(trace*trace-4*determinant)
        cross_ratio=off_product/complex(u[0][0]*u[1][1]).real
        points.append({"distance":distance,"U_complex_re_im":display_complex(u),
                       "trace":trace,"determinant":determinant,
                       "metric_free_mass_split":fingerprint,"diagonal_gauge_free_mixing":mixing,
                       "signed_channel_axis":signed_axis,"channel_cross_ratio":cross_ratio})
    return {"kappa":str(kappa),"masses_slow_fast":[slow,fast],
            "mass_ratio_fast_slow":fast/slow,
            "mixing_angle_radians":atan2(2*sqrt(2*k),1-k)/2,
            "metric_free_mass_split":gap/(5*(1+k)),
            "diagonal_gauge_free_mixing":8*k/(k*k+6*k+1),
            "stationary_limit":stationary_limit(kappa),"kernel_samples":points}


def build_result():
    jet=derive_generator()
    paths=[Path(__file__),*(ROOT/"scripts"/name for name in (
        "p398_anisotropic_cylinder.py","p398_positive_cylinder.py",
        "noncrossing_connectivity_codec.py","p321_homology_trace_certificate.py",
        "p333_generic_q_detach_intertwiner.py","p333_gram_source_intertwiner.py",
        "p333_source_landing_doublet_width4.py","p398_rooted_gr1_completion.py"))]
    return encode({"schema":"matching-one/p398-continuous-kernel/v1",
        "parent":"b35e100a3903c706dceba57c4667386eb4510ac3",
        "limit":"h=epsilon, v=1-kappa*epsilon, s=epsilon*d fixed, kappa>0",
        "exact_corner_jet":jet,
        "mass_matrix":"M=[[3+2*kappa,(1+i)*kappa],[1-i,2+3*kappa]]",
        "masses":"m_slow,fast=[5*(1+kappa)-/+sqrt(kappa^2+6*kappa+1)]/2",
        "kernel":"U(s)=C0^-1 C(s)=exp(-conjugate(M)*s)",
        "hermitian_gauge":"D=diag(1,exp(-i*pi/4)/sqrt(kappa)); D^-1 M D=[[3+2*kappa,sqrt(2*kappa)],[sqrt(2*kappa),2+3*kappa]]",
        "cross_ratio":"X=U12*U21/(U11*U22)=a*sinh(b*s)^2/[1+a*sinh(b*s)^2], a=8*kappa/(kappa^2+6*kappa+1), b=sqrt(kappa^2+6*kappa+1)/2",
        "similarity_and_metric_invariant":"2*acosh(tr(U)/(2*sqrt(det(U))))/[-log(det(U))]=sqrt(kappa^2+6*kappa+1)/[5*(1+kappa)]",
        "diagonal_gauge_and_metric_invariant":"4*U12*U21/[tr(U)^2-4*det(U)]=8*kappa/(kappa^2+6*kappa+1)",
        "joint_invariant_relation":"25*I_mass^2*(2-I_mix)=2",
        "signed_channel_axis":"(U11-U22)/sqrt(tr(U)^2-4*det(U))=(kappa-1)/sqrt(kappa^2+6*kappa+1)",
        "critical_full_kernel":{"C0":"[[6,-4+4i],[-4-4i,6]]/7",
            "channels":"psi_slow=(A-exp(-i*pi/4)*L)/sqrt(2); psi_fast=(A+exp(-i*pi/4)*L)/sqrt(2)",
            "covariance_in_channels":"diag((6+4*sqrt(2))/7*exp(-(5-sqrt(2))*s),(6-4*sqrt(2))/7*exp(-(5+sqrt(2))*s))",
            "normalized_cross_ratio":"tanh(sqrt(2)*s)^2",
            "raw_cross_ratio":"tanh(sqrt(2)*s+acosh(3))^2"},
        "dual_exchange":"H(kappa)=kappa*sigma_x*H(1/kappa)*sigma_x for the Hermitian-gauge mass matrix",
        "samples":[sample(F(n,d)) for n,d in ((1,4),(1,2),(1,1),(2,1),(4,1))],
        "input_sha256":{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},
        "scope":"Exact finite width-four square-bond family. Mass units depend on the common clock rate; kappa is surviving join/detach-rate detuning (only kappa=1 is self-dual). Ratios/invariants depend on kappa and are model-specific, not continuum universal data. Zero MC."})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json",type=Path,default=OUTPUT)
    args=parser.parse_args()
    value=build_result()
    args.json.parent.mkdir(parents=True,exist_ok=True)
    args.json.write_text(json.dumps(value,indent=2)+"\n",encoding="utf-8")
    for row in value["samples"]:
        print(row["kappa"],row["masses_slow_fast"],row["metric_free_mass_split"],
              row["diagonal_gauge_free_mixing"],row["stationary_limit"]["C0_conjugate_M_offdiagonal_skew"])
    print("critical C0:",value["samples"][2]["stationary_limit"]["C0_complex_re_im"])


if __name__=="__main__":
    main()
