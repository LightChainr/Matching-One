#!/usr/bin/env python3
"""Width-eight AP/landing projected memory: deterministic, no fitted kernel."""
from __future__ import annotations

import os
for _key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[_key]="1"

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import linalg
from scipy.integrate import quad
from scipy.optimize import brentq
from threadpoolctl import threadpool_limits

from p398_width8_source_spectrum import complex_display, kreweras, model

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/p398-width8-projected-memory/latest.json"


def decomposition():
    # Reuse the established generator/sector. No repeat of the 93+93 rank test.
    states,mass,f,t2,q,h,source,_,_,pi=model(8)
    index={state:j for j,state in enumerate(states)}
    complement=[index[kreweras(state)] for state in states]
    k=(q.conj().T@q[complement,:]).toarray()
    phase=np.exp(-1j*np.pi/4)
    involution=k/phase
    _,v=linalg.eigh((involution+involution.conj().T)/2)
    weight=(q.conj().T@q.multiply(pi[:,None])).toarray()
    return states,mass,f,t2,q,h,source,pi,weight,v,phase


def exponential_row(matrix,row,column):
    values,vectors=linalg.eig(matrix)
    coefficients=(row@vectors)*linalg.solve(vectors,column)
    return values,coefficients


def series(values,coefficients,t,derivative=0):
    return float(np.real(np.sum(coefficients*(-values)**derivative*np.exp(-values*t))))


def ray_block(sign,data):
    states,mass,f,t2,q,h,source,pi,weight,dual_basis,phase=data
    v=dual_basis[:,:93] if sign<0 else dual_basis[:,93:]
    ray_mass=v.conj().T@h@v
    ray_weight=v.conj().T@weight@v
    psi=(f[:,0]+sign*phase*f[:,1])/np.sqrt(2)
    raw_source=v.conj().T@source@np.array([1,sign*phase])/np.sqrt(2)
    variance=float(np.real(raw_source.conj()@ray_weight@raw_source))
    # W=upper^* upper gives Euclidean coordinates for the physical L2(pi) metric.
    upper=linalg.cholesky(ray_weight,lower=False)
    whitened=upper@ray_mass@linalg.inv(upper)
    visible=upper@raw_source/np.sqrt(variance)
    invisible=linalg.null_space(visible.conj()[None,:])
    unitary=np.column_stack((visible,invisible))
    block=unitary.conj().T@whitened@unitary
    omega=float(block[0,0].real)
    back=block[0,1:]
    force=block[1:,0]
    hidden=block[1:,1:]
    poles,residues=exponential_row(hidden,back,force)
    masses,weights=exponential_row(block,np.eye(93)[0],np.eye(93)[:,0])
    u=lambda t:series(masses,weights,t)
    hazard=lambda t:-series(masses,weights,t,1)/u(t)
    kernel=lambda t:series(poles,residues,t)
    kernel_hat=lambda z:complex(back@linalg.solve(z*np.eye(92)+hidden,force))
    k0=kernel(0)
    assert abs(k0-((block@block)[0,0].real-omega*omega))<1e-10
    variance_force=float(np.vdot(force,force).real)
    left_force_variance=float(np.vdot(back,back).real)
    residual=(mass@psi)-omega*psi
    assert np.max(abs(np.vdot(psi,pi*residual)))<1e-11
    physical_force_variance=float(np.vdot(residual,pi*residual).real/variance)
    assert abs(physical_force_variance-variance_force)<1e-10
    # First explicit hidden observable eta=(M-omega)psi, same physical units.
    hidden_slope=float((np.vdot(force,hidden@force)/np.vdot(force,force)).real)
    hidden_back=float((back@force/np.linalg.norm(force)).real)
    first_hidden_matrix=np.array([[omega,hidden_back],[np.linalg.norm(force),hidden_slope]])
    # Descriptive two-observable Galerkin truncation, not an asserted exact closure.
    vals2,weights2=exponential_row(first_hidden_matrix,np.eye(2)[0],np.eye(2)[:,0])
    forward_unexplained=hidden@force-hidden_slope*force
    return {"sign":sign,"variance":variance,"omega":omega,
        "k0":k0,"k_prime_0":series(poles,residues,0,1),
        "force_variance":variance_force,
        "left_force_variance":left_force_variance,
        "left_right_force_alignment":k0/np.sqrt(variance_force*left_force_variance),
        "stationary_force_formula":"eta=(3-omega)*psi-(R+sign*exp(-i*pi/4)*T2)/sqrt(2)",
        "feedback_vs_force_variance_ratio":k0/variance_force,
        "kernel_integral":float(kernel_hat(0).real),
        "kernel_signed_mean_time":float((back@linalg.solve(hidden,linalg.solve(hidden,force))).real/kernel_hat(0).real),
        "hidden_poles_re_im":complex_display(poles),
        "hidden_residues_re_im":complex_display(residues),
        "minimum_hidden_decay_real":float(poles.real.min()),
        "maximum_hidden_decay_imaginary":float(abs(poles.imag).max()),
        "nonreal_hidden_poles_re_im":complex_display(poles[np.abs(poles.imag)>1e-8]),
        "first_hidden_Galerkin_mass":first_hidden_matrix.tolist(),
        "first_hidden_Galerkin_masses":vals2.real.tolist(),
        "first_hidden_unexplained_forward_norm":float(linalg.norm(forward_unexplained)/linalg.norm(force)),
        "samples":[{"s":t,"u":u(t),"instantaneous_decay":hazard(t),"kernel":kernel(t),
                    "force_autocorrelation":float(np.real(force.conj()@linalg.expm(-hidden*t)@force)),
                    "instant_Markov_u":float(np.exp(-omega*t)),
                    "one_explicit_hidden_u":series(vals2,weights2,t)}
                   for t in (0.,.025,.05,.1,.2,.26565731998357506,.5,1.,2.,4.)]},u,hazard,kernel,kernel_hat


def sign_roots(function,maximum=8):
    grid=np.linspace(0,maximum,401)
    roots=[]
    previous=function(grid[0])
    for lo,hi in zip(grid[:-1],grid[1:]):
        current=function(hi)
        if previous*current<0:
            roots.append(float(brentq(function,lo,hi)))
        previous=current
    return roots


def build_result():
    data=decomposition()
    minus=ray_block(-1,data)
    plus=ray_block(1,data)
    rays=[]
    for row,u,hazard,kernel,khat in (minus,plus):
        row["kernel_sign_changes_0_to_8"]=sign_roots(kernel)
        row["minimum_kernel_on_0_to_8_grid_401"]=min(kernel(t) for t in np.linspace(0,8,401))
        row["zero_frequency_source_integral"]=1/(row["omega"]-row["kernel_integral"])
        # No time stepping is used: check the exact resolvent and a single Volterra evaluation.
        t=.5
        convolution=quad(lambda s:kernel(t-s)*u(s),0,t,epsabs=1e-11)[0]
        row["volterra_at_half"]={"feedback_convolution":convolution,
            "derivative_from_identity":-row["omega"]*u(t)+convolution,
            "derivative_from_spectrum":-hazard(t)*u(t)}
        rays.append(row)
    mrow,um,hm,km,_=minus
    prow,up,hp,kp,_=plus
    crossing=brentq(lambda t:um(t)-up(t),.1,.5)
    rate_crossings=sign_roots(lambda t:hp(t)-hm(t),2)
    dw=prow["omega"]-mrow["omega"]
    dk=prow["k0"]-mrow["k0"]
    approximate=[]
    for row in rays:
        matrix=np.array(row["first_hidden_Galerkin_mass"])
        values,weights=exponential_row(matrix,np.eye(2)[0],np.eye(2)[:,0])
        approximate.append((values,weights))
    approximate_crossing=brentq(lambda t:series(*approximate[0],t)-series(*approximate[1],t),.1,.5)
    return {"schema":"matching-one/p398-width8-projected-memory/v1",
        "parent":"552c45d7595ebcb0d04555cec03b2a5bfd8da44a",
        "projection":"orthogonal L2(pi), separately in the two exact Kreweras rays",
        "mass_block":"M=[[omega,b],[c,D]]; u_prime=-omega*u+integral_0^t b*exp(-D*(t-s))*c*u(s) ds",
        "resolvent":"u_hat(z)=1/(z+omega-b*(zI+D)^(-1)*c)",
        "ray_rows":rays,
        "crossing":{"actual_normalized_ray_crossing":float(crossing),
            "instantaneous_decay_crossings":rate_crossings,
            "delta_initial_decay_plus_minus":dw,"delta_initial_curvature_plus_minus":dk,
            "quadratic_log_ratio_crossing_prediction":2*dw/dk,
            "linear_decay_crossing_prediction":dw/dk,
            "one_geometric_hidden_per_ray_crossing_prediction":float(approximate_crossing),
            "integrated_decay_difference_at_actual_crossing":quad(lambda t:hp(t)-hm(t),0,crossing,epsabs=1e-12)[0]},
        "arithmetic":"single-thread float64 deterministic projection of inherited exact integer generator; no fit, no MC",
        "scope":"finite-state elimination memory, not path-dependent morphism semantics or a Jordan identification"}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json",type=Path,default=OUT)
    args=parser.parse_args()
    with threadpool_limits(limits=1):
        value=build_result()
    value["input_sha256"]={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (Path(__file__),ROOT/"scripts/p398_width8_source_spectrum.py",
                  ROOT/"results/p398-width8-source-spectrum/latest.json")}
    args.json.parent.mkdir(parents=True,exist_ok=True)
    args.json.write_text(json.dumps(value,indent=2)+"\n")
    for row in value["ray_rows"]:
        print({k:row[k] for k in ("sign","omega","k0","force_variance","kernel_integral","kernel_signed_mean_time","minimum_hidden_decay_real","kernel_sign_changes_0_to_8","first_hidden_Galerkin_mass")})
    print(value["crossing"])


if __name__=="__main__":
    main()
