#!/usr/bin/env python3
"""Named, no-fit width-eight cluster geometry dynamics; no lag/basis search."""
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
from scipy.optimize import brentq
import sympy as sp
from threadpoolctl import threadpool_limits

from p398_width8_projected_memory import decomposition,exponential_row,series
from p398_width8_memory_motifs import next_pair_motifs
from p398_width8_source_spectrum import complex_display,kreweras

ROOT=Path(__file__).resolve().parents[1]
PROTOCOL=ROOT/"analysis/p398_width8_geometric_compression_protocol.json"
OUT=ROOT/"results/p398-width8-geometric-compression/latest.json"


def features(states,f,t2):
    pair=next_pair_motifs(states)
    output={"A":f[:,0],"T2":t2,"T3":pair[:,0],"S11":pair[:,1],"B2":pair[:,2]}
    extra=[]
    weight=(1j)**np.arange(8)
    for state in states:
        blocks={label:[j for j in range(8) if state[j]==label] for label in set(state)}
        charge={label:sum(weight[j] for j in sites) for label,sites in blocks.items()}
        sizes={label:len(sites) for label,sites in blocks.items()}
        triplet=[sum(weight[j]*(sizes[state[j]]==3)*
            (sum(state[(j+d)%8]==state[j] for d in (-1,1))==r) for j in range(8)) for r in range(3)]
        t4=sum(charge[label] for label in blocks if sizes[label]==4)
        q3=b3=0j
        for j in range(8):
            a,b=state[j],state[(j+1)%8]
            if a!=b:
                q3+=(charge[a]+charge[b])*(sizes[a]+sizes[b]==3)
                b3+=charge[a]*(sizes[a]==3)+charge[b]*(sizes[b]==3)
        extra.append(triplet+[t4,q3,b3])
    for name,column in zip(("T3_r0","T3_r1","T3_r2","T4","Q3","B3"),np.array(extra).T):
        output[name]=column
    assert np.array_equal(output["T3"],output["T3_r0"]+output["T3_r1"]+output["T3_r2"])
    return output


def exact_columns(matrix,names):
    """A modular nonzero minor plus exact Gaussian-integer null relations."""
    real=np.rint(matrix.real).astype(np.int64)
    imag=np.rint(matrix.imag).astype(np.int64)
    mod=(real+256*imag)%65537
    basis,pivots,rows=[],[],[]
    for j,row in enumerate(mod):
        v=row.copy()
        for pivot,previous in zip(pivots,basis):
            v=(v-v[pivot]*previous)%65537
        nz=np.flatnonzero(v)
        if len(nz):
            pivot=int(nz[0]); v=v*pow(int(v[pivot]),-1,65537)%65537
            pivots.append(pivot); basis.append(v); rows.append(j)
            if len(rows)==len(names):
                break
    exact=sp.Matrix([[int(real[r,c])+sp.I*int(imag[r,c]) for c in range(len(names))] for r in rows])
    _,selected=exact.rref()
    relations=[]
    for vector in exact.nullspace():
        re=[sp.re(x) for x in vector]; im=[sp.im(x) for x in vector]
        denominator=sp.ilcm(*[sp.denom(x) for x in re+im])
        ar=np.array([int(x*denominator) for x in re],dtype=np.int64)
        ai=np.array([int(x*denominator) for x in im],dtype=np.int64)
        assert np.max(abs(real))*(sum(abs(ar))+sum(abs(ai)))<2**60
        assert np.all(real@ar-imag@ai==0) and np.all(real@ai+imag@ar==0)
        relations.append({name:str(vector[j]) for j,name in enumerate(names) if vector[j]})
    assert len(selected)==len(rows) and len(selected)+len(relations)==len(names)
    return {"rank_over_Q_i":len(rows),"selected_columns_in_declaration_order":[names[j] for j in selected],
            "selected_indices":list(selected),"null_relations_equal_zero":relations,
            "modular_prime":65537,"i_image":256,"independent_state_rows":rows}


def response(matrix,source):
    norm=float(np.vdot(source,source).real)
    values,vectors=linalg.eig(matrix)
    coefficients=(source.conj()@vectors)*linalg.solve(vectors,source)/norm
    visible=np.flatnonzero(abs(coefficients)>1e-10*max(abs(coefficients)))
    order=visible[np.argsort(values[visible].real)]
    selected=int(order[0])
    moments=[]
    v=source.copy()
    for _ in range(5):
        moments.append(float((source.conj()@v/norm).real)); v=matrix@v
    return {"mass":float(values[selected].real),"mass_imaginary":float(values[selected].imag),
            "lowest_normalized_residue":complex_display(coefficients[selected]),
            "moments_M0_to_M4":moments,"all_masses_re_im":complex_display(values)},values,coefficients


def build_result(protocol_path=PROTOCOL):
    protocol=json.loads(protocol_path.read_text())
    states,mass,f,t2,q,h,source,pi,weight,_,phase=decomposition()
    index={state:j for j,state in enumerate(states)}
    complement=[index[kreweras(state)] for state in states]
    named=features(states,f,t2)
    assert np.array_equal(-mass@named["T3"],3*(named["T4"]-named["T3"])+named["Q3"]-named["B3"])
    upper=linalg.cholesky(weight,lower=False)
    upper_inv=linalg.inv(upper)
    full=upper@h@upper_inv
    k=upper@(q.conj().T@q[complement,:]).toarray()@upper_inv
    inv=k/phase
    _,dual_basis=linalg.eigh((inv+inv.conj().T)/2)
    rays=[]
    for sign in (-1,1):
        ray_basis=dual_basis[:,:93] if sign<0 else dual_basis[:,93:]
        psi=upper@source@np.array([1,sign*phase])/np.sqrt(2)
        reference=ray_basis.conj().T@full@ray_basis
        raw,values,coeff=response(reference,ray_basis.conj().T@psi)
        ev,left,right=linalg.eig(reference,left=True,right=True)
        target=int(np.argmin(abs(ev-raw["mass"])))
        modes={"left":ray_basis@left[:,target],"right":ray_basis@right[:,target]}
        rays.append({"sign":sign,"basis":ray_basis,"source":psi,"reference":raw,
            "values":values,"coefficients":coeff,"modes":modes})
    result=[]
    final_bases={}
    reference_crossing=float(brentq(lambda t:series(rays[0]["values"],rays[0]["coefficients"],t)
                                   -series(rays[1]["values"],rays[1]["coefficients"],t),
                                   *protocol["crossing_bracket"]))
    for stage in protocol["stages"]:
        names=[]; columns=[]
        for name in stage["seeds"]:
            names.extend((name,"K_"+name))
            columns.extend((named[name],named[name][complement]))
        raw=np.column_stack(columns)
        exact=exact_columns(raw,names)
        selected=raw[:,exact["selected_indices"]]
        z=upper@np.asarray(q.conj().T@selected)
        named_mass=linalg.solve(z.conj().T@z,z.conj().T@full@z)
        stage_rows=[]; functions=[]
        for ray in rays:
            sign=ray["sign"]
            projected=(z+sign*k@z/phase)/2
            vectors,singular,_=linalg.svd(projected,full_matrices=False)
            dimension=exact["rank_over_Q_i"]//2
            basis=vectors[:,:dimension]
            assert singular[dimension-1]>1e-10
            compressed=basis.conj().T@full@basis
            a=basis.conj().T@ray["source"]
            row,values,coeff=response(compressed,a)
            samples=[]
            for t in protocol["distances"]:
                true=series(ray["values"],ray["coefficients"],t)
                prediction=series(values,coeff,t)
                samples.append({"t":t,"u_reference":true,"u_compressed":prediction,
                                "relative_error":prediction/true-1})
            capture={key:float(linalg.norm(basis.conj().T@mode)**2/linalg.norm(mode)**2)
                     for key,mode in ray["modes"].items()}
            row.update({"sign":sign,"dimension":dimension,"mass_relative_error":row["mass"]/ray["reference"]["mass"]-1,
                "source_moment_residuals":[a-b for a,b in zip(row["moments_M0_to_M4"],ray["reference"]["moments_M0_to_M4"])],
                "samples":samples,"lowest_full_mode_capture":capture})
            stage_rows.append(row); functions.append((values,coeff))
            final_bases[sign]=basis
        lo,hi=protocol["crossing_bracket"]
        difference=lambda t:series(*functions[0],t)-series(*functions[1],t)
        crossing=float(brentq(difference,lo,hi)) if difference(lo)*difference(hi)<0 else None
        result.append({"stage":stage["id"],"exact_named_columns":exact,
                       "named_Galerkin_mass_re_im":complex_display(named_mass),"rays":stage_rows,
                       "crossing_in_frozen_bracket":crossing,
                       "crossing_relative_error":crossing/reference_crossing-1 if crossing is not None else None})
        print(stage["id"],"rank",exact["rank_over_Q_i"],"crossing",crossing,
              "masses",[x["mass"] for x in stage_rows],
              "tail4 errors",[x["samples"][-1]["relative_error"] for x in stage_rows],flush=True)
    diagnostics=[]
    for ray in rays:
        basis=final_bases[ray["sign"]]
        slow=ray["modes"]["right"]
        missing_slow=slow-basis@(basis.conj().T@slow)
        left_slow=ray["modes"]["left"]
        missing_left=left_slow-basis@(basis.conj().T@left_slow)
        row={"sign":ray["sign"],"named_next_motifs":{}}
        for name in protocol["missing_direction_diagnostic_only"]:
            x=upper@np.asarray(q.conj().T@named[name])
            x=(x+ray["sign"]*k@x/phase)/2
            residual=x-basis@(basis.conj().T@x)
            row["named_next_motifs"][name]={
                "residual_fraction_of_variance":float(linalg.norm(residual)**2/linalg.norm(x)**2),
                "squared_alignment_with_unresolved_lowest_right_mode":float(abs(np.vdot(residual,missing_slow))**2/(linalg.norm(residual)**2*linalg.norm(missing_slow)**2)),
                "squared_alignment_with_unresolved_lowest_left_mode":float(abs(np.vdot(residual,missing_left))**2/(linalg.norm(residual)**2*linalg.norm(missing_left)**2))}
        diagnostics.append(row)
    return {"schema":"matching-one/p398-width8-geometric-compression/v1",
        "freeze_commit":protocol.get("freeze_commit","588a2fdd"),"protocol_status":protocol["status"],
        "protocol_sha256":hashlib.sha256(protocol_path.read_bytes()).hexdigest(),
        "reference":[{"sign":r["sign"],**r["reference"]} for r in rays],"stages":result,
        "reference_crossing":reference_crossing,
        "next_hierarchy_diagnostic":diagnostics,"visibility_relative_numerical_tolerance":1e-10,
        "scope":"single-thread stationary-L2 Galerkin on fixed named geometry spans; exact column relations, numerical dynamics; no new MC or adjustable fit"}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json",type=Path,default=OUT)
    parser.add_argument("--protocol",type=Path,default=PROTOCOL)
    args=parser.parse_args()
    with threadpool_limits(limits=1):
        value=build_result(args.protocol)
    paths=(Path(__file__),args.protocol,ROOT/"scripts/p398_width8_projected_memory.py",ROOT/"scripts/p398_width8_memory_motifs.py")
    value["input_sha256"]={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
    args.json.parent.mkdir(parents=True,exist_ok=True)
    args.json.write_text(json.dumps(value,indent=2)+"\n")


if __name__=="__main__":
    main()
