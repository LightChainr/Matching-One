#!/usr/bin/env python3
"""P398: weak-source amplification and directed contact/pair memory motifs."""
from __future__ import annotations

import os
for _key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[_key]="1"

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import linalg,sparse
from threadpoolctl import threadpool_limits

from p321_homology_trace_certificate import join_adjacent
from p333_generic_q_detach_intertwiner import detach_state
from p398_width8_projected_memory import decomposition
from p398_width8_source_spectrum import complex_display,kreweras

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/p398-width8-memory-motifs/latest.json"


def directed_events(states):
    index={state:j for j,state in enumerate(states)}
    events=[]
    dual_classes=defaultdict(int)
    grouped={kind:([],[],[]) for kind in ("join","detach")}
    for x,state in enumerate(states):
        dual=kreweras(state)
        for site in range(8):
            for kind,action in (("join",join_adjacent),("detach",detach_state)):
                target=action(state,site)
                if target==state:
                    continue
                y=index[target]
                rows,cols,values=grouped[kind]
                rows.extend((x,x)); cols.extend((x,y)); values.extend((1,-1))
                size=state.count(state[site])
                if kind=="join":
                    other=(site+1)%8
                    sizes=sorted((size,state.count(state[other])))
                    blocks={state[site],state[other]}
                    contacts=sum({state[j],state[(j+1)%8]}==blocks for j in range(8))
                    dual_size=dual.count(dual[site])
                    assert kreweras(target)==detach_state(dual,site)
                    assert (contacts==2)==(dual_size==2)
                    dual_classes[(contacts,dual_size)]+=1
                    key=f"join_size_{sizes[0]}_{sizes[1]}_contacts_{contacts}"
                    coarse=f"join_contacts_{'1' if contacts==1 else 'multiple'}"
                else:
                    contacts=sum(state[(site+offset)%8]==state[site] for offset in (-1,1))
                    key=f"detach_size_{size}_incident_{contacts}"
                    coarse=f"detach_size_{size if size<4 else '4plus'}"
                events.append((x,y,kind,key,coarse))
    parts={kind:sparse.coo_matrix((values,(rows,cols)),shape=(len(states),len(states))).tocsr()
           for kind,(rows,cols,values) in grouped.items()}
    return events,parts,dual_classes


def next_pair_motifs(states):
    weight=(1j)**np.arange(8)
    output=[]
    for state in states:
        size=[state.count(state[j]) for j in range(8)]
        t3=sum(weight[j]*(size[j]==3) for j in range(8))
        s11=sum((weight[j]+weight[(j+1)%8])*(size[j]==size[(j+1)%8]==1) for j in range(8))
        boundary_two=0j
        for label in set(state):
            members=[j for j in range(8) if state[j]==label]
            if len(members)==2:
                charge=sum(weight[j] for j in members)
                boundary=sum((state[j]==label)!=(state[(j+1)%8]==label) for j in range(8))
                boundary_two+=charge*boundary
        output.append((t3,s11,boundary_two))
    return np.array(output)


def build_result():
    states,mass,f,t2,q,h,source,pi,weight,_,phase=decomposition()
    upper=linalg.cholesky(weight,lower=False)
    upper_inv=linalg.inv(upper)
    whitened=upper@h@upper_inv
    sources=upper@source
    covariance=sources.conj().T@sources
    invisible=linalg.null_space(sources.conj().T)
    hidden=invisible.conj().T@whitened@invisible
    right=invisible.conj().T@whitened@sources
    left=sources.conj().T@whitened@invisible
    moments=[right,linalg.solve(hidden,right)]
    moments.append(linalg.solve(hidden,moments[1]))
    source_blocks=[left@value for value in moments]
    events,parts,dual_classes=directed_events(states)
    assert (parts["join"]+parts["detach"]-mass).nnz==0
    compressed={key:upper@(q.conj().T@part@q).toarray()@upper_inv for key,part in parts.items()}
    motifs=next_pair_motifs(states)
    t3,s11,b2=motifs.T
    assert np.array_equal(-mass@t2,2*t3-2*t2+s11-b2)
    # The forward hidden source is purely contact-join or size-two-detachment.
    assert np.max(abs(invisible.conj().T@compressed["detach"]@sources[:,0]))<1e-10
    assert np.max(abs(invisible.conj().T@compressed["join"]@sources[:,1]))<1e-10
    physical_map=q@upper_inv@invisible
    names=("k0","integrated_memory","first_time_moment")
    rows=[]
    for sign in (-1,1):
        vector=np.array([1,sign*phase])/np.sqrt(2)
        psi=f@vector
        variance=float(np.real(vector.conj()@covariance@vector))
        motif_terms=[]
        primitive_terms=[]
        for block,degree in zip(source_blocks,range(3)):
            term=np.conj(vector)[:,None]*block*vector[None,:]/variance
            motif_terms.append({"moment":names[degree],"RR_RT_TR_TT_re_im":complex_display(term),
                "normalized_sum":float(term.sum().real),"raw_sum":float((term.sum()*variance).real)})
            primitive=np.zeros((2,2),dtype=complex)
            for a,ka in enumerate(("join","detach")):
                back=(sources@vector).conj()@compressed[ka]@invisible
                for b,kb in enumerate(("join","detach")):
                    force=invisible.conj().T@compressed[kb]@sources@vector
                    for _ in range(degree):
                        force=linalg.solve(hidden,force)
                    primitive[a,b]=back@force/variance
            assert abs(primitive.sum()-term.sum())<1e-9
            primitive_terms.append({"moment":names[degree],"JJ_JD_DJ_DD_re_im":complex_display(primitive)})
        edge_budgets={}
        for degree,moment in enumerate(moments):
            response=physical_map@moment@vector
            fine=defaultdict(complex)
            coarse=defaultdict(complex)
            for x,y,kind,key,group in events:
                value=pi[x]*np.conj(psi[x])*(response[x]-response[y])/variance
                fine[key]+=value
                coarse[group]+=value
            assert abs(sum(fine.values())-motif_terms[degree]["normalized_sum"])<1e-9
            edge_budgets[names[degree]]={
                "coarse":[{"motif":key,"signed_normalized_re_im":complex_display(val),
                           "signed_raw_re_im":complex_display(val*variance)} for key,val in sorted(coarse.items())],
                "fine":[{"motif":key,"signed_normalized_re_im":complex_display(val),
                         "signed_raw_re_im":complex_display(val*variance)}
                        for key,val in sorted(fine.items(),key=lambda kv:-abs(kv[1]))]}
        rows.append({"sign":sign,"source_variance":variance,"force_motif_blocks":motif_terms,
                     "primitive_generator_blocks":primitive_terms,"directed_edge_budgets":edge_budgets})
    return {"schema":"matching-one/p398-width8-memory-motifs/v1",
        "parent":"39e06607ec3a353b1130acebf770da591acaf340",
        "projection":"global stationary-L2 P onto both A,L; Q removes both; no independent R/T vote inside a protected ray",
        "force_motif_convention":"right=(QMA,QML)=-(QR,QT2); left=(QMdagger A,QMdagger L) are stationary time-reversed dual motifs, not the same right force",
        "source_moment_blocks_re_im":{key:complex_display(value) for key,value in zip(names,source_blocks)},
        "source_covariance_re_im":complex_display(covariance),
        "cancellation_coherence":{"source":float(abs(covariance[0,1])/covariance[0,0].real),
            **{key:float(abs(block[0,1])/block[0,0].real) for key,block in zip(names,source_blocks)}},
        "exact_dual_join_detach_class_counts":[{"join_contacts":key[0],"dual_detach_block_size":key[1],"directed_site_events":count}
            for key,count in sorted(dual_classes.items())],
        "next_pair_identity":{"formula":"G T2=2 T3-2 T2+S11-B2",
            "T3":"sum_j i^j 1{block(j) has size 3}",
            "S11":"sum_adjacent singleton-singleton edges (i^j+i^(j+1))",
            "B2":"sum_size2 blocks (sum_{j in block}i^j)*(number of cut boundary edges)",
            "checked_exactly_on_all_1430_states":True},
        "rays":rows,
        "normalization":{"minus_over_plus_source_variance":rows[0]["source_variance"]/rows[1]["source_variance"],
            "plus_over_minus_raw_k0":rows[1]["force_motif_blocks"][0]["raw_sum"]/rows[0]["force_motif_blocks"][0]["raw_sum"],
            "plus_over_minus_normalized_k0":rows[1]["force_motif_blocks"][0]["normalized_sum"]/rows[0]["force_motif_blocks"][0]["normalized_sum"]},
        "directed_budget_formula":"sum_{x->y in class} pi_x conjugate(psi_x)*(v_x-v_y)/Var(psi), where v=Q eta, Q D^-1 eta, or Q D^-2 eta in physical coordinates",
        "scope":"single-thread deterministic float64 motif decomposition of the existing finite generator; signed response contributions, not independent causal effects or probabilities"}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json",type=Path,default=OUT)
    args=parser.parse_args()
    with threadpool_limits(limits=1):
        value=build_result()
    value["input_sha256"]={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (Path(__file__),ROOT/"scripts/p398_width8_projected_memory.py",ROOT/"scripts/p398_width8_source_spectrum.py")}
    args.json.parent.mkdir(parents=True,exist_ok=True)
    args.json.write_text(json.dumps(value,indent=2)+"\n")
    print(value["normalization"])
    for row in value["rays"]:
        print("sign",row["sign"],"source variance",row["source_variance"])
        print("RT blocks",row["force_motif_blocks"])
        print("primitive blocks",row["primitive_generator_blocks"])
        print("coarse k0",row["directed_edge_budgets"]["k0"]["coarse"])


if __name__=="__main__":
    main()
