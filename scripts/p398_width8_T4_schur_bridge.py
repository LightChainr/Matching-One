#!/usr/bin/env python3
"""One named Schur bridge: the old seven-space plus its T4 residual."""
from __future__ import annotations

import os
for _key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import linalg
from threadpoolctl import threadpool_limits

from p398_width8_geometric_compression import features
from p398_width8_projected_memory import decomposition
from p398_width8_source_spectrum import complex_display, kreweras

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT/"analysis/p398_width8_T4_schur_bridge.json"
GEOMETRY = ROOT/"results/p398-width8-t4-post-reveal/latest.json"
OUT = ROOT/"results/p398-width8-T4-schur-bridge/latest.json"


def scalar(value):
    return complex_display(complex(value))


def bridge_score(h):
    a, b, c, d = h[:7,:7], h[:7,7], h[7,:7], h[7,7]
    e = np.eye(7,dtype=complex)[:,0]
    moments = []
    old, new = e.copy(), np.eye(8,dtype=complex)[:,0]
    for order in range(7):
        moments.append({"order":order,"old7_re_im":scalar(old[0]),"new8_re_im":scalar(new[0]),
                        "increment_re_im":scalar(new[0]-old[0])})
        old=a@old; new=h@new
    eb, ce = b[0], c@e
    eab, cae = (a@b)[0], c@a@e
    eaab, caae = (a@a@b)[0], c@a@a@e
    cubic = {"A_b_c":eab*ce,"b_c_A":eb*cae,"b_d_c":eb*d*ce}
    quartic = {"A2_b_c":eaab*ce,"A_b_c_A":eab*cae,"b_c_A2":eb*caae,
               "A_b_d_c":eab*d*ce,"b_d_c_A":eb*d*cae,"b_d2_c":eb*d*d*ce,
               "b_c_b_c":eb*(c@b)*ce}
    first = next((row["order"] for row in moments
                  if abs(complex(*row["increment_re_im"])) >
                  1e-10*max(1,abs(complex(*row["old7_re_im"])),abs(complex(*row["new8_re_im"])))),None)
    return {"generator_8_re_im":complex_display(h),"left_b_re_im":complex_display(b),
            "right_c_re_im":complex_display(c),"pole_d_re_im":scalar(d),
            "old7_dimension":7,"source_coordinate":0,
            "left_source_w_re_im":scalar(eb),"right_w_source_re_im":scalar(ce),
            "right_w_A_source_re_im":scalar(cae),"right_w_A2_source_re_im":scalar(caae),
            "direct_source_bridge_product_re_im":scalar(eb*ce),
            "moments":moments,"first_numerically_nonzero_increment":first,
            "cubic_path_terms_re_im":{name:scalar(value) for name,value in cubic.items()},
            "quartic_path_terms_re_im":{name:scalar(value) for name,value in quartic.items()},
            "cubic_sum_residual":float(abs(sum(cubic.values())-complex(*moments[3]["increment_re_im"]))),
            "quartic_sum_residual":float(abs(sum(quartic.values())-complex(*moments[4]["increment_re_im"])))}


def build_result():
    protocol=json.loads(PROTOCOL.read_text())
    archived=json.loads(GEOMETRY.read_text())["stages"][0]
    names=archived["exact_named_columns"]["selected_columns_in_declaration_order"]
    states,_,f,t2,q,h,source,_,weight,_,phase=decomposition()
    index={state:j for j,state in enumerate(states)}
    complement=[index[kreweras(state)] for state in states]
    named=features(states,f,t2)
    upper=linalg.cholesky(weight,lower=False)
    inverse=linalg.inv(upper)
    g=-(upper@h@inverse)
    k=upper@(q.conj().T@q[complement,:]).toarray()@inverse
    columns=[named[name[2:]][complement] if name.startswith("K_") else named[name] for name in names]
    z0=upper@np.asarray(q.conj().T@np.column_stack(columns))
    rows=[]
    for sign,name in ((-1,"minus"),(1,"plus")):
        psi=upper@source@np.array([1,sign*phase])/np.sqrt(2)
        e=psi/linalg.norm(psi)
        z=(z0+sign*k@z0/phase)/2
        z_perp=z-e[:,None]*(e.conj()@z)[None,:]
        hidden=linalg.svd(z_perp,full_matrices=False)[0][:,:6]
        basis=np.column_stack((e,hidden))
        t4=upper@np.asarray(q.conj().T@named["T4"])
        t4=(t4+sign*k@t4/phase)/2
        residual=t4-basis@(basis.conj().T@t4)
        w=residual/linalg.norm(residual)
        extended=np.column_stack((basis,w))
        block=extended.conj().T@g@extended
        coordinates=linalg.lstsq(z,basis)[0]
        rows.append({"ray":name,"sign":sign,"selected_named_columns":names,
                     "basis7_in_projected_named_columns_re_im":complex_display(coordinates),
                     "T4_residual_norm":float(linalg.norm(residual)),
                     "original":bridge_score(block),"reversible":bridge_score((block+block.conj().T)/2)})
    return {"schema":protocol["schema"],"protocol":str(PROTOCOL.relative_to(ROOT)),"ray_rows":rows,
            "chronology":"moment order inferred from existing exact hierarchy, then evaluated in only the fixed seven/eight blocks",
            "boundary":"Schur self-energy of a fixed geometric projection; no new physical mode, Markov state chain or fitted pole"}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json",type=Path,default=OUT)
    args=parser.parse_args()
    with threadpool_limits(limits=1):
        result=build_result()
    inputs=(PROTOCOL,Path(__file__),GEOMETRY,ROOT/"scripts/p398_width8_geometric_compression.py",ROOT/"scripts/p398_width8_projected_memory.py")
    result["input_sha256"]={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    args.json.parent.mkdir(parents=True,exist_ok=True)
    args.json.write_text(json.dumps(result,indent=2)+"\n")
    for row in result["ray_rows"]:
        print(row["ray"])
        for process in ("original","reversible"):
            score=row[process]
            print(process,"first",score["first_numerically_nonzero_increment"],"pole",score["pole_d_re_im"])
            print("left,right,A,A2",*[score[key] for key in ("left_source_w_re_im","right_w_source_re_im","right_w_A_source_re_im","right_w_A2_source_re_im")])
            print("moment increments",[r["increment_re_im"] for r in score["moments"]])
            print("quartic",score["quartic_path_terms_re_im"])


if __name__=="__main__":
    main()
