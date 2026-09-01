#!/usr/bin/env python3
"""Exact two-copy closure of the fixed local pair tensor, not a test suite."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import sympy as sp


def equality_patterns(n: int):
    def extend(labels):
        if len(labels) == n:
            yield tuple(labels)
            return
        for label in range(max(labels, default=-1) + 2):
            yield from extend(labels + [label])
    yield from extend([])


def pair_kernel(q, a, b, c, d):
    delta = lambda x, y: sp.Integer(x == y)
    return (1-delta(a,b))*(1-delta(c,d))* (
        delta(a,c)*delta(b,d)+delta(a,d)*delta(b,c)
        -(delta(a,c)+delta(a,d)+delta(b,c)+delta(b,d))/(q-2)
        +2/((q-1)*(q-2))
    )/2


def witness():
    """Construct the four fixed paths and report their actual NN components."""
    paths = [
        {(x,4) for x in range(4,11)},
        {(3,5),(11,5)} | {(x,6) for x in range(3,12)},
        {(3,3),(11,3)} | {(x,2) for x in range(3,12)},
        {(2,4),(12,4)} | {(1,y) for y in range(4,9)}
        | {(13,y) for y in range(4,9)} | {(x,8) for x in range(1,14)},
    ]
    occupied = set().union(*paths)
    labels = {}
    components = []
    for start in sorted(occupied):
        if start in labels:
            continue
        label = len(components)
        labels[start] = label
        stack = [start]
        vertices = []
        while stack:
            x,y = stack.pop()
            vertices.append((x,y))
            for dx,dy in [(1,0),(-1,0),(0,1),(0,-1)]:
                nxt = ((x+dx)%17,(y+dy)%17)
                if nxt in occupied and nxt not in labels:
                    labels[nxt] = label
                    stack.append(nxt)
        components.append(sorted(vertices))
    edges = sum(((x+dx)%17,(y+dy)%17) in occupied
                for x,y in occupied for dx,dy in [(1,0),(0,1)])
    ports = {}
    for name,(x,y) in {"x":(3,4),"y":(11,4)}.items():
        ports[name] = {p:labels[(x+dx,y+dy)]
                       for p,(dx,dy) in {"N":(0,1),"E":(1,0),"S":(0,-1),"W":(-1,0)}.items()}
    b_vac = 2*17**2-4*len(occupied)+edges
    return {
        "L":17, "N":17**2, "holes":[[3,4],[11,4]],
        "occupied_vertices": [list(v) for v in sorted(occupied)],
        "K":len(occupied), "occupied_NN_edges":edges,
        "component_count":len(components),
        "component_sizes":[len(c) for c in components],
        "port_component_labels":ports,
        "vacant_vacant_edges":b_vac,
        "hypergraph_component_count":b_vac+len(components),
        "spectator_colour_blocks":b_vac,
        "rank":0,
        "rank_reason":"All paths, and the graph after adding either or both holes, lie in the contractible rectangle [1,13] x [2,8] without seam edges.",
        "gluing":"xN-yN, xE-yW, xS-yS, xW-yE; reflection of the D4-invariant Kbar",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", default="analysis/local_pair_two_insertion_contract.json")
    parser.add_argument("--output-dir", default="results/local-pair-two-insertion")
    args = parser.parse_args()
    started = time.perf_counter()
    contract_path = Path(args.contract)
    contract = json.loads(contract_path.read_text())
    q = sp.Symbol("Q")
    rows = []
    cross = norm = sp.Integer(0)
    for a,b,c,d in equality_patterns(4):
        blocks = max(a,b,c,d)+1
        multiplicity = sp.prod(q-j for j in range(blocks))
        p = pair_kernel(q,a,b,c,d)
        r = pair_kernel(q,b,c,d,a)
        k = sp.factor((p+r)/2)
        cross += multiplicity*p*r
        norm += multiplicity*k*k
        rows.append({"equality_pattern":[a,b,c,d],"colour_count":str(multiplicity),"Kbar":str(k)})
    cross, norm = sp.factor(cross), sp.factor(norm)
    vx,vy = sp.symbols("v_x v_y", positive=True)
    baseline = q**4+(vx+vy+vx*vy)*q
    susceptibility = sp.factor(norm/baseline)
    result = {
        "schema":"matching-one.local-pair-two-insertion.v1",
        "kernel_source":contract["kernel_source"],
        "symbolic_results":{
            "cross_inner_product":str(cross),
            "G_two_copy_closure":str(norm),
            "Q1_G_residue":str(sp.residue(norm,q,1)),
            "Q1_G_finite_part":str(sp.limit(norm-sp.Rational(1,2)/(q-1),q,1)),
            "relative_fixed_vacant_weight":str(sp.factor(norm/q**4)),
            "conditional_baseline":str(baseline),
            "conditional_connected_mixed_susceptibility":str(susceptibility),
            "Q1_conditional_susceptibility_residue":str(sp.factor(sp.residue(susceptibility,q,1))),
            "sqrt_Qminus1_rescaled_two_copy_limit":str(sp.limit((q-1)*norm,q,1)),
            "sqrt_Qminus1_rescaled_single_response_limit":"0; the original finite nonzero V is not preserved",
        },
        "equality_patterns":rows,
        "physical_witness":witness(),
        "decision":"regular_unrenormalized_two_insertion_Q1_conditional_family_excluded",
        "preserved":"The completed single-insertion original-U response remains finite and nonzero.",
        "not_claimed":["a divergence of the fully summed homogeneous partition", "a continuum field identification", "an independent stochastic result"],
    }
    out = Path(args.output_dir)
    out.mkdir(parents=True,exist_ok=True)
    output = out/"latest.json"
    output.write_text(json.dumps(result,indent=2)+"\n")
    report = out/"REPORT.md"
    report.write_text("""# Single-insertion regularity does not extend to the fixed two-copy local tensor

The exact fixed contraction is

`G(Q)=Q(Q−3)(3Q²−9Q+8)/[8(Q−2)(Q−1)]`.

It has a nonremovable Q1 simple pole with residue **1/2**, although every
one-insertion Bell4 closure is finite. The 17x17 two-hole witness has four
disjoint occupied NN paths, K=52, four components, and rank0. Its physical
reflected gluing equals this Frobenius closure. No occupation ensemble was
enumerated.

After summing the two holes' vacant/occupied states, the conditional
partition, apart from common exterior factors, is

`Q^4+(v_x+v_y+v_x*v_y)Q+epsilon_x epsilon_y G(Q)`.

Both first-insertion terms vanish. The connected mixed log-partition
susceptibility has residue **1/[2(1+v_x)(1+v_y)]**. A common partition
normalizer therefore does not remove this pole. Separate single-site
quadratic counterterms cannot change this mixed derivative.

This excludes an unrenormalized finite-strength Q1 continuation in every
physical exterior for the specified tensor. It does **not** negate the
completed finite linear original-U response, or prove divergence after all
exterior configurations are summed. Rescaling each insertion by
sqrt(Q−1) makes this two-copy limit 1/2 but sends its finite one-insertion
response to zero; that is a different normalization and mechanism.

The 15 exact colour equality patterns, physical coordinates, rational
formulas and provenance are in latest.json; execution and hashes are in
run.json. This is exact algebra plus one graph construction, not Monte
Carlo, a new lattice enumeration, or a scientific test suite.
""")
    sha = lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    receipt = {
        "schema":"matching-one.local-pair-two-insertion.run.v1",
        "definition_commit":subprocess.check_output(["git","rev-parse","HEAD"],text=True).strip(),
        "command":[sys.executable,*sys.argv],
        "created_utc":datetime.now(timezone.utc).isoformat(),
        "python":sys.version,"machine":platform.machine(),"sympy":sp.__version__,
        "elapsed_seconds":time.perf_counter()-started,
        "contract_sha256":sha(contract_path),"script_sha256":sha(__file__),
        "output_sha256":{p.name:sha(p) for p in [output,report]},
        "new_occupation_enumerations":0,"new_random_samples":0,"cloud_jobs":0,"tests_run":0,
        "history":contract["selection_history"],
    }
    (out/"run.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({"results":result["symbolic_results"],"elapsed_seconds":receipt["elapsed_seconds"],"output_dir":str(out)},indent=2))


if __name__ == "__main__":
    main()
