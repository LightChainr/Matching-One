#!/usr/bin/env python3
"""Check PR628 with elementary counterexamples and an optional full L3 census.

The reviewed tree is read-only. Outputs do not replace its artifacts. The full
census compares a one-sign repair of its forest code with the existing weighted
union-find oracle, then separates the rank bug from the dual-transport bug.
"""
import argparse
import inspect
import json
import subprocess
import sys
from fractions import Fraction as F
from pathlib import Path

import numpy as np


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--probe',type=Path,required=True)
    ap.add_argument('--full-bond-census',action='store_true')
    args=ap.parse_args();probe=args.probe.resolve()
    sys.path.insert(0,str(probe/'scripts'))
    sys.path.insert(0,str(probe/'scripts/probe_invariant_shape'))
    import exact_rank_census as c
    import n725_zflow as z
    from square_bond_kappa3 import square_bond_pairs
    from threshold_quantile_lineage import cos_four_theta
    from torus_homology import HomologyUnionFind, _extend_basis

    contractible=[(0,1,1,0),(1,4,0,1),(3,4,1,0),(0,3,0,1)]
    namespace={}
    code=inspect.getsource(c.bond_ambient_rank)
    assert 'sx, sy = ax - dx, ay - dy' in code
    exec(code.replace('sx, sy = ax - dx, ay - dy','sx, sy = ax + dx, ay + dy'),namespace)
    broken=namespace['bond_ambient_rank']
    assert c.bond_ambient_rank(contractible)==0 and broken(contractible)==1
    comps=c.site_rank_pair_components(3);m,f=c.M_and_F_from_components(comps)
    alt={**comps,'P11':[F(3,2)*v for v in comps['P11']]}
    m1,f1=c.M_and_F_from_components(alt)
    ev=lambda a:c.eval_F(a,F(1,2))
    assert ev(f)==F(43,128) and sum(ev(a) for a in alt.values())==F(593,512)
    assert ev(f1)==F(425,1024) and ev(m1)==F(-21,64)
    co=[z.cos_four_theta(*r) for r in [(26,7),(23,14)]]
    correct=[cos_four_theta(*r) for r in [(26,7),(23,14)]]
    weights=lambda cc:[-cc[1]/(cc[0]-cc[1]),cc[0]/(cc[0]-cc[1])]
    report={'reviewed_commit':subprocess.check_output(['git','-C',str(probe),'rev-parse','HEAD'],text=True).strip(),
      'contractible_plaquette':{'edges':contractible,'reviewed_rank':1,'correct_rank':0},
      'float_cdf_at_half_on_site_L3':{'rational':str(ev(f)),'reviewed_float':c.eval_F_float(f,.5)},
      'reweighted_mass_at_half':{'total_mass':str(sum(ev(a) for a in alt.values())),
          'unnormalized_M':str(ev(m1)),'reported_F':str(ev(f1)),
          'contract_F_if_M_kept':str((1+ev(m1))/2),
          'normalized_M':str(ev(m1)/sum(ev(a) for a in alt.values())),
          'normalized_F':str(ev(f1)/sum(ev(a) for a in alt.values()))},
      'n725_orientation':{'reviewed_cos4':co,'actual_cos4':correct,
          'reviewed_weights':weights(co),'actual_weights':weights(correct),
          'actual_spin4_leakage_of_reviewed_weights':float(np.dot(weights(co),correct)),
          'reviewed_levels':list(z.LEVELS),
          'missing_deciles':[u/10 for u in range(1,10) if u/10 not in z.LEVELS]}}
    if args.full_bond_census:
        pairs=square_bond_pairs(3);nb=len(pairs);count=1<<nb
        bad=np.zeros(count,dtype=np.int8);good=np.zeros(count,dtype=np.int8)
        oracle_mismatches=0
        for mask in range(count):
            edges=[p.primal for i,p in enumerate(pairs) if mask>>i&1]
            bad[mask]=broken(edges);good[mask]=c.bond_ambient_rank(edges)
            uf=HomologyUnionFind(9,(3,3))
            for edge in edges:uf.add_edge(*edge)
            basis=[]
            for root in range(9):
                if uf.find(root)[0]==root:
                    for winding in uf.basis[root]:_extend_basis(basis,winding)
            oracle_mismatches+=int(good[mask]!=len(basis))
        primal_index={p.primal:i for i,p in enumerate(pairs)}
        dest=[primal_index[p.dual] for p in pairs]
        dual=np.zeros(count,dtype=np.int64)
        for mask in range(count):
            dual[mask]=sum(1<<dest[i] for i in range(nb) if not(mask>>i&1))
        complement=(count-1)^np.arange(count)
        failures={}
        for name,rank in [('bad_rank',bad),('fixed_rank',good)]:
            for mapping,indices in [('naive_complement',complement),('geometric_dual',dual)]:
                failures[name+'/'+mapping]=int(np.count_nonzero(rank+rank[indices]!=2))
        paircounts={f'{a},{b}':int(np.count_nonzero((good==a)&(good[dual]==b)))
                    for a in range(3) for b in range(3)}
        assert oracle_mismatches==0 and failures['fixed_rank/geometric_dual']==0
        # The broken rank routine is also edge-order dependent. Reproduce the
        # reviewed dual edge order rather than looking up an ascending mask.
        as_written_failures=0
        for mask in range(count):
            dedges=[pairs[dest[i]].primal for i in range(nb)
                    if not(mask>>dest[i]&1)]
            as_written_failures+=int(bad[mask]+broken(dedges)!=2)
        assert as_written_failures==118133
        # The correct float CDF uses raw count coefficients without another binomial factor.
        coefficients=[F(0)]*(nb+1)
        for mask,rank in enumerate(good):coefficients[mask.bit_count()]+=F(int(rank),2)
        assert c.eval_F(coefficients,F(1,2))==F(1,2)
        report['L3_bond_census']={'configurations':count,'identity_failures':failures,
              'reviewed_edge_order_identity_failures':as_written_failures,
              'fixed_forest_vs_union_find_mismatches':oracle_mismatches,
              'correct_rank_pair_counts':paircounts,'correct_F_at_half':'1/2',
              'correct_F_raw_count_coefficients':list(map(str,coefficients))}
    out=Path(__file__).resolve().parents[1]/'results/astra628-counterchecks/latest.json'
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
