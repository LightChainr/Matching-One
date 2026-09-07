#!/usr/bin/env python3
"""Small exact controls for #610: duality is not a sign pattern across scales.

L=2 square bond checks the new readout identities with integer/rational sums.
L=3 square site rejects same-model p -> 1-p symmetry. A three-state transfer
matrix supplies both signs of a quadratic even-observable response.
"""
import json
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from qtangent_scale_decomposition import enumerate_torus, spatial_filtration_levels
from square_bond_duality_exact import geometric_dual_mask
from square_bond_kappa3 import square_bond_pairs
from matched_torus_reference import axis_geometry, cluster_stats


def mean(x):
    return sum(x, F(0))/len(x)


def cov(x,y):
    return mean([a*b for a,b in zip(x,y)])-mean(x)*mean(y)


def main():
    table=enumerate_torus(2); rows=table['rows']; pairs=square_bond_pairs(2)
    dual=[geometric_dual_mask(r['mask'],pairs) for r in rows]
    x=[F(r['x_source']) for r in rows]
    b=[F(r['twice_b_even'],2) for r in rows]
    o=[F(r['wrap_either']) for r in rows]
    cross=[F(r['wrap_cross']) for r in rows]
    assert all(x[d]==-x[i] and b[d]==b[i] for i,d in enumerate(dual))
    assert all(o[i]==1+(x[i]-x[i]**2)/2 for i in range(256))
    assert cov(o,x)==mean(cross)==cov(x,x)/2
    assert cov(o,b)==-cov([t*t for t in x],b)/2
    degree1=[mean([b[m]*(1-2*((m>>i)&1)) for m in range(256)]) for i in range(8)]
    assert not any(degree1)
    # Exact conditional increments for the existing nested spatial filtration.
    previous_o=[mean(o)]*256; previous_x=[mean(x)]*256
    increments=[]
    for subset in spatial_filtration_levels(2):
        selector=sum(1<<i for i in subset)
        def conditional(y):
            buckets={}
            for mask,value in enumerate(y):buckets.setdefault(mask&selector,[]).append(value)
            avg={k:mean(v) for k,v in buckets.items()}
            return [avg[m&selector] for m in range(256)]
        co,cx=conditional(o),conditional(x)
        gamma=mean([(a-b)*(c-d) for a,b,c,d in zip(co,previous_o,cx,previous_x)])
        assert gamma>=0
        increments.append(str(gamma));previous_o,previous_x=co,cx
    assert sum(map(F,increments))==cov(o,x)
    # A genuine square-site torus, not a self-dual bond surrogate.
    geom=axis_geometry(3); total=0
    for mask in range(1<<geom.n):
        black=[bool(mask>>i&1) for i in range(geom.n)]
        _,bp=cluster_stats(black,geom.primal_edges)
        _,wm=cluster_stats([not t for t in black],geom.matching_edges)
        total+=int(bp)-int(wm)
    mhalf=F(total,1<<geom.n)
    assert mhalf!=0
    # Row-stochastic transfer T(eps)=I+(G0+eps H)/4, R exchanges states 1 and 2.
    g=np.array([[-2,1,1],[1,-1,0],[1,0,-1]],dtype=object)
    h=np.array([[0,0,0],[1,-1,0],[-1,0,1]],dtype=object)
    r=np.array([[1,0,0],[0,0,1],[0,1,0]],dtype=object)
    assert np.array_equal(r@g@r,g) and np.array_equal(r@h@r,-h)
    t=np.eye(3,dtype=object)+g*F(1,4); ht=h*F(1,4)
    # Exact polynomial multiplication, four transfer steps, no finite differences.
    coefficients=[np.eye(3,dtype=object)]
    for _ in range(4):
        new=[np.zeros((3,3),dtype=object) for j in range(len(coefficients)+1)]
        for j,c in enumerate(coefficients):new[j]+=c@t;new[j+1]+=c@ht
        coefficients=new
    scalar=[F(c[0,0]) for c in coefficients]
    assert scalar[1]==scalar[3]==0 and scalar[2]!=0
    result={'square_bond_L2':{'configurations':256,'B_even_is_duality_even':True,
               'X_is_duality_odd':True,'cov_wrap_X':str(cov(o,x)),
               'P_rank2':str(mean(cross)),'half_var_X':str(cov(x,x)/2),
               'cov_wrap_B_even':str(cov(o,b)),
               'minus_half_cov_X_squared_B_even':str(-cov([t*t for t in x],b)/2),
               'B_even_degree_one_coefficients':list(map(str,degree1)),
               'spatial_wrap_X_increments':increments},
            'square_site_L3':{'configurations':1<<geom.n,'M_at_half':str(mhalf),
                'F_at_half':str((1+mhalf)/2),'same_model_complement_symmetry':False},
            'three_state_transfer':{'center_probability_coefficients_eps_0_to_4':list(map(str,scalar)),
                 'quadratic_center':str(scalar[2]),'quadratic_complement':str(-scalar[2]),
                 'point':'Both readouts are invariant; quadratic sign is not fixed. A nonzero quadratic coefficient is not forced in general.'}}
    out=Path(__file__).resolve().parents[1]/'results/astra610-independent/symmetry-check.json'
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
