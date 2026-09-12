#!/usr/bin/env python3
"""Width-three square-site torus: exact row-transfer spectrum and its scope.

The variable t used by coefficient recurrences is site fugacity. Coefficients
are the unnormalized Bernstein counts in M(p)=sum c[k] p**k (1-p)**(3*m-k).
No Monte Carlo, fitted exponent or infinite-lattice threshold estimate is used.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from itertools import product
import json
from math import comb
from pathlib import Path


def add(*arrays):
    result = [0] * max(map(len, arrays))
    for a in arrays:
        for i, v in enumerate(a):
            result[i] += v
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def mul(a, b):
    result = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i+j] += x*y
    return add(result)


def scale(a, c):
    return [c*x for x in a]


def power(a, n):
    result = [1]
    for _ in range(n):
        result = mul(result, a)
    return result


def evaluate(a, x):
    result = 0
    for c in reversed(a):
        result = result*x + c
    return result


def block_traces(m):
    """Integer polynomials for tr(B**m), tr(C**m), from Newton identities."""
    if not isinstance(m, int) or m < 0:
        raise ValueError('m must be a nonnegative integer')
    # B eigen-equation L^3 - h L^2 - j L + t^6 = 0.
    h, j, detneg = [0,1,3,1], [0,0,0,1,2], [0,0,0,0,0,0,1]
    b = [[3], h, add(mul(h,h),scale(j,2))]
    for k in range(3, m+1):
        b.append(add(mul(h,b[k-1]),mul(j,b[k-2]),scale(mul(detneg,b[k-3]),-1)))
    # C eigen-equation L^2 - t L - t^3 = 0; multiplicity two in the full transfer.
    c = [[2], [0,1]]
    for k in range(2,m+1):
        c.append(add(mul([0,1],c[k-1]),mul([0,0,0,1],c[k-2])))
    return b[m], c[m]


def bernstein_counts(m):
    if not isinstance(m,int) or m < 2:
        raise ValueError('the torus requires integer m >= 2')
    b,c = block_traces(m)
    out = add(b,scale(c,2),scale(power([1,3,3],m),-1))
    return out + [0]*(3*m+1-len(out))


def power_coefficients(m):
    n = 3*m
    out = [0]*(n+1)
    for k,c in enumerate(bernstein_counts(m)):
        for j in range(n-k+1):
            out[k+j] += c*comb(n-k,j)*(-1)**j
    return add(out)


def ambient_rank(mask, width, length):
    """Independent lift traversal; retains distinct +/- periodic edges.

    Positions are integer lifts of vertices. Every closing edge supplies an
    ambient period. Rank is computed by integer determinants, not row masks.
    """
    if min(width,length)<2:
        raise ValueError('honest axis periods must both be at least two')
    positions = {}
    first = None
    for root in range(width*length):
        if not (mask>>root)&1 or root in positions:
            continue
        positions[root] = (0,0)
        stack = [root]
        while stack:
            current = stack.pop()
            x,y = current%width,current//width
            px,py = positions[current]
            for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
                vertex = ((y+dy)%length)*width+(x+dx)%width
                if not (mask>>vertex)&1:
                    continue
                proposed = (px+dx,py+dy)
                if vertex not in positions:
                    positions[vertex] = proposed
                    stack.append(vertex)
                else:
                    wx,wy = proposed[0]-positions[vertex][0],proposed[1]-positions[vertex][1]
                    if wx%width or wy%length:
                        raise AssertionError('cycle displacement is not a torus period')
                    wx,wy=wx//width,wy//length
                    if wx or wy:
                        if first is None:
                            first=(wx,wy)
                        elif first[0]*wy-first[1]*wx:
                            return 2
    return int(first is not None)


def row_events(mask, width, length):
    rows=[(mask>>(width*j))&((1<<width)-1) for j in range(length)]
    full=any(row==(1<<width)-1 for row in rows)
    overlaps=all(rows[j]&rows[(j+1)%length] for j in range(length))
    return full,overlaps


def check_census(m):
    counts=[0]*(3*m+1)
    failures=0
    for mask in range(1<<(3*m)):
        rank=ambient_rank(mask,3,m)
        full,overlap=row_events(mask,3,m)
        failures+=rank != int(full)+int(overlap)
        counts[mask.bit_count()]+=rank-1
    if failures or counts != bernstein_counts(m):
        raise AssertionError(f'width 3, length {m}: row classification or coefficient failure')
    return {'configurations':1<<(3*m),'classification_failures':failures,
            'bernstein_counts':counts,'power_coefficients':power_coefficients(m)}


def overlap_trace_half(width,length):
    """Independent integer-matrix trace at p=1/2, not the block recurrence."""
    states=list(range(1,1<<width)); d=len(states)
    K=[[int(bool(s&t)) for t in states] for s in states]
    P=[[int(i==j) for j in range(d)] for i in range(d)]
    for _ in range(length):
        P=[[sum(P[i][k]*K[k][j] for k in range(d)) for j in range(d)] for i in range(d)]
    return F(sum(P[i][i] for i in range(d))-((1<<width)-1)**length, 1<<(width*length))


def bisect_sign(poly,lo=F(0),hi=F(1),steps=150):
    """An exact sign bracket, allowing either crossing orientation."""
    a,b=evaluate(poly,lo),evaluate(poly,hi)
    if not a*b < 0:
        raise ValueError('need opposite endpoint signs')
    for _ in range(steps):
        mid=(lo+hi)/2; value=evaluate(poly,mid)
        if value==0:
            return [str(mid),str(mid)]
        if a*value>0:
            lo,a=mid,value
        else:
            hi,b=mid,value
    return [str(lo),str(hi)]


def width4_witnesses():
    """Shows both missing within-row connectivity and unrecorded transverse loops."""
    answer=[]
    # All overlaps, no full row; the two sites of 0101 cannot connect.
    # Alternating pair of equal 1101/0111 rows makes a transverse winding loop
    # without any full row: it is a second, distinct width-four failure.
    for rows in ([1,5,4,5],[11,14]):
        width=4;m=len(rows)
        mask=sum(row<<(width*i) for i,row in enumerate(rows))
        full,overlap=row_events(mask,width,m)
        r=ambient_rank(mask,width,m)
        answer.append({'width':width,'length':m,'rows_as_bitmasks':rows,
                       'occupied_vertices':[i for i in range(width*m) if mask>>i&1],
                       'full_row':full,'every_interface_overlaps':overlap,
                       'actual_ambient_rank':r,'naive_rank':int(full)+int(overlap)})
    return answer



def interval_certificate():
    """Exact rational sign enclosures for q and all five distinct eigenvalues.

    Three disjoint sign brackets for a cubic and two for a quadratic exhaust
    their roots. Coarse fixed bands suffice for the all-m positive-defect proof.
    """
    sextic=[1,1,0,-4,-5,-3,1]
    # A wider fixed rational interval makes this certificate short and readable.
    lo,hi=F(58888069991785,10**14),F(58888069991786,10**14)
    assert evaluate(sextic,lo)>0>evaluate(sextic,hi)
    def point(x): return (F(x),F(x))
    def plus(a,b):return (a[0]+b[0],a[1]+b[1])
    def minus(a,b):return (a[0]-b[1],a[1]-b[0])
    def times(a,b):
        v=[x*y for x in a for y in b]
        return (min(v),max(v))
    def val(poly,x):
        result=point(0)
        for c in reversed(poly):result=plus(times(result,x),point(c))
        return result
    def sign(a):
        if a[0]>0:return 1
        if a[1]<0:return -1
        raise AssertionError('interval sign unresolved')
    q=(lo,hi);q2=times(q,q);q3=times(q,q2);qm=minus(point(1),q)
    U=times(q,times(qm,qm));V=times(q2,qm);W=q3
    H=plus(plus(U,times(point(3),V)),W)
    J=plus(times(U,V),times(point(2),times(U,W)))
    det=times(times(U,V),W)
    def charb(x):return plus(minus(minus(times(x,times(x,x)),times(H,times(x,x))),times(J,x)),det)
    def charc(x):return minus(minus(times(x,x),times(U,x)),times(U,V))
    bands={'B_negative':(charb,F(-102,1000),F(-99,1000)),
           'B_small_positive':(charb,F(35,1000),F(38,1000)),
           'B_perron':(charb,F(794,1000),F(797,1000)),
           'C_negative':(charc,F(-81,1000),F(-78,1000)),
           'C_positive':(charc,F(178,1000),F(180,1000))}
    roots={}
    for name,(f,a,b) in bands.items():
        signs=[sign(f(point(a))),sign(f(point(b)))]
        assert signs[0]*signs[1]<0
        roots[name]={'rational_interval':[str(a),str(b)],'endpoint_signs':signs}
    # M(q)>0: |beta_-|/mu+ < 3/5, |mu-|/mu+ < 1/2.
    assert F(102,178)<F(3,5) and F(81,178)<F(1,2)
    # h'(q)=-(q-1)^2 P'(q)/(F_lambda(q,lambda_c)*lambda_c)>0.
    lc=minus(point(1),q3)
    derivative=val([i*sextic[i] for i in range(1,len(sextic))],q)
    fl=minus(minus(times(point(3),times(lc,lc)),times(point(2),times(H,lc))),J)
    assert sign(derivative)==-1 and sign(fl)==1 and sign(lc)==1
    return {'q_rational_interval':[str(lo),str(hi)],
            'q_polynomial_endpoint_signs':[1,-1], 'eigenvalue_intervals':roots,
            'negative_beta_over_positive_mu_upper':'3/5',
            'negative_mu_over_positive_mu_upper':'1/2',
            'P_prime_at_q_sign':-1,'F_lambda_at_perron_sign':1,
            'hprime_positive':True,'all_m_at_least_2_defect_positive':True,
            'arithmetic':'Fraction interval arithmetic; every endpoint sign certified'}


def continuation_memory_witness():
    # Both open prefixes have frontier 0101 and no winding. Only the second
    # prefix connects its two frontier vertices. A common suffix exposes this.
    prefixes=[[0,5],[7,5]];suffix=[13,0];ranks=[]
    for prefix in prefixes:
        rows=prefix+suffix
        mask=sum(row<<(4*j) for j,row in enumerate(rows))
        ranks.append(ambient_rank(mask,4,len(rows)))
    assert ranks==[0,1]
    return {'width':4,'prefixes':prefixes,'common_frontier_mask':5,
            'common_suffix':suffix,'completed_torus_ranks':ranks,
            'completed_matching_D':[-1,0],
            'conclusion':'Any deterministic continuation state identifying these prefixes cannot preserve the rank readout.'}



def open_prefix_signature(rows, width=4):
    """Exact open-cylinder connectivity with horizontal lifted displacements."""
    height=len(rows)
    occupied={width*y+x for y,row in enumerate(rows) for x in range(width) if row>>x&1}
    positions={};components={};winding=False;component=0
    for root in sorted(occupied):
        if root in positions:continue
        positions[root]=(0,0);components[root]=component;stack=[root]
        while stack:
            vertex=stack.pop();x,y=vertex%width,vertex//width
            for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
                if not 0<=y+dy<height:continue
                other=width*(y+dy)+(x+dx)%width
                if other not in occupied:continue
                proposed=(positions[vertex][0]+dx,positions[vertex][1]+dy)
                if other not in positions:
                    positions[other]=proposed;components[other]=component;stack.append(other)
                elif proposed!=positions[other]:
                    assert (proposed[0]-positions[other][0])%width==0
                    assert proposed[1]==positions[other][1]
                    winding=True
        component+=1
    last=width*(height-1)
    partition={}
    for x in range(width):
        if last+x in occupied:partition.setdefault(components[last+x],[]).append(x)
    groups=sorted(partition.values())
    marked={}
    for group in groups:
        for x in group[1:]:
            marked[f'{group[0]}->{x}']=positions[last+x][0]-positions[last+group[0]][0]
    return {'occupied_count':len(occupied),'frontier_mask':rows[-1],
            'ordinary_frontier_partition':groups,'existing_transverse_winding':winding,
            'lifted_frontier_displacements':marked}


def topology_memory_witness():
    prefixes=[[13,5],[7,5]];suffix=[13,0]
    signatures=[open_prefix_signature(rows) for rows in prefixes]
    ordinary_keys=('occupied_count','frontier_mask','ordinary_frontier_partition',
                   'existing_transverse_winding')
    assert all(signatures[0][k]==signatures[1][k] for k in ordinary_keys)
    assert signatures[0]['lifted_frontier_displacements']=={'0->2':-2}
    assert signatures[1]['lifted_frontier_displacements']=={'0->2':2}
    ranks=[]
    for prefix in prefixes:
        rows=prefix+suffix;mask=sum(row<<(4*i) for i,row in enumerate(rows))
        ranks.append(ambient_rank(mask,4,len(rows)))
    assert ranks==[0,1]
    return {'width':4,'prefixes':prefixes,'common_suffix':suffix,
            'prefix_signatures':signatures,'completed_torus_ranks':ranks,
            'completed_matching_D':[-1,0],
            'exact_common_prefix_probability':'p^5*(1-p)^3',
            'conclusion':'Even frontier mask, ordinary set partition, occupation count, and current winding flag do not preserve continuation; lifted/annular path information is needed.'}


def report(max_check=5):
    import mpmath as mp
    if not 2<=max_check<=6:
        raise ValueError('local exact check restricted to 2 <= max_check <= 6')
    checks={str(m):check_census(m) for m in range(2,max_check+1)}
    for m in range(2,11):
        assert overlap_trace_half(3,m)==evaluate(power_coefficients(m),F(1,2))
    # Independent 3x3 table from #668/#684 (the existing square-torus rung).
    assert bernstein_counts(3)==[-1,-9,-36,-78,-90,-36,36,36,9,1]
    witness=width4_witnesses()
    assert witness[0]['actual_ambient_rank']==0 and witness[0]['naive_rank']==1
    sextic=[1,1,0,-4,-5,-3,1]
    with mp.workdps(80):
        q=mp.findroot(lambda p:evaluate(sextic,p),(.58,.60))
        def blocks(p):
            u=p*(1-p)**2;v=p*p*(1-p);w=p**3
            return mp.matrix([[u,2*v,w],[2*u,3*v,w],[3*u,3*v,w]]),u,v
        def leading(p):
            b,_,_=blocks(p)
            vals=mp.eig(b,left=False,right=False)
            return max(mp.re(v) for v in vals)
        b,u,v=blocks(q)
        betas=sorted([mp.re(x) for x in mp.eig(b,left=False,right=False)],reverse=True)
        mu=(u+mp.sqrt(u*u+4*u*v))/2
        nu=(u-mp.sqrt(u*u+4*u*v))/2
        lc=1-q**3
        hprime=mp.diff(lambda p:mp.log(leading(p)/(1-p**3)),q)
        r=mu/lc
        roots={}
        for m in (2,3,4,5,6,8,12,20):
            bracket=bisect_sign(power_coefficients(m))
            lf,hf=map(F,bracket)
            root=(mp.mpf(lf.numerator)/lf.denominator+mp.mpf(hf.numerator)/hf.denominator)/2
            shift=-2*r**m/(m*hprime)
            roots[str(m)]={'exact_rational_root_bracket':bracket,
                          'root_diagnostic':mp.nstr(root,40),
                          'root_minus_cylinder_diagnostic':mp.nstr(root-q,25),
                          'leading_shift_diagnostic':mp.nstr(shift,25),
                          'shift_ratio_diagnostic':mp.nstr((root-q)/shift,25)}
        s=lambda x:mp.nstr(x,50)
        return {'schema':'matching-one.width3-cylinder-exact.v1',
            'scope':'square-site width 3, length m>=2; NOT infinite-lattice pc or an all-width pTL intertwiner',
            'cylinder_sextic_low_first':sextic,
            'cylinder_exact_rational_bracket':bisect_sign(sextic),
            'cylinder_q_diagnostic':s(q),
            'published_n3_value':'0.5888806999178529980514426957517049337221',
            'source':'Jacobsen 2015 arXiv:1507.03027v1 section 6.1 Table 2, primary HTML read 2026-09-12',
            'published_match_below_1e_minus_40':abs(q-mp.mpf('0.5888806999178529980514426957517049337221'))<mp.mpf('1e-40'),
            'block_B_eigenvalues_diagnostic':[s(x) for x in betas],
            'block_C_eigenvalues_diagnostic':[s(mu),s(nu)],
            'block_C_eigenvalue_multiplicity':2,
            'closed_rate_diagnostic':s(lc),
            'rho_diagnostic':s(r),'hprime_diagnostic':s(hprime),
            'all_leading_coefficients':1,'dominant_subleading_coefficient':2,
            'exact_enumeration_checks':checks,'integer_matrix_trace_checks_m_2_to_10':True,
            'finite_roots':roots,'width4_obstructions':witness,
            'exact_interval_certificate':interval_certificate(),
            'width4_continuation_memory_witness':continuation_memory_witness(),
            'width4_topology_memory_witness':topology_memory_witness(),
            'boundary':'Rows remember occupancy only; width 4 requires within-row/accumulated connectivity. Width-three proof does not justify mask-only transfer at width 4.'}


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path)
    ap.add_argument('--max-check',type=int,default=5)
    args=ap.parse_args()
    text=json.dumps(report(args.max_check),indent=2,allow_nan=False)+'\n'
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        with args.out.open('x',encoding='utf-8') as f:f.write(text)
    else:
        print(text,end='')
