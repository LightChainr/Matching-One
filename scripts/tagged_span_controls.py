#!/usr/bin/env python3
"""Independent physical controls and deterministic report for tagged span law."""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
from time import perf_counter
import tagged_winding_span as T


def components(width, rows, matching=False):
    occ={(x,y) for y,mask in enumerate(rows) for x in range(width) if mask>>x&1}
    steps=[(1,0),(-1,0),(0,1),(0,-1)]
    if matching:steps += [(a,b) for a in (-1,1) for b in (-1,1)]
    unseen=set(occ); out=[]
    while unseen:
        root=min(unseen); potentials={root:0}; stack=[root]; unseen.remove(root); winding=False
        while stack:
            x,y=stack.pop()
            for dx,dy in steps:
                v=((x+dx)%width,y+dy)
                if v not in occ:continue
                q=potentials[(x,y)]+dx
                if v in potentials:
                    if potentials[v]!=q:
                        assert (potentials[v]-q)%width==0
                        winding=True
                else:
                    potentials[v]=q;unseen.discard(v);stack.append(v)
        out.append((frozenset(potentials),winding))
    return out


def external_boundary(width, C, matching=False):
    steps=[(1,0),(-1,0),(0,1),(0,-1)]
    if matching:steps += [(a,b) for a in (-1,1) for b in (-1,1)]
    return {((x+dx)%width,y+dy) for x,y in C for dx,dy in steps}-set(C)


def shape_activity_counts(width,height,matching=False):
    counts=Counter(); checked=0
    for mask in range(1,1<<(width*height)):
        checked+=1
        rows=[(mask>>(width*y))&((1<<width)-1) for y in range(height)]
        if not rows[0] or not rows[-1]:continue
        cs=components(width,rows,matching)
        if len(cs)!=1 or not cs[0][1]:continue
        C=cs[0][0]
        counts[(len(C),len(external_boundary(width,C,matching)))]+=1
    return counts,checked


def activity_value(counts,p):
    return sum(c*p**k*(1-p)**b for (k,b),c in counts.items())


def enumerate_tag_paths(width,height,matching=False,p=F(1,2)):
    """Check all source/row words with retirement at given height against BFS."""
    count=0; accepted_mass=F(0); physical_mass=F(0)
    sources={}
    for prev,mask,anchor,s in T.source_paths(width,matching):
        sources.setdefault((prev,mask),[]).append((anchor,s))
    for rows in product(range(1<<width),repeat=height+2):
        count+=1;occ=sum(m.bit_count() for m in rows)
        wt=p**occ*(1-p)**(width*(height+2)-occ)
        actual=[]
        for C,wind in components(width,rows,matching):
            if wind and min(y for x,y in C)==1 and max(y for x,y in C)==height:
                actual.append(min(x for x,y in C if y==1))
        predicted=[]
        for anchor,s in sources.get((rows[0],rows[1]),[]):
            terminal=0
            for time,mask in enumerate(rows[2:],1):
                s,terminal=T.step(s,mask,matching)
                if terminal:
                    if terminal==1 and time==height:predicted.append(anchor)
                    break
        assert sorted(actual)==sorted(predicted),(width,height,matching,rows,actual,predicted)
        accepted_mass+=wt*len(predicted);physical_mass+=wt*len(actual)
    assert accepted_mass==physical_mass
    return count,physical_mass


def verify_sample(width,matching,sample):
    rows=sample['rows'];anchor=sample['anchor']
    selected=[(C,wind) for C,wind in components(width,rows,matching) if (anchor,1) in C]
    assert len(selected)==1
    C,wind=selected[0]
    assert wind and min(y for x,y in C)==1 and max(y for x,y in C)==len(rows)-2
    assert max(y for x,y in C)-min(y for x,y in C)+1==sample['span']
    assert anchor==min(x for x,y in C if y==1)
    return {'sites':sorted([list(v) for v in C]),'occupation':len(C),'span':sample['span']}


def symbolic_width_two():
    import sympy as sp
    p,z=sp.symbols('p z');out=[]
    for g in (False,True):
        states,tr,src=T.build(2,g);tr,src,_=T.lump(tr,src);n=len(tr)
        R=sp.zeros(n);a=sp.zeros(1,n);b=sp.zeros(n,1)
        for i,row in enumerate(tr):
            for m,(j,_) in enumerate(row):
                wt=p**m.bit_count()*(1-p)**(2-m.bit_count())
                if j>=0:R[i,j]+=wt
                elif j==-1:b[i]+=wt
        for (i,k),c in src.items():a[i]+=c*p**k*(1-p)**(4-k)
        pgf=sp.factor((z*a*(sp.eye(n)-z*R).inv()*b)[0])
        # First coefficient is the exact isolated full-row activity.
        assert sp.simplify(sp.diff(pgf,z).subs(z,0)-p**2*(1-p)**4)==0
        out.append({'matching':g,'density_pgf':str(pgf),'p_half':str(sp.factor(pgf.subs(p,sp.Rational(1,2))))})
    return out


def laplace_diagnostic(tr,src,p,mean,scales=(1,2,4)):
    import numpy as np
    from scipy.linalg import solve
    alpha,R,b,delta=T.numeric_system(tr,src,p,False)
    R=np.array(R);alpha=np.array(alpha);b=np.array(b);I=np.eye(len(b))
    nu=float(alpha@solve(I-R,b));out={}
    for s in scales:
        z=float(np.exp(-s/mean))
        out[str(s)]=float(z*alpha@solve(I-z*R,b)/nu)
    return out


def brownian_laplace(scales=(1,2,4)):
    import mpmath as mp
    mp.mp.dps=50
    mean=mp.sqrt(mp.pi/2)
    def cdf(x):
        if x<=0:return mp.mpf(0)
        if x<mp.mpf('1.0'):
            return mp.sqrt(2*mp.pi)*mp.pi**2/x**3*sum(n*n*mp.exp(-mp.pi**2*n*n/(2*x*x)) for n in range(1,40))
        return 1+2*sum((1-4*n*n*x*x)*mp.exp(-2*n*n*x*x) for n in range(1,30))
    out={}
    for s in scales:
        rate=mp.mpf(s)/mean
        val=rate*mp.quad(lambda x:mp.exp(-rate*x)*cdf(x),[0,mp.mpf('.5'),1,2,4,8,mp.inf])
        out[str(s)]=mp.nstr(val,26)
    return out


def activity_coefficients(tr,src,height):
    cur=Counter({(i,0,b):n for (i,b),n in src.items()});outputs=[]
    for _ in range(height):
        nxt=Counter();coeff=Counter()
        for (i,kk,bb),n in cur.items():
            for j,k,b in tr[i]:
                if j<0:coeff[(kk+k,bb+b)]+=n
                else:nxt[(j,kk+k,bb+b)]+=n
        outputs.append(coeff);cur=nxt
    return outputs


def report(max_width=8):
    if not 4<=max_width<=8: raise ValueError('report supports max width 4..8')
    out={'scope':'new one-tag construction; no age cutoff; finite state closures checked only at recorded widths',
         'source_commit_reviewed':'7226a2c6d099535f34486eccb1bd996f7affda13',
         'systems':[],'physical_activity_controls':[],'finite_word_controls':[], 'samples':[], 'direct_activity_joint':[]}
    cache={}
    for g in (False,True):
        ps=[F(1,8),F(1,4)] if not g else [F(1,16),F(1,8)]
        for w in range(2,max_width+1):
            t=perf_counter();states,tr,src=T.build(w,g);tr,src,blocks=T.lump(tr,src)
            cache[(w,g)]=(tr,src)
            sysrec={'width':w,'matching':g,'tagged_states':len(states),'exit_lumps':len(tr),
                    'mask_transitions':len(states)*(1<<w),'runs':[]}
            for p in ps+([F(1,2)] if w<=4 else []):
                if w<=4:
                    values=T.moments(tr,src,p,bins=8,exact=True)
                    rec={'p':str(p),'mode':'exact Fraction, all heights','moments':T.jsonable(values)}
                    mean=float(values['mean']);nu=values['nu']
                else:
                    values=T.certified_moments(tr,src,p)
                    rec={'p':str(p),'mode':'rational centres + exact residual enclosures, all heights',
                         'moments':T.jsonable(values)}
                    mean=float(values['centres']['mean'])
                rec['laplace_mean_scaled_float_diagnostic']=laplace_diagnostic(tr,src,p,mean)
                sysrec['runs'].append(rec)
            sysrec['elapsed_seconds']=perf_counter()-t;out['systems'].append(sysrec)
            print('system',w,g,'states',len(states),'blocks',len(tr),'seconds',round(sysrec['elapsed_seconds'],3),flush=True)
    total=0
    for g in (False,True):
        for w,H in [(2,1),(2,2),(2,3),(2,4),(3,1),(3,2),(3,3),(3,4),(4,1),(4,2),(4,3)]:
            counts,checked=shape_activity_counts(w,H,g);total+=checked
            tr,src=cache[(w,g)]
            for p in (F(1,4),F(1,2),F(3,4)):
                val=activity_value(counts,p);vals=T.moments(tr,src,p,bins=H,exact=True)
                assert vals['d_h'][H-1]==val
                out['physical_activity_controls'].append({'width':w,'height':H,'matching':g,'p':str(p),'density':str(val)})
    out['nonempty_shape_masks_checked']=total
    for g in (False,True):
        for w in range(2,min(5,max_width)+1):
            st,tt,src=T.activity_transfer(w,g)
            p=F(1,8) if g else F(1,4)
            mm=T.activity_joint_moments(tt,src,p)
            tagged_tr,tagged_src=cache[(w,g)]
            reference=T.moments(tagged_tr,tagged_src,p,bins=3,exact=True)
            assert mm['nu']==reference['nu'] and mm['mean_span']==reference['mean'] and mm['variance_span']==reference['variance']
            out['direct_activity_joint'].append({'width':w,'matching':g,'states':len(st),'blocks':len(tt),'p':str(p),'moments':T.jsonable(mm)})
            if w<=4:
                cf=activity_coefficients(tt,src,3)
                for H,co in enumerate(cf,1):
                    physical,_=shape_activity_counts(w,H,g)
                    assert co==physical,(w,g,H,co,physical)
    out['joint_activity_coefficient_checks']=18
    total=0
    for g in (False,True):
        for w,H in ((2,1),(2,2),(2,3),(3,1),(3,2)):
            count,value=enumerate_tag_paths(w,H,g);total+=count
            out['finite_word_controls'].append({'width':w,'height':H,'matching':g,'words':count,'density':str(value)})
    out['full_source_word_controls']=total
    # Known complete-component controls from the prior physical activity calculation.
    for g in (False,True):
        tr,src=cache[(4,g)];vals=T.moments(tr,src,F(1,2),bins=3,exact=True)
        assert vals['nu']==F(323849,5576960)
        if not g:assert sum(vals['d_h'])==F(9087,1048576)
    for g,p in ((False,F(1,4)),(True,F(1,8))):
        for seed in range(4):
            sample=T.sample_component(4,g,p,seed)
            actual=verify_sample(4,g,sample)
            out['samples'].append({'matching':g,'p':str(p),'seed':seed,'sample':T.jsonable(sample),'physical_check':actual})
    out['width_two_symbolic_pgf']=symbolic_width_two()
    out['brownian_range_laplace_mean_scaled']=brownian_laplace()
    out['note']='Exact certificates bound arithmetic on the newly constructed finite operator. No fixed-p Brownian claim or global repository CI is inferred.'
    return out

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--max-width',type=int,default=8)
    ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    r=report(args.max_width);args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(T.jsonable(r),indent=2)+'\n')
