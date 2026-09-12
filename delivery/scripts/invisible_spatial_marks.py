#!/usr/bin/env python3
"""Strictly positive spatial marks invisible to all uniform rank/K data.

Finite countermodels, not alternative claims about the actual Bernoulli ensemble.
The marks are D4/translation invariant and independent of p.
"""
from __future__ import annotations
import argparse,json
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import time
from mpmath import mp
from verify_source_hessian import read_cert,all_ranks,pattern,group_counts,bivariate

def orbit(mask):
    out=set()
    points=[(i%4,i//4) for i in range(16) if (mask>>i)&1]
    for reflect in (False,True):
        for rotation in range(4):
            rotated=[]
            for x,y in points:
                if reflect:x=-x
                for _ in range(rotation):x,y=-y,x
                rotated.append((x,y))
            for dx in range(4):
                for dy in range(4):out.add(sum(1<<(((y+dy)%4)*4+(x+dx)%4) for x,y in rotated))
    return sorted(out)

def weighted_field_sum(mark,field,p,e):
    plus=sum(1<<i for i,h in enumerate(field) if h==1)
    ans=0
    for mask,f in mark.items():
        a=(mask&plus).bit_count();b=mask.bit_count()-a
        ans+=f*(p+e)**a*(1-p-e)**(8-a)*(p-e)**b*(1-p+e)**(8-b)
    return ans

def pevalbern(c,p):return sum(a*p**i*(1-p)**(len(c)-1-i) for i,a in enumerate(c))

def report(certpath,hessianpath,dps=70):
    start=time.perf_counter();ranks=all_ranks(read_cert(certpath));h=json.loads(hessianpath.read_text())
    square=3|(3<<4)
    Lshape=3|(1<<4)|(1<<8)
    Tshape=7|(2<<4)
    oa=orbit(square);rows=[]
    with mp.workdps(dps):
        root=mp.findroot(lambda p:pevalbern(h['M_coefficients'],p),('.58','.60'))
        mp1=pevalbern(h['Mprime_coefficients'],root)
        for label,mask in [('L',Lshape),('T',Tshape)]:
            ob=orbit(mask);g=gcd(len(oa),len(ob));fa,fb=len(ob)//g,-len(oa)//g
            assert not set(oa)&set(ob)
            mark={x:fa for x in oa};mark.update({x:fb for x in ob})
            strength=F(1,2*max(abs(fa),abs(fb)))
            assert sum(mark.values())==0
            assert all(ranks[x]==0 and x.bit_count()==4 for x in mark)
            assert min(1+strength*f for f in mark.values())>=F(1,2)
            assert min(1-strength*f for f in mark.values())>=F(1,2)
            modes=[]
            for kx,ky in ((1,0),(2,0),(1,1),(2,1),(2,2)):
                field=pattern(kx,ky)
                contrast=sum(f*sum(field[i] for i in range(16) if (mask>>i)&1)**2 for mask,f in mark.items())
                eps=mp.mpf(strength.numerator)/strength.denominator
                correction=eps*root**2*(1-root)**10*contrast/mp1
                modes.append({'wavevector':[kx,ky],'integer_spatial_variance_contrast':contrast,
                              'plus_mark_additive_root_curvature_change':mp.nstr(correction,40)})
            finite=[]
            for kx,ky in ((2,0),(2,2)):
                field=pattern(kx,ky);base=group_counts(ranks,field)
                contrast=next(v['integer_spatial_variance_contrast'] for v in modes if v['wavevector']==[kx,ky])
                expected=2*eps*root**2*(1-root)**10*contrast/mp1
                for den in (64,128):
                    amp=mp.mpf(1)/den
                    # X=-1 on support of mark; normalizer stays positive, so zeros use numerator.
                    roots=[]
                    for sign in (1,-1):
                        roots.append(mp.findroot(lambda p:bivariate(base,p,amp)-sign*eps*weighted_field_sum(mark,field,p,amp),('.58','.60')))
                    estimate=2*(roots[0]-roots[1])/amp**2
                    finite.append({'wavevector':[kx,ky],'field_amplitude':f'1/{den}',
                                   'plus_root':mp.nstr(roots[0],40),'minus_root':mp.nstr(roots[1],40),
                                   'difference_of_root_curvatures_estimate':mp.nstr(estimate,30),
                                   'predicted_difference':mp.nstr(expected,30)})
            rows.append({'shape':label,'reference_square_mask':square,'marked_shape_mask':mask,
                         'reference_orbit':oa,'comparison_orbit':ob,'reference_orbit_size':len(oa),'comparison_orbit_size':len(ob),
                         'reference_mark_value':fa,'comparison_mark_value':fb,'strength':str(strength),
                         'all_support_in_rank0_K4':True,'all_rank_K_coefficient_changes_zero':True,
                         'weight_range_plus':[str(min([F(1)]+[1+strength*f for f in mark.values()])),str(max([F(1)]+[1+strength*f for f in mark.values()]))],
                         'weight_range_minus':[str(min([F(1)]+[1-strength*f for f in mark.values()])),str(max([F(1)]+[1-strength*f for f in mark.values()]))],
                         'source_mode_contrasts':modes,'finite_amplitude_checks':finite})
    matrix=[[next(v['integer_spatial_variance_contrast'] for v in row['source_mode_contrasts'] if v['wavevector']==list(k)) for row in rows] for k in ((2,0),(2,2))]
    det=matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0]
    assert det!=0
    return {'schema':'matching-one.invisible-spatial-marks.v1','scope':'4x4 finite positive marked product measures, not the unmodified site ensemble',
            'marks':rows,'two_source_two_mark_integer_matrix':matrix,'determinant':det,
            'conclusion':'Identical entire uniform-p joint laws of (r,K), but independent nonzero spatial-source distinctions.',
            'elapsed_seconds':time.perf_counter()-start}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--certificate',type=Path,required=True)
    p.add_argument('--hessian',type=Path,required=True);p.add_argument('--out',type=Path);p.add_argument('--dps',type=int,default=70)
    a=p.parse_args();text=json.dumps(report(a.certificate,a.hessian,a.dps),indent=2)+'\n'
    if a.out:
        with a.out.open('x') as f:f.write(text)
    else:print(text,end='')
if __name__=='__main__':main()
