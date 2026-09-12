#!/usr/bin/env python3
"""The intrinsic rank source exposes modes canceled by Matching One.

Reads the already certified width-four automaton and the preceding parametric
trace definition. Exact integer arithmetic only. No new state exploration.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/research-control-20260912'


def pvalue(coefficients,t):
    out=0
    for c in reversed(coefficients):out=out*t+c
    return out


def polynomial_multiply(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out


def trace_moments(descending_coefficients,maximum):
    d=len(descending_coefficients)-1
    if descending_coefficients[0]!=1:raise ValueError('monic polynomial required')
    result=[d]
    for m in range(1,maximum+1):
        if m<=d:
            result.append(-sum(descending_coefficients[i]*result[m-i]
                               for i in range(1,m))-m*descending_coefficients[m])
        else:
            result.append(-sum(descending_coefficients[i]*result[m-i]
                               for i in range(1,d+1)))
    return result


def det_bareiss(matrix):
    a=[list(row) for row in matrix];n=len(a);last=1;sign=1
    for k in range(n-1):
        if not a[k][k]:
            j=next((j for j in range(k+1,n) if a[j][k]),None)
            if j is None:return 0
            a[k],a[j]=a[j],a[k];sign=-sign
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=pivot*a[i][j]-a[i][k]*a[k][j]
                q,r=divmod(numerator,last)
                if r:raise AssertionError('Bareiss division was not exact')
                a[i][j]=q
            a[i][k]=0
        last=pivot
    return sign*a[-1][-1]


def det_modular(matrix,prime):
    """Independent finite-field Gaussian elimination of the original minor."""
    a=[[x%prime for x in row] for row in matrix];n=len(a);det=1
    for k in range(n):
        j=next((j for j in range(k,n) if a[j][k]),None)
        if j is None:return 0
        if j!=k:a[j],a[k]=a[k],a[j];det=-det
        pivot=a[k][k];det=det*pivot%prime;inverse=pow(pivot,-1,prime)
        for i in range(k+1,n):
            multiplier=a[i][k]*inverse%prime
            for j in range(k+1,n):a[i][j]=(a[i][j]-multiplier*a[k][j])%prime
    return det%prime


def sector_weights(definition):
    output={}
    for label in ('P0','P2'):
        weights=Counter(dict(definition['blocks']['5' if label=='P0' else '15']['factorization']))
        for name,n in definition['blocks']['16']['factorization']:weights[name]-=n
        weights['positive_square']+=2
        output[label]=weights
    # Do not use Counter subtraction: it would silently discard negative weights.
    combinations={'M':(-1,1),'E':(1,1),'twice_tilted_mean_numerator':(-1,4),
                  'twice_Z_at_exp_s_2':(-1,2)}
    for label,(a,b) in combinations.items():
        output[label]=Counter({name:a*output['P0'][name]+b*output['P2'][name]
                               for name in sorted(set(output['P0'])|set(output['P2']))})
    output['twice_Z_at_exp_s_2']['total_row_weight']+=2
    return {label:{name:c for name,c in values.items() if c}
            for label,values in output.items()}


def exact_sequences(source,t,maximum):
    rows=source['quotient_transitions'];outputs=source['quotient_rank_output']
    initial=source['quotient_initial'];n=len(rows)
    v=[0]*n
    for mask,state in enumerate(initial):v[state]+=t**mask.bit_count()
    transitions=[list(Counter((dest,mask.bit_count()) for mask,dest in enumerate(row)).items())
                 for row in rows]
    result={name:[] for name in ('P0','P2','M','E','twice_tilted_mean_numerator','twice_Z_at_exp_s_2')}
    for m in range(1,maximum+1):
        p0=sum(v[i] for i,r in enumerate(outputs) if r==0)
        p2=sum(v[i] for i,r in enumerate(outputs) if r==2)
        total=sum(v)
        if total!=(1+t)**(4*m):raise AssertionError('occupation weight not conserved')
        values={'P0':p0,'P2':p2,'M':p2-p0,'E':p0+p2,
                'twice_tilted_mean_numerator':4*p2-p0,
                'twice_Z_at_exp_s_2':2*total-p0+2*p2}
        for name,value in values.items():result[name].append(value)
        new=[0]*n
        for i,row in enumerate(transitions):
            for (dest,k),count in row:new[dest]+=v[i]*count*t**k
        v=new
    return result


def verify(source_path,definition_path):
    started=time.perf_counter()
    raw=Path(source_path).read_bytes();source=json.loads(raw)
    defraw=Path(definition_path).read_bytes();definition=json.loads(defraw)
    actual_blob=hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
    if actual_blob!=definition['source_blob']:raise ValueError('definition refers to a different automaton')
    factors={name:v['coefficients_desc_x_ascending_t'] for name,v in definition['factors'].items()}
    factors['positive_square']=[[1],[0,0,-1]]
    factors['total_row_weight']=[[1],[-1,-4,-6,-4,-1]]
    weights=sector_weights(definition)
    # t=3 avoids the already known specializations t=1 and t=2.
    t=3;maximum=60;sequences=exact_sequences(source,t,maximum)
    traces={name:trace_moments([pvalue(row,t) for row in co],maximum)
            for name,co in factors.items()}
    for label,terms in weights.items():
        assembled=[sum(c*traces[name][m] for name,c in terms.items())
                   for m in range(1,maximum+1)]
        if assembled!=sequences[label]:raise AssertionError(f'{label}: raw state DP and trace formula disagree')
    channels=[]
    for label in ('M','E','twice_tilted_mean_numerator','twice_Z_at_exp_s_2'):
        order=sum(len(factors[name])-1 for name in weights[label])
        sequence=sequences[label]
        # Physical tail starts at m=2, not the convenient one-row initializer.
        minor=[[sequence[1+i+j] for j in range(order)] for i in range(order)]
        determinant=det_bareiss(minor)
        if determinant==0:raise AssertionError('candidate generic upper bound is not witnessed')
        modchecks={}
        for prime in (1000003,1000033):
            residue=det_modular(minor,prime)
            if residue!=determinant%prime:raise AssertionError('independent determinant check failed')
            modchecks[str(prime)]=residue
        channels.append({'channel':label,'generic_minimum_scalar_order':order,
                         'fugacity_witness':t,'probability_witness':'3/4',
                         'upper_bound':'Product of the listed distinct trace factors over Q(t), using preceding all-length sector identities.',
                         'Hankel_first_length':2,'Hankel_dimension':order,
                         'exact_Hankel_minor':minor,'exact_Hankel_determinant':str(determinant),
                         'independent_modular_determinant_checks':modchecks,
                         'trace_weights':weights[label]})
    return {'schema':'matching-one.width4-topological-source-spectrum.v1',
            'source_pr':708,'source_head':definition['source_head'],'source_blob':actual_blob,
            'parametric_definition_sha256':hashlib.sha256(defraw).hexdigest(),
            'standing':'Exact finite-width source response; no new state exploration or Monte Carlo.',
            'identities':{'Z':'1 + M*sinh(s) + E*(cosh(s)-1)',
                          'source_mean':'[M*cosh(s)+E*sinh(s)]/[1+M*sinh(s)+E*(cosh(s)-1)]',
                          'source_susceptibility_at_zero':'E-M^2',
                          'shared_block_in_tilted_mean_numerator':'-2*sinh(s)*tr(B16^m)/(1+t)^(4m)'},
            'all_length_basis':'P0/P2 trace identities proved in preceding parametric delivery; linear combinations here need no new interpolation.',
            'raw_state_DP_vs_trace_lengths_checked':maximum,
            'newly_visible_factors_relative_to_M':sorted(name for name in weights['E'] if name not in weights['M']),
            'channels':channels,
            'first_60_exact_sequences_at_t3':sequences,'elapsed_seconds':time.perf_counter()-started,
            'limits':['Orders concern unnormalized fugacity sequences; physical normalization rescales eigenvalues and not their generic orders.',
                      'The normalized tilted mean is a ratio of sequences; no order-28 linear recurrence for that ratio is asserted.',
                      'The finite nonzero source example is exp(s)=2; no claim that every exceptional source or p has the same order.',
                      'Three-dimensional single-configuration source algebra is not a three-state dynamical closure.',
                      'No local continuum operator or original-U candidate has been identified.']}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--certificate',type=Path,default=BASE/'width4-rank-closure-certificate.json')
    parser.add_argument('--definition',type=Path,default=BASE/'width4-parametric-definition.json')
    parser.add_argument('--out',type=Path)
    args=parser.parse_args()
    for path in (args.certificate,args.definition):
        if not path.exists():parser.error('missing preceding delivery input: '+str(path))
    report=verify(args.certificate,args.definition)
    text=json.dumps(report,indent=2,allow_nan=False)+'\n'
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        with args.out.open('x',encoding='utf-8') as f:f.write(text)
    else:print(text,end='')

if __name__=='__main__':main()
