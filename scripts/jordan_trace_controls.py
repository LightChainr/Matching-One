#!/usr/bin/env python3
"""Exact controls: a linear eigenvalue split does not certify semisimplicity.

Build positive stochastic / irreducible Markov families with identical ordinary
trace data at every parameter and length, but different Jordan structure and
ordinary probability readouts. No percolation or LCFT identification is made.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import time
import sympy as sp


def matrices():
    e=sp.ones(3,1);u=sp.Matrix([1,-1,0]);v=sp.Matrix([1,1,-2])
    P=e*e.T/3
    Q=sp.eye(3)-P
    D=u*u.T/2-v*v.T/6
    N=u*v.T/6
    return P,Q,D,N


def markov_pair(t):
    P,Q,D,N=matrices()
    GD=-Q+t*D
    GJ=GD+N/2
    return GD,GJ


def encode_matrix(a):return [[str(x) for x in row] for row in a.tolist()]


def run_checks():
    start=time.perf_counter();t,x,z,s=sp.symbols('t x z s');P,Q,D,N=matrices()
    O=sp.zeros(3)
    assert P*P==P and Q*Q==Q and P*Q==O
    assert N*N==O and N!=O and P*N==O and N*P==O
    assert D*D==Q and D*N==N and N*D==-N
    GD,GJ=markov_pair(t)
    assert GD*sp.ones(3,1)==sp.zeros(3,1)==GJ*sp.ones(3,1)
    assert sp.ones(1,3)*GD==sp.zeros(1,3)==sp.ones(1,3)*GJ
    char=z*((z+1)**2-t*t)
    assert sp.expand(GD.charpoly(z).as_expr()-char)==0
    assert sp.expand(GJ.charpoly(z).as_expr()-char)==0
    minima=[]
    for G in (GD,GJ):
        vals=[G[i,j].subs(t,a) for i in range(3) for j in range(3) if i!=j
              for a in (-sp.Rational(1,4),sp.Rational(1,4))]
        minima.append(min(vals))
    assert minima==[sp.Rational(1,6),sp.Rational(1,12)]
    TD=sp.eye(3)+GD/2;TJ=sp.eye(3)+GJ/2
    for T in (TD,TJ):
        assert min(T.subs(t,a)[i,j] for i in range(3) for j in range(3)
                   for a in (-sp.Rational(1,4),sp.Rational(1,4)))>0
    G0D=GD.subs(t,0);G0J=GJ.subs(t,0)
    assert (G0D+sp.eye(3)).nullspace().__len__()==2
    assert len((G0J+sp.eye(3)).nullspace())==1
    assert G0D*(G0D+sp.eye(3))==O
    assert G0J*(G0J+sp.eye(3))!=O
    assert G0J*(G0J+sp.eye(3))**2==O

    # Polynomial trace equality on a finite grid of lengths is a regression
    # check. The all-length proof follows from triangular blocks / spectrum.
    AD=sp.eye(3);AJ=sp.eye(3);thermal_checks=0
    traces=[]
    for m in range(13):
        a=sp.expand(sp.trace(AD));b=sp.expand(sp.trace(AJ))
        expected=1+((1+t)/2)**m+((1-t)/2)**m
        assert sp.expand(a-expected)==0==sp.expand(a-b)
        for order in range(5):
            assert sp.diff(a-b,t,order)==0;thermal_checks+=1
        traces.append(str(sp.factor(a)))
        AD=(AD*TD).applyfunc(sp.expand);AJ=(AJ*TJ).applyfunc(sp.expand)

    # A completely ordinary probability response observes the nilpotent.
    yd=[];yj=[];AD=sp.eye(3);AJ=sp.eye(3)
    for m in range(12):
        yd.append(AD[0,2]);yj.append(AJ[0,2])
        expectedD=(1-sp.Rational(1,2)**m)/3
        expectedJ=expectedD-sp.Rational(m,6)*sp.Rational(1,2)**m
        assert yd[-1]==expectedD and yj[-1]==expectedJ
        assert 0<=yd[-1]<=1 and 0<=yj[-1]<=1
        AD=AD*TD.subs(t,0);AJ=AJ*TJ.subs(t,0)
    HD=sp.Matrix(2,2,lambda i,j:yd[i+j]);HJ=sp.Matrix(3,3,lambda i,j:yj[i+j])
    detD=sp.factor(HD.det());detJ=sp.factor(HJ.det())
    assert detD!=0 and detJ!=0
    # Resolve exactly the fixed-operator minimal polynomials on the readout.
    for values,poly in [(yd,(x-1)*(x-sp.Rational(1,2))),
                        (yj,(x-1)*(x-sp.Rational(1,2))**2)]:
        cc=list(reversed(sp.Poly(poly,x).all_coeffs()))
        assert all(sum(cc[k]*values[i+k] for k in range(len(cc)))==0
                   for i in range(len(values)-len(cc)+1))

    # Counterexample 1: analytic linear split THROUGH a genuine Jordan block.
    C1=sp.Matrix([[1+t,1],[0,1-t]])
    assert sp.expand(C1.charpoly(x).as_expr()-(x-1-t)*(x-1+t))==0
    assert len((C1.subs(t,0)-sp.eye(2)).nullspace())==1
    # Counterexample 2: semisimple at zero and pointwise diagonalizable, yet
    # eigenbranches +/-t^(3/2) fail holomorphicity at zero.
    C2=sp.Matrix([[0,t],[t*t,0]])
    assert C2.subs(t,0)==sp.zeros(2)
    assert sp.expand(C2.charpoly(x).as_expr()-(x*x-t**3))==0
    # Counterexample 3: the same issue occurs even for an AFFINE pencil.
    C3=sp.Matrix([[0,t,0],[0,0,t],[t,0,1]])
    cp=sp.expand(C3.charpoly(x).as_expr())
    assert sp.expand(cp-(x*x*(x-1)-t**3))==0
    disc=sp.factor(sp.discriminant(cp,x))
    assert disc==-t**3*(4+27*t**3)
    assert len(C3.subs(t,0).nullspace())==2

    # The parameter jet is a different matrix from the physical matrix.
    A=sp.diag(1+t,1-t);A0=A.subs(t,0);A1=A.diff(t)
    jet=A0.row_join(A1).col_join(sp.zeros(2).row_join(A0))
    assert (jet-sp.eye(4))**2==sp.zeros(4) and jet!=sp.eye(4)
    for m in range(1,10):
        assert (jet**m)[:2,2:]==(A**m).diff(t).subs(t,0)

    return {'schema':'matching-one.jordan-trace-controls.v1',
        'scope':'exact finite analytic matrix controls; NOT a percolation or LCFT candidate identification',
        'domain':'-1/4 <= t <= 1/4','sympy_version':sp.__version__,
        'P':encode_matrix(P),'Q':encode_matrix(Q),'D':encode_matrix(D),'N':encode_matrix(N),
        'semisimple_generator':encode_matrix(GD),'defective_at_zero_generator':encode_matrix(GJ),
        'generator_offdiagonal_lower_bounds':list(map(str,minima)),
        'stationary_distribution':['1/3']*3,
        'shared_generator_characteristic_polynomial':str(sp.factor(char)),
        'shared_stochastic_power_trace':'1+((1+t)/2)^m+((1-t)/2)^m',
        'shared_semigroup_trace':'1+2*exp(-s)*cosh(s*t)',
        'nilpotent_coefficient':'1/2','geometric_multiplicities_at_minus_one':[2,1],
        'minimal_generator_polynomials':['z*(z+1)','z*(z+1)^2'],
        'thermal_trace_regressions':thermal_checks,'regression_max_length':12,
        'all_length_proof':'ordinary traces depend only on eigenvalues with algebraic multiplicities; both analytic families have identical characteristic polynomial',
        'probability_readout':{'source':'delta_state_1','readout':'indicator_state_3',
          'semisimple':'(1-2^(-m))/3','Jordan':'(1-2^(-m))/3-m*2^(-m)/6',
          'semisimple_values':list(map(str,yd)),'Jordan_values':list(map(str,yj)),
          'minimal_recurrence_orders':[2,3],'hankel_determinants':[str(detD),str(detJ)],
          'continuous_time_difference':'-s*exp(-s)/6'},
        'counterexamples':{
          'linear_split_defective':encode_matrix(C1),
          'nonlinear_pencil_semisimple_Puiseux':encode_matrix(C2),
          'affine_pencil_semisimple_Puiseux':encode_matrix(C3),
          'affine_pencil_characteristic':str(cp),'affine_pencil_discriminant':str(disc),
          'analytic_branch_obstruction':'an analytic root with x(0)=0 would have integer order k; x^2*(x-1)=t^3 forces 2k=3'},
        'jet_lift':encode_matrix(jet),'jet_lift_nontrivial_Jordan':True,
        'not_claimed':['No new general Jordan or realization theorem.',
          'No numerical eigenvalue fit certifies exact Jordan structure.',
          'Noncommuting marks / specified matrix elements can distinguish these examples.',
          'The Markov pair is not a model of square-site percolation.'],
        'elapsed_seconds':time.perf_counter()-start}

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path)
    args=ap.parse_args();text=json.dumps(run_checks(),indent=2,allow_nan=False)+'\n'
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        with args.out.open('x') as f:f.write(text)
    else:print(text,end='')
