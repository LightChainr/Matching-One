#!/usr/bin/env python3
"""Exact two-rotation H4/H8 identifiability contract for a REAL C3 readout.
This is a conditional forward-design calculation, not a fitted field identity.
No empirical vector, MC, new geometry or lattice enumeration is generated.
"""
from __future__ import annotations
import argparse,json,sys
from fractions import Fraction as F
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
def rot(x):return s.Matrix([[s.cos(x),-s.sin(x)],[s.sin(x),s.cos(x)]])
def design(spin,delta):return s.eye(2).col_join(rot(spin*delta))/s.sqrt(2)
def simp(M):return M.applyfunc(s.simplify)
def gp(a,b,n):
    z=(1,0)
    for _ in range(n):z=(z[0]*a-z[1]*b,z[0]*b+z[1]*a)
    return z

def build():
    checks=0;table=[]
    for numerator in (0,1,2,3,4,6,8,12):
        delta=s.pi*s.Rational(numerator,24)
        X=design(4,delta);Y=design(-8,delta)
        cross=simp(X.T*Y)
        assert simp(X.T*X)==s.eye(2);checks+=1
        assert simp(cross.T*cross)==s.eye(2)*s.simplify(s.cos(6*delta)**2);checks+=1
        residual=simp(X.T*(s.eye(4)-Y*Y.T)*X)
        expected=s.simplify(s.sin(6*delta)**2)
        assert residual==s.eye(2)*expected;checks+=1
        joint_rank=X.row_join(Y).rank()
        assert joint_rank==(2 if s.simplify(s.sin(6*delta))==0 else 4);checks+=1
        signed_alias=s.simplify(s.sin(12*delta))==0
        positive_alias=s.simplify(s.cos(12*delta))==1
        if signed_alias:assert s.simplify(s.cos(12*delta)) in (-1,1)
        checks+=1
        table.append({'angle_degrees':str(s.Rational(180*numerator,24)),
                      'shared_complex_amplitude_joint_rank':joint_rank,
                      'isotropic_balanced_wrong_model_residual_energy_fraction':str(expected),
                      'independent_signed_real_gain_alias':bool(signed_alias),
                      'independent_positive_real_gain_alias':bool(positive_alias),
                      'signed_real_gain_phase_test_sin2':str(s.simplify(s.sin(12*delta)**2))})
    # Direct real C3 readout -> first Fourier coordinate, independently of the 4x2 reduction.
    direct_checks=0
    for numerator in (0,1,2,3,4,6):
        theta=s.pi*s.Rational(numerator,24)
        for spin in (4,8):
            A=1+2*s.I
            vals=[s.re(s.expand_complex(A*s.exp(s.I*spin*(theta+2*s.pi*j/3)))) for j in range(3)]
            z=s.Rational(2,3)*sum(vals[j]*s.exp(-2*s.pi*s.I*j/3) for j in range(3))
            target=A*s.exp(4*s.I*theta) if spin==4 else s.conjugate(A)*s.exp(-8*s.I*theta)
            assert s.simplify(s.expand_complex(z-target))==0
            direct_checks+=1
    checks+=direct_checks
    # Exact counterexample: a 15-degree second rotation, unknown signed gain.
    d=s.pi/12;v=s.Matrix([1,2]);z4=rot(4*d)*v;z8=-rot(-8*d)*v
    assert simp(z4-z8)==s.zeros(2,1);checks+=1
    # At 7.5 degrees any nonzero real gains have distinct relative phases.
    d=s.pi/24
    for gain1,gain2 in [(1,1),(2,-3),(-2,3),(-2,-3)]:
        prod=s.expand_complex(gain2*(1+2*s.I)*s.exp(4*s.I*d)*s.conjugate(gain1*(1+2*s.I)))
        h4=s.simplify(s.im(s.expand_complex(prod*s.exp(-4*s.I*d))))
        h8=s.simplify(s.im(s.expand_complex(prod*s.exp(8*s.I*d))))
        assert h4==0 and h8!=0;checks+=1
    gaussian=[]
    for a,b in [(1,1),(2,1),(8,1),(7,4)]:
        re,im=gp(a,b,12);norm=(a*a+b*b)**6
        assert re*re+im*im==norm*norm;checks+=1
        gaussian.append({'multiplier':[a,b],'exp_i12theta_real':str(F(re,norm)),
                         'exp_i12theta_imag':str(F(im,norm)),
                         'fixed_gain_residual_fraction':str(F(norm-re,2*norm)),
                         'signed_gain_phase_separation_sin2':str(F(im*im,norm*norm)),
                         'boundary':'conditional phase arithmetic only; multiplying a torus changes size and does not by itself establish common amplitude/phase transport'})
    # Full-rank covariance whitening changes separation, not identifiability.
    covariance=s.Matrix([[2,1,0,0],[1,3,1,0],[0,1,4,1],[0,0,1,2]])
    assert covariance.is_positive_definite;checks+=1
    W=covariance.inv()
    gls=[]
    for num in (0,1,2,4):
        d=s.pi*s.Rational(num,24);X=design(4,d);Y=design(-8,d)
        Fmat=simp(X.T*(W-W*Y*(Y.T*W*Y).inv()*Y.T*W)*X)
        if num in (0,4):assert Fmat==s.zeros(2)
        else:assert s.simplify(Fmat.det())>0 and s.simplify(Fmat[0,0])>0
        checks+=1
        gls.append({'angle_degrees':str(s.Rational(180*num,24)),
                    'profile_information_matrix':[[str(s.simplify(x)) for x in row] for row in Fmat.tolist()]})
    return {'scope':'conditional two-rotation prediction vectors for the same phase-calibrated real C3 observable; no empirical field decision',
            'total_exact_checks':checks,'direct_C3_checks':direct_checks,'angles':table,
            'gaussian_arithmetic':gaussian,'gls_controls':gls,
            'counterexample_15_degree':{'amplitude':'1+2i','H4_relative_gain':1,'H8_relative_gain':-1,'same_two_readouts':True},
            'shared_amplitude_prediction':'z2=e^(i4delta)z1 or z2=e^(-i8delta)z1',
            'signed_gain_prediction':'Im[z2*conj(z1)*e^(-i4delta)]=0 or Im[z2*conj(z1)*e^(i8delta)]=0; both amplitudes nonzero',
            'arbitrary_complex_gain':'unidentifiable at every angle: each model spans all two-complex-coordinate data',
            'nonclaims':['no existing C3 data reanalysis','no choice of a new lattice geometry','no claim that raw trace or normalized source response obeys a pure spin law','isotropic angle optimum is not optimal for arbitrary covariance','conditional transport must be justified before using these vectors']}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    data=build();a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'exact_checks':data['total_exact_checks'],'output':str(a.out)}))
if __name__=='__main__':main()
