#!/usr/bin/env python3
"""Exact level-4 quasiprimary check in the c=0, h=1/3 V_<1,3> module.

The diagonal Potts singlet V_<1,3> has x=2/3.  Its level-4 chiral descendant
would give a bulk spin-4 field of x=14/3, so it is a serious *formal* competitor
below the thermal-family x=21/4 field unless matching/interchiral parity excludes
it from the matching-odd channel.
"""
from fractions import Fraction as F

# PBW level-4 basis:
# [L_-4, L_-3 L_-1, L_-2^2, L_-2 L_-1^2, L_-1^4]
# L1 -> [L_-3, L_-2 L_-1, L_-1^3] at h=1/3.
H=F(1,3)
L1=[
    [F(5),2*H,F(3),F(0),F(0)],
    [F(0),F(4),F(6),4*H+2,F(0)],
    [F(0),F(0),F(0),F(3),8*H+12],
]

# Level-3 singular vector:
# chi3=(4 L_-3 -24 L_-2 L_-1 +9 L_-1^3)|h>.
CHI3=[F(4),F(-24),F(9)]
# L_-1 chi3, reordered into the level-4 PBW basis.
NULL4=[F(8),F(-20),F(0),F(-24),F(9)]

# One exact L1 quasiprimary representative.
Q4=[F(-4),F(-15),F(10),F(0),F(0)]

def matvec(matrix,vector):
    return [sum(row[j]*vector[j] for j in range(len(vector))) for row in matrix]

def proportional(a,b)->bool:
    ratio=None
    for x,y in zip(a,b):
        if y==0:
            if x!=0: return False
            continue
        q=x/y
        if ratio is None: ratio=q
        elif q!=ratio: return False
    return True

def partitions(n:int)->int:
    if n<0:return 0
    dp=[0]*(n+1);dp[0]=1
    for k in range(1,n+1):
        for j in range(k,n+1):dp[j]+=dp[j-k]
    return dp[n]

def main()->int:
    assert matvec(L1,Q4)==[0,0,0]
    assert not proportional(Q4,NULL4)
    # The first singular level is 3; the next Kac solution with h=1/3 is level 9,
    # so through level 4 the quotient dimensions are p(n)-p(n-3).
    d3=partitions(3)-partitions(0)
    d4=partitions(4)-partitions(1)
    q4=d4-d3
    assert (d3,d4,q4)==(2,4,2)
    x_primary=F(2,3)
    x_spin4=x_primary+4
    assert x_spin4==F(14,3)
    residual_N=(x_spin4-2)/2
    root_L=x_spin4-F(5,4)
    assert residual_N==F(4,3)
    assert root_L==F(41,12)
    print("c=0, V_<1,3>: h=hbar=1/3, x_primary=2/3")
    print("chi3 = 4 L_-3 -24 L_-2 L_-1 +9 L_-1^3")
    print("Q4 = -4 L_-4 -15 L_-3 L_-1 +10 L_-2^2")
    print("L1 Q4=0; Q4 is nonzero modulo the level-3 null descendant")
    print(f"bulk spin4 x={x_spin4}; residual~N^-{residual_N}; root~L^-{root_L}")
    print("This field must be excluded from matching-odd D by parity/interchiral structure, not by nullness.")
    return 0

if __name__=="__main__":raise SystemExit(main())
