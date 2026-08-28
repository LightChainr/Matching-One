#!/usr/bin/env python3
"""Lowest formal non-diagonal spin-4 Potts primaries and singlet filter.

At percolation c=0, for W(r,s) with weights (h_{r,s},h_{r,-s}),
spin=-r*s.  Setting r*s=4 gives formal spin -4 primaries.

The conformal data alone are insufficient: each W(r,s) carries a Potts S_Q
representation Xi(r,s), periodic under s->s+1.  From the explicit generic-Q
space-of-states decomposition of Jacobsen--Ribault--Saleur (2023):

  Xi(2,0)   = [2]
  Xi(3,1/3) = [21]
  Xi(4,0)   contains [] (the trivial representation)

Hence W(2,2) and W(3,4/3) are not ordinary generic-Q singlet perturbations,
whereas W(4,1) is the first member of this formal non-diagonal spin-4 sequence
whose internal representation contains the singlet.
"""
from __future__ import annotations
from fractions import Fraction as F


def h(r:F,s:F)->F:
    return ((2*r-3*s)**2-1)/24

def row(r:int)->dict[str,object]:
    s=F(4,r)
    hl=h(F(r),s); hr=h(F(r),-s); x=hl+hr
    assert hl-hr==-4
    return {
        "r":r,"s":s,"s_mod_1":s%1,"h":hl,"hbar":hr,"x":x,
        "residual_N_power":(x-2)/2,
        "root_L_power":x-F(5,4),
    }

def main()->int:
    rows=[row(r) for r in (2,3,4)]
    assert rows[0]["x"]==F(17,4)
    assert rows[1]["x"]==F(17,4)
    assert rows[2]["x"]==F(6)
    labels={2:"Xi(2,2)=Xi(2,0)=[2] : nontrivial",3:"Xi(3,4/3)=Xi(3,1/3)=[21] : nontrivial",4:"Xi(4,1)=Xi(4,0) contains [] : singlet allowed"}
    for r in (2,3,4):
        q=rows[r-2]
        print(f"r={r} s={q['s']} x={q['x']} residual~N^-{q['residual_N_power']} root~L^-{q['root_L_power']}")
        print("  "+labels[r])
    print("\nConclusion at generic Q: the first two x=17/4 primary spin-4 fields are not ordinary S_Q-singlet lattice perturbations; the first non-diagonal primary in this sequence whose Xi contains [] is x=6.")
    print("Q->1 logarithmic mixing is a separate caveat and must be tested rather than erased.")
    return 0

if __name__=="__main__": raise SystemExit(main())
