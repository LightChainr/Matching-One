#!/usr/bin/env python3
"""Enumerate formal c=0 non-diagonal spin-4 Kac competitors.

At percolation on the critical Potts branch,
h_{r,s}=((2r-3s)^2-1)/24,
a non-diagonal primary W(r,s) has weights (h_{r,s},h_{r,-s}).  Its
conformal spin is -r*s and its scaling dimension is

    x = (4 r^2 + 9 s^2 - 1)/12.

The Potts space of states further attaches an S_Q representation Xi(r,s), so
formal Kac admissibility alone is NOT enough to make a field an allowed
S_Q-singlet lattice perturbation.

This utility prints the lowest r>=2 spin-4 formal primaries and the induced
finite-size powers for a dimensionless torus observable:

    residual ~ L^(2-x) = N^(-(x-2)/2)
    root bias ~ residual/L^(3/4) = L^(-(x-5/4)).
"""
from __future__ import annotations

from fractions import Fraction as F


def h(r: F, s: F) -> F:
    return ((2*r-3*s)**2-1)/24


def record(r: int) -> dict[str, object]:
    s=F(4,r)
    hl=h(F(r),s); hr=h(F(r),-s)
    spin=hl-hr; x=hl+hr
    assert spin==-4
    residual_L=x-2
    residual_N=residual_L/2
    root_L=x-F(5,4)
    return {
        "r":r,"s":s,"h_left":hl,"h_right":hr,"x":x,"spin":spin,
        "residual_L_power":residual_L,"residual_N_power":residual_N,
        "root_L_power":root_L,
    }


def main()->int:
    rows=[record(r) for r in range(2,9)]
    assert rows[0]["x"]==F(17,4)
    assert rows[0]["residual_N_power"]==F(9,8)
    assert rows[0]["root_L_power"]==3
    print("formal c=0 non-diagonal spin -4 primaries (r>=2, s=4/r)")
    for row in rows:
        print(
            f"r={row['r']} s={row['s']} x={row['x']} "
            f"residual~N^-{row['residual_N_power']} root~L^-{row['root_L_power']}"
        )
    print("\nRepresentation-theory warning:")
    print("Potts W(r,s) also carries Xi(r,s). Formal Kac states are not automatically S_Q singlets.")
    print("For r=2, Xi(2,2)=Xi(2,0)=[2], so x=17/4 is a two-cluster/non-singlet sector at generic Q.")
    return 0


if __name__=="__main__": raise SystemExit(main())
