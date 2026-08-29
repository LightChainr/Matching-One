# Exact geometry closure and no-go rules for a pure g2 root bias

The master law in Issue #161 becomes sharply falsifiable once square, rectangular, hexagonal, and degree-2 child geometries are put in one frame.

## 1. The typed root observable

The clean object is the root or residual-to-slope response in which the common thermal-primary block cancels.  At fixed physical area define

```text
W(tau)=A^2 delta p_root(tau).
```

The pure Weierstrass law is

```text
W(tau)=C E4hat(tau),
E4hat(tau)=Im(tau)^2 E4(tau),                  (1)
```

with one common microscopic coupling `C`.  A raw residual need not satisfy (1) because it may retain a modulus-dependent primary block.

## 2. A three-modulus amplitude vector

The exact CM/Hecke identities give

```text
E4hat(2i)/E4hat(i)=11/4,
E4hat(rho)=0,
rho=exp(i pi/3).
```

Hence

```text
(W(i),W(2i),W(rho))=A(4,11,0).                (2)
```

The two no-fit residuals are

```text
4W(2i)-11W(i)=0,
W(rho)=0.                                      (3)
```

One geometry enhances the signal and the other annihilates it.  A new amplitude per modulus would destroy the content of this test.

## 3. Degree-2 hexagonal closure

On the three children

```text
(2rho,rho/2,(rho+1)/2),
```

the full chiral response must be

```text
Y=B(1,zeta,zeta^2).                            (4)
```

Thus

```text
Y1-zeta Y0=0,
Y2-zeta^2 Y0=0,
Y0+Y1+Y2=0.
```

For a real aligned coupling, (4) reduces to

```text
1:-1/2:-1/2.
```

The complex relation is stronger.  From the #164 ring grading, E6 contamination occupies child DFT character `r=0`, while E4-squared contamination occupies `r=2`; pure g2 occupies `r=1`.

## 4. A stronger C4/hexagonal selection rule

A local square-lattice coupling invariant under microscopic C4 can contain only spins divisible by four.  In the ordinary ring `E4^aE6^b`, this means the total weight is divisible by four.

At exact rho, every term with `a>0` vanishes.  A surviving pure `E6^b` term is C4-allowed only when `b` is even.  Therefore the first ordinary holomorphic survivor after the g2 zero is

```text
E6^2, weight 12,
delta p_root ~ L^-12=N^-6.                     (5)
```

Neither E6/weight 6 nor E4E6/weight 10 is allowed by the microscopic C4 coupling.  E4-squared/weight 8 is allowed but also vanishes at rho.

This is an exact-shape statement.  On an integer-period Pell sequence with

```text
tau_L-rho=O(L^-2),
```

the simple E4 zero leaks the leading term as

```text
L^-4 E4(tau_L)=O(L^-6)=O(N^-3),                (6)
```

which dominates the intrinsic exact-rho `N^-6` survivor.  Confusing (5) and (6) would assign the wrong correction spectrum.

## 5. When pure g2 is impossible

The one-coupling law is rejected as the sole leading mechanism if any of these persist in a consistently typed observable:

1. a nonzero leading `N^-2` coefficient at exact rho;
2. failure of the `11/4` rectangular ratio after the same slope normalization;
3. leading hex-child DFT support outside `r=1`;
4. failure of the cube-root phase closure despite an acceptable real `cos(4theta)` fit;
5. a claimed ordinary local-g2 term with a half-integer/non-ring correction exponent.

These failures point to a logarithmic/quasimodular tangent, another spin family, a vector-valued sector response, or quotient/topological dependence.

## 6. What passing does not prove

The modular space `M4` is one-dimensional.  Therefore every ordinary scalar/chiral weight-4 modular one-point has E4 shape.  Passing (2)--(4) supports a local weight-4 law but cannot by itself distinguish the thermal Q4 module from another weight-4 module.

In particular, the sector-valued vacuum-KdV response of #231 can be nonzero at rho by occupying a nontrivial sector character.  That behavior is not a scalar g2 root law and must not be averaged into one.

## Frozen score

Reconstruct intrinsic centers, H4 projections, slopes, root responses, area normalization, and all residuals within the same synchronized replicas.  Score (3) jointly, then score the complex child residuals in (4).  Fit only the one common amplitude `C` if an amplitude representation is desired.

## Boundary

All modular ratios and selection rules are exact.  Their application to matching roots remains conditional on a common typed lattice-to-CFT bridge and microscopic frame.
