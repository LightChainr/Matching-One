# Executable triangular-site E6 control from the center score

The triangular self-matching identity kills the center value and pins the root, but it does **not** kill the center derivative.  That distinction supplies a cheap nondegenerate E6 instrument.

## 1. Correct self-matching parity

Complement/matching symmetry gives

```text
M_L(1-p)=-M_L(p).
```

Therefore

```text
M_L(1/2)=0,
K_L:=M'_L(1/2) is complement-even and may be nonzero.
```

The committed exact triangular tori already give

```text
L=2: K=3,
L=3: K=261/64,
L=4: K=5147/1024.
```

Thus no threshold search is needed.  The replacement for the degenerate root is the same-N sextic orientation component of `K`, divided by its scalar mean.

## 2. Ward-active level-6 line

After quotienting the `c=0,h=5/8` level-2 null, the level-6 quasiprimary quotient is two-dimensional, but the torus Ward map has rank one.  One active representative is

```text
Q6=-25L_-6+28L_-5L_-1-56L_-4L_-2+35L_-3^2,
```

with

```text
<Q6 phi>/<phi>=-3975 g3/224
                 =-1325 pi^6 E6/252.
```

The second quotient direction is Ward-null.  Therefore a failed lattice overlap is meaningful: it may select the null line rather than refute the one-dimensional `M6` image.

## 3. Minimal same-N orientation projector

Use `omega=exp(i pi/3)` and `N(a+bomega)=a^2+ab+b^2`.  The smallest convenient reflection-even pair with distinct sextic cosines is

```text
g1=1+9omega,
g2=5+6omega,
N(g1)=N(g2)=91.
```

Their exact normalized sixth phases have

```text
cos(6theta1)= 644221/753571,
cos(6theta2)=-716579/753571,
Delta cos6=1360800/753571.
```

Define

```text
D6=[K(g1)-K(g2)]/Delta cos6,
Kbar=[K(g1)+K(g2)]/2,
R6=D6/Kbar,
A6=N^3 R6.                                      (1)
```

Same `N` cancels the scalar size dependence.  Dividing by `Kbar` cancels the thermal-primary block, analogously to the residual-to-slope construction in #161.  The bridge assumption is that the remaining direction overlaps the Ward-active Q6 line.

## 4. Exact Eisenstein multiplier discriminator

For

```text
m=1+omega,
N(m)=3,
M_m=[[1,-1],[1,2]],
```

one has exactly

```text
m^-6=-1/27.
```

Apply `m` to both N91 orientations, producing an N273 child pair.  The two equivalent predictions are

```text
R6(child)/R6(parent)=-1/27,
A6(child)/A6(parent)=-1.                         (2)
```

A spin-12 alias has scaled factor `+1`, so (2) is a sign-level discriminator rather than merely a radial decay test.

## 5. Independent modulus control

At the natural hexagonal modulus, `E6hat(rho)` is nonzero.  On its degree-2 Hecke children,

```text
E6hat(child_j)/E6hat(rho)=11/4
```

for all three children, while

```text
E6hat(i)=0.
```

Thus a shape-normalized E6 coefficient obeys

```text
W6(child1)=W6(child2)=W6(child3),
4W6(child)-11W6(parent)=0,
W6(i)=0.                                          (3)
```

Equations (2) and (3) are complementary: (2) is an Eisenstein similarity/phase test; (3) is a modulus/Hecke test.  A generic scalar contaminant may mimic equal children but cannot simultaneously satisfy the signed multiplier and square zero after the same center-slope normalization.

## 6. Minimal run

At exact `p=1/2`, estimate

```text
K=E[M * 4(k-N/2)]
```

with common random fields for the two N91 orientations and then their N273 children.  No off-critical scan or root estimation is needed.  Score both forms of (2), retaining their covariance; assert `M(1/2)=0` and `Kbar!=0` as runtime controls.

## Boundary

The Ward coefficient, Eisenstein phases, multiplier factors, Hecke ratios, and tiny center derivatives are exact.  The lattice overlap of (1) with the active Q6 line is the one conditional step.  Raw `K` is not itself an E6 amplitude, and self-matching does not imply `M(p)` vanishes away from `p=1/2`.
