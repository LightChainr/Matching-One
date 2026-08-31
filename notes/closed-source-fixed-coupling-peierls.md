# A uniform fixed-coupling winding bound from the local colour gas

**Theorem.** For the original projected closed-source law on the axis
LxL torus, let 4 divide L, L>=16, and fix an integer m=exp(t)>=256.
For every h>0, including an L-dependent matching root,

```text
P_star(r=1) <= min{1, [4m^2 L^2/(3(1-rho_m))] rho_m^L},
rho_m = 3(2m^(-15/16))^(1/4) < 1.                  (1)
```

This is genuinely fixed t followed by L growing, uniform in the activity
h. It does not require h=1, a pure-empty/pure-full approximation, or a
previous phase-coexistence theorem. The constants are deliberately coarse;
256 is a sufficient integer coupling, not a transition estimate.

The proof starts from `85fd4923` and the exact
[local colour gas](closed-source-local-colour-gas.md). No new numerical
score, simulation, enumeration or coupling point enters it. The final
section identifies precisely what (1) does and does not settle for the
two-geometry original U.

## 1. Reflection positivity of the local law, before the rank projection

Write Qc=m^2 and let a spin take values 0,1,...,Qc. Put occupation n=0
for spin0 and n=1 otherwise. An equivalent edge gauge for the existing
local colour law is

```text
w_col(c) = h^K product_(NN edges uv) R(c_u,c_v),
R00=Raa=1,   R0a=Ra0=1/m,
Rab=0 for distinct active a,b.
```

Indeed the occupation marginal is `h^K m^(-Bmix+2C_B)`. The original
law has the additional factor m^(-r), exactly as before. If e_a are the
standard basis of R^Qc, then

```text
v_a=e_a,    v_0=(1/m) sum_a e_a,    R(c,d)=<v_c,v_d>.  (2)
```

Here ||v0||=1; R is positive semidefinite with rank Qc. Expand the Gram
factor on each edge crossing a bond-reflection plane. The unnormalized
expectation of F times its reflected copy becomes a sum of squares of
half-torus partition functionals. Identical site activities on reflected
vertices preserve this factorization. Thus the local law is reflection
positive across every horizontal or vertical bond plane on the axis
torus, at any h>0. Hard zeros in R cause no problem: the argument is a
finite Gram expansion, and Cauchy--Schwarz also holds on its null-space
quotient.

The resulting chessboard inequality applies to disjoint 2x1 dominoes.
There are M=N/2 dominoes per tiling. For events D_j on distinct dominoes,
their joint probability is at most the product of their fully reflected
dissemination probabilities to the power 1/M. This is the coordinatewise
rectangular-block version of the reflection/Cauchy--Schwarz proof in
[Biskup, Section5.3, Theorem5.8](https://www.math.ucla.edu/~biskup/PDFs/private/2008-2010/Biskup-LNM08.pdf).
That proof uses an even number of blocks along each coordinate, not a
power-of-two assumption; here these numbers are L/2 and L. Its proof is
unchanged for disjoint blocks with reflection planes through bonds.

We do **not** assert reflection positivity of the rank-projected law.
Instead its exact Radon--Nikodym factor gives, for every event D,

```text
P_star(D)=E_col[1_D m^(-r)]/E_col[m^(-r)]
         <= m^2 P_col(D).                           (3)
```

The projection is retained, as a fixed-volume-independent comparison
constant at fixed m. It is not discarded from phase weights or roots.

## 2. An exactly evaluated mixed-domino dissemination

On one horizontal domino, prescribe the occupation pattern 01, leaving
the active colour unrestricted. Reflection dissemination makes alternating
vacant and occupied vertical stripes of width two: the column pattern
is 0110, repeated. All rows have the same occupation pattern. Therefore

```text
K=N/2,   Bmix=N/2,   C_B=L/4,
Z_dis = Qc^(L/4) h^(N/2) m^(-N/2).
```

Different occupied stripes may choose their colours independently. The
factor Qc^(L/4) is included; dropping it would give a false finite-size
bound. The unconstrained local partition sum satisfies
`Z_col>=max(1,h^N)`, from the empty state and one full active-colour state.
Thus the chessboard norm of this ordered event obeys

```text
z(01)=(Z_dis/Z_col)^(2/N)
     <= m^(-1+1/L) min(h,h^(-1))
     <= beta_L,      beta_L=m^(-1+1/L).              (4)
```

The same calculation holds for 10, either horizontal dimer tiling, and
both vertical dimer tilings. Summing over the 2^k possible ordered patterns
on k distinct dominoes in a single tiling gives

```text
P_col(all k specified dimers are mixed) <= (2 beta_L)^k.  (5)
```

No independence assumption is used. This is a multi-edge estimate, not
a one-site finite-energy estimate. Most importantly, (4)--(5) are uniform
in h; the entropy displacement of the matching root never has to be
set to zero in this argument.

## 3. From mixed dimers to an essential interface

Resolve alternating cut-dual vertices around occupied corners as in
the [winding-barrier proof](closed-source-winding-barrier.md). A rank-one
configuration has an essential boundary curve. Its unsmoothed dual walk
has length n>=L, uses every cut edge at most once, and never immediately
backtracks. All n crossed primal edges are mixed.

Partition the primal edge set into four domino tilings: horizontal edges
with even or odd starting x, and vertical edges with even or odd starting
y. At least ceil(n/4) of a specified contour's edges belong to one
tiling. They are distinct dominoes. If 2 beta_L<1, equation(5) implies

```text
P_col(the specified contour is mixed) <= (2 beta_L)^(n/4).
```

There are at most `4N 3^(n-1)` oriented nonbacktracking dual walks of
length n with a chosen starting vertex, hence no more possible boundary
curves. Overcounting closures and repeated vertices is harmless. Summing
this bound over n>=L and using (3) proves

```text
P_star(r=1) <= [4m^2 N/(3(1-rho_L))] rho_L^L,
rho_L=3(2m^(-1+1/L))^(1/4),             rho_L<1.     (6)
```

For L>=16, rho_L<=rho_m in (1). At m=256,
`rho_m=3*2^(-13/8)<1`, since `3^8<2^13`; increasing m only improves
the bound. This completes the uniform theorem. Multiple winding
components, contractible holes and dual corner contacts do not change
the counting argument.

## 4. A moving single-geometry root is already covered

Let h_L be the unique matching root of this one axis geometry. No value
of h_L was assumed in (1). At that root, exactly

```text
P_star(r=0)=P_star(r=2)=(1-P_star(r=1))/2.             (7)
```

Thus (1) gives exponentially accurate coexistence of the two **rank
sectors** at a fixed finite coupling. It says nothing about concentration
on the individual empty and full configurations. Local defects are
allowed throughout the proof; the previously necessary condition
N exp(-2t)->0 for a literal pure-state mixture is neither imposed nor
silently reused.

There is also a useful normalized derivative consequence. Since q=r-1
is increasing and its one-site increment is at most2, `2K-q` is increasing.
The original law's FKG property gives

```text
Cov_star(q,K) >= Var_star(q)/2.
```

At its own root, `Var(q)=1-P1`, while E=1-1_(r=1) gives
`|Cov(E,K)|<=N P1`. Because h derivatives are covariances with K divided
by h,

```text
|E_h/Q_h|_(h=h_L) <= 2N P1/(1-P1).                 (8)
```

So even the single-geometry thermal slope ratio is exponentially small
in L, up to a polynomial prefactor, once the right side of (1) is below1.
This is not just an unnormalized probability statement.

## 5. The exact remaining boundary for the original two-geometry U

The present reflection argument is for axis square tori with 4|L.
An oblique Gaussian period lattice need not be invariant under the
coordinate bond reflections. Its Gram edge matrix alone does not make
that quotient reflection positive, so (1) has not been proved for the
original tilted companion by this calculation.

Even if winding bounds for both geometries are supplied, their **pooled**
root requires a separate denominator distinction. At the equal-weight
pooled root write

```text
Q_f=a,   Q_s=-a,
Pbar=(P1_f+P1_s)/2,
kappa=(Var_f(q)+Var_s(q))/2=1-Pbar-a^2.
```

These are existing rank means and variances, not an added model or source.
Applying the same FKG inequality within each geometry and retaining the
original angular denominator Delta gives the rigorous sufficient bound

```text
|U/A_N| <= 4N Pbar/(|Delta| kappa),
A_N=N^(13/8)/2.                                     (9)
```

Small Pbar alone does not bound this ratio if a tends to1: the pooled
zero could balance a nearly rank-two first geometry against a nearly
rank-zero second geometry while both within-geometry slopes are small.
Replacing the within-geometry covariance by the covariance of a mixture
of geometries would incorrectly erase this issue. A subexponential lower
bound on kappa, together with exponential winding bounds for both
geometries, would close the original fixed-t U conclusion through (9).
Neither bulk-pressure equality nor the bounded factor m^(-r) supplies
that lower bound by itself.

The new result is therefore a closed, explicit fixed-t Peierls theorem
for the projected axis law and its own moving-root slope ratio. It does
not assert an unproved full phase diagram, a critical activity exactly1,
or completion of the oblique pooled-U thermodynamic problem. The
unresolved object is the existing ratio in (9), with its within-geometry
normalization retained.
