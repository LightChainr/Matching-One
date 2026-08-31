# Oblique winding: an exact order-25 twist comparison

**New reduction.** The Gaussian torus with periods `(4k,3k),(-3k,4k)`
is exactly a `W x H = 25k x k` rectangle with a horizontal shift by 7k
at its vertical seam. The shift has order25, independent of k. Cutting
at two suitably chosen rows and applying Gram Cauchy--Schwarz removes
the twist from a constrained numerator. This gives a quantitative
comparison with the already controlled axis rectangular law:

```text
P_star,oblique(r=1)
 <= [4m^2 N/(3(1-rho_k))] exp(Delta_k) rho_k^(7k),
rho_k=3(2m^(-1+1/k))^((1-2/k)/4),
Delta_k=log Z_rect-log Z_twist >=0,       N=25k^2.    (1)
```

Equation(1) holds when `4|k`, `k>=16`, `rho_k<1`, and for every h>0.
For example any fixed integer m>=1024 makes rho_k uniformly less than1
over this range. It does not yet prove exponential suppression in the
oblique law: the exact remaining requirement is a bound on Delta_k.
Below we identify that quantity spectrally and show why PSD alone cannot
bound it as needed. This is a new surface-scale comparison, not an
assertion that an oblique quotient is reflection positive.

The note starts at `e17b286b`. It neither recomputes the axis theorem nor
uses new sampling, fixed-volume coupling expansions or a phase diagram.

## 1. Geometry and the finite transfer operator

For alpha=(4k,3k) and beta=(-3k,4k),

```text
4 alpha-3 beta=(25k,0),     -alpha+beta=(-7k,k).
```

Consequently `(x,y+k)` is identified with `(x+7k,y)`, while the horizontal
period is W=25k. Let P rotate a row configuration by 7k sites, with the
convention that the seam joins the top row at x to the bottom row at
x+7k. Then P^25=I and no smaller positive power is the identity.

Use the same local colour gas, Qc=m^2, and the same Gram edge matrix R
as in the [fixed-coupling proof](closed-source-fixed-coupling-peierls.md).
For a row sigma of W spins, put

```text
D_sigma=h^(K(sigma)/2) product_x R(sigma_x,sigma_(x+1))^(1/2),
V_(sigma,tau)=product_x R(sigma_x,tau_x),
T=D V D.
```

This is a finite real symmetric positive-semidefinite, entrywise
nonnegative matrix. Zero-weight horizontal rows are simply retained as
zero rows of T. Row translations commute with T. With a consistent
permutation-matrix convention for P, the two local partition functions
are exactly

```text
Z_twist=Tr(T^k P),           Z_rect=Tr(T^k).         (2)
```

Changing the convention replaces P by P^(-1) and leaves the trace the
same. Equation(2) is a spatial identification, not a spin-colour twist.
Both partition functions use h and m unchanged and have the same number
N of sites. Tr(T^k P)>0 follows also directly from the positive empty
configuration. Since T^k is PSD and P is unitary,
`Z_twist<=Z_rect`, proving the sign of Delta_k in (1).

These are **local colour** partition functions. The original occupation
law still contains m^(-r); it is applied in Section3, not silently
included in the unmarked row transfer matrix.

## 2. Cauchy--Schwarz removes the seam in a constrained numerator

Cut the H=k rows into two halves along opposite bond seams. For events
F and G determined entirely by the spins and edges inside the respective
halves, expand each crossing R in its Gram channels. Each half is a real
matrix A_F or A_G on the seam-channel indices. The oblique constrained
partition sum has the form

```text
Z_twist(F,G)=Tr(A_F A_G P_channels).
```

The last factor is a permutation of the seam channels and hence unitary.
Hilbert--Schmidt Cauchy--Schwarz gives

```text
Z_twist(F,G)
 <= [Tr(A_F A_F^T) Tr(A_G^T A_G)]^(1/2)
  = [Z_rect(F reflected F) Z_rect(G reflected G)]^(1/2).   (3)
```

Each trace glues a half to its reflected copy with no shift; both use
the same W x H rectangle. Site weights occur once per reflected site.
This step requires no reflection symmetry of the original oblique torus.
It uses its local Gram factorization and pays the normalization ratio

```text
P_twist(F,G)
 <= exp(Delta_k)
    [P_rect(F reflected F) P_rect(G reflected G)]^(1/2).   (4)
```

The rectangle has coordinate reflection symmetry. The existing domino
calculation extends directly to its two side lengths: the horizontal
ordered-mixed-domino norm is at most m^(-1+1/H), and the vertical one
at most m^(-1+1/W), uniformly in h. The disseminated stripe counts are
respectively C_B=W/4 and H/4; the colour multiplicities are retained.
Thus beta=m^(-1+1/k) bounds both norms. For any collection of b distinct
specified edges in that rectangle, at least b/4 belong to one of its
four domino tilings, so

```text
P_rect(all b specified edges mixed) <= (2 beta)^(b/4),
                                               2 beta<1.          (5)
```

## 3. A contour loses only a 2/k fraction of constraints at the cuts

Fix an edge-simple resolved cut contour of length n in the oblique
torus. Among the k choices of an opposite pair of horizontal cut seams,
each of its vertical primal edges is cut exactly twice. Some pair
therefore cuts at most `2n/k` contour edges. Translate the cylindrical
fundamental domain so this pair is used in (3); the partition ratio in
(4) is unchanged.

Drop only those crossing-edge requirements. The remaining mixed edges
lie inside the two halves, with counts b1,b2 satisfying
`b1+b2>=n(1-2/k)`. Reflection doubles each set to 2b_i distinct edges
of an axis rectangle. Equations(4)--(5) consequently yield

```text
P_twist(the specified contour is mixed)
 <= exp(Delta_k) (2 beta)^((b1+b2)/4)
 <= exp(Delta_k) (2 beta)^(n(1-2/k)/4).              (6)
```

No additional combinatorial factor for the cut location is necessary:
for each specified contour choose one minimizing pair deterministically.

A rank-one occupied configuration has an essential resolved boundary,
whose length is at least `ell1((4k+3ki)Z[i])=7k`. The same nonbacktracking
walk count `4N 3^(n-1)` as before applies locally. Summing (6) over n>=7k
gives the local-colour version of (1). Finally

```text
P_star,oblique(D)
 =E_twist[1_D m^(-r)]/E_twist[m^(-r)]
 <= m^2 P_twist(D)
```

proves (1) for the original source. For k>=16 one can use the uniform
constant `rho=3(2m^(-15/16))^(7/32)`; at m=1024 it is already less than1.
All estimates so far hold for arbitrary h, including an unknown moving
pooled-root activity. There is no assumption that h=1.

## 4. The remaining quantity is an explicit 25-sector cancellation

Let Pi_j be the spectral projection of P with eigenvalue exp(2 pi i j/25),
and define

```text
w_j=Tr(Pi_j T^k)/Tr(T^k),     j=0,...,24.
```

Commutation and positivity imply `w_j>=0`, `sum w_j=1`, and reality gives
`w_j=w_(25-j)`. The normalization penalty is exactly

```text
exp(-Delta_k)=sum_j w_j cos(2 pi j/25)
            =1-sum_j w_j[1-cos(2 pi j/25)] >0.       (7)
```

In particular `w_0>=1/2+epsilon` would imply
`Delta_k<=-log(2epsilon)`. More generally, a subexponential lower bound
on (7) suffices. Quantitatively, put
`tau_infty=-log[3(2/m)^(1/4)]>0`. Since rho_k tends to this limiting
rate's exponential, (1) is exponentially small whenever

```text
limsup_(k->infinity,4|k) Delta_k(h_k,m)/k < 7 tau_infty.  (8)
```

To make the conclusion uniform in h, the same inequality must hold after
taking the supremum over h. The factor in (7) is a translation-sector
partition ratio, not an observable/source fit or a field identification.

There is a simple exact but inadequate comparison. Force the row at the
twisted seam to be all vacant. The seam shift then has no effect, so

```text
Z_twist(row0)=Z_rect(row0),
Delta_k <= -log P_rect(row0)
        <= W log(1+h m^4).                          (9)
```

For the last step the local colour conditional probability of a vacant
site is at least `1/(1+h m^4)`: with no active neighbour the total active
odds are h/m^2; with one active neighbour colour they are at most h m^4;
with conflicting active neighbour colours they vanish. Iterated
conditioning proves the row bound. This costs only O(k), rather than
the O(k^2/m^2) volume term in a raw contour-subset relaxation. But near
h=1 its coefficient approaches `100 log m`, whereas the available
7 tau_infty approaches `(7/4)log m`. Equation(9) therefore cannot close (8).

Nor does PSD plus a fixed-order twist automatically fix that loss. To
see the algebraic obstruction, let D=25M+1 and let P fix one state and
cycle the other states in M cycles of length25. The positive PSD matrix
`T=I+D^(-2) 11^T` commutes with P, yet

```text
Tr(T^k P)/Tr(T^k)
 = (1+1/D)^k / [D-1+(1+1/D)^k].                    (10)
```

When D grows exponentially in k, this ratio can be exponentially small.
This is not a counterexample to the local colour model; it proves that
the Gram/commutation/order25 premises alone cannot remove Delta_k.

## Boundary of the result

The new exact product is (1), together with the fully specified transfer
ratio (7). It reduces the oblique fixed-t problem to a surface-scale
twist stability estimate and identifies an explicit threshold (8) for
success. No fixed-t oblique suppression theorem is claimed without it.
The bounded m^(-r) projection is retained throughout. Nothing here uses
literal pure-state concentration or duplicates the separate dilute
double-scaling calculation. The pooled-U within-geometry denominator
remains a subsequent issue even after (8) is established.
