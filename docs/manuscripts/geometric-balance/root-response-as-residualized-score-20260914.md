# Root response as an exact residualized-score quotient

Date: 2026-09-14

Status: exact finite probability calculus for a normalized source family.  It reorganizes #802 after the source-normalization erratum and gives a gauge-invariant definition of the nonthermal source direction.  No continuum field identification is assumed.

## 1. Setup

Let `P_{p,g}` be a **normalized** finite-torus probability law with baseline `g=0` equal to the Bernoulli site model.  Let

```text
X(omega)=r(omega)-1 in {-1,0,+1},
M(p,g)=E_{p,g}[X].
```

The finite matching root `p_*(g)` is defined by

```text
M(p_*(g),g)=0.
```

Write the normalized likelihood scores at `g=0` as

```text
S_p = partial_p log P_{p,0}(omega),
H   = partial_g log P_{p,g}(omega)|_(g=0).
```

Normalization implies

```text
E S_p=0,
E H=0.
```

For iid Bernoulli sites on `N` vertices,

```text
S_p=(K-Np)/[p(1-p)].
```

## 2. Root tangent is an exact covariance ratio

Differentiate `M`:

```text
M_p = Cov(X,S_p),
M_g = Cov(X,H).
```

At the root `E X=0`, these are simply `E[X S_p]` and `E[XH]`.

Implicit differentiation gives

```text
boxed:
T_g := dp_*/dg
     = - Cov(X,H) / Cov(X,S_p).                       (2.1)
```

The denominator is positive for the monotone rank observable in the interior regime.

Equation (2.1) is the normalized finite-torus analogue of the safe-Perron Clapeyron formula.  It makes no reference to a continuum field or a chosen metric on source space.

## 3. Residualize the source against the thermal nuisance

Define

```text
beta_H
 = Cov(X,H)/Cov(X,S_p)
 = -T_g.
```

Now set

```text
boxed:
H_perp = H - beta_H S_p.                              (3.1)
```

Then exactly

```text
Cov(X,H_perp)=0.                                      (3.2)
```

So `H_perp` is the source score after subtracting precisely the amount of thermal motion required to keep the root fixed to first order.

This is the intrinsic content of the “fixed-b normal direction” in #802.  No Euclidean/Fisher orthogonality language is needed: it is a nuisance-score quotient defined by the root estimating equation itself.

## 4. Shape response of any observable is a partial covariance

Let `A(omega)` be any finite observable with no explicit `p,g` dependence.  Along the moving root,

```text
d/dg E_{p_*(g),g}[A]|_0
 = Cov(A,H) + T_g Cov(A,S_p)
 = Cov(A,H_perp).                                      (4.1)
```

Thus the pair

```text
T_g = -beta_H,
N_A(g)=Cov(A,H_perp)
```

is an exact tangent/shape decomposition.

For a smooth functional of the rank law rather than a single `A`, replace `A` by its finite influence function; the same residual-score formula holds.

This is the probability-theory form of the source quotient already used in #773, now tied directly to the root.

## 5. Gauge invariance explains the #802 erratum

Suppose a second source score differs by

```text
H' = H + c S_p.
```

Then

```text
beta_H' = beta_H + c,
T_g'     = T_g - c,
```

but

```text
H'_perp
 = H + c S_p - (beta_H+c) S_p
 = H_perp.                                             (5.1)
```

Hence every fixed-root/shape response is identical.

A constant added to an **unnormalized** source weight disappears after centering the likelihood score, so the full gauge class is

```text
H ~ H + constant + c * thermal_score.                 (5.2)
```

This is exactly why two raw source formulas can give different root tangents but no independent shape information.

## 6. Black-pair / white-pair alias becomes one line

For a periodic row, the exact source identity is

```text
H_W = N - 2K + H_B.
```

After normalization/centering,

```text
H_W^c - H_B^c
 = -2(K-Np)
 = -2 p(1-p) S_p.                                     (6.1)
```

Therefore the two sources satisfy

```text
H_W,perp = H_B,perp                                    (6.2)
```

configurationwise as source classes.

Their root tangents differ only by the thermal gauge term:

```text
T_W-T_B=2p(1-p),                                      (6.3)
```

which is precisely the independent #802 audit result.

The earlier apparent “two-source discrimination” was therefore impossible in principle after thermal profiling.

## 7. Relation to Bernoulli-chaos grading

The centered score quotient should be performed **before** interpreting the source by chaos degree.

For example, black and white pair counts have the same degree-2 Hoeffding/Walsh component and differ only by degree 0/1 pieces.  Equation (6.2) is the score-space expression of that exact fact.

A genuinely new local source direction must have

```text
H_perp != 0
```

and be linearly independent of previously tested residualized scores.

Odd chaos degree three is a natural first exact complement-odd nonthermal candidate, but its continuum RG image remains an empirical/theoretical question.

## 8. Rank-law normal coordinate as an influence function

For the canonical rank probabilities `P0,P1,P2`, define

```text
b = (1/2) log(P0/P2),
d = log[P1/sqrt(P0 P2)].
```

At any interior point their source derivatives are linear functionals of the score `H`.  In particular one can write

```text
b_g = Cov(A_b,H),
d_g = Cov(A_d,H)
```

for rank-measurable influence functions `A_b,A_d`.

Then the fixed-b shape response

```text
N_g=d_g-(d_p/b_p)b_g
```

is exactly

```text
boxed:
N_g = Cov(A_d,H_perp_b),                               (8.1)
```

where

```text
H_perp_b=H-[b_g/b_p] S_p.
```

This recovers #802's chain-rule formula but exposes its statistical structure: it is a partial covariance after profiling the thermal nuisance.

## 9. Safe-transfer version

For a semi-infinite safe transfer, the normalized physical one-row source score should be used.  If an implementation starts from unnormalized row weights `W_tilde`, with

```text
f(g)=log Z_row(g),
I(g)=f(g)-log lambda_tilde(g),
```

the physical derivative is

```text
I_g=f_g-E_Q[H_raw].
```

One should first convert all source columns to physical normalized free-energy derivatives, then perform the same thermal residualization.

This prevents the omitted-normalizer failure found in the original #802 delivery.

## 10. A new source experiment with real information gain

Instead of comparing two motif sources before quotienting them, pre-register a small normalized score basis:

```text
S_p : thermal control;
H_2 : one genuine degree-2 residualized source;
H_3 : one local degree-3 residualized source with exact complement sign;
```

and compute the matrix

```text
Cov(X,H_a),
Cov(A_shape,H_a,perp),
```

or the safe-phase Perron analogues.

The first row tells which sources move the root; the second tells which change shape after thermal motion is removed.

This is a direct lattice realization of the #61 RG-tangent programme without prematurely naming local CFT parity.

## 11. Relation to original-U / Rao--Blackwell work

The same algebra is the standard nuisance projection of an estimating equation.  For an original-U influence function, replacing raw source content by its residual after projection onto the root/thermal nuisance is exactly the kind of object that conditional integration should target.

Thus the current topological source programme and #578's Rao--Blackwell programme share a mathematical quotient, even though their physical observables are different.

This does not change #275's frozen source contract.

## 12. Claim boundary

Exact finite probability identities:

```text
T_g=-Cov(X,H)/Cov(X,S_p),
H_perp=H-beta_H S_p,
Cov(X,H_perp)=0,
d/dg E[A] along root = Cov(A,H_perp),
source gauge invariance under H->H+c S_p.
```

Programme:

- use residualized microscopic source basis to estimate an RG tangent map;
- connect particular score blocks to continuum fields only after angular/radial and representation information is available.

The central correction is methodological: **normalize first, quotient the thermal nuisance second, and only then interpret a source physically.**