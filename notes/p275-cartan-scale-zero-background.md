# P275 scale-zero Cartan background and the next independent readout

Status: post-reveal mechanism analysis of the frozen P275 nine-geometry block.
The frozen selector and its `selected=none` conclusion are not changed.

## What 8bef10b closes, and what is new here

Commit `8bef10b` already proves the rank-gate contact identity

```text
Cov(q,J_D) = E[J_S]/2 + (p-1/2-E[q]) E[J_D].
```

This note does not claim that proof again. The exact identity is used as a
regression constraint on the nine revealed geometries. The new work is (i) a
full-covariance decomposition of the scale-zero response into its thermal and
relative-source pieces and (ii) a frozen way to leave the closed rank gate.

## Full-covariance decomposition

After the frozen frame transport and birth-mass normalization, write

```text
Gamma_T = phase * E[J_S]/(2B),
Gamma_R = phase * (p-1/2-E[q]) E[J_D]/B,
Gamma   = Gamma_T + Gamma_R.
```

The script recomputes the finite matching root inside every delete-one batch.
It stores the complete 36 by 36 covariance of `(Gamma_T,Gamma_R)`, with all
three common-field moduli coupled within a size and different seeds independent
across sizes. The identity residual is at most `2.78e-16` over all delete-one
replicates.

The decomposition is decisive in magnitude. `Gamma_T` supplies more than
99.8% of `|Gamma|` in every geometry; the largest
`|Gamma_R|/|Gamma|` is `0.159%`. Thus the order-one profile is specifically
the thermal/source-contact term, not an unexplained additive constant.

The post-reveal GLS is deliberately unforgiving:

| descriptive model for Gamma | chi2 / dof |
|---|---:|
| constant by modulus | 9275.92 / 12 |
| constants + frozen-Q4-shaped N^-13/8 tail | 4184.73 / 10 |
| constants + free tail at each modulus | 74.745 / 6 |

Even the fully shape-free `N^-13/8` tail is rejected. Therefore the data do
not support parsing the same statistic into a fitted scale-zero constant plus
an independent Q4 remainder. The visually stable N170 values are useful shape
coordinates, but not established continuum constants at this precision.

## Background-annihilated future observable

The non-tautological escape is to change the observer, not rearrange the same
`q J_D` algebra. A complete priority path already has two essential-birth
times `K1,K2`. Freeze

```text
O_ext = (K1+K2-(N+1))/(2(N+1)).
```

It is complement-odd but is not the instantaneous rank gate `q_k`. Hence
`Cov(O_ext,J_D)` is allowed while the three-state contact identity does not
determine it from source means.

Thermal source mixing is removed at field level, rather than with batch-total
regression. Retain configuration-level complex products and cross-fit

```text
alpha  = Cov(J_D,J_S)/Var(J_S),
J_perp = J_D-alpha J_S,
Gamma_ext_perp = Cov(O_ext,J_perp)/B.
```

Within each outer delete-one replicate, train `alpha` on even retained batches
and score odd retained batches, swap folds, and average with fixed equal
weights. This makes each fitted source coefficient independent of the fold it
scores. The frozen first holdout is the next common cyclic size `N=250` at
`tau=i,2i,5i/2`; its six-real GLS tests the fixed Q4 modulus vector with four
residual degrees of freedom. It does not yet test a size exponent. Only a
passing nonzero modulus test authorizes a second size for the already-fixed
`N^-13/8` transfer.

The manifest is not production-authorized. Tiny exact complement/nonclosure
oracles and a nested cross-fit scorer test are required first.

## Scientific card

1. MECHANISM: the revealed scale-zero response is the thermal half of an exact rank-gate contact Ward identity.
2. EXACT BOUNDARY: contact subtraction of the same `q J_D` statistic is identically zero, not a hidden Q4 estimator.
3. DISCOVERY: full-covariance constant-plus-`N^-13/8` decompositions all fail, including a free tail at every modulus.
4. NEXT OBSERVABLE: centered essential-birth clock `O_ext` coupled to a cross-fitted field-level `J_D` residual.
5. HOLDOUT: one independent common-field N250 three-modulus block tests Q4 shape before any second-size scaling claim.
