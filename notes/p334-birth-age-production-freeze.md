# P334/P337 birth-age production freeze

## Scope

This is a zero-new-Monte-Carlo, retrospective score of the two P267
two-observer production archives.  It does not repeat the exact `1/57`
witness, the six finite collision censuses, or the conditional six-arm
scaling argument in `fee33287`.

The two raw sparse tables contain the full `(K1,K2,ell)` law by batch and
orientation.  They therefore answer the next production question directly:

> Among paths that are rank one at the observed intrinsic-center layer, does
> the probability of exiting rank one on the next insertion still depend on
> when rank one was born?

## Locked layers and inputs

The intrinsic centers were estimated by the already completed P267 score and
are not re-estimated from the birth-age outcome:

| size | intrinsic center | `k0=round(N p0)` | complement layer |
|---:|---:|---:|---:|
| 325 | 0.5927647401191154 | 193 | 132 |
| 425 | 0.5927586490022030 | 252 | 173 |

Each size has 2,000,000 paths per orientation in 100 paired batches.  The two
orientations share their permutation stream and remain one covariance block.
The size streams have disjoint seeds and counter domains.

## One-degree birth-age statistic

At fixed size and orientation retain rank-one survivors

```text
K1 <= k0 < K2.
```

The next-step outcome and density-normalized age are

```text
y = 1[K2=k0+1],
x = (k0-K1)/N.
```

Within every primitive projective line `ell`, center `x` by its survivor-count
weighted mean.  With `n_lj` survivor counts and `y_lj` next-step exit counts,
the one-free-degree line-fixed-effect risk slope is

```text
beta_age
 = sum_lj (x_j-xbar_l)y_lj
   / sum_lj n_lj(x_j-xbar_l)^2.
```

Positive beta means older rank-one plateaux have larger next-step exit risk;
negative beta means the hazard is concentrated among younger plateaux.  The
null `beta_age=0` is independence of the next-step hazard from `K1` after
conditioning on the current layer and primitive line.

All site, landing and local-H4 mark fields are summed over.  They are not extra
strata or evidence rows.

## Collision mass and covariance

In the same pass compute

```text
D_N = P(K1=K2).
```

For each size, delete one common batch index from both orientations and
recompute the vector

```text
(beta_first, D_first, beta_first_complement,
 beta_second, D_second, beta_second_complement).
```

The 100 delete-one values supply one `6 x 6` jackknife covariance.  This keeps
the birth-age/collision covariance and the shared-randomness orientation
covariance.  Orientation slopes receive one-degree Student-t tests; a two-
degree joint slope-zero Wald test is reported once per size/archive.  The
frozen alpha is `0.01`.  `D_N` receives a 99% interval but no size-law fit.

## Complement view

Use the exact stored mapping

```text
K1c=N+1-K2,
K2c=N+1-K1,
k0c=N-k0,
ell_c=ell.
```

The identical age estimator is applied after this transform.  The complement
score is a paired view of the same paths, not independent confirmation and
not assumed equal to the forward slope.  The scorer must also verify that all
stored complement-audit failure sums are zero and that direct-birth mass is
unchanged by the mapping.

## Claim boundary

A resolved slope would establish production-scale finite-size predictive
memory beyond `(k,rank,ell)`.  It would not prove that memory survives a
near-critical scaling limit, identify a Jordan/CFT state, or imply that a
finite-dimensional spatial transfer description fails.  Collision mass is a
separate channel and is not allowed to explain away or subtract the slope.
