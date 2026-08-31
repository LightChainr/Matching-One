# Where an Euler-invisible control changes the joint birth clock

**The dominant response is transport of prefix mean clocks, and most of that
transport is between the means of the fixed rank cells.** Rank-cell masses
do not change. A smaller, resolved response remains between prefixes *within
the same rank cell*, so the response cannot be reduced to rank-cell identity
alone. Changes in covariance inside a complete fixed prefix are smaller.

The exact-census estimator now uses all eligible score differences, removing
the old same-class pair mask. It targets the same derivative on the same
saved data. The original prefix law is unchanged by the common-next-label
policy; the following layers locate where its joint fluctuation response lives.

## Main numerical decomposition

All entries below are in units of **1e-7**, with original20-batch standard
errors. Products are computed in each physical geometry before S/D.

| N / source -> output | Total covariance response | Between rank-cell means | Between prefixes within a rank cell | Within a complete prefix |
|---|---:|---:|---:|---:|
| 325 plus -> S | +4.4578 +/- .3567 | +3.9947 +/- .2939 | +.3699 +/- .0935 | +.0932 +/- .0641 |
| 425 plus -> S | +4.2191 +/- .3398 | +3.7760 +/- .2693 | +.3192 +/- .0833 | +.1239 +/- .0822 |
| 325 minus -> D | -13.0906 +/- .9209 | -11.2591 +/- .7222 | -1.3344 +/- .2759 | -.4971 +/- .2313 |
| 425 minus -> D | -9.7331 +/- .5824 | -8.2790 +/- .5297 | -1.0643 +/- .2125 | -.3897 +/- .1887 |

The last three columns sum to the first, with their joint covariance retained.
Between-prefix transport accounts for about96--98% of the total point
response; between-rank-cell mean transport alone accounts for about85--90%.
These are **signed response contributions**, not fractions of baseline
variance explained or independent evidence counts.

The [same-batch coordinator](https://github.com/LightChainr/Matching-One/blob/ce20158a5928e55b67324cba7ed3a18a5c163b39/notes/p334-birth-covariance-hierarchy-joint.md)
reports between-rank-cell shares of89.612+/-2.176% and89.498+/-2.759%
for plus -> S, and86.009+/-2.458% and85.061+/-2.172% for minus -> D
(errors in percentage points). Old matched-mask versus new exact-score
total differences are at most1.95 paired SE; those same-source changes
are not a new physical effect.

The same-rank-cell prefix term is positive for plus -> S and negative for
minus -> D at both sizes, about3.8--5.0 SE from zero. If each conditional
mean response depended only on the rank pair G, this term would vanish.
Its survival therefore localizes a response dependence on information inside
the rank cell. It does not by itself identify which microscopic mark carries
that dependence, or show two response directions at each individual prefix.

The fixed-prefix covariance response in minus -> D is only about2.1 SE at
each size; plus -> S is weaker. Its split into suffix selection and label-mean
dispersion is:

| N / source -> output | Suffix selection | Label-mean dispersion |
|---|---:|---:|
| 325 plus -> S | +1.3710e-8 +/- 1.1148e-8 | -4.3895e-9 +/- 8.7960e-9 |
| 425 plus -> S | +7.9699e-9 +/- 1.1241e-8 | +4.4217e-9 +/- 7.5385e-9 |
| 325 minus -> D | -6.6384e-8 +/- 2.8568e-8 | +1.6670e-8 +/- 2.1030e-8 |
| 425 minus -> D | -4.1659e-8 +/- 1.9199e-8 | +2.6855e-9 +/- 1.2252e-8 |

The large total covariance effect should consequently not be described as
equally strong evidence for a change in the joint law inside a fixed prefix.

## Lifetime variation can cancel after rank-cell pooling

For plus -> S the total variance response of W=Y-X is weak, but this sum
contains opposing components:

| Variance-of-W response | N325 +/- SE | N425 +/- SE |
|---|---:|---:|
| within rank cells, including prefix and suffix variation | +5.5053e-8 +/- 2.8110e-8 | +7.6234e-8 +/- 2.6345e-8 |
| between the rank-cell mean lifetimes | -6.3045e-8 +/- 2.0183e-8 | -6.5041e-8 +/- 1.7952e-8 |
| total | -7.9924e-9 +/- 2.6770e-8 | +1.1193e-8 +/- 3.1780e-8 |

The policy changes the alignment of rank-cell mean lifetimes while their
masses stay fixed. Positive within-cell variation partly offsets that
compression. A globally weak lifetime-variance response therefore does not
justify assuming that lifetime is preserved path by path. Such preservation
would leave each of these fixed-cell lifetime distributions unchanged.
The within-cell positive result is about2.0 SE atN325 and2.9 SE atN425;
this is an exploratory mechanism readout, not a claim of a completed
microscopic transport identification.

## Exact hierarchy and estimators

Use X=K1/(N+1), Y=K2/(N+1). These are normalized birth ranks, without the
extra uniform-order-statistic clock. Z denotes the full original ordered
prefix, U the next label, and R the remaining suffix. With
mu_i(Z)=E[i|Z] and m_i(Z,U)=E[i|Z,U], the exact three-level identity is

```
Cov(X,Y)
 = E_Z E_U Cov_R(X,Y|Z,U)
 + E_Z Cov_U(m_X,m_Y|Z)
 + Cov_Z(mu_X,mu_Y).
```

The complete vacant-label census gives a centered score
`s(Z,U)=pi_a*(g(U)-mean_a(g))` in each joint-safe degree class a, zero outside.
Hence E_U[s|Z]=0 exactly. The policy changes neither the law of Z nor the
conditional suffix law at fixed (Z,U). The three derivative terms are

```
suffix selection = E_Z E_U[s Cov_R(X,Y|Z,U)]
label dispersion = E_Z E_U[s (m_X-mu_X)(m_Y-mu_Y)]
prefix transport = Cov_Z(H mu_X,mu_Y)+Cov_Z(mu_X,H mu_Y).
```

The first is selection among labels with different residual joint noise;
it is not a change of the suffix dynamics at a given label. The sum of the
first two is the average within-prefix covariance derivative. If this is
nonzero, a model that merely translates both birth distributions by constants
depending on Z is insufficient. Label- or suffix-dependent motion remains
possible, as does a change of conditional spread or correlation.

For G=(first_rank,second_rank), the final prefix term has the further
conditional-covariance split

```
prefix transport
 = E_G[Cov_Z(mu_X,h_Y|G)+Cov_Z(mu_Y,h_X|G)]
 + Cov_G(E_Z[mu_X|G],E_Z[h_Y|G])
 + Cov_G(E_Z[mu_Y|G],E_Z[h_X|G]).
```

The first line is the measured same-cell prefix residual; the last two
terms give transport among rank-cell means. The observed cell masses are
retained. Within-cell mean products use distinct-prefix U-products at the
cell's own retained prefix count, then return to the full population weight.

## The saved experiment already contains the necessary products

Each prefix has eight independent quartets, each comprising independent U,V
and two independent suffixes at each label. Set f_U to a two-suffix mean for
X,Y,X^2,XY,Y^2, and use

```
b_f=(f_U+f_V)/2,
h_f=(s_U-s_V)(f_U-f_V)/2.
```

The same-label cross-suffix product
`c_XY=(X_0 Y_1+X_1 Y_0)/2` has conditional expectation m_X*m_Y.
Thus h_XY-h_cXY estimates the suffix-selection term. For Q=8,

```
U_bh(X,Y)
 = sum_(q != r) [b_X,q h_Y,r+b_Y,q h_X,r] / [Q(Q-1)]
```

is unbiased for mu_X*Hmu_Y+mu_Y*Hmu_X at a fixed prefix. Then the average
label-dispersion term is `mean_Z(h_cXY-U_bh)`. Using the same quartet for
the product would mix conditional signal with Monte Carlo covariance.

For global products the scorer uses distinct *prefixes*, not products of
separately reported errors. For P sampled prefixes, the symmetric global
product is

```
[P*(mean(b_X)*mean(h_Y)+mean(b_Y)*mean(h_X))
 -mean(b_X*h_Y+b_Y*h_X)]/(P-1).
```

This is the unbiased global mean-product estimator under the original
prefix sampling law. All quantities and retained P are recomputed inside
each original-batch deletion. Cell contributions keep the full population
denominator and global centering; they are not separately renormalized
conditional effects. Compute each physical geometry first, then S/D.

The reader retains XX and YY as well, so the same decomposition supplies
center variance, lifetime variance and endpoint variance imbalance. No
new prefix, quartet, suffix, cloud job, model fit or matrix inverse is needed.

## Result lifecycle

The scorer definition is `f34bcd6f` at
`scripts/p334_birth_covariance_hierarchy.py`; its complete new result is
`44dc9e3396e39105cae85a29d04b39d0afc82d84`,
`results/p334-birth-covariance-hierarchy/score.json`. Rank-cell transport and
cell-centered lifetime readouts are at
`2bc3529468fbcba589182acaf98fa4855eb0a85e`,
`results/p334-rankcell-covariance-transport/score.json`.

Both consume the [once-extracted exact-score archive](https://github.com/LightChainr/Matching-One/tree/375cd3a12b2b7a87d79148a59f62b95898f9e471/results/p334-exact-score-quartet-moments),
with baseline/tangent8-moment blocks, integer birth clocks, exact score
numerators and the original prefix/batch IDs. The source pass took5.04s;
subsequent decompositions consume its arrays or pooled batch moments, not
fresh trajectories. The [nested-response identities](https://github.com/LightChainr/Matching-One/blob/03603388e6c0bee5889a64229a124d3f5e89790b/notes/p334-nested-covariance-response.md)
give the exact translation predictions and the distinction between
source-population and fixed-empirical-mixture U-products.

The old matched-mask result and new exact-score result share an estimand
and the original e32a8593/959a7fa2 block; they are not separate confirmations.
All new raw batch vectors and derived LOO are now appended to the common
factor at `ce20158a5928e55b67324cba7ed3a18a5c163b39`, under
`results/p334-birth-covariance-hierarchy-joint/`; it includes old/new paired
estimator contrasts and the response-share uncertainty. Results are branch deliveries, not an assertion of
integration into main. No other team's local-rank or finite-q_t calculation
was repeated, and no cloud machine was used.

The next mechanism question is which within-rank-cell prefix feature predicts
the conditional mean response. Rank-cell identity explains the dominant
organization but leaves a measurable residual; its population share need not
change to generate either part. Repository notes and Issue/PR handoffs carry
this result; routine cross-task status messages remain suppressed.
