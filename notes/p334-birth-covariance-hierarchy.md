# Where an Euler-invisible control changes the joint birth clock

The new question concerns the *location* of the joint fluctuation response.
The original prefix law is unchanged by the common-next-label policy. Is its
covariance response carried by differences between prefixes, by differences
between next labels inside one prefix, or by the residual uncertainty after
the next label is already fixed?

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
`scripts/p334_birth_covariance_hierarchy.py`. It consumes exact-score
prefix/quartet moments; the old matched-mask result is a same-estimand,
same-data precursor, not a separate confirmation. Quantitative results
will be attached here once the new source archive is committed.
