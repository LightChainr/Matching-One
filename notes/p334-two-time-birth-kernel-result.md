# P334 two-time homology-rank kernel

## Result

The two production birth archives determine the same-path two-time rank kernel
on seven frozen near-critical layers.  A single separable temporal amplitude is
not remotely adequate:

| archive | joint adjacent-minor chi-square / df | decision |
|---|---:|---|
| N325 | `5.2921e7 / 12` | reject rank-one separability |
| N425 | `5.0319e7 / 12` | reject rank-one separability |

The effect is structural rather than a marginal significance artifact.  The
six normalized adjacent principal-minor defects are between `0.874` and
`0.979` in every orientation; rank one would require all six to be zero.

The kernel nevertheless has a reproducible low-effective-rank shape.  The
first three eigenmodes carry `96.57%--96.63%` of the trace in all four
size/orientation rows.  Their fractions are approximately

```text
(0.633--0.638, 0.178--0.179, 0.150--0.154).
```

After one descriptive optimal amplitude rescaling, the N325 and N425
equal-orientation mean kernels differ by `3.67%` in Frobenius norm; the two
correlation matrices differ by `2.03%`.  The frozen exact-equality score still
rejects (`chi2=6395.39/28`) because two million paths per size resolve these
small finite-size/shape differences.  Thus the useful observation is a stable
three-mode hierarchy, not exact finite-size collapse or an exact rank-three
claim.

## Exact reconstruction gate

For every frozen pair `p<=q`, both orientations and both sizes satisfy

```text
P(K1<=p,K2<=q) = E[r(p)r(q)] - F1(p) - 2 F2(p)
```

to at most `4.45e-16`.  The summaries therefore retain a genuine same-path
joint birth kernel; they are not products of independent single-time samples.

## Scientific card

- Mechanism changed: one scalar temporal amplitude cannot represent the
  production ambient-homology birth process.
- Not proved: a CFT state count, Jordan structure, intrinsic temporal memory
  after conditioning on complete geometry, exact rank three, or universal
  kernel collapse.
- Observer/source/geometry: same-permutation `r(k)` process in the paired N325
  and N425 norm-five quotient orientations.
- Dependency groups: one paired-orientation archive per size; N325 and N425
  use disjoint seeds and counter domains.
- Next lift: project the newly acquired current-`k0` geometry covariates onto
  the second and third temporal eigenmodes and test which mode survives.

