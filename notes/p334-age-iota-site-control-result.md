# P334 result: exact index and birth-site Smith class do not absorb birth age

## Result

The line-only production slope from `2e99533` replays exactly at the point
level; the largest replay covariance difference is `1.36e-20`.  The two new
fixed-effect layers then give a sharp negative mechanism result.

### Saturation indices are an exact no-op

Every rank-one risk-set path has

```text
iota01 = iota12 = 1.
```

This holds for all `1,824,619` N325 and `1,800,118` N425 survivors across both
orientations.  Consequently `(ell,iota01,iota12)` is exactly the same partition
as `ell`; every point and delete-one slope is identical.  This is production-
scale support for the integral-saturation theorem, not an estimated null index
effect.

### Exact site-pair Smith class changes essentially nothing

The translation-invariant site class

```text
g = gcd((site12-site01) mod N, N)
```

expands the fixed-effect partition from 8 line strata to 38--40 populated
line/site strata while retaining more than `99.997%` of the within-stratum age
information.

| archive | orientation | primary beta | site-controlled beta | magnitude retained | paired difference p |
|---|---|---:|---:|---:|---:|
| N325 | first | `-0.067253` | `-0.067287` | `100.051%` | `0.581` |
| N325 | second | `-0.081399` | `-0.081340` | `99.927%` | `0.333` |
| N425 | first | `-0.069438` | `-0.069316` | `99.825%` | `0.0364` |
| N425 | second | `-0.066918` | `-0.066841` | `99.885%` | `0.175` |

The two-orientation controlled associations remain overwhelming:

```text
N325: chi2=127.566/2, p=1.99e-28
N425: chi2=151.474/2, p=1.28e-33.
```

The paired control-induced changes are not resolved at the frozen
`alpha=0.01`:

```text
N325: p=0.540
N425: p=0.0375.
```

Thus neither the integral index nor the exact relative Smith class of the two
birth sites explains the age coefficient.

## Mechanism update

This rules out a concrete omitted-embedding explanation: the strong negative
age slope is not produced by mixing different saturation indices or cyclic
orders of the two birth-site displacement.  Together with the earlier
landing/H4 control, the recorded birth-local and topological embedding fields
now explain at most a few parts in a thousand of the slope magnitude.

The conclusion should still stop short of intrinsic temporal memory.  The
archive does not contain the full microscopic configuration at the observation
layer `k0`.  Age may remain a proxy for an unrecorded current cluster shape,
bottleneck, cycle-rank profile, or contact geometry.

There is also a strict temporal caveat: `site12/iota12` are future relative to
the k0 prediction.  Their failure to attenuate is a useful retrospective
negative control, but a future field cannot be promoted into a predictive
state or causal mediator.

The next most discriminating acquisition is therefore small and explicit:
record a compact current-geometry summary at `k0` itself.  Repeating more
birth-local marks or exact site classes is now low information.

## Provenance

No new simulation was run.  Both 2M raw hashes and line counts were verified.
The first two scorer calls stopped before emitting any score: identical iota
columns produced an algebraic zero difference variance that needed an explicit
no-op representation.  Both failure logs and the two narrow interface-fix
commits are retained.  The sole completed score uses `aebd329`.

The 200-row batch artifact contains one common-batch deletion across both
orientations per size and reconstructs both complete `6 x 6` covariances.

## Scientific card

- Mechanism space changed: saturation-index and cyclic site-pair embedding
  mixtures are removed as explanations of the age slope.
- Not proved: intrinsic temporal memory, causal history, or scaling-limit
  non-Markovianity; unrecorded current geometry remains live.
- Observer/source/geometry: next-step K2 hazard at frozen k0; N325/N425 paired
  orientations; line, exact iota pair and cyclic birth-site Smith class.
- Dependency groups: one raw P267 block per size; every control is a nested view
  of the same paths.
- Next lift: a compact geometry covariate measured at k0, not another future or
  birth-local classification.
