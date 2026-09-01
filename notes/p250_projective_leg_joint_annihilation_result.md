# P250 result: the radius-five R2 direction does not satisfy the joint null

## Regression gate

The mandatory identity-linear replay passed before any other candidate was
reported.  It reproduces the `a770ac9` shared-block rank-five calculation with
the same pivot and 38 resolved covariance modes.  The Schur residual is
identical, the maximum covariance difference is `3.71e-22`, and the finite
batch p-value differs by `1.39e-19`:

```text
reference p = 2.8401996780954676e-5
observed  p = 2.8401996780954537e-5
```

This establishes that the new scorer is a strict mapped extension of the old
rank calculation rather than another null-line direction comparison.

## Joint result

Every individual transformed hand remains compatible with rank at most five
at the frozen `alpha=0.01` gate.  Jointly imposing one candidate-mapped q is
more restrictive:

| map | radius-four joint p | decision |
|---|---:|---|
| identity linear | `2.8402e-5` | rejected |
| identity + conjugation | `1.8595e-4` | rejected |
| Alexander R0 + conjugation | `0.0041203` | rejected |
| Alexander R1 + conjugation | `0.0040736` | rejected |
| Alexander R2 + conjugation | `0.0058241` | rejected |
| Alexander R3 + conjugation | `0.0774864` | survives |

Among all sixteen frozen maps, only orientation-preserving R3 plus
conjugation (`p=0.0126626`) and Alexander R3 plus conjugation survive.  No map
is selected by its largest p-value; all scores are correlated views of the
same 80k/400-batch archive.

## New cross-round tension

The independent degree-five acquisition at `11130ae` left only Alexander R2
plus conjugation above its line-direction threshold (`p=0.01338`), with R3
immediately below (`p=0.009686`).  The stricter radius-four joint null instead
rejects R2 and retains R3.  The two frozen survivor sets therefore have empty
intersection.

This does not justify multiplying or combining the two p-values.  Although the
fresh 1.2M stream has an independent seed, the published radius-five map score
also reuses the old 80k moments when re-extracting its annihilator lines.  The
result is best read as a mechanism-level tension: direction-only extension and
joint annihilation do not currently identify one fixed parameter-free map.

The numerical margin also matters.  R2 is moderately below the old joint
threshold, while the radius-five R2/R3 distinction straddles `0.01` narrowly.
The result excludes the frozen degree-two joint R2 chart at the declared gate;
it does not exclude modulus-dependent intertwiners, higher-degree/context
states, a general `5+5` direct sum, or other physical sector maps.

## Next discriminator

Do not collect another shell merely to repeat the map vote.  The next useful
zero-sample object is an augmented joint operator that uses the old radius-four
Hankel rows and the already acquired independent degree-five shift rows in one
candidate-constrained rank calculation, with old and fresh jackknife influence
covariances kept separate.  Its primary question is whether any fixed map,
especially the R2/R3 pair, satisfies both annihilation and extension in the
same model.  Until then, retain the empty survivor intersection as an explicit
tension rather than promoting either convention to a physical intertwiner.
