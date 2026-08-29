# Same-parent conjugate norm-5 coalescence tomography

Status: **prospective exact-control design before the new C-node data are read**.

Related: #57, #180, #200, #205, PRs #199, #202 and #204.

## Why this is now the next high-information target

The N145->290 full-curve result changes the mechanism question. The frozen
finite-size center-slope correction predicts well, `P4[D']` transfers almost as
an eigenchannel, but the three-level matching-odd curve has a resolved shape
failure and `P4[S']` still drifts. The post-N290 two-generator model therefore
introduces a conjugation-even radial/Jordan direction and a tentative
conjugation-odd spin-4 quadrature.

The weak point is that the norm-5 evidence for the conjugation-odd direction is
confounded: the `2-i` branch was measured on the N65 lineage and the `2+i`
branch on the different N85 lineage. A difference between those residuals may
therefore be multiplier conjugation, parent geometry, parent size, or finite
quotient arithmetic.

There is an unusually cheap exact way to separate these possibilities before
spending the frozen norm-4 production budget: apply the missing conjugate
multiplier to the **same two parents**. In both lineages the two missing paths
coalesce, up to D4, to one new noncyclic child.

```text
N65 parents: 8+i, 7+4i
  observed 2-i -> A=(17,6), B=(18,1)
  missing  2+i -> C=(15,10) from both parents

N85 parents: 9+2i, 7+6i
  observed 2+i -> A=(16,13), B=(19,8)
  missing  2-i -> C=(20,5) from both parents
```

The new C quotients have Smith invariants `(5,65)` and `(5,85)`, while A and B
are cyclic. One new geometric node per size therefore tests conjugation and
quotient class at the same time.

## Bold conjecture: coalescence is a kernel test of the transfer representation

Suppose the leading reflection-even angular response at fixed N is a scalar
background plus one square harmonic,

```text
M(theta) = S_N + A_N cos(s theta).
```

Three same-N angles must then lie in a two-dimensional affine subspace. The C
node is not fitted: its value is fixed by A and B. For H4 the exact nulls are

```text
N325:  5 M_C - 11 M_A +  6 M_B = 0
N425: 20 M_C + 13 M_A - 33 M_B = 0.
```

These equations use neither `N^-13/8`, `p_c`, a root, a thermal metric nor a
parent amplitude. The coefficient sums vanish, so an arbitrary scalar H0
background cancels exactly. This makes the experiment a representation/kernel
test rather than another exponent fit.

The equivalent conjugate ratios are

```text
N325: (M_C-M_A)/(M_C-M_B) = 6/11
N425: (M_C-M_A)/(M_C-M_B) = 33/13.
```

The fixed H8 and H12 affine alternatives are frozen in
`predictions/norm5_conjugate_coalescence_20260829.yaml`; no harmonic is selected
from the new C values.

## Why this can beat norm-4 in information per CPU

PR #204 shows that a roughly three-sigma norm-4 q2/Jordan decision requires a
large balanced source extension plus N260/N340 production, with projected
incremental compute around 260k CPU seconds. That run is valuable, but its
interpretation depends on treating the noncyclic quotient as a clean radial
semigroup test.

The coalescence control attacks that assumption directly at N325/N425 while
also unconfounding the proposed conjugation-odd quadrature. A 10M common-field
stage at each size is small compared with the frozen norm-4 campaign and can
change what a later norm-4 failure means.

This is a priority argument, not a permission gate. Norm-4 may proceed in
parallel when compute is already allocated.

## Common-field construction

For each size run two arbitrary-period pairs with exactly the same seed and
replica-counter interval, always placing C first:

```text
N325: C=(15,10) vs A=(17,6)
      C=(15,10) vs B=(18,1)

N425: C=(20,5)  vs A=(16,13)
      C=(20,5)  vs B=(19,8).
```

The priority permutation is geometry-independent at fixed N. Consequently the
first-orientation C threshold stream should be exactly duplicated across the
two pair jobs. Check this before covariance construction. Then form the H4
linear residual inside synchronized batches. Across N325/N425 use independent
seeds.

The first stage is 10M pair replicas with 100 batches. It is a fixed-sample
prospective score, not a blinded variance-only artifact. If it is underpowered,
its covariance may set a later sample count, but the harmonic residuals remain
unchanged.

## Decision map

- **H4 passes at both sizes while H8/H12 are worse:** the simplest complex H4
  transfer survives same-parent conjugation and a quotient-class change. The
  noncyclic norm-4 experiment becomes much cleaner to interpret.
- **H4 fails but one frozen H8/H12 affine law passes:** reopen the angular
  representation before launching norm-10 tomography.
- **all fixed harmonics fail mainly through C:** prioritize quotient arithmetic,
  nonlocal topology, or a quotient-sensitive readout before adding continuum
  fields.
- **all models are simply noisy:** accumulate more samples with the same
  residuals. Do not add a fitted exponent or target-chosen harmonic.

## Claim boundary

A pass would not prove a unique H4/LCFT operator. It would establish something
more operational and immediately useful: the same-N angular transfer law is
stable under the first deliberately noncyclic conjugate-coalescence control.
A failure would be equally valuable because it would prevent a quotient effect
from being misnamed as q2/Jordan or higher-rank RG mixing.
