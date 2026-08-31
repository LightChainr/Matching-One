# Population clock loading, persistent shape leverage, and the T4 slow pole

This round uses the tools on complete archived populations and an already
completed independent N900 source. Its main new observation is that exact
suffix averaging removes about half the canonical R1-contrast noise but
less than one percent of its integrated noise. The source of fluctuations,
not merely the accuracy of a conditional clock, now becomes the question.

## 1. The population readout changes the conditional-noise story

[0d1e586d](https://github.com/LightChainr/Matching-One/commit/0d1e586dafbade5e7d1f9bfc598170d0c881e337)
contains all40 original batches,20000 shared counters per size, including
every checkpoint stratum in the denominator. The target is explicitly
`1_R1 * BinomTail(K2;N,p_ref)` and its p integral; zero outside R1 is part
of this target, not an assertion that missing rows are rank two.

| H4-normalized R1 orientation contribution | N325 | N425 |
|---|---:|---:|
| Canonical mean +/- batch SE | .00060998 +/- .00156001 | .00113588 +/- .00118417 |
| Integrated mean +/- batch SE | -.00058818 +/- .00315435 | .00349118 +/- .00220881 |
| Estimated canonical suffix-noise fraction removed | 49.15% | 50.03% |
| Estimated integrated suffix-noise fraction removed | 0.816% | 0.681% |

The mean differences are not yet resolved sharply. The noise allocation is
nevertheless informative: the earlier selected147-prefix mixture cannot
stand in for a full-population contribution. Its83.95% integrated noise
fraction described a different mixture.

For fixed k0 and d=N-k0, the integrated conditional mean within R1 is

```
m_i = (d+1-E[T_i | R1])/(N+1).
```

The stratum indicator contributes an approximately order-one step, while
the waiting-time fluctuation is divided byN+1. The paired prevalence/clock
decomposition and four-state risk variance therefore directly test the
natural explanation for the very small integrated noise reduction.

The completed [paired-risk decomposition](https://github.com/LightChainr/Matching-One/blob/3d760b86/notes/p334-r1-prevalence-clock-loading-result.md)
now explains **99.84895% and
99.86380%** of the hybrid integrated-contrast variance by the four states
`(R_first,R_second)`. For the canonical contrast the corresponding fractions
are84.2390% and83.0724%. The within-R1 clock is therefore not the dominant
remaining source of noise in this gated integrated observable.

This suggests a specific next bridge to the global observable: retain all
checkpoint strata jointly. Their contributions can cancel both mean and
indicator noise. The present result does not assign99.85% of full A_top's
variance to rank-one prevalence.

Mean loading has a different decomposition: `D=C+L`, with C the symmetric
prevalence term and L the conditional-clock term. N325 has canonical
`L=.00096906 +/- .00039612`, against `C=-.00035908 +/- .00150262`;
N425 has `L=-.00055355 +/- .00043189`, against
`C=.00168943 +/- .00102780`. The mean clock contribution can oppose the
prevalence contribution even when prevalence dominates the variance. The
two terms remain in one covariance block.

All required orientations in a pair were conditionally replaced together,
or both original values were kept. There were47/164 whole-pair fallbacks,
not discarded difficult prefixes. This is a completed archive-derived R1
contribution, not full F2/A_top or a runtime-speedup claim.

### The micro-source partition is now available on the same population

[32270fa2](https://github.com/LightChainr/Matching-One/commit/32270fa2f8c5dfb19bf534b364fde26e2ac117f6)
uses monotonicity to obtain the original-singleton-gate birth law directly:

```
P(T=j, final site in original H2 gates) = S(j-1)*H2/(d-j+1).
```

Subtracting this from the full birth mass gives collective completion;
it includes sites that become singleton triggers only after subsequent
safe insertions. All unsolved whole pairs retain their original value in
an explicit unclassified channel. The three channels add back to Y.

Direct gates carry about82% of classified positive canonical loading and
66% of integrated loading. Yet the smaller collective source can determine
an orientation difference: N325's integrated H4 contributions are
`direct +.00101609`, `collective -.00156328`, and
`unclassified -.00004098`. N425's classified integrated contributions
are both positive, `+.00167889` and `+.00172133`. These are dependent
archive point estimates, not established population signs. Source signs
must be read with common batch covariance and the unresolved-source
allocation envelopes; a small net unclassified contrast need not imply
small uncertainty in its individual-source allocation.

There is also a new [exact source-competition identity](https://github.com/LightChainr/Matching-One/blob/e41b8e0014e747854000b2512ed36736ab8a98ef/notes/p334-integrated-source-competition-identity.md).
For each solved prefix, h=H2 and mu=E[T] give
`I_direct=h*mu/(N+1)` and
`I_collective=[d+1-(h+1)*mu]/(N+1)`.
Equivalently, the collective integral is
`(h+1)*E[T_original_gates_only-T]/(N+1)` on the same permutation.
Collective completion can therefore raise total loading while pre-empting
and reducing direct-source loading. Its integrated source information beyond
prevalence and mean waiting time is exactly the H2-weighted first moment;
the fixed-p readout, in contrast, still probes the full clock shape.

The [shared source-by-C/L crosswalk](https://github.com/LightChainr/Matching-One/blob/6133b39d20a198a61857dd30ddcf676e2b0a1a65/notes/p334-source-loading-crosswalk.md)
locates N325's canonical conditional-clock
term: direct `+.00114380 +/- .00050986`, collective
`-.00020310 +/- .00018303`, and unclassified
`+.00002836 +/- .00006657`. For its integrated L, direct
`+.00152158 +/- .00074444` and collective
`-.00131605 +/- .00070585` nearly cancel, leaving a much smaller net clock
contribution. This is exactly why separate source error bars must not be
combined as independent evidence: all45 derived coordinates share the same
20 original batches and a rank-at-most19 covariance matrix.

## 2. A thinning shoulder retains shape leverage

[3bacf19a](https://github.com/LightChainr/Matching-One/commit/3bacf19a)
turns the maximum-Gaussian three-center description into an affine- and
common-Gaussian-blur-invariant polynomial gap to two centers. Delta2 is
.20968 +/-.01044, .17327 +/-.02874 and .17272 +/-.04112 atN100/400/900.
The early weight falls18.06% ->6.54% ->3.20%, but its relative leverage
rises1.16 ->2.65 ->5.40. Small mass does not imply loss of shape influence.

The separate broader-class calculation
[b6db7ba5](https://github.com/LightChainr/Matching-One/commit/b6db7ba57c3c5bcb6e25558b5274f08aeef1ce63)
applies the original common-positive-symmetric-kernel two-translation
condition directly to N900. It requires kernel sixth moment
**-2.142812 +/-.481800**, and moment determinant **-.260440 +/-.059149**.
Both violate positivity at the empirical moment estimate. All800 aligned
delete-one rows retain the unique admissible reconstruction branch; both
margins stay negative. These overlapping rows are not800 independent tests.

The Gaussian-center gap alone would not rule out arbitrary symmetric
kernels. The second calculation is essential and uses the same N900 source,
not independent corroboration. Unequal kernels, asymmetric lobes, additional
components and signed cancellation remain different mechanistic directions.
The propagated SEs are exploratory, not a calibrated boundary certificate.

## 3. T4 repairs the tail by moving the slow pole

[1f19fc1a](https://github.com/LightChainr/Matching-One/commit/1f19fc1a2d9fc59dce650e95268c716762725985)
uses the saved7-to8-dimensional Schur blocks. The plus slow mass increases
1.931369756 ->1.947928395; its residue also increases
.461292288 ->.466944056. At t=4 the exact endpoint decomposition is:

```
pole movement       -6.967914 percentage points
residue restoration +1.281542 percentage points
other modes          +.000753 percentage points
net tail repair     -5.685620 percentage points
```

The residue increase partly opposes the repair. T4 is a fast feedback
coordinate that changes slow propagation, not an added slow component or
a removal of source weight. This extends the already identified fourth-order
Schur bridge without recomputing the full1430-state generator.

## 4. Clock slope and spatial noise sensitivity separate on real prefixes

[795908fb](https://github.com/LightChainr/Matching-One/commit/795908fbc9a781a0cda704864c237deaf0327f37)
evaluates the two original realN425 witnesses at their own S(u)=1/2.
The first positive noise energy factorizes as

```
E1 = u(1-u) * (sum I_v)^2 * [sum I_v^2/(sum I_v)^2].
```

WitnessB has a smaller clock prefactor (ratio .72854), but more concentrated
pivotal weights (ratio1.79257), giving stronger noise sensitivity
(E1 ratio1.30595). Effective pivotal counts are37.43 versus20.88, despite
83 versus127 sites with positive influence. This is a concrete microscopic
separation of clock rate from spatial sensitivity, not a population-frequency
claim. All site polynomials were already available; no new DP or MC was run.

## Lifecycle and dependencies

| Result | Observer / sector / geometry | Shared dependency group | Lifecycle |
|---|---|---|---|
| Paired population loading | R1-weighted F2; canonical/integral; N325/N425 orientation pairs | original e81dd59 counter blocks; includes the old selected witnesses | whole archived pool -> exact/hybrid means -> original20-batch covariance |
| Shoulder leverage and two-lobe obstruction | signed odd D_A rank-profile; 4i versus2i; N100/400/900 | existing shape archives; N90032M/800 shared batches5f30397c | completed source -> saved moment coordinates -> aligned LOO |
| T4 pole mechanism | charged physical source rays; fixed width-eight cylinder | same deterministic generator/Schur blocks074a5f53 | saved blocks -> endpoint poles/residues -> additive tail split |
| Real-prefix noise slope | local birth event at equal survival; two fixed N425 prefixes | original marked-site laws1c06230b and clocks6358ba49 | saved exact polynomials -> fixed S=1/2 -> concentration/clock factorization |

No new GPU or Huawei workload was launched. ZyTrST and TV2N0X are recorded
locally as user-reported Ready/off CPU capacity. No task was closed, merged
or reprioritized by this scientific delivery.
