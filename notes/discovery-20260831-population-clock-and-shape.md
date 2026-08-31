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

All required orientations in a pair were conditionally replaced together,
or both original values were kept. There were47/164 whole-pair fallbacks,
not discarded difficult prefixes. This is a completed archive-derived R1
contribution, not full F2/A_top or a runtime-speedup claim.

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
