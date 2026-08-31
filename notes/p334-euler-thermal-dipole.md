# Euler-invisible controls reshape the joint birth clock

The subsequent [nested covariance decomposition](p334-birth-covariance-hierarchy.md)
now locates this response: prefix mean transport dominates, mostly between
fixed rank-cell means, with a resolved within-rank-cell residual. Its complete
census score improves the estimator on the same data; the earlier numbers
below retain their original matched-mask provenance.

The common-label perturbation changes the intrinsic covariance of the two
birth ranks at both sizes, even after the shared uniform-order-statistic
timing contribution is removed. The common source increases the even
covariance response; the difference source has a negative geometry-contrast
response. A weak lifetime integral therefore misses a substantial part of
this perturbation's effect on the birth process.

Two complementary readouts explain the shape. The normalized rank-one
occupation moves later and broadens under the common source. In N425's
source-minus/geometry-difference channel, center displacement and unequal
endpoint broadening partly cancel, leaving a weak net thermal dipole.
All are projections of the same saved source block and its twenty batches,
not independent replications.

## The new joint direction survives the timing correction

For continuous birth times tau1<=tau2, mu1=E(tau1), and m12=E(tau1*tau2),

```
Cov(tau1,tau2) = Cov(K1,K2)/(N+1)^2 + (mu1-m12)/(N+1).
H_intrinsic = H[Cov(tau1,tau2)] - (H[mu1]-H[m12])/(N+1).
```

The second term integrates the common order-statistic clock exactly. The
intrinsic response, still expressed in continuous-clock units, is:

| Source/output | N325 +/- original-batch SE | N425 +/- original-batch SE |
|---|---:|---:|
| plus -> S intrinsic rank covariance | +4.59145e-7 +/- 7.51861e-8 | +5.34794e-7 +/- 6.33625e-8 |
| minus -> D intrinsic rank covariance | -1.19492e-6 +/- 3.25239e-7 | -1.08880e-6 +/- 1.61897e-7 |
| plus -> S continuous Pearson correlation | +6.71241e-5 +/- 2.42816e-5 | +1.08497e-4 +/- 2.04602e-5 |
| minus -> D continuous Pearson correlation | -2.22375e-4 +/- 7.95354e-5 | -2.38442e-4 +/- 5.62496e-5 |

At N425 minus -> D, the total covariance derivative is
`-1.08148e-6 +/- 1.60878e-7`; the timing term is only
`+7.32127e-9 +/- 3.03203e-9`, with the opposite sign. This response is not
created by uniform-priority smoothing.

These are derivatives of the covariance **between the first and second
birth within each geometry**, followed by S/D. They are distinct from the
previously weak cross-geometry source-to-observer mean susceptibilities.
A negative D derivative does not mean negative unconditional covariance,
and a Pearson coordinate does not identify the complete copula or a field.
The [joint coordinator's result](https://github.com/LightChainr/Matching-One/blob/e2ef9983f426890a299f5a6e1a2eba8b6d072855/notes/p334-euler-dipole-connected-clock.md)
retains these moments, the preceding full factor, and the new plateau LOO
columns in one shared covariance representation.

## Fixed readouts, no small-integral denominator

Let H_E(p) be the derivative under the previously defined common-label
Euler/rank-preserving policy, and let p0=.59274605079. Read

```
I_m = integral_0^1 p^m H_E(p) dp,
D1 = I1-p0 I0,
D2 = I2-2p0 I1+p0^2 I0.
```

Every coordinate is a linear functional of the same signed integer birth
histograms. No response centroid is defined by dividing I1 by a weak I0.

| Source/output | N325: I0 +/- SE | N425: I0 +/- SE |
|---|---:|---:|
| plus -> S | -1.45945e-6 +/- 1.24142e-6 | -1.53866e-6 +/- 1.19255e-6 |
| minus -> D | +5.68783e-6 +/- 3.98966e-6 | -5.25797e-7 +/- 2.21120e-6 |

| Source/output | N325: D1 +/- SE | N425: D1 +/- SE |
|---|---:|---:|
| plus -> S | -2.04556e-7 +/- 8.26059e-8 | -2.00956e-7 +/- 7.01655e-8 |
| minus -> D | +4.88315e-7 +/- 2.76289e-7 | +1.91398e-7 +/- 1.17965e-7 |

The plus -> S second centered moments are
`-1.81862e-8 +/- 8.03261e-9` and `-1.83638e-8 +/- 6.03633e-9`.
These describe the same thermal profile, not new independent confirmations.
The full score also retains plus -> D, minus -> S, A, the five original
rank-cell contributions and all twenty aligned batch vectors.

## Why these are birth center/lifetime observables

Conditioned on the permutation birth indices K1<=K2, their continuous
priority values tau1,tau2 are uniform order statistics. With rising factorials,

```
E[tau_j^q | K_j] = (K_j)_(q)/(N+1)_(q).
```

For C=(tau1+tau2)/2 and W=tau2-tau1, exact integration gives

```
I0 = -H[W],
I1 = -H[C W],
I2 = -H[C^2 W + W^3/12].
```

Thus, at each individual orientation's baseline means muC,muW,

```
D1 = -muW H[C] -(muC-p0) H[W] - H[Cov(C,W)].
```

The nonlinear centering must be performed orientation by orientation before
forming S=(first+second)/2 and D=(first-second)/delta_cos4. Multiplying paired
contrast means would not give the same estimand.

The identity `Cov(C,W)=(Var(tau2)-Var(tau1))/2` shows what that connected
term measures: unequal broadening of the two endpoints. In contrast,
`Cov(tau1,tau2)=Var(C)-Var(W)/4` needs their joint second moment. This is the
reason for the additional joint-moment reader rather than reinterpreting
marginal histograms as new pair-dependence information.

At N425 minus -> D the dipole decomposes as follows:

| Contribution | Response +/- original-batch SE |
|---|---:|
| center displacement, -muW H[C] | +5.42584e-7 +/- 7.55604e-8 |
| mean lifetime, -(muC-p0) H[W] | -1.64931e-9 +/- 1.88437e-9 |
| endpoint spread imbalance, -H[Cov(C,W)] | -3.49536e-7 +/- 9.69285e-8 |
| total dipole | +1.91398e-7 +/- 1.17965e-7 |

The small total conceals a resolved cancellation. At N325 the center term
is positive too, but its opposing spread term is not resolved. Errors for
the sum are propagated jointly from the original batches. N425 minus -> D
also has `H[Var(C)]=-1.05259e-6 +/- 1.42045e-7`, while `H[Var(W)]` remains
weak. This rejects a description consisting only of a rigid constant
translation; prefix-dependent shifts can still change these moments.

## The lifetime-normalized plateau also broadens

The saved continuous moments now supply a positive, unperturbed denominator
muW, about .044--.048, for the rank-one occupation measure. Normalize
`P(tau1<=p<tau2) dp` by muW and define its centroid eta and variance V:

```
eta=E(CW)/muW,
V=E(C^2 W+W^3/12)/muW-eta^2.
```

This is not the ratio of two small signed responses. Differentiate the
positive-baseline ratios separately in each orientation before S/D. For the
common plus -> S source the results are:

| Plateau response | N325 +/- batch SE | N425 +/- batch SE |
|---|---:|---:|
| centroid H_eta | +4.34757e-6 +/- 1.77826e-6 | +4.67072e-6 +/- 1.66208e-6 |
| variance H_V | +3.17944e-7 +/- 1.28394e-7 | +3.62163e-7 +/- 1.04955e-7 |
| width H_sqrt(V) | +3.09664e-6 +/- 1.25383e-6 | +3.88703e-6 +/- 1.12769e-6 |

The common perturbation moves the normalized rank-one occupation later and
broadens it. A rigid translation would leave V unchanged; a broader
path-dependent center motion is not excluded. V includes both the spread
of lifetime-weighted centers and the widths of individual plateaux, so it
is not identified with either component alone.

N425 minus -> D has variance response
`-4.41008e-7 +/- 1.69308e-7`. Its centroid response is
`-4.35046e-6 +/- 2.80147e-6`, while the difference between the
lifetime-weighted and ordinary-center responses is
`+8.08586e-6 +/- 2.25570e-6`. The latter is resolved even though the
centroid derivative itself is weak. This is a difference of complete
within-orientation derivatives, not the variance of a paired difference
or a claim that both individual geometries narrow.

The [exact hierarchy and shape maps](https://github.com/LightChainr/Matching-One/blob/06abeeaefe7063365fee36c6399dbab0d1c06b9b/notes/p334-thermal-dipole-clock-moments.md)
give all finite-N corrections and the rigid-translation conditions.
`scripts/p334_plateau_shape_tangent.py` produced these new ratio derivatives
in about .08s from f4682eb3, with complete20-batch LOO at
`results/p334-plateau-shape-tangent/score.json` (1e8549b5). No curve was
re-fitted or re-selected.

## Source and computation

`scripts/p334_euler_thermal_moments.py` consumes only the complete signed
histograms at4db356e1b026853468f94d59d938895a2367ceb7. Outputs at
9059776d866287e0cdb95e0a5f079843905cbeb9 are in
`results/p334-euler-thermal-moments/score.json`. The readout took about0.11s;
there were no new trajectories, root search, DP, package install or cloud
operation. All uncertainty uses the same original20 batches per size, with
the existing e32a8593/959a7fa2 dependency. The connected decomposition is now
delivered at e2ef9983, using all supplied original-batch means and re-forming
products after each batch deletion. It adds no independent sampling block.

The [new joint source](https://github.com/LightChainr/Matching-One/blob/d179bdf6fdef4b29f7d43f3c1a60d842fb35faf5/notes/p334-continuous-center-lifetime-moments.md)
reads all old forks once and stores baseline/tangent moments. Its `all`
baseline includes all nine rank cells, whereas the five R0-containing
cells exhaust only the nonzero tangent. Continuous-order-statistic
variance is integrated analytically rather than sampled again.

## Scientific handoff

- Changed mechanism space: instantaneous Euler/rank invisibility coexists
  with a resolved response of joint birth fluctuations and plateau shape.
  Weak integral/dipole means do not imply weak distributional response.
- Scope: finite N325/N425 common-next-label policy; plus/minus sources;
  per-geometry birth moments before S/D; original e32a8593/959a7fa2 block.
- Lifecycle: scripts, raw moment archive, exact identities, derived scores
  and shared covariance are pushed on their named analysis/theory branches.
  This synthesis is branch-only, not a main-branch integration claim.
- Next discriminant: a conditional or geometry-resolved response that
  separates prefix-dependent clock displacement from a change in the
  conditional joint birth law. The other team's finite-q_t and local-prefix
  rank work is complementary and should not be repeated here.
- Coordination: consolidate scientific results and necessary handoffs in
  repository notes, Issues and PRs. Routine partial-result/status messages
  between the three teams are no longer needed.
