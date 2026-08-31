# A near-zero lifetime integral does not remove the thermal moment response

The common-label perturbation now has a fixed-moment readout, independent
of selecting any particular peak or numerical crossing. Its orientation-even
E dipole is negative at both sizes, although its zeroth moment is weak.
The source-odd geometry-difference dipole is weaker at the present precision.

The scientific next distinction is exact: the first thermal moment is a
center-times-lifetime response. After the two mean changes are removed, its
remainder is half the difference of the endpoint variance responses. It is
not independent evidence about the joint birth-time copula. Actual joint
fluctuation requires the extra K1*K2 information from the same saved tails.

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

## Source and computation

`scripts/p334_euler_thermal_moments.py` consumes only the complete signed
histograms at4db356e1b026853468f94d59d938895a2367ceb7. Outputs at
9059776d866287e0cdb95e0a5f079843905cbeb9 are in
`results/p334-euler-thermal-moments/score.json`. The readout took about0.11s;
there were no new trajectories, root search, DP, package install or cloud
operation. All uncertainty uses the same original20 batches per size, with
the existing e32a8593/959a7fa2 dependency. A covariance-aware connected
decomposition is a separate downstream result, not assumed from this table.
