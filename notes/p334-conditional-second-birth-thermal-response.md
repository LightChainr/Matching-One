# The complete birth clock now gives an exact conditional thermal readout

The full physical T distributions at
[`6358ba49`](https://github.com/LightChainr/Matching-One/commit/6358ba49ef390c10a3f501b589ba7ba1d4e05b09)
already determine the entire canonical second-birth profile for each real
N425 prefix. No continuation simulation or reliability-DP rerun is needed.
They also determine exactly which **suffix noise** conditional averaging
removes from a production binomial-tail estimator.

## The estimand is a binomial tail, not a Bernoulli event

For a fixed prefix X, let `q_X(t)=Pr(T=t|X)`, obtained from the exact
survival coefficients as `q_X(t)=S_X(t-1)-S_X(t)`. Here N=425, k0=252,
there are 173 uninserted sites, and `K2=k0+T`. Define

\[
g_t(p)=\Pr[\operatorname{Bin}(425,p)\ge252+t],\qquad
f_{2,X}(p)=\sum_t q_X(t)g_t(p).
\]

An original fresh-permutation continuation produces the **already
canonicalized real number** g_T(p). Its conditional covariance is

\[
V_X(p,q)=\sum_tq_X(t)g_t(p)g_t(q)-f_{2,X}(p)f_{2,X}(q).
\tag{1}
\]

In particular its variance is `E g^2-(E g)^2`, not `f2*(1-f2)`. Exact
conditional averaging replaces g_T by f2_X and removes (1) entirely for
that fixed prefix. The returned profile can be reused for every p.

The integrated readout also has a closed form:

\[
\int_0^1g_t(p)\,dp=\frac{426-252-t}{426},\quad
E[\text{integrated }g\mid X]=\frac{174-E[T\mid X]}{426},\quad
\operatorname{Var}(\text{integrated }g\mid X)
=\frac{\operatorname{Var}(T\mid X)}{426^2}.
\tag{2}
\]

These integrated quantities are saved as exact rational numbers. The
thermal tails and covariance are deterministic double-precision evaluations
of exact rational T weights, not new estimates with Monte Carlo error.

## Actual noise removed on the two prefixes

The common reference is `p_ref=0.59274605079`. Both records use seed
`20260831430425`, but their marginal T distributions do not specify a joint
coupling of continuations across A and B.

| Conditional readout | A: counter 43042514269 | B: counter 43042505280 |
|---|---:|---:|
| f2(p_ref) | 0.1016281483 | 0.08617025757 |
| **suffix variance of g(p_ref)** | **0.01185405546** | **0.01171741881** |
| integrated f2 | 0.3668254042 | 0.3596742285 |
| **suffix variance of integrated g** | **0.0005374238001** | **0.0008102606562** |
| Cov[g(p_ref), integrated g] | 0.002048430905 | 0.002369455486 |

The mistaken Bernoulli baselines would be 0.09129987 and 0.07874494,
substantially overestimating the noise of the actual production readout.
Those numbers are recorded only as a semantic contrast and are not used in
any variance gain calculation.

The exact physical delay in B produces
`f2_B(p_ref)-f2_A(p_ref)=-0.01545789078` and an integrated-clock difference
`-0.007151175715`. The inherited full stochastic ordering of T implies
`f2_B(p)<=f2_A(p)` at every p, not just the saved grid. This is a comparison
of two fixed conditional laws, not population replication or a new p-value.

## Connection to the ordinary A/E readout

In `A=F1+F2-1` and `E=1+F2-F1`, F2 enters both fields with coefficient +1.
Thus this *isolated second-birth contribution* has A/E conditional covariance

\[
V_X(p,q)\begin{pmatrix}1&1\\1&1\end{pmatrix}.
\]

Conditional averaging removes that suffix component along the common A+E
direction. The source record does not explicitly give K1, so no K1 is
guessed from its age label, and no full A/E mean or covariance is claimed.
No orientation-difference covariance is inferred from these two marginal
prefix distributions either.

## What this does and does not imply for production

The total-covariance identity is

\[
\operatorname{Cov}(G)=
\operatorname{Cov}(E[G\mid X])+E[\operatorname{Cov}(G\mid X)].
\]

The present calculation gives the second term's **conditional value** at
each of two saved X, and exact averaging sets it to zero there. It does not
estimate the production average over X, the remaining between-prefix
variation, the fraction of prefixes with cheap two-port reductions, or the
cost of obtaining/solving their networks. Consequently no total production
speedup or ensemble variance-reduction factor is reported.

This is nevertheless a usable finite estimator: once a prefix's full T
law is available, the complete second-birth canonical curve and integral
can be read without drawing more suffix permutations.

## Artifacts and science card

Run `python3 scripts/p334_conditional_thermal_averaging.py`. The output
`results/p334-conditional-thermal-averaging/score.json` stores each exact T
law, 18 p values plus integrated clock, their full 19x19 conditional
covariance, source hash, exact clock moments, and explicit K1/coupling
boundaries. A focused analytic two-point test distinguishes tail variance
from Bernoulli variance. No new prefix or suffix is generated.

- **Mechanism changed:** full physical reliability now yields the canonical
  K2 profile and its removable continuation noise, not merely mean waiting
  time or a finite-trigger surrogate.
- **Not established:** no population law, overall speedup, independent
  replication, guessed K1, or cross-prefix/common-random-number covariance.
- **Observer / source:** ordinary K2 binomial-tail thermal readout,
  N425/k0=252 fixed prefixes, full physical archive 6358ba49.
- **Dependency:** the same two selected real prefixes and their exact
  conditional distributions; zero new random samples.
- **Next production quantity:** the frequency, cost, and conditional-noise
  weight of exactly solvable prefixes in an independently collected prefix
  ensemble. That—not these two conditional examples—determines total gain.
