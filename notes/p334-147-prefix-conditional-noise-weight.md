# P334: the conditional clock removes 84–86% of noise in a fixed 147-prefix mixture

This is the first multi-prefix noise budget for the exact continuation clock,
not another clock solver. It consumes the complete supplied 147-member eligible
old-source set in commit `9cca7bc60e26db5ec47b5e00fbc5d98532447c29`,
`results/p334-all147-prefix-clocks/full_clocks.json`. Every member has
`N=425`, orientation `second`, `k0=252`, age `10`, and primitive line
`(12,-19)`. The original source exclusions are retained. This is not a draw
from all possible prefixes, and the reused twelve clocks are not an independent
block.

## The two fixed readouts

Let `X` be uniform on these 147 named prefixes and, conditional on `X`, let
`T` have its committed exact continuation law. No new continuation is drawn.
For each possible `T`, the ordinary production readout is the *real-valued*
canonical binomial tail, not a Bernoulli event:

\[
g_T(p)=\Pr\{\operatorname{Bin}(425,p)\ge252+T\},\qquad
I_T=\int_0^1g_T(p)\,dp=\frac{174-T}{426}.
\]

The endpoints are fixed at `p_ref=0.59274605079` and the integrated clock.
Write `G=(g_T(p_ref),I_T)`, `m_X=E[G|X]`, and
`V_X=Cov(G|X)`. The exact finite-mixture decomposition is

\[
\Sigma_{\rm fresh}=\frac1{147}\sum_X V_X+
\frac1{147}\sum_X(m_X-\bar m)(m_X-\bar m)^T
=\Sigma_{\rm suffix}+\Sigma_{\rm prefix}.
\]

Replacing one fresh suffix readout by its exact conditional average `m_X`
removes `Sigma_suffix`, and leaves `Sigma_prefix`. Denominator 147 is deliberate:
this is a fully specified finite empirical distribution, not an estimated
population variance needing an `n-1` correction.

| Readout | Mixture mean | Suffix variance | Prefix variance | Total variance | Removed fraction |
|---|---:|---:|---:|---:|---:|
| `g_T(p_ref)` | 0.1706204126 | 0.02092217261 | 0.003543006813 | 0.02446517943 | 85.5182% |
| integrated `I_T` | 0.3744119018 | 0.0005775498650 | 0.0001104170921 | 0.0006879669571 | 83.9502% |

The off-diagonal covariance is `0.00283480067916` within prefixes and
`0.000607719826059` between prefixes, totaling `0.00344252050522`.
Both full matrices and all 147 per-prefix means/covariances are in the score.
Integrated-clock quantities are also stored as exact rational numbers.
Canonical quantities use exact probability weights evaluated numerically
against the binomial tails. No sampling SE or population CI is assigned to
these deterministic functionals of the declared mixture.

## What changes scientifically

The exact clock is not merely capable of averaging isolated witnesses: in
this actual multi-prefix conditional stratum, most of the variance in both
the near-critical second-birth thermal readout and its integrated clock is
continuation noise. Nontrivial prefix structure remains after it is removed.
This supplies a concrete variance target for reusing the solver in production.

The claim is about the specified experiment “choose one of these 147 prefixes
uniformly, then choose a fresh uniform suffix.” If instead one averages the
known conditional means over *all* 147 prefixes, that finite-mixture mean is
already a deterministic number. Neither observation estimates the cost of
generating prefixes, the frequency of this stratum, or solver success/cost in
other strata, so neither is a global MC speedup claim. It also does not specify
cross-orientation coupling. `K1` is not inferred from age: the output is the
second-birth contribution alone, not full `A_top`/`E_top` covariance.

## Scientific card

- **Changed mechanism space:** exact conditional averaging has an explicit
  83.95–85.52% removable noise share in the fixed empirical stratum; continuation
  noise dominates but is not the entire signal variance.
- **Not proved:** a global production speedup, population effect size, or
  sufficiency of any reduced prefix state.
- **Observer / sector / source / geometry:** canonical second-birth CDF and
  its integral; fixed second-orientation `N425`, `k0=252`, age10,
  `ell=(12,-19)`; source `9cca7bc6`.
- **Dependency group:** the same 147 old-source prefix set used by the full-clock
  and crossing analyses; no independent new block.
- **Next quantity that would raise its scope:** stratum-weighted solver coverage
  and cost, together with the covariance of the actual coupled production
  readouts. This note does not launch that collection.

Reproduce using only the committed input:

```sh
python3 scripts/p334_147_prefix_noise_mixture.py \
  --source-commit 9cca7bc60e26db5ec47b5e00fbc5d98532447c29
```

No network solve, new Monte Carlo, new prefix, p-grid expansion, or extra test
run was used for this readout.
