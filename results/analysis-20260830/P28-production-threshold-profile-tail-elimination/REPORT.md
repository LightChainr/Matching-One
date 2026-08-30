# P28 production threshold-profile tail elimination

## Decision

All three frozen pure tail families are **eliminated** on the four held-out sizes.

| fixed tail law | exponent | chi-square / effective df | decision |
|---|---:|---:|---|
| Gaussian | `2` | `3,452,161.03 / 48` | eliminated |
| stretched `4/3` | `4/3` | `6,280,338.87 / 48` | eliminated |
| exponential | `1` | `18,513,710.28 / 48` | eliminated |

The survival probabilities underflow ordinary double precision.  Their high-precision base-10
logs are approximately `-749506`, `-1363631`, and `-4020063`, respectively.  Every individual
held-out size (`N=265,290,325,425`) rejects every family, and every nominal covariance direction
is retained (`48/48`).

Gaussian is the least-wrong of these three one-term approximations; it does **not** survive the
frozen test.  The result therefore does not select a replacement asymptotic law.

## Gates and covariance

The source sizes `N=130,145,170,185` fixed mean/standard-deviation standardization, the window
`|z|=2.5..3.5`, and the endpoint count gate before scoring.  All source and held-out gates pass.
The weakest source endpoint contains an expected `15,787` thresholds (`157.9` per batch), and the
weakest held-out endpoint contains `38,514` (`385.1` per batch).

For each held-out size, deleting a batch is synchronized across the two orientations, both tails,
and every point on the same reconstructed curve.  The resulting 12-dimensional delete-one
covariance block retains shared-permutation and same-curve correlations.  The four held-out RNG
domains are disjoint and are combined as independent blocks.

## Scientific interpretation

On the accessible standardized window, `log rho` has resolvable curvature beyond
`a-c|z|^alpha` for all three fixed exponents, even after fitting a separate intercept and decay
constant for every orientation and side.  The pure `4/3` proposal is therefore not an adequate
finite-size production model here.  Plausible unresolved structures include subleading tail
terms, finite-`N` order-statistic smoothing/bounded support, and the mixture of onset/completion
channels.  Those are hypotheses for a separately frozen reanalysis of the existing archives;
this score supplies no reason to generate additional Monte Carlo.

This is a Level-S production statistical model-elimination result.  It is not an exact no-go
theorem for an eventual `4/3` asymptotic tail and not a cross-microscopic universality proof.

## Post-reveal rejection decomposition

The frozen fit already gives every size, orientation, and tail side its own intercept and decay
constant.  Consequently the rejection cannot be produced by a left/right amplitude mismatch,
geometry amplitude, or cross-size coefficient drift.  It is a within-tail curvature failure on
the five frozen `z` points.

The covariance-aware marginal side scores are `2,758,737.89` on the left and `3,525,748.88` on
the right (these correlated marginal scores are not additive).  Both sides reject enormously.
The two orientation totals, `2,996,561.11` and `3,158,304.45`, are close, so no single geometry
drives the result.  P43, P50, and P57 all reject independently.

The lowest curvature mode dominates the marginal diagnostic (`1,439,918.77`), versus
`13,613.10` and `4,260.50` for the next two modes.  Signed GLS attribution alternates most strongly
over `z=2.75,3.0,3.25`, rather than accumulating at the sparse outer endpoint.  This identifies a
smooth missing correction, not an extreme-count accident.

Secant-based local effective exponents also make the failure concrete.  Across orientations they
range from `1.658..1.701` (left) and `1.945..2.099` (right) at `N=265`, drifting to
`1.609..1.639` and `1.832..1.894` at `N=425`.  Every single-side constant-effective-exponent score
rejects (`p<2e-21`), so neither side is locally a pure power on this window.  The side split and
size drift are real descriptive features, but the main rejection remains present inside each
individual size/orientation/side.

One explicitly post-reveal nested description,
`a-c z^(4/3)+d z^(2/3)`, reduces chi-square from `6,280,338.87 / 48` to
`49,893.28 / 32`, a factor of about 126.  It still fails overwhelmingly
(`log10 p approximately -10780`) and therefore does not close or receive a new model vote.
It only shows that one smooth low-order correction captures most of the visible curvature.

## Post-reveal K1/K2 mechanism decomposition

The marginal histograms identify `K1=K_minus` as the first ambient-H1 birth and `K2=K_plus` as
the second.  Their separately standardized `4/3` scores are even larger than the composite:

| birth clock | chi-square / df |
|---|---:|
| K1 | `14,918,477.15 / 48` |
| K2 | `15,965,901.67 / 48` |

Thus the composite curvature is not created by superposing two individually simple clocks.
Both clocks already have strong internal curvature.  The directional split is complementary:
K1 is dominated by its right tail (`11,950,834.22` versus `1,809,197.88` on the left), while K2
is dominated by its left tail (`12,517,817.38` versus `2,730,295.64` on the right).  These are
correlated marginal scores and are not additive model evidence.

On the composite coordinates the exact pointwise identity

```text
log rho_mix = w1 log rho_K1 + w2 log rho_K2 + H(w)-log 2
```

closes for every archived curve (maximum CDF/density/log reconstruction error below `1.6e-15`).
With the frozen composite covariance, the responsibility-weighted component shape contributes
`+9,110,751.10` to the composite quadratic form, while the clock-separation entropy contributes
`-2,830,412.23`; the latter cancels about 31% of the former and leaves the observed
`6,280,338.87`.  Mixing therefore masks part of the clock curvature rather than causing it.

Effective exponents drift downward with size in all four clock/side channels.  From `N=265` to
`N=425`, K1-left moves from `1.618–1.657` to `1.568–1.599`, K1-right from `2.195–2.238` to
`2.116–2.118`, K2-left from `1.917–1.939` to `1.880–1.914`, and K2-right from `1.873–1.997` to
`1.767–1.814`.  All component constant-beta diagnostics reject at `p<0.01`.

One boundary is explicit: the self-standardized `N=290` K1 right tails contain only
`6.6–6.8k` events (`66–68` per batch), below the original `10k/100` composite gate.  They are
diagnostic only and cannot independently carry elimination.  K1 curvature is nevertheless
resolved with ample counts in the independent `N=265,325,425` blocks; all K2 component gates pass.

## Provenance

- base: `939f7ecc26f18c68977aed626821767207a46c89`;
- scientific freeze: `9c2c595`;
- SHA-only fail-closed correction: `ba94295`;
- first invocation stopped on the mistyped source SHA before loading held-out archives;
- no new Monte Carlo was run.
