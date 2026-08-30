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

## Provenance

- base: `939f7ecc26f18c68977aed626821767207a46c89`;
- scientific freeze: `9c2c595`;
- SHA-only fail-closed correction: `ba94295`;
- first invocation stopped on the mistyped source SHA before loading held-out archives;
- no new Monte Carlo was run.
