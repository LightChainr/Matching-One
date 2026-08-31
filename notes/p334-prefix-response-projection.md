# Named prefix features after fixing the rank pair

## Result: contact structure organizes the remaining clock response

The response differences left after fixing the rank pair now have a small,
physical descriptor family. Exact source energy alone captures about half
of the **own-source, own-geometry within-rank-cell covariance loading**.
The four declared contact descriptors capture 80–99% at the point estimates.
This is a concrete organization of the response beyond rank-cell identity,
using existing data rather than another sampling run.

Original8 result `90226598`, with original20-batch delete-one standard errors:

| N / receiving geometry | full own-source loading | source-energy captured share | four-contact captured share | remaining loading |
|---|---:|---:|---:|---:|
| 325 / first | `8.9633e-8 ± 1.3470e-8` | `51.68 ± 8.05 pp` | `90.01 ± 9.29 pp` | `8.9578e-9 ± 9.2003e-9` |
| 325 / second | `8.6211e-8 ± 1.1830e-8` | `49.09 ± 9.25 pp` | `80.28 ± 8.08 pp` | `1.6999e-8 ± 7.7924e-9` |
| 425 / first | `8.2643e-8 ± 8.2704e-9` | `51.86 ± 7.63 pp` | `97.02 ± 6.97 pp` | `2.4642e-9 ± 5.9149e-9` |
| 425 / second | `7.6224e-8 ± 1.1376e-8` | `53.10 ± 7.60 pp` | `99.00 ± 10.05 pp` | `7.6293e-10 ± 7.7432e-9` |

Shares are signed loading ratios, with uncertainty in percentage points;
they are not variance fractions. All rows use the full20000-prefix mass,
not a rank-zero-renormalized population. The two orientations within a size
share the original prefixes and suffixes. Cross-source response coordinates
are retained in the result, not pooled into these own-source summaries.

### Contact information remaining after baseline clock means

The declared joint-safe mass has positive partial covariance with own-source
center response after removing the two latent baseline clock means:

| N / geometry | `Cov(T_safe, H_C | linear mu_C,mu_W; within G)` |
|---|---:|
| 325 / first | `9.5478e-8 ± 3.9203e-8` |
| 325 / second | `8.2904e-8 ± 3.3880e-8` |
| 425 / first | `9.0137e-8 ± 3.7228e-8` |
| 425 / second | `1.07475e-7 ± 2.3633e-8` |

This is a non-tautological center-response direction shared in sign by all
four readouts. The analogous lifetime partial covariances are weaker. The
four-feature incremental projected center-response variance is
`3.43e-11 ± 2.21e-11`, `4.01e-11 ± 2.68e-11`,
`2.68e-11 ± 1.91e-11`, `4.81e-11 ± 1.74e-11` in the same row order.
The underlying deconvolved response-variance moments are still noisy,
including negative estimates; no response R-squared is identified here.

The mechanism lead is **contact-regulated clock susceptibility**: source
energy measures available perturbation strength, while contact structure
adds organization of the center response that rank pair and the two mean
clocks do not exhaust linearly. The fitted four-feature block does not
identify a unique causal contact feature or assert that the complete
response law is linear. A high share of this signed loading need not be a
high share of the full response heterogeneity.

The paired contact-minus-energy loading gains are `3.435e-8 ± 8.664e-9`,
`2.689e-8 ± 8.869e-9`, `3.732e-8 ± 7.759e-9`,
`3.499e-8 ± 4.950e-9` in the same row order. Thus the contact block's added
loading is measured directly with shared-batch covariance, rather than
inferred by comparing two separate error bars.

### Relation to the new local-two-dimensional result

The separately delivered
[fixed-prefix local response study](https://github.com/LightChainr/Matching-One/blob/8ad30617b0a3076a5c01a208eb213096d8879b32/experiments/p334-mechanism-response-20260831/REPORT.md)
finds positive mean local determinants for the two birth-center responses
after adding64 quartets on each original00 prefix. That result answers
whether two source directions survive locally. This readout answers which
prefix structure organizes their signed response loading. Local response
rank and population loading organization can coexist; neither is a count
of continuum fields. The original8 projection below remains frozen,
separate from the new64 stream-reuse readout.

## Same00 new64 reuse: about80% captured, about20% remains

The [standalone readout `323de7d5`](https://github.com/LightChainr/Matching-One/blob/323de7d5ee4a980b3c77e1a972cb6c812a9f88e5/notes/p334-new64-feature-loading.md)
now isolates the original00 subset:1502 prefixes
at N325 and1551 at N425. It pairs the original8 estimates of baseline
clock means with the already collected independent new64 estimates of
the same conditional responses. No trajectory was regenerated. The
four-feature definitions stay fixed, with separate descriptive coefficients
for this population and stream; this is not a held-out prediction test.

| N / geometry | source-energy share | four-contact share | remaining loading |
|---|---:|---:|---:|
| 325 / first | `52.89% ± 5.05 pp` | `80.36% ± 5.66 pp` | `4.976e-9 ± 1.727e-9` |
| 325 / second | `46.35% ± 7.43 pp` | `78.20% ± 5.84 pp` | `5.190e-9 ± 1.851e-9` |
| 425 / first | `55.40% ± 5.16 pp` | `79.49% ± 6.24 pp` | `4.651e-9 ± 1.676e-9` |
| 425 / second | `49.71% ± 6.77 pp` | `79.66% ± 6.17 pp` | `3.885e-9 ± 1.351e-9` |

The two-source-supporting population therefore has a consistent point
pattern: about half from source energy, about four fifths aligned with
the full contact block, and about one fifth left outside that linear
description. The four residuals are about2.8–2.9 original-batch SE each;
these paired geometries are correlated, not four independent replications.
This makes the residual an explicit next target without declaring the
four-descriptor description a sufficient state.

The safe-mass/own-center partial cross-moments after both mean clocks are
`(2.837 ± 0.668, 2.967 ± 0.913, 3.005 ± 0.659, 3.396 ± 0.639)e-8`
in table order. The cross-stream product is conditionally unbiased because
old8 clocks and new64 tangents use independent suffix streams at the same
prefix. The clock Gram still uses old8 distinct-quartet products. Both
streams retain the same prefix population and original20 batch deletions.

This is a different population from the receiver-R0 three-cell pool above;
its approximately80% shares should not be subtracted from the pooled
80–99% shares as a stream effect. Same00 old8/new64 paired contrasts are
saved in the standalone result. All00 covariances remain weighted by
their prevalence among the original20000 prefixes; the source population
has not been renormalized or independently resampled.

## Fixed, low-dimensional readout

For each receiving geometry, consider original prefixes where that geometry
has rank zero. Retain the three paired-rank cells separately for centering,
then pool their within-cell moments with their original population weights.
There is one slope vector, with a separate intercept for each rank cell.

The four exact prefix predictors are joint-safe label mass, own physical
score energy, safe mean contact degree, and safe mean R0 loop count. All safe
means have the original vacant-label denominator, not the safe-subset
denominator. The physical score is

```
s_o(u)=pi_a*(L_o(u)-mean_a L_o),
source_energy_o=E_uniform[s_o^2],
```

and is zero outside jointly safe labels. The other two predictors are the
latent baseline means mu_C and mu_W. Responses are conditional mean-clock
derivatives under either physical source:

```
H_Lfirst=H_plus+H_minus,
H_Lsecond=H_plus-H_minus.
```

The source-strength-only, four-contact-feature, two-clock, and combined
six-feature projections are fixed in `ffad2d75`, before descriptor results
were consumed. No feature or regularization search is performed.

## Latent products instead of noisy mean squares

Each prefix has eight independent quartets. With the exact census feature
repeated across quartets, append the baseline and tangent clock estimates to
one vector z_q. The ordered off-diagonal product

```
[(sum_q z_q)(sum_q z_q)^T-sum_q z_q z_q^T]/[Q(Q-1)]
```

estimates products of its latent prefix means. Distinct-prefix products
center those moments within each rank cell. This removes shared-label/tail
noise from the predictor and response covariance moments; it does not
create new population replicates. Original20-batch deletions carry the
uncertainty of the full projection, including the latent Gram matrix.

For K=Cov(F,F) and v=Cov(F,r), the slope is K^{-1}v. Matrices are rescaled
to correlation units for the solve; an unidentified or indefinite predictor
Gram is not repaired with a ridge. The true projected response variance is
a lower bound on latent response variance. Its finite-sample estimate and
the unbiased response-variance moment remain estimates; no R-squared is
formed from an unidentified or negative response-variance estimate.

## Interpretation

The main physical summary is the own-source contribution

```
2 Cov(mu_C,H_C)-Cov(mu_W,H_W)/2
```

to the within-rank-cell covariance response of the normalized birth ranks.
For the source-strength and contact projections, compare the predicted
part and its residual. A signed loading share can exceed one or be negative;
it is not a probability or a fraction of baseline noise explained.

A clock-only predictor necessarily reproduces clock-response cross moments.
That algebraic identity is not evidence that two means determine the full
conditional law. The separate partial covariance of exact contact features
with the response after projecting out mu_C,mu_W asks a non-tautological
question. It remains a descriptive conditional association, not causal
mediation or an out-of-sample result.

This analysis reuses original source `375cd3a1` and the complete label census;
it adds no prefix, quartet, suffix, cloud run, local-rank determinant test,
finite-policy run or independent evidence block.

## Sources and lifecycle

- Reader `011f50e3`; result `90226598`,
  [`score.json`](../results/p334-prefix-response-projection/score.json),
  including726 raw coordinates, all derived estimates and original20 LOO
  vectors per size. The numerical reader completed in under a second.
- Exact descriptors
  [`1cfa4ae8`](https://github.com/LightChainr/Matching-One/blob/1cfa4ae892a2f7f4168e9a71690efd7a5560d4cd/notes/p334-exact-prefix-structure.md),
  source identity `375cd3a1`, census `ac5761ce`.
- The complete original-block covariance join and paired contact-minus-energy
  contrasts are delivered in
  [`172fbeb1`](https://github.com/LightChainr/Matching-One/blob/172fbeb1ed28019b04f132859da3ea247942bb1c/notes/p334-prefix-response-projection-joint.md).
  It reuses the producer's LOO without replay or refitting and retains all
  prior same-block coordinates. The new64 descriptive readout has its own
  aligned20-batch LOO, rather than being relabeled as original8 data.
- The new64 readout `323de7d5` uses reader `5c4c9b45`, consumes only
  existing committed prefix-level NPZ statistics, and completed in0.43s.
  Its same00 old8/new64 loading differences are within one paired SE for
  all four cases; no changed stream estimand or independent replication
  is claimed. No response-variance ratio is computed from new64 mean-only
  coordinates.
- Projection identities and population/empirical-mixture distinction:
  [`a940e65d`](https://github.com/LightChainr/Matching-One/blob/a940e65dd73a08815b627660cc4ccd529840afed/notes/p334-rankcell-response-projection.md).
- Lifecycle: exploratory finite-size mechanism result / normalized birth
  clocks `C,W` / exact physical `L_first,L_second` sources / paired N325,
  N425 geometries / original `9c495ab1` prefix population and `e32a8593`
  conditional suffix block / branch-delivered, not independently replicated.
- Next useful upgrade: account for the remaining approximately20% loading
  on original00 while preserving the center-response direction. The PR509
  team already owns held-out `J=B G`/contact prediction and conditional
  birth-shape analysis; these descriptive moments are inputs, not competing
  repetitions of those tasks. No new production or local-rank rerun is needed.
