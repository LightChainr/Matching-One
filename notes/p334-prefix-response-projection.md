# Named prefix features after fixing the rank pair

The preceding covariance hierarchy found a resolved response difference among
prefixes with the same rank pair. This readout asks whether it can be traced
to the amount of available source, to elementary contact structure, or to
information beyond the two baseline birth-clock means.

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
