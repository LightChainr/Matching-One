# P406 archived-data spatial Fourier-cone result

## Decision

The positive finite-group Fourier cone is not rejected in any of the four
archived same-hand/charge autocorrelation channels.

| channel | cone distance squared | bootstrap p | 99% bootstrap quantile |
|---|---:|---:|---:|
| plus, charge 1 | 70.8780 | 0.3147 | 98.8022 |
| plus, charge 2 | 75.6920 | 0.2151 | 95.1019 |
| minus, charge 1 | 76.7658 | 0.1633 | 100.9268 |
| minus, charge 2 | 68.0672 | 0.3785 | 97.6075 |

For every channel the nonnegative-cone optimum equals the unconstrained
weighted least-squares optimum to about `1.3e-11` in squared distance.  Thus
the archived radius-four, radius-five, and radius-six endpoint data show no
resolved violation of the exact positive spatial-autocorrelation completion.

## What changed scientifically

The earlier endpoint Hankel ladder remains a correct finite-window model
elimination: ranks below eight fail for the frozen radius-six moment block.
It should no longer be read as evidence for eight physical fields or eight
transfer states.  The actual observable is a spatial autocorrelation on
`H = Z/101`, and ordinary Fourier content can inflate additive-translation
Hankel rank even for one underlying random field.

The observed support resolves an unconstrained design rank of 69 in each
channel.  One descriptive NNLS completion has 69 positive weights and inverse
participation counts between 44.11 and 48.58.  These are approximation
summaries, not unique spectra or exact mode counts: the measured displacement
window underdetermines the 101-frequency completion.  In particular, this fit
does not replace the exact global lower bound in Issue #406 or show that all of
those exact modes are statistically resolved.

The radius-six coordinates also contain three exact endpoint alias classes:

- residue 6: `(-4,-1)` and `(6,0)`;
- residue 41: `(1,-4)` and `(0,6)`;
- residue 96: `(-5,0)` and `(5,1)`.

They share Fourier-design rows and were retained as repeated noisy estimators,
not counted as distinct spatial vertices.

## Statistical contract

The score used all archived batches from independent radius-four, radius-five,
and radius-six blocks.  Each block supplied its full plug-in batch covariance;
only singular modes below the frozen relative cutoff `1e-10` were removed.
The nonnegative least-squares cone distance was calibrated by 250 plug-in
Gaussian bootstrap replicates using the frozen seed.  No new simulation was
run and no prior rank score was modified.

## Scientific card

- Mechanism space changed: ordinary finite-group spatial spectral content is a
  sufficient explanation of the observed endpoint rank growth.
- Not proved: a unique spectrum, a small physical state count, absence of a
  Jordan dilation sector, or absence of intervention-sensitive memory.
- Observer/sector/source/geometry: same-field projective-leg spatial
  autocorrelation; plus/minus hands, charges 1/2; N505 children with parent
  translation group of order 101.
- Dependency group: archived radius-four 80k, radius-five 1.2M, and radius-six
  1.2M batch blocks.
- Next lift: use held-out mesoscopic/log-radius endpoints or an explicitly
  intervention-sensitive ordered observable.  Another endpoint-only rank vote
  will mostly measure how much spatial spectrum the larger window resolves.
