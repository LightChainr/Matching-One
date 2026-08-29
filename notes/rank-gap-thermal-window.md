# Paired threshold-rank gap as a thermal-window observable

The threshold-rank production files contain more information than the matching
curve.  The curve uses the separate marginal distributions of `K_minus` and
`K_plus`; the retained joint moments also determine the paired rank gap

```text
G = K_plus - K_minus >= 0.
```

This pairing is lost if the two threshold histograms are analyzed separately.
It therefore supplies a small but genuinely new observable bridge rather than
another estimator of the matching root.

## Frozen scaling diagnostic

For a system with linear scale `L`, the thermal window scales as
`delta_p ~ L^(-1/nu)`.  A rank interval contains `N delta_p` occupation steps.
Using `N ~ L^2` and the percolation value `nu=4/3` gives

```text
E[G] ~ N^(1-1/(2 nu)) = N^(5/8).
```

The committed score fixes `5/8`; it does not estimate an exponent from the six
available sizes.  It reports `N^(-5/8) E[G]`, a one-amplitude GLS lack-of-fit
score, and the two aligned norm-two transfer ratios.  First and second
orientations are pooled with equal event weight before the statistic is formed,
reducing the leading orientation contrast without changing the topology.

This is a retrospective source-data analysis.  The exponent is theory-fixed,
but the observable and score were developed after these runs existed and the
result is not prospective evidence.  Its value is to expose a new joint
observable and a precise correction-flow target for a future freeze.

Every uncertainty is obtained by deleting the same batch from both
orientations.  N=65/85/130/170 share a counter interval and retain their full
cross-size jackknife covariance.  N=185 and N=265 use disjoint intervals and
are independent groups even though their numeric batch labels agree.

## Secondary joint-shape views

The same delete-one calculation reports

- `Var(G)` and `N^(-5/4) Var(G)`;
- the coefficient of variation of `G`;
- the event-level correlation of `K_minus` and `K_plus`.

These are descriptive views of the limiting thermal-window shape.  They do not
change the fixed-exponent score and are not promoted to separate evidence
blocks.  In particular, mode/root statistics reconstructed from the marginal
histograms must not be counted as independent of this result merely because
they use a different coordinate.

## Reproduction

```bash
python3 scripts/analyze_rank_gap_thermal_window.py \
  --run 65:results/server-20260828/P45-root-amplitude/n65.moments.csv:results/server-20260828/P45-root-amplitude/n65.metadata.json \
  --run 85:results/server-20260828/P45-root-amplitude/n85.moments.csv:results/server-20260828/P45-root-amplitude/n85.metadata.json \
  --run 130:results/server-20260828/P49-fullcurve-doubling-100m/raw/n130.moments.csv:results/server-20260828/P49-fullcurve-doubling-100m/raw/n130.metadata.json \
  --run 170:results/server-20260828/P49-fullcurve-doubling-100m/raw/n170.moments.csv:results/server-20260828/P49-fullcurve-doubling-100m/raw/n170.metadata.json \
  --run 185:results/server-20260828/P43-heldout-fullcurve-500m/raw/n185.moments.csv:results/server-20260828/P43-heldout-fullcurve-500m/raw/n185.metadata.json \
  --run 265:results/server-20260828/P43-heldout-fullcurve-500m/raw/n265.moments.csv:results/server-20260828/P43-heldout-fullcurve-500m/raw/n265.metadata.json \
  --output results/server-20260828/rank-gap-thermal-window/score.json
```
