# P250 charged three-point 20k variance smoke

- Status: model-development variance smoke; not production evidence.
- Exact map: 65 translations and 1,950 root labels passed; maximum DFT
  conjugacy residual `1.71e-17`.
- Runtime: 194.69 seconds on eight local workers for 20,000 replicas.
- Frozen-score plumbing: closure p `0.594`; nonneutral joint-zero p `0.953`;
  all three H4/H8/H12 templates ran without fitting a channel-specific phase.
- Production consequence: freeze one million fresh-seed replicas, giving a
  covariance-only projected worst primary component SE near `3.02e-7`.

The smoke point estimates are neither selected nor combined with production.
