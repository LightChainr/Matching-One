# RNG domain policy for new production

The production default is intentionally asymmetric:

- orientations compared at the same `N` share one common field;
- distinct sizes use domains derived from `(experiment_tag, N, base_seed)`;
- cross-size coupling is an explicit, named experimental design rather than an
  accidental consequence of reusing a seed.

`scripts/rng_domain_policy.py` derives a stable 64-bit engine seed and emits the
metadata record that belongs beside a production command.  For example:

```bash
python scripts/rng_domain_policy.py \
  --base-seed 2026105001 --experiment-tag P50-fullcurve --size 290
```

The same record covers every orientation in one same-size paired job.  A job
orchestrator must not derive a separate seed for each orientation because that
would discard the useful common-random-number covariance.

Intentional cross-size coupling must instead use
`--mode intentional_cross_size_coupling --coupled-residual NAME`.  It keeps the
base seed across sizes and records that aligned batches and full cross-size
covariance are mandatory.  The named residual and the decision to couple must
be frozen before production.  Historical streams are never rewritten.

The P45 root-amplitude campaign is the clean replay required by Issue #39: it
records clean source and binary hashes, commands and counter interval, verifies
one-thread/four-thread equality, and uses aligned full covariance for its
deliberate coupled-size design.
