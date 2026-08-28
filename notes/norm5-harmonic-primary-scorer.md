# Frozen norm-5 primary scorer

Issue #57 already freezes the two norm-5 Gaussian lineages, the fixed
probability, the radial exponent and the H4/H12/H8 model order.  This change
turns that protocol into the executable scorer used when the N=325 and N=425
production histograms arrive.  It does not add or fit a post-target model.

For each lineage and fixed model it evaluates

```text
R = DeltaM_child - r_model DeltaM_parent.
```

The two residuals are scored jointly with two degrees of freedom.  The H4,
H12 and H8 ratios are recomputed from the exact rational angular factors and
`5^(-13/8)` in the original frozen prediction artifact.  The zero model is
reported fourth.

## Counter-aware covariance

The N=65 and N=85 parent runs deliberately share the same seed and counter
interval, so their delete-one replicates supply a real cross-size covariance.
The Huawei N=325 and N=425 runs use disjoint counter intervals.  They are
therefore independent groups even though both files label their batches
`0..99`.  The scorer infers groups from the exact `(seed, first counter, last
counter)` tuple and rejects partial overlaps whose covariance is not defined by
the protocol.

Each model has a different multiplier ratio, so its residual covariance is
propagated separately.  A diagonal-only chi-square is also reported as the
frozen sensitivity check.

## Production invocation

```bash
python scripts/score_norm5_harmonic_primary.py \
  --run 65:parent/n65.hist.csv:parent/n65.metadata.json \
  --run 85:parent/n85.hist.csv:parent/n85.metadata.json \
  --run 325:child/n325.hist.csv:child/n325.metadata.json \
  --run 425:child/n425.hist.csv:child/n425.metadata.json \
  --output results/issue-57/primary_score.json
```

The output records SHA-256 hashes for every raw input and the prediction
artifact.  No result is emitted until all four frozen sizes are present.
