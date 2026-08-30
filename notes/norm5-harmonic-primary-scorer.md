# Frozen norm-5 primary scorer

Issue #57 freezes the two norm-5 Gaussian lineages, fixed probability, radial exponent, and H4/H12/H8/zero scoring order before the N=325/N=425 target reveal. The original 2026-08-28 scorer contract and computational kernel are retained unchanged.

## Canonical entrypoint

After the Issue #43 cross/either erratum and the type-safe channel layer in PR #170, the canonical production entrypoint is

```text
scripts/score_norm5_harmonic_primary_typed.py
```

The underlying frozen kernel remains

```text
scripts/score_norm5_harmonic_primary.py
```

with its original implementation contract in

```text
predictions/norm5_harmonic_primary_scorer_20260828.yaml.
```

The added semantic gate is

```text
predictions/norm5_harmonic_semantic_gate_20260829.yaml.
```

This separation preserves chronology: the numerical score is not redefined after the protocol correction; only the source/target observable semantics are made explicit before target scoring.

## Primary residual

For each lineage and fixed model the kernel evaluates

```text
R = DeltaM_child - r_model DeltaM_parent.
```

The two residuals are scored jointly with two degrees of freedom. H4, H12, and H8 ratios are recomputed from the exact rational angular factors and `5^(-13/8)` in the original frozen prediction artifact. The zero model remains fourth.

No target exponent or amplitude is fitted.

## Exact observable map

The frozen raw matching-function contrast is matching-odd. The historical source convention is represented by P31 `either/odd`, while the threshold-rank `K_minus/K_plus` target machinery reconstructs `cross/odd`.

For the declared orientation order, complementary torus topology gives exactly

```text
D_either = D_cross.
```

Therefore the semantic map is

```text
source: either / odd / p / first-minus-second / raw
 target: cross / odd / p / first-minus-second / raw
 scale:  +1
 offset: 0
```

The typed wrapper validates this registered map before invoking the frozen kernel. If the map is unavailable or changes from unit scale, scoring fails. After a successful score, the wrapper appends the source descriptor, target descriptor, exact transform, and semantic manifest to the output JSON.

This is materially different from the matching-even Issue #43 correction, where `either/even -> cross/even` has scale `-1` for orientation contrasts.

## Counter-aware covariance

The N=65 and N=85 parent runs deliberately share the same seed and counter interval, so their delete-one replicates supply a real cross-size covariance. The N=325 and N=425 target runs use disjoint counter intervals and are independent groups even if both files label batches `0..99`.

The frozen kernel infers covariance groups from the exact `(seed, first counter, last counter)` tuple and rejects partial overlaps whose covariance is not defined by the protocol. Each model has its own multiplier ratio, so residual covariance is propagated separately. A diagonal-only chi-square is also reported as the frozen sensitivity check.

## Production invocation

```bash
python scripts/score_norm5_harmonic_primary_typed.py \
  --run 65:parent/n65.hist.csv:parent/n65.metadata.json \
  --run 85:parent/n85.hist.csv:parent/n85.metadata.json \
  --run 325:child/n325.hist.csv:child/n325.metadata.json \
  --run 425:child/n425.hist.csv:child/n425.metadata.json \
  --output results/issue-57/primary_score.json
```

The output records SHA-256 hashes for the prediction artifact and every raw input, plus the typed observable semantics. No result is emitted until all four frozen sizes are present.

## Claim boundary

The semantic gate does not alter:

- the frozen child geometries;
- the H4/H12/H8/zero order;
- `13/8`;
- any multiplier ratio;
- covariance rules;
- target sample data;
- the zero target-refit parameter count.

It only makes the exact channel identity part of the executable score contract, as required by Issue #146.
