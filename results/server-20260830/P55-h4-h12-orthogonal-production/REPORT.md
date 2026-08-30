# Issue 55 H4/H12 opposite-alias production

## Integrity and coverage

- Frozen acquisition/scorer commit: `0b8bb6300554e2fbd6e6e2712c2a882dab5f6417`.
- Each size has three 200M-replica shards and 300 retained batches, giving
  600M replicas per design.
- The exact counter coverage is contiguous and disjoint:
  `[15551000000,15751000000)`, `[15751000000,15951000000)`, and
  `[15951000000,16151000000)`.
- N305 uses seed `10276216240036424455`; N325 uses
  `11868444676518621238` on every shard.  Period matrices, signed orientation
  order and reconstruction coordinate are unchanged.
- All three hosts report the same ARM64 binary hash.  `REMOTE_SHA256SUMS.txt`
  verifies every copied file byte for byte.  The exact design gate and ten
  focused tests pass.
- Shard 1 ran on Zy rather than its nominal XP assignment because Zy was the
  first Target1 host released.  The replacement is recorded in
  `shard1/HOST_REPLACEMENT_PROVENANCE.txt`; its seed, counters, matrices,
  binary and source commit are the frozen shard-1 values.

## Frozen score order

### 1. H4-only target

The observed signed matching differences are

| size | observed DeltaM | sampling SE | frozen H4-only mean |
|---:|---:|---:|---:|
| 305 | `6.1632e-6` | `2.9671e-5` | `4.9321e-5` |
| 325 | `3.9094e-5` | `3.0269e-5` | `4.9863e-5` |

With the frozen fully correlated source-amplitude uncertainty, the H4-only
target gives chi-square `2.224/2`, p `0.329`.  It survives.

### 2. Zero effect

The sampling-covariance zero-effect score is chi-square `1.711/2`, p `0.425`.
Zero also survives.  Thus this campaign does not independently establish the
H4 amplitude even though it remains compatible with it.

### 3. Orthogonal H4/H12 coordinates

The exact opposite-alias transform gives

```text
A4  = 0.3510 +/- 0.3370
A12 = 0.1474 +/- 0.1911
A12 z = 0.771, two-sided p = 0.441
corr(A4,A12) = -0.0192
```

The H12 coordinate is unresolved.  The central A12 is also far below the
preregistered equal-amplitude benchmark `A12=0.7885` that set the production
budget, but the primary inferential statement remains the frozen zero
contrast: there is no resolved H12 contribution.

### 4. Held-out third alias

Not run.  The protocol permits it only after A12 resolves.

## Interpretation

The exact opposite-alias construction worked as intended: it converted the
historical harmonic ambiguity into a separately estimable A12 coordinate
without fitting per-size coefficients.  The result is a clean null for H12,
not an H4 identification.  Both the H4-only prediction and zero remain
compatible because the achieved sampling errors are of the same order as the
target DeltaM.

The next high-information move should not be another H4/H12 row or the
held-out alias.  It should first increase the source amplitude of the same
signed matching-odd channel through an independently motivated observable or
geometry, while preserving the exact opposite-alias coordinate.

## Five-line scientific card

1. **Question:** does the historical H4-compatible correction hide an H12
   completion with the opposite exact alias sign?
2. **Acquisition:** 600M paired replicas per N, six exact disjoint shards,
   fixed seeds/matrices and full batch covariance.
3. **Result:** H4-only p `0.329`; zero p `0.425`; A12 `0.147 +/- 0.191`,
   z `0.771`, p `0.441`.
4. **Meaning:** H12 is not resolved; H4 remains compatible but is not newly
   established because zero also survives.
5. **Stop rule:** the third alias remains closed; seek a stronger independent
   matching-odd source before another harmonic discrimination campaign.

