# Issue #55 H4/H12 acquisition standby

Status: preregistered and target-blind; production is not started. The three
Huawei machines remain reserved for Issue #250 until that campaign releases
them and a separate launch authorization is given.

## What was actually missing

PR #325 froze two signed Gaussian rows and the H4/H12 model, but intentionally
contained no target acquisition. The existing
`threshold_rank_integer_period_mc.cpp` already accepts arbitrary period
matrices, shares each counter permutation within a same-N pair, and retains
equal-batch threshold-rank histograms. No new Monte Carlo engine or target
model was needed. The missing layer was the acquisition coordinate, distinct
RNG domains, exact period-matrix commands, three disjoint shards, target-blind
variance pilot, and covariance-aware score.

The reconstruction coordinate is fixed at
`p_ref=0.592746050790`, inherited from the clean P31 source-amplitude campaign.
This is an acquisition/scorer coordinate, not a change to `A4=0.7885+/-0.0352`,
the exponent `13/8`, or either angular column.

## Why this is an H4/H12 discriminator

For each frozen row define

\[
y_N=N^{13/8}\Delta M_N/\Delta c_{4,N}=A_4+r_NA_{12}.
\]

The exact aliases are `r305=-1.7136858893...` and
`r325=+1.8130844968...`. Therefore the two identifiable coordinates are

\[
\widehat A_{12}={y_{325}-y_{305}\over r_{325}-r_{305}},\qquad
\widehat A_4={r_{325}y_{305}-r_{305}y_{325}\over r_{325}-r_{305}}.
\]

This is not another harmonic vote: the common row component estimates A4,
while the opposite-alias component estimates A12. With only two rows the
two-column fit is saturated (`df=0`), so the scorer reports the A12 contrast
and its covariance rather than inventing an omnibus goodness-of-fit claim.

## Target-blind 20k variance pilot and frozen count

The smoke used 20,000 paired replicas per design, 20 batches, fresh
domain-separated seeds, and counters `[15550000000,15550020000)`. The scorer
withheld both observed target means. Centered sampling SEs were `0.00550622`
at N305 and `0.00543458` at N325.

The primary budget rule is the first preregistered grid point where the full
two-row covariance gives conditional Mahalanobis separation at least 3 between
H4-only and the equal-amplitude alternative `A12=A4=0.7885`, profiling the
shared A4 column. The H12 shifts are truly opposite:

```text
N305  -8.45203e-5
N325  +9.04051e-5
```

| replicas/design | H4 vs equal-A12 distance | zero-effect H4 distance | zero-effect df2 power |
|---:|---:|---:|---:|
| 300M | 2.768 | 1.570 | 0.270 |
| 600M | 3.915 | 2.221 | 0.498 |
| 1.2B | 5.537 | 3.141 | 0.810 |
| 2.4B | 7.830 | 4.442 | 0.984 |

Thus the frozen production count is **600,000,000 paired replicas per
design**. At that count the projected `SE(A12)=0.2014`; two-sided power is
0.975 for `|A12/A4|=1`, 0.499 for `0.5`, 0.165 for `0.25`, and 0.068 for
`0.1`. Zero-effect power is secondary and did not choose the count.

## Frozen production domains and one-command launch

Production uses fresh domain-separated effective seeds:

```text
N305  10276216240036424455
N325  11868444676518621238
```

Each machine runs both sizes concurrently with eight threads per size. One
machine receives exactly one shard:

```text
shard 0  [15551000000,15751000000)  200M/design
shard 1  [15751000000,15951000000)  200M/design
shard 2  [15951000000,16151000000)  200M/design
```

After checking out the pushed branch on an idle machine, the full build,
exact self-test, ARM64 binary hash, and acquisition are one command:

```bash
scripts/launch_issue55_h4_h12_huawei_shard.sh 0
```

Use shard indices 0, 1, and 2 once each. The launcher refuses tracked local
changes and existing output files. The local single-core smoke took 0.60--0.62
seconds per 20k row; linear/OpenMP extrapolation is about 13 minutes per 200M
row at eight ideal cores. Allow approximately 15--25 minutes wall time per
machine for the two concurrent rows until an ARM64 calibration is observed.

## Frozen final scoring order

1. H4-only target with fully correlated source-A4 uncertainty plus target
   sampling covariance.
2. Zero effect with target sampling covariance.
3. The declared two-column H4/H12 coordinates above, including `A12_z`.
4. A third alias row only if A12 is resolved; it is not part of this campaign.

### Five-line scientific card

Exact: the signed Gaussian rows have opposite H12/H4 aliases, so their normalized difference identifies A12 while their weighted common component identifies A4.

Acquisition: 20k/design target-blind variance smoke froze 600M/design on the first grid point with conditional H4-vs-equal-H12 Mahalanobis distance at least 3.

Power: projected distance is 3.915 and equal-amplitude A12 power is 0.975; zero-effect H4 power is only 0.498 and is explicitly secondary.

Boundary: two rows saturate the two-column model, so a resolved A12 is a contrast result, not a goodness-of-fit result or an H12-only identification.

Standby: all seeds, counters, matrices, shards, and score order are frozen, but no Huawei production was started while Issue #250 occupies the machines.
