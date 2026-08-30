# P50-A: prospective third Gaussian doubling lineage

## Decision

The frozen `145 -> 290` Gaussian-lineage prediction **passes**.  In lineage
order `(13,11) - (17,1)`, the fresh target run gives

```text
Delta M_290 = -0.000160648 +/- 0.000040542
frozen target = -0.0001376564 +/- 0.000024997  (source only)
residual = -0.000022992 +/- 0.000047628
z = -0.4827, chi-square = 0.2330 / 1 df
```

The sampling-only zero benchmark is `z=-3.9625`
(`chi-square=15.7018/1`).  The observed point ratio to the frozen P31 parent is
`-0.37836`, compatible with the fixed `-2^(-13/8)=-0.3242099` prediction.
This is the third independent exact Gaussian-doubling lineage and the first
target generated after its lineage-specific prediction was committed.

## Frozen design and power

The pre-production pilot used 1,000,000 replicas in counters
`[6000000000,6001000000)` and measured a target standard error of
`8.73e-4`.  Scaling that variance before viewing the production result showed
that 500,000,000 replicas would give about `2.97` expected combined standard
errors against zero after including the frozen source uncertainty.  The
production count was therefore fixed at 500,000,000, with 100 batches and the
disjoint counter interval `[6100000000,6600000000)`.

## Provenance

- Huawei environment: `DevEnvC_ZyTrST`, Kunpeng/AArch64, 16 vCPU, 32 GiB.
- Clean detached source commit:
  `4bbbe90180a9960dfae613b47548897ae1defa8b`.
- Source SHA-256:
  `9507452e451b050a3c47d5563faf992b7a491640222702f15aef46720968ce37`.
- Binary SHA-256:
  `3028b67759af20bef7cfbdf6a0ce08f7957c2f50a91601de1f93d1f28f4462fe`.
- Compiler/flags: GCC 10.3.1, `-O3 -std=c++17 -fopenmp`.
- Seed: `2026105001`; 16 OpenMP threads.
- Production wall time: 1259.71 seconds.
- Exact exhaustive self-test passed before the pilot and production run.
- Stderr is empty; retained files are covered by `checksums.sha256`.

## Scope

This completes the fixed-`p` P47-A component of Issue #50.  It does not yet
complete the full-curve slope/root triptych, the logarithmic residual fit, or
the N=1105 commuting square.  The aligned N=130/170 threshold-rank production
was launched immediately afterward to advance the first two full-curve
doubling lineages.
