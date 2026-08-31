# Fresh regular-pair spatial kernel: completed raw batches

## Completed parent readout

The parent has now run the frozen score once: the L64/r16 zero-transmission
null is rejected, with C64=6.85546875e-6 and a 99% Monte Carlo interval
[5.20339728e-6,8.50754022e-6]. The fixed C64/C32 ratio is 0.1873498799
(99% Fieller interval [0.1401616411,0.2381648820]). See the
[main result](../../notes/regular-pair-spatial-transmission-result.md),
[score report](score/REPORT.md), and [joint numerical output](score/score.json).
The original producer handoff below is retained as a raw-data receipt;
its statement that no scoring had occurred describes that earlier handoff.

## Original raw-data handoff

One frozen production block per size is complete. L32 (`r=8`, seed `2026083123593201`) took **2.134110625 s**; L64 (`r=16`, seed `2026083123596401`) took **7.974634458 s**. Each has **200 independent batches × 1000 iid configurations**. Two local workers ran concurrently after one compilation; both exited zero. No benchmark, L16 run, restart, top-up, cloud computation or scoring was performed.

The source is the prescribed two-site regular-completion Q-jet kernel `g_pi`, including its signed values. It is not a covariance of the one-site activation mark. Kernel source `32ff99fa5361ba0fe435fac835be2dbb206e0a6c`, theorem `7f60e92d5cdb58e7542db06cd49547a4451ba022`, contract `3210aeb338ca7bb52c799d1de9048232f50ab921`, producer `9f6ff44dd41764fc34a251b202494172e62228b6`, and scorer `0096e79469fc9b8de00ebdafe52226345c65c364` all preceded the data.

`L32.csv` and `L64.csv` retain the full batch sums, their shared-component-count 0–4 contributions, eligible-pair counts and nonzero-support counts. The 16 fixed anchors and two fixed directions produce **32 correlated pairs per configuration**, not 32 independent samples. Divide every `sum_g16` or shared-component sum by **`16*32*1000`** for a batch mean. Both-endpoint vacancy gates the contribution, never the denominator. Missing sparse-lookup keys are exact zero; negative values are retained.

The two sizes use separate fresh RNG streams with the fixed full-uint64 Bernoulli rule `word < 10934234699625173385`. No old exact or Monte Carlo archive is pooled. Per-size metadata records the exact probability and RNG semantics. `run.json` records actual commands, completion times, source hashes and CSV/metadata hashes. The frozen scorer is left to the parent task; this delivery makes no numerical transmission or field claim.
