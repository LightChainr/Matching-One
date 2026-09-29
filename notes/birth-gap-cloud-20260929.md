# Independent fine-resolution birth-gap block

2026-09-29. This block follows the local pilot; its design is informed by that
pilot, not retrospectively described as an untouched research programme.

## Question and design chosen before cloud outputs

The local L=16--128 pilot has a rapidly falling direct atom and soft
near-diagonal mass approximately proportional to the examined resolution.
The useful extension is finer resolution and larger L, not a rerun of the
small exact census or a free-exponent fit.

- Square and triangular site-percolation tori, paired births in one uniform
  site-permutation filtration.
- Planned L=64,128,256,512, 14 independent batches of 10,000 samples per cell.
  A largest-cell benchmark decides whether this fits the bounded run. A
  benchmark-only stop is not production.
- New seed namespace `cloud-independent-20260929`; no local-pilot samples
  are pooled into its estimates.
- Primary descriptive grid: delta=.00625,.0125,.025,.05 for width-normalised
  label-gap CDF and soft mass Z. Also retain .1,.2 and the prescribed L^(3/4)
  reference scale, with the same covariance.
- Report Z(delta)/delta and resolution-halving behaviour without fitting an
  intercept or exponent. A visible plateau challenges a continuous-density
  picture on the resolved scales; a roughly proportional decrease does not
  prove there is no unresolved atom or logarithmically slow crossover.
- No numerical hypothesis-rejection threshold was selected. This is an
  independent exploratory block, not model-elimination certification.

## Execution context

Only DevEnvC_TV2N0X, UUID `4a8d1d443419434889e49148ed0a7ba6`, is used.
The cloud agreement initially blocked listing. The owner restored the original
account through the web interface; authenticated listing then showed Ready.
This round powered on that instance, performed one successful SSH key reset,
checked key mode 600 and connected through its dedicated port 10026.
No other instance was started, reset or stopped.

Measured: aarch64, 16 visible CPUs, CPU quota 14.5 cores, memory quota 25 GiB,
Python 3.9.9, about 14 GiB free on the 40 GiB /workspace disk. No research job
was running in the process snapshot after this startup. The base image lacked
a compiler; GCC C++ 10.3.1 and its dependencies were installed from the
configured official Huawei Cloud EulerOS repositories.

Task files stay in the new directory
`/workspace/Matching-One-TV2N0X/birth-gap-20260929`.
The prior local L=3 exact control is reused, not repeated for each worker.
Actual run outcome and interpretation follow below.

The first largest-cell benchmark took 3.188 s (square) and 3.300 s
(triangular) per 200 L=512 samples. It predicted about 431 s at 14 workers
for the originally planned block, exceeding the 300 s acquisition budget;
the runner exited before producing any research samples. Its benchmark-only
manifest is retained. Before seeing any cloud production outputs, the sample
count was reduced to **5,000 per batch**, preserving all four sizes, both
lattices, 14 independent batches and the same finer-resolution grid:
**560,000 planned production filtrations**. This changes precision, not the
readout or hypothesis after seeing a result.

## Completed result

All **560,000** cloud filtrations completed: 112 batches, 70,000 samples per
lattice/size. Acquisition including build and benchmark took **329.71 s**;
the preflight 300 s limit was a benchmark-estimate bound, not a process kill
timer. Actual concurrent runtime exceeded the extrapolation. The subsequent
standard-library analysis took **17.02 s**. No production was repeated and no
exact control was rerun on the server.

With W the pooled birth-mixture IQR, G=(N+1)(T2-T1)/W:

| Lattice | L | Direct atoms / 70,000 | Z(.00625) | Z(.0125) | Z(.025) | Z(.05) |
|---|---:|---:|---:|---:|---:|---:|
| square | 128 | 8 | .00348 ± .00016 | .00695 ± .00028 | .01390 ± .00041 | .02769 ± .00049 |
| square | 256 | 4 | .00356 ± .00015 | .00707 ± .00022 | .01421 ± .00033 | .02807 ± .00049 |
| square | 512 | 1 | .00336 ± .00015 | .00683 ± .00024 | .01399 ± .00036 | .02792 ± .00047 |
| triangular | 128 | 5 | .00373 ± .00014 | .00734 ± .00022 | .01443 ± .00036 | .02839 ± .00060 |
| triangular | 256 | 2 | .00391 ± .00015 | .00763 ± .00022 | .01479 ± .00026 | .02883 ± .00037 |
| triangular | 512 | 0 | .00378 ± .00014 | .00742 ± .00022 | .01454 ± .00038 | .02843 ± .00068 |

Errors are complete-batch delete-one SE, with W recalculated in each deletion.
CDF, Z, moments and widths from a cell are correlated views, not independent
confirmations. The full covariance and every delete-one estimate are saved.

At L=512, Z(delta)/delta on the four columns above is approximately
(.538,.546,.560,.558) for square and (.605,.594,.582,.569) for triangular.
These are descriptive ratios of the prespecified readouts, **not** fitted
slopes or independent tests. There is no resolved resolution-independent
plateau down to delta=.00625; soft mass continues approximately proportionally
to resolution. This extends the local pilot by a factor four in both L and
resolved delta, using new streams rather than a new view of its old samples.

The statement is finite: it does not rule out a smaller unresolved atom, a
slow crossover, or a different scaling limit. At zero direct counts, the
one-sided 95% bound is only 1-0.05^(1/70000), about 4.28e-5.
Square and triangular primitive tori have different conformal moduli; their
nearby numbers after width normalisation are not a universality test.

## What should follow, rather than another automatic size ladder

The evidence now favours spending the next main effort on a **two-time
anti-concentration bound**, not on L6's single-site slope or merely larger
Monte Carlo. A useful sufficient target is

    Z_L(delta) <= C delta^alpha + e_L,  alpha>0, e_L->0,

uniformly for fixed small scaled delta, together with tight birth locations.
The observed approximately linear resolution behaviour makes alpha=1 a
working conjecture to investigate, not a law established by this experiment;
any positive alpha would settle noncoalescence under the stated tightness.
The two-time pivotal/arm estimate must be derived in its actual source and
torus geometry. Triangular inputs and square universality assumptions stay
separate. Neither further static rank precision nor the direct atom supplies
this missing bound.

## Reproduction and retained artifacts

Source commit: `abc78c6a1c69e83c9a2daf4d92911833bd6f3937`.
Only engine.cpp, run_pilot.py and analyze.py were transferred from that commit;
the runner records their source hashes, GCC version, Python version, seeds,
individual commands and runtime. No full repository or credentials were sent.

Executed in the task directory:

```bash
TMPDIR=/workspace/Matching-One-TV2N0X/birth-gap-20260929/tmp \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 run_pilot.py \
  --compiler g++ --role cloud-independent-20260929 \
  --sizes 64 128 256 512 --batches 14 --samples 5000 --workers 14 \
  --budget-seconds 300 --skip-exact-control \
  --exact-control-reference abc78c6a:analysis/birth-gap-20260929/data/run.json \
  --git-base-commit abc78c6a1c69e83c9a2daf4d92911833bd6f3937 \
  --output cloud-data-5k
python3 analyze.py --data cloud-data-5k \
  --epsilons .00625 .0125 .025 .05 .1 .2 \
  --output-json cloud-summary.json --output-md cloud-SUMMARY.md
```

Use a fresh output directory to reproduce; use a fresh role to obtain an
independent block. The same seed namespace reproduces, not independently
confirms, this block.

- [Cloud tables](../analysis/birth-gap-20260929/cloud-SUMMARY.md)
- [Cloud JSON, including covariance](../analysis/birth-gap-20260929/cloud-summary.json)
- [Production manifest](../analysis/birth-gap-20260929/cloud-data-5k/run.json)
- [Initial benchmark-only stop](../analysis/birth-gap-20260929/cloud-data/run.json)
- [Analysis timing](../analysis/birth-gap-20260929/cloud-analysis-time.txt)

Raw per-batch paired histograms are the .json.gz files beside the production
manifest (about 2.2 MiB total). The initial and production benchmarks use the
same excluded benchmark streams; they are not additional independent data.

Returned archive SHA256:
`f775b25b884b26a9841d98af9ecde87bb85bba0a86a06da46cc5a1538fbd16a2`;
local transfer check matched. Cloud JSON SHA256:
`6a4b1c1e0f2e223436b822cab44225c831dc7d78b2b16b66afccbffa002bc4a0`.
The archive itself is not duplicated in Git because its contents are committed.

After transfer, a process snapshot showed no remaining research process.
The task's managed tunnel was stopped and TV2N0X was powered off; authenticated
listing confirmed **Ready**. Other instances were not operated. Raw task files
remain on the instance's data disk as well as in this research branch.
