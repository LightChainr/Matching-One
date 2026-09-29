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
Actual run outcome and interpretation are appended below after acquisition.
