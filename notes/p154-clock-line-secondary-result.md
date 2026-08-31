# P154: weak fresh transmission cannot identify a surviving clock

## Scientific decision

The completed independent N85/N340 lag1 experiment activates its original
stop rule: stop prioritizing this particular source as the main explanation
of the H4 response at the declared resolution. Both entry-selective and
completion-selective numerical templates fail; all four channel intervals
lie inside the weak-channel band. No lag change, additional sample or
replacement numerical template is attached to this result.

The pre-output clock map establishes **sensitivity**, not source occupancy:
a large coefficient for a possible angular clock does not imply that the
actual centered cluster-count source excites that direction. A null or
weak response can lie on several amplitude-free lines simultaneously.
Consequently a secondary line that is not excluded cannot rescue a failed
primary template, demonstrate transmission, or identify the source's clock.

The primary estimates below are copied from the completed official result,
not recomputed here. Intervals are its original six-coordinate simultaneous
95% normal intervals, using the frozen critical value 2.6382572735.

| N | Entry response | Completion response | Net response | Net simultaneous interval |
|---:|---:|---:|---:|---:|
|85| .044492 | -.001022 | .043470 |[-.071640,.158581]|
|340| -.053571 | .114247 | .060675 |[-.157394,.278745]|

This concerns the fixed lag1 rank-centered CB+CW conditional policy, in bulk
source units, for the two declared sizes. It does not prove exact zero,
reject other sources, identify continuum fields or establish an all-size
scaling law.

## Secondary registration and source boundary

The two pure locally flat clock restrictions M10 and M11 were derived in
`c2828e3430fe1ac7e02fbe0e5ddc0e6a24c99847`. Their four-residual family,
gain-line normalization, joint omission covariance and rejection rule were
frozen in `83f3eba88d7f1290704f82610c28669dc5e12f3c` before outcome inspection.
Implementation `9c6e353cf1325401069cecd08a1a9ddd0339c951`, integrated as
`4129d47d`, was committed and pushed before root read the fresh decisions.

This is **outcome-blind secondary registration during/after production**,
not a second pre-production primary protocol. The root first encountered
the primary decision in the local fleet registry around20:05 CST, after
implementation had been pushed, and then read the final report. The score
author did not read fresh outcomes while implementing the fixed rule.

Only unmarked q/E baseline jets are newly reduced from the nine completed
raw shards; source response values and the corresponding200 omissions per
N are reused from the official result. The same omission recomputes the
baseline root and line gains. The saved joint10x10 covariance keeps all six
official coordinates and four secondary residuals in one dependency block
per N. They are not ten independent observations or a revised ten-test
primary family. No source alpha, mixture amplitude or new observer is fit.

Formal secondary intervals are Bonferroni-adjusted normal intervals over
four residuals, not finite-sample certificates. Both lines contain zero;
failure to reject is not equivalence. An arbitrary two-amplitude mixture
would generally saturate the two-channel response plane, so it is not
offered as a post-result explanation.

## Completed secondary score

All four simultaneous intervals contain zero. Both M10 and M11 have the
frozen status **not_excluded**; neither is identified or promoted. The
primary stop rule remains active.

| N | Pure flat restriction | Normalized line residual | SE | Secondary simultaneous95% interval |
|---:|---|---:|---:|---:|
|85|M10, angular/common-birth|-.03462443|.03171942|[-.11385020,.04460134]|
|85|M11, angular/relative-birth|-.03330387|.03124154|[-.11133604,.04472829]|
|340|M10, angular/common-birth|.11639950|.06634273|[-.04930510,.28210410]|
|340|M11, angular/relative-birth|-.03698585|.05853971|[-.18320081,.10922911]|

The normal critical value is2.4977054744 for these four secondary residuals.
The residual is `(C_entry*v_completion-C_completion*v_entry)/hypot(C_entry,C_completion)`.
Fresh baseline gains `(entry,completion)` are `(27.921618,32.996813)` and
`(-27.921701,32.996730)` atN85; `(117.001939,129.186168)` and
`(-117.001950,129.186157)` atN340. These are unmarked-baseline sensitivities;
no source amplitudes are inferred from them. Gain uncertainty and its
covariance with the response are included through the paired omissions.

**What changed:** the independent primary experiment bounds this source's
finite-scale global transmission, and the secondary comparison supplies no
clock identity with which to override that result. The exact high-gain map
remains useful algebra; this source is not demonstrated to occupy its
high-gain direction. A future revival requires a separately justified
experiment, not further fitting of this block.

## Immutable delivery

- Official raw and result delivery: [f4999e29612da16a3650f24d124fb59137f053d7](https://github.com/LightChainr/Matching-One/blob/f4999e29612da16a3650f24d124fb59137f053d7/experiments/p154-prospective-transmission-20260831/REPORT.md),
  produced under the original primary freeze0820b8d2; confirmed pushed.
- Official JSON SHA256:
  `2c02bf3214ba4c9b31e8ad7ae65addf6d1ff0a82882518d1151556e7d0ce6821`.
- Secondary reader checkout: `b7b40ab3bf2cbe484a59295363ec97bf186543a8`;
  output [score.json](../results/p154-clock-line-secondary/score.json).
  All nine raw/receipt identities and all joint factors are saved there.
- The official source scorer was not run again. One local secondary pass
  used the managed research-py311 environment with BLAS/OMP threads1;
  no installation, simulation, source refit or additional test suite.
- According to the producing team's20:02:29 CST final receipt, its five
  hosts were returned Ready and all ten hosts were observed Ready.
  This team made no cloud lifecycle change in the secondary analysis.

Reproduce the secondary result from the fixed package and implementation,
using a new output path rather than overwriting the delivered result:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/p154_clock_line_secondary.py \
  --official-package /path/to/f4999e29/experiments/p154-prospective-transmission-20260831 \
  --official-result /path/to/f4999e29/experiments/p154-prospective-transmission-20260831/PROSPECTIVE_RESULT.json \
  --raw-source-commit f4999e29612da16a3650f24d124fb59137f053d7 \
  --official-result-commit f4999e29612da16a3650f24d124fb59137f053d7 \
  --output /path/to/new-clock-line-score.json
```
