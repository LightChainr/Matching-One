# Absolute cluster fugacity: a new measured source, not an identified energy field

## The E_top response remains unresolved; the source is not a linear clock alias

The missing **absolute paired-cluster source** has now been measured on the
original N65/N130 Phase-E configurations. Its directional E_top response is
unresolved, both before and after the declared linear clock/Euler controls.
After those controls the source retains about **65–70% of its variance**.
Thus an inconclusive response is not evidence that the microscopic source is
just the already measured clock or local Euler span.

N130 has a **2.47-standard-error matching-response hint** after the full
controls; N65 does not reproduce it (−.45 standard errors). This is an
exploratory, auxiliary finite-size observation, not a new field, a norm-4
scaling law or a confirmed directional mechanism.

## What was measured

For each configuration, count black NN components and white matching
components, including only sites of the corresponding colour. Define
`S=(c_black+c_white)/N`, `q=I2−I0`, and `E_top=q²=I2+I0`.
The finite positive measure proportional to `Bernoulli(p) exp(lambda S)` has
derivatives `J_A=Cov(q,S)` and `J_E=Cov(E_top,S)` at lambda=0.
This is paired site-cluster fugacity, not an assumed one-colour bond-FK Q
derivative. Coupling to the extensive cluster count multiplies each response
by N; it does not introduce a fitted scaling exponent.

All values below are `(first−second)/Delta cos(4theta)` at
`p=.59274605079`, with the same orientation order as the original pilot.
Uncertainties are aligned delete-one-batch standard errors. The two sizes
use independent original seeds; directions and all derived views within a
size remain correlated.

### N65: no directional response is resolved

- Absolute source: `J_A=+0.000132952 ±0.000128874`;
  `J_E=−0.0000825911 ±0.0000814201`.
- After linear K/K(K−1) control: `J_E=−0.0000612535 ±0.0000782448`.
- After linear clock/Euler control:
  `J_A=−0.0000361966 ±0.0000805613` and
  `J_E=−0.0000508189 ±0.0000685958`.
- Source variance retained after full controls: first direction
  `65.5032% ±0.6859 percentage points`; second
  `64.5744% ±0.6484 percentage points`.

### N130: an auxiliary matching hint, but still no E_top response

- Absolute source: `J_A=+0.000120002 ±0.0000944509`;
  `J_E=−0.0000134106 ±0.0000529152`.
- After linear K/K(K−1) control: `J_E=−0.00000864990 ±0.0000503164`.
- After linear clock/Euler control:
  `J_A=+0.000147413 ±0.0000596019` (z=2.4733) and
  `J_E=+0.0000254518 ±0.0000454270` (z=.5603).
- Source variance retained after full controls: first direction
  `68.9498% ±0.7265 percentage points`; second
  `69.7531% ±0.7342 percentage points`.

The nominal two-size zero-response statistics for J_E are `1.09320/2`
(p=.578915) for the absolute source, `.642399/2` (p=.725279) after the
clock controls, and `.862766/2` (p=.649610) after all controls. These are
three correlated summaries, not independent votes. The estimated
delete-one covariance is not an exact finite-sample confidence certificate.
The interval widths, not a claim of zero, describe the result.

## The controls allocate variation; they do not define thermal orthogonality

The two declared auxiliary spans are `(K,K(K−1))` and
`(K,K(K−1),T_NN,chi_local)`, including an intercept through centering.
Here `chi_local=K−T_NN+F4`; F4 counts fully black elementary square faces.
This Euler count is computed independently of cluster counts. On all
80,000 paired-geometry configurations it agrees with `c_black−c_white−q`.

For each direction the source projection is fitted separately, and its
coefficients are re-estimated after every delete-one batch. This is an
empirical **linear** projection, not the unknown exact conditional mean
`E[S|K]`, not exact fixed-K centering, and not a renormalization-group
thermal subtraction. The controls contain neither q nor E_top.

A separate diagnostic records the fraction of source variance outside
`span(1,q,q²)`: about94–96%. It is not used to manufacture a residual
q/E response. Indeed, any microscopic source is seen by these two readouts
only through `E[S|q]=a+bq+cq²`; the saved sector means reproduce both
responses. New microscopic variation and three-state readout closure
coexist. Neither the residual source variance nor the 2.47-SE hint counts
continuum fields.

## Existing data were used, not a new simulation or another test campaign

The declaration was committed at `3e7eb742` before extracting these new
moments. The analysis is nevertheless retrospective: the parent pilot was
already public. It replays exactly N65 counters
`15466000000..15466019999` and N130 counters
`15466200000..15466219999`, 20,000 each. Each block has100 aligned batches
of200 configurations. No new random sample, suffix, geometry or probability
was generated, and no server was used.

All400 batch/orientation rows reproduce the original sample counts,
K1/K2 and I0/I1/I2 sums, plus every field from the completed Q/R/H edge
replay. The new information is the cluster mark and its cross moments,
including the previously missing **q² × absolute cluster count**. This
is not a repeated B or Q/R/H score.

The source is not a newly invented observable: P34's18-dimensional and
P40's20-dimensional Euler Gram archives already contain black/white
cluster counts and q cross moments. They do not store this E_top mixed
moment, and their N85 comparison is not the N130 norm-4 child.

The historical metadata declares full SHA `0578105d92d3822cb48f5c421bd23ff339295cc6`;
the resolvable original runner is instead
`05781051b76001f2b18560d7b0914f2481412584`. We preserve that discrepancy
rather than rewrite old metadata. Provenance is additionally anchored by
the unchanged backend blob, explicit seeds/counters and archived file hashes.

## What the team should do next

Use this as a **completed physical-source response**, not a new pending
pilot. Keep attention on the original norm-4 mechanism rather than repeat
the first E_top detection, synthetic certificate or these20k rows.

The inexpensive next use is to extract the already available q × absolute
cluster response from the P34/P40 million-sample Gram archives. It can
clarify matching coupling at its own N65/N85 geometries, but cannot stand
in for E_top or the missing norm-4 lineage. For the latter, the actual
missing data are source × rank moments along additional dyadic generations
and/or across the thermal coordinate—not another scalar source label.
The present seven-coordinate controls cannot reconstruct
`d_p Cov(O,S)=kappa(O,S,K)/[p(1−p)]` without the third mixed moment.

The unresolved E response does not rule out this source; it tells us what
the inherited pilot has and has not measured. A future source/geometry
prediction should be explicit enough to distinguish a clock effect from a
different coupling, while parallel connectivity and continuation work
remain open without an approval sequence.

## Reproduction and files

```bash
c++ -O3 -std=c++17 src/p154_absolute_cluster_replay.cpp -o /tmp/p154-absolute-replay
mkdir -p results/p154-absolute-cluster/raw
/tmp/p154-absolute-replay 65 results/p154-absolute-cluster/raw/n65.csv
/tmp/p154-absolute-replay 130 results/p154-absolute-cluster/raw/n130.csv
python3 scripts/analyze_p154_absolute_cluster.py
```

Use a fresh output directory when replaying an already populated checkout;
both executables refuse to overwrite saved artifacts. The scorer accepts
`--output-dir` and requires the corresponding `raw/n65.csv` and `raw/n130.csv`.
The [JSON](latest.json) preserves39 ordered outputs, the complete39×39
covariance for each N, every aligned delete-one vector and source-projection
coefficients. No inverse of that deliberately redundant full matrix is used.
The [declaration](../../analysis/p154_absolute_cluster_replay.json),
[design note](../../notes/p154-absolute-cluster-source-design.md) and
[run receipt](run.json) preserve definitions, source and output hashes.
Native ARM64 replay took .626796s and1.20048s, respectively. No repeated
scientific test suite was run; the stream/Euler checks were part of the
single analysis execution. A compact numerical readout is used here rather
than a chart: two finite-size estimates do not establish a scaling curve.
