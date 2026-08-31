# Zero-source Hessian: resolved own-source curvature, weak mixed response

The fixed two-source family now has its complete `ff/fs/ss` second-response
tensor for physical receivers and `A_ref,E_ref,C,W`. The separately retained
New64 stream resolves negative own-source center curvature in both directions
and sizes. Mixed `fs` response does not show a common cross-size signal at the
present precision. This is a new nonlinear-response readout of existing paths,
not another fit of the contact-loading residual.

## Direct New64 readout

All entries retain the original20000-prefix denominator, though New64 contains
only the original1502/1551 double-R0 prefixes. Errors are one original20-batch
SE, without treating source/receiver directions as independent observations.

|N / receiver|Own A_ref curvature ×10⁷|Own C curvature ×10⁸|Mixed A_ref curvature ×10⁷|Mixed C curvature ×10⁸|
|---|---:|---:|---:|---:|
|325 first|5.404 ±1.248|−6.222 ±1.208|−0.023 ±0.949|−0.311 ±0.970|
|325 second|5.876 ±1.294|−5.985 ±1.148|1.711 ±0.972|−2.229 ±0.924|
|425 first|10.232 ±1.386|−8.493 ±1.014|0.308 ±1.192|0.833 ±0.896|
|425 second|6.771 ±0.973|−6.520 ±0.836|0.185 ±1.187|−0.418 ±1.015|

Here own means `ff` for first receiver and `ss` for second. The single N325
second-receiver mixed-C entry is about2.4 SE from zero; its counterpart at425
does not resolve the same effect. The mixed E_ref and W entries also have no
common resolved pattern. Thus individual-source curvature is visible while a
joint-source nonadditivity claim awaits the full coordinated finite-vs-zero
comparison. Weak mixed entries are not evidence of exact additivity.

All old8 results, including every rank cell, are retained. They are not pooled
with New64 and do not become independent evidence from the same prefixes.
No pure/mixed ratio, inverse-covariance test, source search or sign selection is
introduced here.

## Exact density derivative and paired estimator

For each joint-safe contact-degree class, use its full count `n`, loop sums
`S_i`, product sums `Q_ij`, and `d=N-k0`. Both physical marks `L_i` are R0-only.
The first density score is `(n L_i-S_i)/d`; the second is

\[
t_{ij}(u)=\frac{(nL_i-S_i)(nL_j-S_j)-(nQ_{ij}-S_iS_j)}{d^2}.
\]

All numerator arithmetic is integer. The class score sums to zero, and the
score is zero outside joint safety. Its normalizer correction is essential:
using the score product alone would not be the second density derivative.
The mixed score is exactly zero outside original00; pure scores retain their
respective one-R0 cells.

The saved response averages, over8 or64 quartets separately,

\[
\frac12\,[t_{ij}(U)-t_{ij}(V)]\,[F_U-F_V],
\]

where each `F_label` is the mean of its two saved suffixes. Cross-class U,V and
safe/unsafe pairs remain valid under this exact score. `fs` denotes the actual
mixed derivative with **no factor2**; diagonal entries likewise have no Taylor
factor1/2. In a Taylor expansion the mixed term is `t_f t_s H_fs F`.

The fixed observables use `p_ref=.59274605079` and
`F_i=Pr(Binomial(N,p_ref)>=K_i)`:

- `A_ref=F1+F2-1`, `E_ref=1-F1+F2`;
- `C=(K1+K2)/(2(N+1))`, `W=(K2-K1)/(N+1)`.

The parameters are the specified commuting exponential-family source
coordinates. A nonzero mixed derivative would refute additivity in those
coordinates, not establish path memory, noncommuting perturbations or a field
identity. No finite-source weight or rectangle is evaluated by this reader.

## Immutable products and covariance interface

Sources are old8 NPZ `375cd3a12b2b7a87d79148a59f62b95898f9e471`, exact class
products `1cfa4ae892a2f7f4168e9a71690efd7a5560d4cd`, and New64 extension gzip
`8ad30617b0a3076a5c01a208eb213096d8879b32`. Allocation is `6bace935`; reader
`efcb015f` is `scripts/p334_mixed_source_curvature.py`.

`results/p334-mixed-source-curvature/` contains four prefix NPZs:
`old8_N{N}.npz`, `new64_N{N}.npz`. Each retains
`counter,batch,rankcell,old_rank,k0,d`, `mean_observable[prefix,receiver,4]`, and
`mean_response2[prefix,receiver,4,3]`. The last axes are
`A_ref,E_ref,C,W` and `ff,fs,ss`; values are means across8/64 quartets, not local
products or local determinants. These prefix means support reuse without
another raw pass.

`score.json` contains every physical tensor entry,264 named columns per size:
old8 `all` plus9 cells, and New64 `00` only. Every cell contribution sums over
its prefixes and divides by original1000 per batch, without a cell-mass
renormalization. Original20batch rows, estimates, SE, covariance and factors
are preserved. The sign conventions are explicit:

- `factor=(raw_batch-mean)/sqrt(20*19)`;
- `LOO=(20*mean-raw_batch)/19`, `LOO_factor=-factor`;
- `mean_covariance=factor.T@factor`.

Use the matching factor sign when joining the finite-source result; the two
sign conventions give the same individual covariance but opposite cross-block
products if mixed incorrectly. The coordinator receives all physical source,
receiver, endpoint and stream directions together.

The reader consumed no old fork gzip and each of the40 existing New64 gzip
files exactly once. No new sampling, tail replay, DP or prediction/shape test
was run. This extraction took1.77 seconds in the managed local Python research
environment; input and prefix-output hashes are recorded in `score.json`.
