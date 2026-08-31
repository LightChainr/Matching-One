# Exact-score quartet moments: the three levels of birth-clock covariance

The existing archive now supplies eight independent quartet estimates per
original prefix, including the same-label independent-suffix cross-products.
This is enough to separate covariance changes within a fixed next label,
between next-label means within a prefix, and between prefix means. No new
simulation, tail replay, DP, finite-policy scan, determinant, or covariance
significance calculation is performed here.

## Source and exact score

- Forks: `e32a85939279b8574278024d647b56d2d1485247`.
- Contacts: `959a7fa26677c416b874d272f1ba66523fb38f73`.
- Full vacant-label census and score definition: PR509 head
  `ac5761ce504c3cd170fa42c86c17d6fb87f0375b`,
  `experiments/p334-finite-source-20260831/analyze_finite_source.py`.
- Reader: `ad83bb6c`, `scripts/p334_exact_score_quartet_moments.py`.

Within a prefix let `d=N-k0`. A joint-safe label belongs to its contact-degree
class `a=(e_first,e_second)`, with census count `n_a`. For either common
half-sum/difference loop mark `g`, the exact centered score is

\[
H(u)={n_a\over d}\bigl[g(u)-\overline g_a\bigr].
\]

It is zero outside joint safety. The implementation preserves the integer
numerator `2*d*H`; class masses and mark means are obtained from the entire
census, not estimated from the two sampled labels. For two iid labels U,V,

\[
h_f={1\over2}[H(U)-H(V)][f_U-f_V]
\]

is the requested response estimator. It targets the same tangent as the old
matched-class pair mask, but does **not** reimpose that mask on U,V: cross-class
and safe/unsafe pairs can contribute to the exact-score estimator. Reimposing
the old pair mask would change this estimator and discard valid terms.

## Eight features, with no extra averaging of the last three

Set `x=K1/(N+1)` and `y=K2/(N+1)` for each saved suffix. The feature order is

`x, y, xx, xy, yy, cross_xx, cross_xy, cross_yy`.

For the first five, `f_U` averages the two suffix values `x,y,x*x,x*y,y*y`.
For the final three it is directly

\[
x_0x_1,\qquad (x_0y_1+x_1y_0)/2,\qquad y_0y_1.
\]

The baseline `b=(f_U+f_V)/2` thus averages all four tails for the first five,
and averages the two label cross-products for the last three. Given one label,
the cross-products estimate products of its conditional suffix means.
Consequently `b_xy-b_cross_xy` and `h_xy-h_cross_xy` expose the same-label suffix
covariance and its response without estimating a label-specific mean. To remove
prefix mean products, use **different** quartets; do not multiply a quartet's
own baseline and response and call it a conditional covariance derivative.

All nine prefix rank cells and their baselines are retained. The two scores
vanish outside the five cells `00,01,02,10,20`. No cell prevalence is divided
out. The original twenty batches remain the dependence blocks.

## Continuous order-statistic recovery remains exact

The stored `xx,xy,yy` follow this task's normalized integer-clock convention.
For the continuous common uniform order statistics, either baseline or tangent
transforms without another source pass as

\[
xx_{\tau}={(N+1)xx+x\over N+2},\quad
xy_{\tau}={(N+1)xy+x\over N+2},\quad
yy_{\tau}={(N+1)yy+y\over N+2}.
\]

The cross-products of independent suffixes remain unchanged. This preserves
the option to include the Beta/order-statistic sampling variance rather than
silently conflating the two conventions. All actual integer birth clocks are
also retained compactly, allowing later coordinate changes without CSV replay.

## Data interface

`results/p334-exact-score-quartet-moments/metadata.json` gives the full schema,
input/output SHA256, exact source commits, formulas, dtypes, and signed H4
orientation factors. Each `N325.npz` / `N425.npz` has:

- `b`: float64 `[20000,8,2,8]` = prefix,quartet,orientation,feature.
- `h`: float64 `[20000,8,2,2,8]` = prefix,quartet,mark,orientation,feature.
- `score`: float64 `[20000,8,2,2]` = prefix,quartet,label,mark, plus its exact
  integer `score_numerator` with denominator `2*d`.
- `batch,counter,k0,d,rankcell,old_rank` and `next_label`.
- `birth_k`: uint16 `[20000,8,2,2,2,2]` =
  prefix,quartet,label,suffix,orientation,birth.
- `joint_safe,contact_e,contact_c,r0_loop,rank_after` and full per-prefix
  `census_class_count,census_class_loop_sum` for score/geometry reuse.

Each source CSV was consumed once. Compression produced 18.284 MB (N325) and
19.584 MB (N425); both sizes together took 5.04 seconds locally. The managed
Python research environment was used with single-thread BLAS. The files add no
independent evidence to their common raw block and perform no new population
or model comparison.
