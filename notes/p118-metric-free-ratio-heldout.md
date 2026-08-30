# Issue 118: held-out metric-free ratio score

## Protocol chronology

The ratio definitions were frozen at `920c6393f6db7927887df0905dfaed81838ca062`
on 2026-08-29 13:36:57 +08:00.  The first commit containing the target N145/N290
full-curve result was `9675bce5b406247e15c03bca20abef954f26a3a2` at 14:47:18
+08:00, 1 hour 10 minutes 21 seconds later.

The execution manifest was locked in `29daa19`.  A single mistyped nibble in
the recorded N145 metadata digest was corrected in `76d1dd4`, before any P118
score was run.  All four raw/metadata digests then matched.  The scorer and its
synthetic/contract tests were pushed in `e162374` before the formal reveal.

## Frozen estimand

At each size and inside every delete-one replicate,

\[
R_I=\frac{P_4[D']}{P_4[S]\,\overline M'},\qquad
R_T=\frac{P_4[S']}{P_4[D]\,\overline M'}.
\]

The primary null was the two-vector equality
`(R_I,R_T)_N290 - (R_I,R_T)_N145 = 0`.  The two production streams are
independent, so their size-local nonlinear-jackknife covariance matrices were
added.  The joint score was evaluated before either marginal.

## Result

| coordinate | N145 | N290 | difference | SE of difference | marginal p |
|---|---:|---:|---:|---:|---:|
| R_I | 1.33555 | 1.87466 | 0.53911 | 1.20372 | 0.65425 |
| R_T | 3.23632 | 2.68309 | -0.55323 | 0.66851 | 0.40792 |

The frozen joint statistic is chi-square 0.857475 on 2 degrees of freedom,
with p=0.651331.  The constant-response null therefore **survives this block at
alpha=0.01**.  The observed child/parent ratios are 1.40366 for R_I and 0.82906
for R_T, but the full covariance does not resolve either displacement.

No effect-size or precision gate was preregistered, so the result is not
relabelled post reveal as “underpowered.”  Survival is not evidence that the
ratios are universal.  This score reuses the P50 raw and is therefore a
correlated held-out diagnostic rather than an additional independent primary
experiment.

