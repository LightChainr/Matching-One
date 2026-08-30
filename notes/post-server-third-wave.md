# Post-server third-wave interpretation

Status: research update after Huawei server commit `107ddc948c8ae8e5c8a30549c9b123ac9be091b7` on PR #21.

This note separates what is now empirically resolved from what remains a conjecture, and converts the new signal into a frozen follow-up program.

## 1. The orientation program has crossed a qualitative threshold

The first large-Pell fixed-`p` scan was underpowered. The follow-up results are materially different.

### Small Pell calibration

For the Pell pair `(a,d)=(7,5)`, two independent 5,000,000-replica runs gave

```text
+4.6267e-4 +/- 5.06e-5
+4.0417e-4 +/- 5.10e-5
```

for the diamond-minus-axis matching-root gap. The inverse-variance pool is

```text
+4.3365e-4 +/- 3.59e-5,
```

about twelve standard errors from zero. Halving the finite-difference step changed the result by less than `1e-6`, far below the sampling error.

Together with the exact tiny `(3,2)` pair, the two-point effective root-gap exponent is about `4.23`, and the `(7,5)` gap multiplied by the fourth power of the common physical scale is about `1.06`. This is not an asymptotic fit, but it is a successful positive power calibration for an `L^-4` orientation gap.

### Same-N Gaussian tomography

For exact square tori with periods

\[
(a,b),\qquad(-b,a),\qquad N=a^2+b^2,
\]

the 30,000,000-replica production run measured the matching-function orientation difference at one frozen `p_ref`:

| `N` | representations | `Delta M` | SE | z | `N^(13/8) Delta M / Delta cos(4 theta)` |
|---:|---|---:|---:|---:|---:|
| 65 | `(8,1)` vs `(7,4)` | `+1.00377e-3` | `1.68448e-4` | 5.96 | `0.65016 +/- 0.10911` |
| 85 | `(9,2)` vs `(7,6)` | `+7.60333e-4` | `1.58350e-4` | 4.80 | `0.65117 +/- 0.13562` |
| 145 | `(12,1)` vs `(9,8)` | `+3.28667e-4` | `1.95415e-4` | 1.68 | `0.55743 +/- 0.33143` |

All signs agree with `Delta cos(4 theta)`. The first two scaled amplitudes agree to about `0.2%`, while their Monte Carlo errors are of order `20%`. A fixed-amplitude inverse-variance combination of all three points gives, as a descriptive statistic only,

\[
A_4^{\rm MC}=0.645\pm0.082,
\]

with a very small internal chi-square. The pairwise exponent inferred from `N=65,85` after dividing out `Delta cos(4 theta)` is

\[
\alpha_{65,85}=1.619,
\]

close to the preregistered `13/8=1.625` power for `M(p_c)`. Neither observation is yet an asymptotic determination, but the joint agreement is too structured to treat as the previous null result.

The appropriate current statement is:

> A same-area, same-shape, threshold-insensitive orientation signal has been resolved in the square-site matching function at `N=65` and `N=85`, with the sign and first two size amplitudes predicted by a leading `cos(4 theta) N^-13/8` term.

## 2. Sharpened working hypothesis

For fixed square-torus modulus and microscopic orientation `theta`, write

\[
M_N(p_c;\theta)
=N^{-13/8}
\left[A_0+A_4\cos(4\theta)
+A_8\cos(8\theta)+\cdots\right]
+\text{subleading terms}.
\]

For two orientations at identical `N`, the scalar term cancels:

\[
\Delta M_N
=M_N(\theta_1)-M_N(\theta_2).
\]

The primary frozen model for confirmation is

\[
\boxed{
\Delta M_N
=A_4\,\Delta\cos(4\theta)\,N^{-13/8}
\left[1+bN^{-\omega/2}+\cdots\right].
}
\]

A logarithmic alternative must also remain live:

\[
\Delta M_N
=\Delta\cos(4\theta)\,N^{-13/8}
\left[A_4+B_4\log N+\cdots\right].
\]

The root difference should then scale as

\[
\Delta p_N^*\sim N^{-2},
\]

because `M'_N(p_c) ~ N^(3/8)`.

A natural operator-level candidate remains a matching-odd, spin-4 field of total dimension

\[
x=21/4,
\]

possibly a level-four thermal-family descendant/quasiprimary in the `c=0` logarithmic theory. This remains a conjecture until both the angular harmonic and radial exponent are independently resolved.

## 3. A previous gate is now rejected

The production matching-even orientation differences did not dominate the matching-odd sector. Their z-scores were approximately `4.19`, `0.18`, and `1.66`, with no clean amplitude collapse.

Therefore the earlier planning claim

```text
resolve a large matching-even L^-2 spin-4 sector first,
then spend statistics on a small matching-odd residue
```

is not supported as a required gate. The matching-even sector remains scientifically useful, but a null or unstable even sector must not invalidate the directly resolved matching-function difference.

## 4. Important dependence and topology cautions

### `either` and `cross` are not independent confirmations of `M`

For the complementary primal/matching construction on the validated periodic quotients, the `either` and `cross` matching differences are configuration-level identical. Their matching-function outputs therefore duplicate the same random variable. They are useful topology regressions, not independent replications.

Independent evidence must come from another seed/counter range, another size/orientation, another topology implementation, or another microscopic model.

### Primitive versus nonprimitive Gaussian representations

The cyclic fast path

\[
\mathbb Z^2/\langle(a,b),(-b,a)\rangle\cong\mathbb Z_N
\]

requires `gcd(a,b)=1`. A nonprimitive axis representation such as `(13,0)` at `N=169` has quotient group `Z_13 x Z_13`, not `Z_169`. It is still a valid exact square torus for a general period-matrix engine, but it cannot share the same canonical cyclic labeling with primitive `(12,5)`.

Consequently, `N=169` is a useful magic-orientation comparison, but not a zero-caveat same-cyclic-field pair. Report the group/Smith-normal-form difference explicitly.

### High-precision finite-size scorer

`scripts/summarize_finite_size_grid.py` currently parses decimal observations before setting `mp.mp.dps` to the selected precision. This silently rounds the input at the previous mpmath precision. The scientific effect on the existing `1e-11` tail score is tiny, but the implementation must set precision before reading the CSV and add a regression test using more than 30 significant digits.

## 5. Highest-information next sizes

The original confirmation list `(205,425)` is not optimal for discovery. Under the observed amplitude `A_4 ~ 0.645`, candidate same-N pairs have approximate expected signals:

| `N` | pair | `Delta cos4` | expected `Delta M` |
|---:|---|---:|---:|
| 130 | `(11,3)` vs `(9,7)` | 1.363 | `3.2e-4` |
| 145 | `(12,1)` vs `(9,8)` | 1.918 | `3.8e-4` |
| 170 | `(13,1)` vs `(11,7)` | 1.594 | `2.4e-4` |
| 185 | `(13,4)` vs `(11,8)` | 1.178 | `1.6e-4` |
| 221 | `(14,5)` vs `(11,10)` | 1.179 | `1.2e-4` |
| 265 | `(16,3)` vs `(12,11)` | 1.722 | `1.3e-4` |
| 205 | `(14,3)` vs `(13,6)` | 0.822 | `9.3e-5` |
| 425 | `(19,8)` vs `(16,13)` | 0.893 | `3.1e-5` |

Run `N=130` and `N=170` before `N=205` or `N=425`. They extend the radial range while retaining substantially larger angular leverage.

## 6. Frozen confirmation protocol

### Stage A: independent replication

Run a new production seed/counter range, with no pilot retuning, for

```text
N = 65, 85, 130, 145, 170
```

using the same `p_ref`, orientation order, cyclic labeling, and topology definitions. At least `100,000,000` replicas per listed pair are justified by the current throughput; retain 100 or more equal batches.

Report each seed separately before pooling. The primary statistics are

\[
\Delta M_N,
\qquad
A_4(N)=
\frac{N^{13/8}\Delta M_N}{\Delta\cos4\theta}.
\]

Do not select or drop sizes by observed z-score.

### Stage B: radial model comparison

Freeze the following models before viewing Stage A confirmation data:

1. `alpha=13/8`, constant amplitude;
2. `alpha=13/8` plus one declared inverse-power correction;
3. `alpha=13/8` plus a logarithmic amplitude term;
4. one free-`alpha` model trained only on the smaller sizes.

Score held-out size prediction, signed residuals, conditioning, and amplitude drift. A visual straight line is insufficient.

### Stage C: thermal-coordinate check

The current result is at one `p_ref`. Use the same random fields or, preferably, threshold-rank Newman-Ziff data to measure the orientation difference at symmetric nearby coordinates. Confirm that the central amplitude is not produced by an orientation-dependent slope multiplied by a small `p_ref-p_c` error.

### Stage D: multi-angle confirmation

Only after the two-angle radial test succeeds should the expensive `N=1105` four-orientation regression be attempted. It must use leave-one-orientation-out prediction between `constant`, `constant+cos4`, and `constant+cos4+cos8` models.

## 7. Kappa3 interpretation

The square-bond exact-threshold sequence is a valid method control:

```text
L=4   -1.56230 +/- 0.01483
L=6   -1.60727 +/- 0.02169
L=8   -1.57106 +/- 0.02584
L=12  -1.65012 +/- 0.03406
L=16  -1.68356 +/- 0.03909
```

A fixed `L^-3/2` model trained through `L=12` predicts `L=16` within about `1.55` standard errors. The all-size intercept near `-1.649 +/- 0.023` is compatible with `-5/3`, but this is not a sharp test because the higher-size score estimator is noisy and only one microscopic realization has been implemented.

The next kappa3 resources should go to:

1. a second independent seed with sample counts powered per size;
2. triangular-site and union-jack exact-threshold controls;
3. covariance-aware derivative estimation or microcanonical reconstruction;
4. a coordinate-free parametric `M(U)` comparison.

Do not spend the next large run merely extending square-bond `L` while holding one million samples fixed.

## 8. Control-variate program must be rewritten

The five wrapping-difference channels are configuration-level identical on the tested symmetric quotients, so wrapping-only GLS cannot reduce variance. Ridge regularization may produce numerical weights, but not new information.

Use the exact Euler decomposition instead. With occupied-site count `V`, occupied NN edges `E`, occupied elementary faces `F_0`, and topological variable `q`:

\[
D_N=q+(V-Np)-(E-2Np^2)+(F_0-Np^4).
\]

The centered local quantities are exact zero-mean controls. Extend them with prespecified local motif counts and, in the microcanonical ensemble, exact hypergeometric conditional expectations.

The new target is an Euler/motif control-variate hierarchy evaluated by pilot-frozen or cross-fitted weights, not a combination of duplicate wrapping indicators.

## 9. Decision

The orientation hypothesis is promoted from exploratory to **confirmation-stage P0**.

The decisive object is no longer an axis/diamond sign anecdote. It is the jointly constrained law

\[
\boxed{
\Delta M_N
\propto
\Delta\cos(4\theta)\,N^{-13/8}
}
\]

with independent-seed, multi-size, thermal-coordinate, and eventually multi-angle tests.

A successful confirmation would explain the empirical `L^-4` matching-root bias through a rotation-odd matching sector, not merely improve a decimal estimate of `p_c`.