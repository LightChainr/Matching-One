# Minimal two-sector correction for the full-curve doubling slope

Status: prospective model for the unrevealed full-curve `145 -> 290` lineage. It uses only the already committed clean P49 `65 -> 130` and `85 -> 170` source data.

## 1. The raw exponent is not the main problem

The clean 100M P49 score resolves

```text
65 -> 130 slope ratio = 1.2939835119594856
85 -> 170 slope ratio = 1.294377573998004
```

against the asymptotic thermal target

\[
r_\infty=2^{3/8}=1.2968395546510096\ldots
\]

with very high precision. The point deficit is only about two parts in one thousand, so the natural question is the leading finite-size correction to the thermal metric, not a free replacement of `y_t=3/4`.

## 2. One scalar `1/N` term is too simple

Writing

\[
\bar M'_N=B_\infty N^{3/8}(1+c/N+\cdots)
\]

and solving each lineage independently gives effective values near `c=0.288` and `c=0.324`. At the committed P49 precision a single shared `c` is not adequate.

This is not surprising because P48/P49 also resolve an orientation-dependent `D'` sector. Its absolute scaling is `N^-5/8`, exactly one factor `1/N` relative to the universal slope `N^(3/8)`.

## 3. Minimal scalar + H4 model

Let

\[
\bar c_4=\frac{\cos4\theta_1+\cos4\theta_2}{2}.
\]

Use

\[
\boxed{
\bar M'_N(\bar c_4)
=B_\infty N^{3/8}
\left[1+\frac{a+b\bar c_4}{N}+\cdots\right].
}
\]

The `a` term is the smallest ordinary scalar `L^-2=N^-1` correction. The `b` term is the smallest correction compatible with the already observed H4 derivative sector.

Gaussian multiplication by `1+i` sends

\[
N\to2N,\qquad \bar c_4\to-\bar c_4.
\]

Therefore

\[
\frac{\bar M'_{2N}}{\bar M'_N}
=2^{3/8}
\frac{1+(a-b\bar c_4)/(2N)}
     {1+(a+b\bar c_4)/N}.
\]

For a measured ratio `r`, this gives a linear equation in `(a,b)` after dividing by `2^(3/8)`:

\[
(q-1/2)a+(q+1/2)\bar c_4 b=N(1-q),
\qquad q=r/2^{3/8}.
\]

The two clean source lineages thus determine `(a,b)` without a nonlinear fit.

Exact parent mean harmonics are

```text
N=65:   cbar4 =  833/4225
N=85:   cbar4 = -1127/7225
N=145:  cbar4 = -287/21025
```

The central source values give approximately

\[
a=0.30789201,\qquad b=-0.03426210.
\]

## 4. Frozen third-lineage prediction

For the unrevealed full-curve lineage

```text
N=145: (12,1)-(9,8)
        |
        | multiply by 1+i
        v
N=290: (13,11)-(17,1)  [Gaussian lineage order]
```

the frozen slope target is

\[
\boxed{
\bar M'_{290}/\bar M'_{145}=1.2954593652984558.
}
\]

A first-order source-covariance propagation gives an approximate source SE `1.81e-5`. The final scorer must recompute the full delete-one transformation rather than treat this approximation as publication-grade uncertainty.

If the primary matching residual continues to obey

\[
\Delta M_{290}/\Delta M_{145}=-2^{-13/8},
\]

then the same model predicts the root-gap ratio

\[
\boxed{
\Delta p^*_{290}/\Delta p^*_{145}
=-0.2502663513402128.
}
\]

## 5. Interpretation

A pass would support a compact decomposition of the resolved slope drift into:

- the universal thermal eigenvalue `y_t=3/4`;
- a relative scalar `L^-2` correction;
- an H4 derivative correction already visible in the P48 parity spectrum.

It would **not** uniquely identify the scalar correction operator. A failure should trigger a model with additional radial powers/harmonics using a new training lever, not a free exponent fit on the same two source lineages.
