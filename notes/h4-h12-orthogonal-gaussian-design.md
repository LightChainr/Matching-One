# Low-cost Gaussian design for H4 versus H12 discrimination

Status: prospective angular-design note.  This is a supplement to issue #36 and is intended to avoid using the expensive N=1105 four-angle run as the first H4/H12 discriminator.

## 1. Why doubling is not enough

The successful `theta -> theta+pi/4` Gaussian-doubling test proves the odd-`m` harmonic phase class

\[
\cos(4m(\theta+\pi/4))=(-1)^m\cos(4m\theta).
\]

It therefore distinguishes H4/H12/H20/... from H8/H16/... but does **not** by itself identify H4 over H12.

For a same-N pair define

\[
\Delta c_4=\cos4\theta_1-\cos4\theta_2,
\qquad
\Delta c_{12}=\cos12\theta_1-\cos12\theta_2.
\]

Since

\[
\cos12\theta=4\cos^3 4\theta-3\cos4\theta,
\]

the pair-dependent ratio

\[
r_{12}=\Delta c_{12}/\Delta c_4
\]

controls how an H12 admixture aliases into an H4-only fit.  We want moderate-N pairs with comparable H4 signal sizes but very different, preferably opposite, `r12`.

## 2. Selected pair A: N=305

Use orientation order

```text
first  = (17,4)
second = (16,7)
```

Exact harmonics:

\[
\Delta c_4=\frac{12672}{18605}
=0.6811072292394518\ldots,
\]

\[
\Delta c_{12}
=-\frac{187921184989824}{161001169878125}
=-1.1672038478482918\ldots,
\]

so

\[
\boxed{r_{12}=-\frac{14829638967}{8653650625}
=-1.713685889300621\ldots}
\]

## 3. Selected pair B: N=325

Use orientation order

```text
first  = (18,1)
second = (17,6)
```

Exact harmonics:

\[
\Delta c_4=\frac{16128}{21125}
=0.7634556213017751\ldots,
\]

\[
\Delta c_{12}
=\frac{25095083827968}{18129541015625}
=1.3842095509389745\ldots,
\]

so

\[
\boxed{r_{12}=\frac{1555994781}{858203125}
=1.813084496750114\ldots}
\]

The two H12/H4 alias ratios have nearly equal magnitude and opposite sign.

## 4. Frozen H4-only prediction

Using the already measured source coefficient

\[
A_4=0.7885\pm0.0352
\]

in

\[
\Delta M=A_4\Delta c_4 N^{-13/8},
\]

the no-H12 predictions are

```text
N=305: DeltaM = +4.9320781401e-5
       source-coefficient-only SE = 2.20176e-6

N=325: DeltaM = +4.9862612964e-5
       source-coefficient-only SE = 2.22595e-6
```

These source uncertainties do not include the future target Monte Carlo error.

The important design feature is not that the two predicted values are close; it is that an H12 term shifts them in opposite relative directions.

## 5. Two-column prospective model

After the H4-only frozen score, fit only the declared two-column model

\[
\Delta M_N=N^{-13/8}
\left(A_4\Delta c_{4,N}+A_{12}\Delta c_{12,N}\right).
\]

Because the two selected rows have opposite `r12`, the H4 and H12 columns are substantially better conditioned than for many arbitrary same-N pairs.

Do not select the pair after observing percolation outputs.  Do not add H8/H20 or radial powers until the frozen H4/H12 score is reported.

## 6. Statistics/hardware

The expected signal is only about `5e-5`, so fixed-p production will likely require several hundred million paired replicas per size for a decisive result.  Use a bounded pilot to freeze sample count from variance only.

A threshold-rank implementation is preferable if #49/#48 clean production is already available because it preserves full thermal information and lets the same run score derivative/homology models.

## 7. Acceptance

The H4 interpretation is strengthened if both target means agree with the frozen H4-only predictions and a shared H12 coefficient is unnecessary.

If the two residuals move coherently in the opposite directions predicted by their `r12` signs, fit one shared `A12` and test it on a later held-out pair with a third distinct `r12`.  Do not rescue a failure by letting every N have its own harmonic mixture.
