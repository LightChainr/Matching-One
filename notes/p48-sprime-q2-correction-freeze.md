# P48 `P4[S']` correction freeze: choose the minimal `q=2` RG term

Status: **retrospective model-development note**. This note is written after the P48/P33 retrospective spectrum is known and before fresh `N=185,265` target statistics are to be scored. Nothing below is confirmatory evidence.

Frozen machine-readable artifact: `predictions/p48_sprime_q2_correction_20260828.yaml`.

## 1. What failed in the pure-power package

P48 tests the matching-parity derivative spectrum. For the matching-odd thermal spin-4 field (`x=21/4`, `eta=-1`) the first matching-even derivative is predicted to scale as

```text
P4[S'] ~ N^-5/4.
```

On the retrospective P33 threshold-rank data, the original training split `N=65,85,130` predicted held-out `N=145,170` poorly in this channel:

```text
pure N^-5/4 held-out chi2 = 10.1908 / 2.
```

The scaled values `Y_N=N^(5/4) P4[S']` visibly rise over the available range:

```text
N=65   1.82726 +/- 0.09526
N=85   2.05521 +/- 0.11784
N=130  2.70212 +/- 0.31739
N=145  2.45593 +/- 0.26371
N=170  3.09606 +/- 0.43984
```

The alternating parity signal pattern remains useful, but the four pure powers cannot be declared a conjunction pass.

## 2. Minimal symmetry-allowed correction

Write the leading matching-odd spin-4 contribution as `T4`. A matching-even scalar irrelevant coupling `S0` with length correction exponent `omega=2` preserves both the matching parity and H4 angular character of `T4`:

```text
T4 * S0 : eta = (-1)*(+1) = -1,
           H4 * H0 -> H4,
           relative correction L^-2 = N^-1.
```

The same `omega=2` local-RG scale is independently singled out by the percolation subleading thermal eigenvalue `y_t2=-2` and by the post-annihilator bookkeeping in #47/#58. Therefore the smallest ordinary correction to test is

```text
P4[S'] = N^-5/4 (A + B/N).
```

This is not a new exponent fit. The leading `5/4` stays fixed and the correction exponent is fixed to the already-motivated `q=2` length power.

## 3. Retrospective GLS diagnostic

Using all five P48 retrospective sizes and the committed synchronized delete-one-batch covariance gives

```text
A =  3.203310807356976
B = -90.59560328584558
chi2 = 1.89235 / 3.
```

For comparison, a constant scaled amplitude (`B=0`) gives

```text
chi2 = 17.0871 / 4.
```

More importantly for the old development split, fitting `A+B/N` only on `N=65,85,130` and predicting `145,170` gives

```text
held-out chi2 = 1.12494 / 2,
```

versus `10.19081/2` for the pure power.

These are **selection diagnostics on already seen data**. They justify freezing a candidate; they do not count as a prospective success.

## 4. Why not freeze the logarithm as the one secondary model?

A same-power logarithmic form

```text
Y_N = A + B log N
```

also describes these five retrospective points well (and is numerically slightly better on this tiny sample). That is expected: over `N=65..170`, `1/N`, `1/sqrt(N)`, and `log N` are highly confounded once two amplitudes are allowed.

The choice of `A+B/N` is therefore **mechanistic, not data-ranked**:

1. `q=2` has a concrete matching-even scalar mixing mechanism that preserves the required H4/matching-odd quantum numbers.
2. The three prospective fixed-p Gaussian-doubling lineages currently do not require a coherent logarithmic/Jordan residual in the central matching-odd H4 channel: the first two pass the pure-power residual jointly, and the third N145->290 lineage also passes separately. This does not exclude a derivative-specific logarithm, but it gives no independent reason to make the logarithm the first correction.
3. Freezing exactly one ordinary correction now prevents choosing among `1/N`, `log N`, or a free exponent after seeing `N=185,265`.

If the fresh derivative data reject this `q=2` model, the next step should be a separately frozen logarithmic/Jordan alternative, not a retuned correction exponent.

## 5. Frozen fresh-target predictions

With the full retrospective GLS coefficients frozen above:

```text
N=185:
  Y = N^(5/4) P4[S'] = 2.71360484365 +/- 0.19668994
  P4[S']             = 0.003977241344 +/- 0.000288281975

N=265:
  Y = N^(5/4) P4[S'] = 2.86144060628 +/- 0.23233037
  P4[S']             = 0.002676254395 +/- 0.000217294455
```

The quoted uncertainties are source-fit uncertainty only. Fresh target sampling covariance must be added during scoring.

## 6. Scoring order

Do not rewrite the original P48 primary after this retrospective discovery. On fresh `N=185,265` data:

1. score the original four pure-power P48 laws unchanged;
2. report the zero-effect benchmark;
3. score this frozen `q=2` correction for `P4[S']` without refitting `A` or `B`;
4. only then consider a separately preregistered log/Jordan or free-exponent model.

A q=2 pass would support an ordinary correction to the thermal spin-4 sector, not uniquely identify its CFT operator. A failure of both the original and q=2 laws is evidence to enlarge the operator/logarithmic model, not permission to retune this artifact.
